%% ---------------- Evaluate Objectives and Constraints -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: evaluateObjectivesAndConstraints
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       The core interface for the optimizer. It calculates the Cost Function
%       (Objective) and the Non-Linear Constraints (c, ceq) for a given
%       design vector 'x'.
%
%   Features:
%       - Memoization: Caches results to avoid re-calculating physics 
%         separately for Objective and Constraints calls by the optimizer.
%       - Mode Handling: Switches between 'GA' (Relaxed constraints) and 
%         'SQP' (Strict constraints) logic.
%       - Penalties: Handles failed BEMT convergence with soft penalties.
%
%   Inputs:
%       x       - Design vector (12x1 double)
%       mission, prop, motor, env, opt - Configuration structs
%       mode    - 'GA' or 'SQP'
%
%   Outputs:
%       J       - Cost function (Weighted negative total efficiency)
%       c       - Inequality constraints vector
%       ceq     - Equality constraints vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J, c, ceq] = evaluateObjectivesAndConstraints(x_norm, mission, prop, motor, env, opt, mode)

    % 0. MEMOIZATION (CACHE)
    % For FitnessFunction, ObjectivesAndConstraints is directly used, but
    % for Constraints, ConstraintsInterface is used.
    % That means two consecutive calls to ObjectivesAndConstraints for each
    % x the optimizer tries.
    % By using persistent variables, we allow the function to "remember" 
    % the last calculation. If the optimizer calls it again with the same 
    % x, it simply returns the stored result instantly instead of 
    % recalculating.
    % This reduces optimization time to half.

    % No noticeable increase in speed was observed by using persistent
    % variables. That is because the optimization uses Parallel Pool, but
    % different workers can not share Cache.
    %
    % % Persistent variables stay in memory between function calls
    % persistent last_x last_mode last_J last_c last_ceq
    % 
    % % Check if we have seen this x AND this mode before
    % % (We check 'mode' to ensure we don't return GA constraints to SQP)
    % if ~isempty(last_x) && isequal(x, last_x) && strcmpi(mode, last_mode)
    %     J   = last_J;
    %     c   = last_c;
    %     ceq = last_ceq;
    %     return; % <--- EXIT EARLY (Skip all physics)
    % end

    % 1. Undo the normalization of the design variables:
    x = x_norm .* (opt.ub - opt.lb) + opt.lb;

    % 2. Input Handling
    if nargin < 6
        mode = 'SQP'; % Default to strict mode if not specified
    end

    % 3. Generate Geometry from Design Vector
    % x(1:7) define the twist and chord distributions
    % x(8) is the dimensional radius (NOT normalized): R = x(8);
    [R, c_dist, theta_deg] = getPropGeometry(x, prop);

    % 4. Analyze Mission Segments
    % We loop through the 4 flight conditions defined in the Mission Profile.
    % x(9:12) are the Normalized RPMs for these 4 segments.
    num_segments = numel(mission.v_ref);
    geometry_indices = numel(x)-num_segments;
    objs = zeros(1, num_segments);
    results_for_con = cell(1, num_segments);
    success_flags = false(1, num_segments);

    % Indices in 'x' for RPMs: Slow(9), Cruise(10), Climb(11), Fast(12)
    for i = 1:num_segments
        v_current = mission.v_ref(i);
        rpm_current = x(geometry_indices + i);
        
        % Call the single-point evaluator
        [~, ~, objs(i), results_for_con{i}, success_flags(i)] = evaluateModels(v_current, rpm_current, R, ...
                                                          c_dist, theta_deg, prop, motor, env, opt);
    end

    % 5. Fail-Safe Check
    % If BEMT failed for any segment (complex numbers, NaN), set as
    % unfeasible solution and with massive penalty.
    if ~all(success_flags)
        J_nominal = 10; % Big penalty, so that the optimizer never prefers an unfeasible rather than a bad feasible.
                        % The objective is negative and normalized (eta_mp in this case)

        penalty = 0;
        penalty = penalty + sum(~success_flags) * 5;
        J = J_nominal + penalty; % Optimizer now sees directional improvement towards realistic solutions

        % We count how many flags are set to '1' in the settings
        
        % 1. Count Physical Constraints (Inequalities)
        % (Summing the logic vectors gives the total number of active checks)
        n_phys = sum(opt.constr_flags.current) + ...
                 sum(opt.constr_flags.voltage) + ...
                 sum(opt.constr_flags.M_tip) + ...
                 sum(opt.constr_flags.stress);
        
        % 2. Count Thrust Constraints
        n_thrust = sum(opt.constr_flags.thrust);

        % --- ASSIGN PENALTY VECTORS ---
        if strcmpi(mode, 'GA')
            % GA Mode: 
            % - Thrust is an Inequality Band (Lower + Upper) -> 2 inequality constr per active thrust
            % - Physical are Inequalities
            total_ineq = (n_thrust * 2) + n_phys;
            c   = ones(1, total_ineq)*10; % Setting c to >0 is unfeasible
            ceq = [];
        else
            % SQP Mode:
            % - Thrust is an Equality -> 1 equality constr per active thrust
            % - Physical are Inequalities
            c   = ones(1, n_phys)*10;   % Physical penalty
            ceq = ones(1, n_thrust)*10; % Thrust penalty
        end
        
        % % Save to cache (implemented memoization) before returning
        % last_x = x; last_mode = mode; last_J = J; last_c = c; last_ceq = ceq;
        
        return;
    end


    % 6. Calculate Objective Function
    % Maximize Weighted Average Efficiency (Minimize Negative)
    % Note: objs from EvaluateModels are already normalized and negative.
    J = objs * mission.weights(:); % Dot product


    % 7. Calculate Constraints
    % Initialize empty arrays
    c_phys = [];      % Physical Inequalities
    c_thrust_ga = []; % Thrust Inequalities (Band for GA)
    ceq_thrust = [];  % Thrust Equalities (Strict for SQP)
    
    % Loop through each flight segment (1 to 4)
    for i = 1:num_segments
        res = results_for_con{i};
        
        % --- A. PHYSICAL CONSTRAINTS (Inequalities) ---
        
        % 1. Motor Current
        if opt.constr_flags.current(i)
            c_phys = [c_phys, (res.I_mot / motor.I_max) - 1]; 
        end
        
        % 2. Battery Voltage
        if opt.constr_flags.voltage(i)
            c_phys = [c_phys, (res.V_req / motor.V_batt) - 1];
        end
        
        % 3. Mach Number
        if opt.constr_flags.M_tip(i)
            c_phys = [c_phys, (res.M_tip / prop.M_tip_max) - 1];
        end
        
        % 4. Structural Stress
        if opt.constr_flags.stress(i)
            c_phys = [c_phys, (opt.constr.stress_sf * res.max_sigma / prop.sigma_p) - 1];
        end
        
        % --- B. THRUST CONSTRAINTS (Mode Dependent) ---
        
        if opt.constr_flags.thrust(i)
            T_calc = res.T;
            T_targ = mission.T_ref(i);
            
            if strcmpi(mode, 'GA')
                % GA Mode: Band Constraint (+/- (tol*100) %)
                tol = opt.constr.GA_thrust_tol;
                % Lower Bound
                c_thrust_ga = [c_thrust_ga, ((T_targ*(1-tol)) - T_calc)/T_targ];
                % Upper Bound
                c_thrust_ga = [c_thrust_ga, (T_calc - (T_targ*(1+tol)))/T_targ];
            else
                % SQP Mode: Strict Equality
                ceq_thrust = [ceq_thrust, (T_calc - T_targ)/T_targ];
            end
        end
    end

    % 8. Final Assembly
    if strcmpi(mode, 'GA')
        thrust_error_ratio = 0;
        for i = 1:num_segments
            thrust_error_ratio = thrust_error_ratio + abs(results_for_con{i}.T - mission.T_ref(i)) / mission.T_ref(i);
            % Note that max_thrust_error_ratio = tol + tol + tol + tol = 4*tol
        end
        J = J + 0.1 * thrust_error_ratio; % Drive GA toward thrust targets naturally (not the whole tol band)
        % Combine all inequalities
        c = [c_thrust_ga, c_phys];
        ceq = [];
    else
        % SQP: Physical are Inequalities, Thrust are Equalities
        c = c_phys;
        ceq = ceq_thrust;
    end

end