%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script: run_optimization.m
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Executes the two-stage optimization process defined by the settings
%       in 'main_runner.m'.
%
%   Workflow:
%       1. Pre-calculation: Prepares aerodynamic interpolants (Cl, Cd vs 
%          Alpha/Re) and structural properties to speed up the loop.
%       2. Global Search (Genetic Algorithm): Finds a global optimum region 
%          using relaxed constraints (Thrust bands).
%       3. Local Refinement (SQP - fmincon): Refines the GA solution using 
%          strict equality constraints for target thrust.
%
%   Outputs:
%       - Saves 'Results/GA_opt_result.mat'
%       - Saves 'Results/SQP_final_result.mat'
%
%   Dependencies:
%       - Global structs: mission, prop, motor, env, opt
%       - evaluateObjectivesAndConstraints.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           OPTIMIZATION                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Expects 'mission', 'prop', 'motor', 'env', 'opt' structs to exist.

if ~exist('mission','var')||~exist('prop','var')||~exist('motor','var')||~exist('env','var')||~exist('opt','var')
    error('Configuration variables not found. Please run main_runner.m before.'); 
end

% Ensure Results folder exists
if ~exist('results', 'dir')
    mkdir('results');
end


%% ========================================================================
%  1. PRE-CALCULATE AERODYNAMICS AND INERTIAS TO MAKE SABEMMT RUN FASTER
%  ========================================================================

% Aerodynamics (Load & Interpolate)
airfoil_data = prop.airfoil_data; % Unpack for easier access

% 'linear' ensures smooth transitions between Reynolds numbers
Fc = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cl_matrix, 'linear', 'linear');                      
Fd = griddedInterpolant({airfoil_data.alpha_deg*(pi/180), airfoil_data.Re}, ...
                        airfoil_data.cd_matrix, 'linear', 'linear');

% Linear Aerodynamics (for SABEMMT init)
nom_idx = find(airfoil_data.Re >= 0.6 * max(airfoil_data.Re), 1);
alpha_lin = airfoil_data.alpha_deg * (pi/180);
cl_lin = airfoil_data.cl_matrix(:, nom_idx);
lin_region = alpha_lin >= deg2rad(-1) & alpha_lin <= deg2rad(6);
p_poly = polyfit(alpha_lin(lin_region), cl_lin(lin_region), 1);
cla_data = p_poly(1); cl0_data = p_poly(2);

% Airfoil geometry and inertias (Structural Properties)
geo = getAirfoilInertias(airfoil_data.x_u, airfoil_data.y_u, ...
                             airfoil_data.x_l, airfoil_data.y_l);

% Pack values in the prop struct
prop.airfoil.Fc = Fc; prop.airfoil.Fd = Fd;
prop.airfoil.cla_data = cla_data; prop.airfoil.cl0_data = cl0_data;
prop.airfoil.alpha_min = min(alpha_lin); prop.airfoil.alpha_max = max(alpha_lin);
prop.airfoil.cd_max_val = max(airfoil_data.cd_matrix(:));
prop.structural = geo;



%  ========================================================================
%  2. OPTIMIZATION
%  ========================================================================
fprintf('\n================ INITIALIZING OPTIMIZATION ================\n');

% The current problem is identified by:
signature = struct();
signature.mission = mission;
signature.prop    = prop;
signature.motor   = motor;
signature.env     = env;
signature.opt     = opt;


%% ========================================================================
%  2.1. STEP 1: GENETIC ALGORITHM (Global Search)
%  ========================================================================
%  Uses 'GA' mode: Relaxed constraints (Thrust within +/- tol %) to find
%  a feasible region without getting stuck on strict equality constraints.

% Check for existing GA result of the SAME problem to save time
GA_file = 'results/GA_opt_result.mat';
if exist(GA_file, 'file') == 2
    tmp = load(GA_file, 'x_norm_opt_GA', 'fopt_GA', 'signature');
    if isfield(tmp, 'signature') && isequal(tmp.signature, signature)
        fprintf('\n  > GA results file found and with SAME problem values → skipping GA\n');
        pause(3); % To let the user read the message
        run_GA = false;
        x_norm_opt_GA = tmp.x_norm_opt_GA; % Assign the previous value
        fopt_GA  = tmp.fopt_GA;
    else
        fprintf('\n  > GA results file found but PROBLEM VALUES HAVE CHANGED\n');
        fprintf('\n    You can run GA for the new problem and then run SQP from there,\n');
        fprintf('\n    or you can run SQP for the new problem but from the saved GA result.\n');
        reply = input('\n  > Run GA for the new problem? (y/n): ','s');
        if strcmpi(reply,'y')
            run_GA = true;
        else
            fprintf('  > Moving on to SQP starting from the previous GA solution saved\n');
            pause(2); % To let the user read the message
            run_GA = false;
            x_norm_opt_GA = tmp.x_norm_opt_GA; % Assign the previous value
            fopt_GA  = tmp.fopt_GA;
        end
    end
else
    fprintf('\n  > GA results file NOT found → running GA \n');
    pause(2); % To let the user read the message
    run_GA = true;
end


if run_GA
    fprintf('\n------- Starting Genetic Algorithm (Global Search) --------\n');
    pause(2);

    % Cost Function (returns J only)
    FitFcn_GA = @(x_norm) evaluateObjectivesAndConstraints(x_norm, mission, prop, motor, env, opt, 'GA');
    
    % Non-Linear Constraints (returns [c, ceq])
    % We use a helper 'constraintsInterface' to discard J and return c, ceq
    Con_GA = @(x_norm) constraintsInterface(x_norm, mission, prop, motor, env, opt, 'GA');

    % --- GA Options ---
    % Bigger GA ≠ better GA for constrained aero problems.
    options_GA = optimoptions('ga', ...
        'PopulationSize', 700, ...        % Big PopulationSize for better search space coverage
        'MaxGenerations', 2500, ...       
        'MaxStallGenerations', 200, ...   % Wait 200 stall gens (default is 50) before stopping if no improvement
        'FunctionTolerance', 1e-6, ...    % Strict tolerance to force refinement
        'MaxTime', 4*3600, ...
        'UseParallel', true, ...
        'Display', 'iter', ...
        'PlotFcn', {@gaplotbestf, ...     % Plots best and mean fitness
                    @gaplotdistance, ...  % Plots average distance (diversity) - Critical to see convergence
                    @gaplotmaxconstr, ... % Plots max constraint violation - Critical for feasibility
                    @gaplotrange}, ...    % Plots range of fitness values
        'ConstraintTolerance', 1e-4);
    
    % Conditionally use initial population
    % Use initial guess if it is non-empty AND useInitialGuess == true:
    if opt.use_x_norm_init_GA && ~isempty(opt.x_norm_init_GA)
        fprintf('\n  > Using Initial Guess for GA ...\n');
        pause(1);
        options_GA = optimoptions(options_GA, 'InitialPopulationMatrix', opt.x_norm_init_GA);
    else
        fprintf('\n  > NOT using Initial Guess for GA ...\n');
        pause(1);
    end

    rng default % For reproducibility
    tic;
    [x_norm_opt_GA, fopt_GA, exitflag, output] = ga(FitFcn_GA, 12, [],[],[],[], ...
                                                    opt.lb_norm, opt.ub_norm, Con_GA, options_GA);
    tGA = toc;

    fprintf('\n  GA Finished. Runtime: %02dh %02dm %05.2fs\n', floor(tGA/3600), floor(mod(tGA,3600)/60), mod(tGA,60));
    save('results/GA_opt_result.mat', 'x_norm_opt_GA', 'fopt_GA', 'signature');
    fprintf('  GA result saved. GA objective result is: %.6f', fopt_GA);
    pause(3);
end


%% ========================================================================
%  2.2. STEP 2: SQP REFINEMENT (Local Optimization)
%  ========================================================================
%  Uses 'SQP' mode: Strict Equality constraints (Thrust == Target)
%  Starts from the best GA point to refine the solution.

fprintf('\n--------- Starting SQP Refinement (Local Search) ----------\n');
pause(1);

% Wrappers for SQP:
FitFcn_SQP = @(x_norm) evaluateObjectivesAndConstraints(x_norm, mission, prop, motor, env, opt, 'SQP');
Con_SQP = @(x_norm) constraintsInterface(x_norm, mission, prop, motor, env, opt, 'SQP');

options_sqp = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 5000, ...
    'MaxFunctionEvaluations', 20000, ...
    'StepTolerance', 1e-9, ...
    'ConstraintTolerance', 1e-2, ... % 1% tolerance
    'UseParallel', true);

tic;
[x_norm_final_SQP, fopt_final_SQP, exitflag_sqp, output_sqp] = fmincon(FitFcn_SQP, x_norm_opt_GA, ...
                                                  [],[],[],[], opt.lb_norm, opt.ub_norm, Con_SQP, options_sqp);
tSQP = toc;

fprintf('\n  SQP Finished. Runtime: %02dh %02dm %05.2fs\n', floor(tSQP/3600), floor(mod(tSQP,3600)/60), mod(tSQP,60))
fprintf('\n  SQP improvement over GA: %.6f -> %.6f\n', fopt_GA, fopt_final_SQP);
save('results/SQP_final_result.mat', 'x_norm_final_SQP', 'fopt_final_SQP', 'signature');
fprintf('\n  Final Results saved to results/SQP_final_result.mat\n');
pause(3);


%% HELPER FUNCTION FOR SOLVERS
function [c, ceq] = constraintsInterface(x_norm, mission, prop, motor, env, opt, mode)
    % Memoization is used inside ObjectivesAndConstraints to NOT execute 
    % ObjectivesAndConstraints twice for each x (One for the
    % FitnessFunction and other for the Constraints)
    [~, c, ceq] = evaluateObjectivesAndConstraints(x_norm, mission, prop, motor, env, opt, mode);
end