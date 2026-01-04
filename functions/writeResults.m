%% ----------- Write numerical results in the Command Window --------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: writeResults
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Outputs a textual summary of the optimization results to the 
%       Command Window.
%
%   Displays:
%       - Final Geometric parameters (Bezier coefficients).
%       - Detailed performance breakdown for each Mission Segment:
%         (RPM, Thrust, Torque, Electrical Power, Efficiency).
%
%   Inputs:
%       x_final, mission, prop, motor, env, opt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeResults(x_norm_final, mission, prop, motor, env, opt)
    
    x_final = x_norm_final .* (opt.ub - opt.lb) + opt.lb;
    
    num_segments = numel(mission.v_ref);
    geometry_indices = numel(x_final)-num_segments;

    fprintf('\n==================== NUMERICAL RESULTS ====================\n');
    
    % 1. Reconstruct Geometry
    [R_opt, c_dist_opt, theta_deg_opt] = getPropGeometry(x_final, prop);
    fprintf('\n--- Geometric characteristics (numbers) ---\n');
    fprintf('Optimized Radius: %.4f cm\n', R_opt*100);
    fprintf('Bezier cubic coefficients for the pitch distribution are:')
    fprintf('a0 = %g, a1 = %g, a2 = %g, a3 = %g\n', x_final(1:4));
    fprintf('Bezier cubic coefficients for the chord distribution (NOT NORMALIZED) are:')
    fprintf('c0 = %g, c1 = %g, c2 = %g, c3 = %g\n', R_opt*prop.c0_norm, R_opt*x_final(5:7));
    
    % 2. Re-Evaluate all segments to get full details
    out_prop  = cell(1,num_segments);
    out_motor = cell(1,num_segments);
    res = cell(1,num_segments);
    segments = {'Slow Cruise', 'Standard Cruise', 'Climb', 'Fast Cruise'};
    for i = 1:num_segments
        v_curr = mission.v_ref(i);

        % rpm_indices = [9, 10, 11, 12];
        rpm_curr = x_final(geometry_indices + i);
        
        % Use the user-defined EvaluateModels function
        [out_prop{i}, out_motor{i}, ~, res{i}, ~] = evaluateModels(v_curr, rpm_curr, R_opt, ...
                                                    c_dist_opt, theta_deg_opt, prop, motor, env, opt);
        
        fprintf('\n--- %s ---\n', segments{i});
        fprintf('  Velocity:   %.1f m/s\n', v_curr);
        fprintf('  RPM:        %.0f\n', rpm_curr);
        fprintf('  Thrust:     %.3f N (Target: %.3f)\n', out_prop{i}.T, mission.T_ref(i));
        fprintf('  Torque:     %.3f Nm\n', out_prop{i}.Q);
        fprintf('  eta_mp:  %.2f%% (eta_p: %.1f%%, eta_m: %.1f%%)\n', ...
            res{i}.eta_mp*100, out_prop{i}.eta_p*100, out_motor{i}.eta_m*100);
        fprintf('  Elec Power: %.2f W (%.2f V * %.2f A)\n', ...
            out_motor{i}.V_req * out_motor{i}.I_mot, out_motor{i}.V_req, out_motor{i}.I_mot);
    end

end