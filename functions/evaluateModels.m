%% ---------------- Evaluate Propeller and Motor models -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: evaluateModels
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Couples the Aerodynamic Propeller Model (SABEMMT) with the Electro-
%       Mechanical Motor Model. It evaluates the performance of a specific
%       propeller geometry at a specific flight condition.
%
%   Inputs:
%       v           - Flight velocity [m/s]
%       rpm         - RPM
%       R           - Propeller Radius [m]
%       c_dist      - Chord distribution vector [m]
%       theta_deg   - Twist distribution vector [deg]
%       prop, motor, env, opt - Configuration structs
%
%   Outputs:
%       out_prop    - Propeller physics output (Thrust, Torque, Efficiency)
%       out_motor   - Motor electric output (Current, Voltage, Efficiency)
%       obj         - Normalized single-point objective
%       res         - Summary struct for constraints (Sigma, M_tip, etc.)
%       success     - Boolean flag for valid BEMT convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_prop, out_motor, obj, res, success] = evaluateModels(v, rpm, R, ...
                                                    c_dist, theta_deg, prop, motor, env, opt)
    
    % 1. Run BEMT to get Propeller results
    out_prop = runSABEMMT(v, rpm, prop.Nb, R, c_dist, theta_deg, prop, env);

    % For the stress constraint, we only want to consider now the part of
    % the blade that is pure aerodynamics (without structural
    % reinforcement):
    % The sections closer to the root will have structural reinforcement,
    % so we won't worry about them breaking.
    blade_stresses = out_prop.sigma_total_max(prop.mesh > prop.x_finish_reinforcement);
    max_sigma = max(blade_stresses);

    
    % Validation of BEMT
    check = [out_prop.eta_p, out_prop.T]; % Do NOT include out_prop.sigma_total_max
    if any(isnan(check) | isinf(check) | ~isreal(check))
        success = false; 
        obj = 10; 
        res = []; 
        out_motor = [];
        return;
    end
    
    success = true;

    % 2. Calculate Motor Performance
    out_motor = runMotorModel(out_prop.Q, rpm, motor);
    
    % 3. Calculate Combined Efficiency
    eta_mp = out_motor.eta_m * out_prop.eta_p;

    % Normalize Objective and make it negative: We want to Maximize Eta_mp
    obj = -(eta_mp - opt.eta_mp_min) / (opt.eta_mp_max - opt.eta_mp_min);
    
    % 4. Pack Results
    res.T = out_prop.T;
    res.max_sigma = max_sigma;
    res.M_tip = out_prop.M_tip;
    res.I_mot = out_motor.I_mot;
    res.V_req = out_motor.V_req;
    % We also output eta_mp as a result in case we want it (the obj output
    % was eta_mp normalized and with changed sign):
    res.eta_mp = eta_mp;

end

