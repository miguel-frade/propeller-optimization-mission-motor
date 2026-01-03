%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIMAL PROPELLER DESIGN FOR MISSION PROFILE WITH MOTOR RESTRICTIONS  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script: main_runner.m
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       This is the main script for the optimization framework.
%       It defines the mission profile, propeller constants, material 
%       properties, and brushless motor specifications. 
%       
%       The script orchestrates the following workflow:
%       1. parameter setup (Physics, Mission, Motor, Optimization bounds).
%       2. Execution of the optimization logic (GA + SQP).
%       3. Audit of final constraints.
%       4. Post-processing and visualization of results.
%
%   Design Variables (12 Total):
%       x(1:4): Twist distribution [rad] control parameters (Bezier cubic)
%       x(5:7): Normalized chord distribution control params (Bezier cubic)
%               c_norm is normalized with the actual radius: c_norm = c/R
%               Root chord normalized (c0_norm) is a constraint.
%       x(8): Radius normalized with R_ref ( R_norm = x(8) = R/R_ref )
%       x(9:12): RPM Scaling for each flight segment ( x(i) = rpm/rpm_ref )
%       
%       x = [a0, a1, a2, a3, c1n, c2n, c3n, Rn, rpm1n,rpm2n,rpm3n,rpm4n]
%
%       Design variables should all have the same order of magnitude (exact
%       normalization is not needed)
%
%   Dependencies & Data Flow Map:
%   
%   main_runner.m (User inputs)
%    │
%    ├── [1. Initialization]
%    │    └── getAirfoilDataE214_MultiRe_cleaned.m (Loads Airfoil Database)
%    │
%    ├── [2. Optimization Process]
%    │    └── run_optimization.m                  (GA + SQP Solvers)
%    │         └── evaluateObjectivesAndConstraints.m
%    │
%    ├── [3. Evaluation]
%    │    └── evaluateObjectivesAndConstraints.m  (Cost & Constraints Calc)
%    │         ├── getPropGeometry.m              (Builds Geometry from x)
%    │         └── evaluateModels.m               (Couples Prop + Motor)
%    │              ├── runSABEMMT.m              (BEMT Aerodynamics)
%    │              └── runMotorModel.m           (Electric Motor)
%    │
%    └── [4. Post-Process Results]
%         ├── writeResults.m                      (Text Output)
%         │    ├── getPropGeometry.m
%         │    └── evaluateModels.m
%         │
%         └── plotGeometryAndPerformances.m       (Visualization)
%              ├── getPropGeometry.m
%              └── runSABEMMT.m
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clear evaluateObjectivesAndConstraints; clc; close all;
% Clear the persistent variables in evaluateObjectivesAndConstraints.
addpath('functions');


%% ========================================================================
%  SECTION 0: LOGGING SETUP (SAVE EVERYTHING THAT APPEARS IN THE COMMAND
%             WINDOW AS .TXT)
%  ========================================================================
% 1. Create directory if needed
if ~exist('results', 'dir')
    mkdir('results');
end
% 2. Define log file path
logFile = fullfile('results', 'optimization_results.txt');
% 3. Delete existing log so we start fresh (otherwise diary appends)
if exist(logFile, 'file')
    delete(logFile);
end
% 4. Start recording to file
diary(logFile);


%% ========================================================================
%  SECTION 1: MISSION PROFILE, PHYSICS CONSTANTS & OPT SETTINGS (USER INPUTS)
%  ========================================================================
mission = struct();
prop = struct();
motor = struct();
env = struct();
opt = struct();

% ----------------------- 1.1 Mission Profile -----------------------------
%  Define the 4 flight segments
%  For our specific example (hand thrown UAV), we don't care about takeoff

% Define the names for the mission segments (that will appear above the plots)
mission.names = {'Slow cruise', 'Cruise', 'Climb', 'Fast cruise'};

% Segments:       [   Slow,   Cruise,  Climb,     Fast]
mission.v_ref   = [ 7.8401,     10.0,   10.0,    16.67]; % Velocity
mission.T_ref   = [0.39578,  0.34343,   2.00,  0.51298]; % Thrust
mission.weights = [   0.10,     0.70,   0.15,     0.05]; % Cost weights
% For v_ref(3)=20: T_ref(3)=0.6895

fprintf('\n============= MISSION PROFILE SPECIFICATIONS ==============\n');
fprintf('  1. Slow Cruise (%.0f%%): %.1f m/s, %.2f N\n', mission.weights(1)*100, mission.v_ref(1), mission.T_ref(1));
fprintf('  2. Standard Cruise (%.0f%%): %.1f m/s, %.2f N\n', mission.weights(2)*100, mission.v_ref(2), mission.T_ref(2));
fprintf('  3. Exigent Climb  (%.0f%%): %.1f m/s, %.2f N\n', mission.weights(3)*100, mission.v_ref(3), mission.T_ref(3));
fprintf('  4. Fast Cruise  (%.0f%%): %.1f m/s, %.2f N\n', mission.weights(4)*100, mission.v_ref(4), mission.T_ref(4));
pause(2); % To let the user read the message


% --------------- 1.2 Propeller Constants & Normalization -----------------
% The same airfoil is used all along the radius (only one airfoil is used)
prop.airfoil_data = getAirfoilDataMultiReE214();  % Airfoil selection
prop.Nb = 2;     % Number of blades
prop.x0 = 0.05;  % Section at which the aeroydnamic blade starts (normalized radius, x0=r0/R)
                
% prop.x0 = 0.18 is a realistic value based on real propeller data.
% However, since we are using Bezier cubic shape functions, the
% optimization will always try to make a0 very big in order for the pitch
% distribution to be DECREASING MONOTONIC (the optimization will always
% want negative slope of theta(r) at r0).
% If we set prop.x0 = 0.18, the optimization will asign almost the max
% theta it can to that x0, and it is irrealistic to have very big theta at
% 20% of the radius.
% This problem can be solved by just setting a very small x0 (for example,
% prop.x0 = 0.05). It is obviously irreal that the aerodynamics start at 5%
% of the radius, but what is more important is that that allows for a
% realistic decreasing monotonic theta(r). In the CAD design of the
% propeller after this optimization process, the region of the propeller
% from x0=0.05 to x=0.18 can just be removed and designed focused on
% structural reinforcement.
% IMPORTANT: In order for the structural constraints to be realistic even
% when x0 is very small, we need to only consider the stresses starting at:

prop.x_finish_reinforcement = 0.25; % This is a realistic value (do not
                                    % make it smaller).
                                    % For smaller x, we can thicken the
                                    % airfoils while keeping c constant.

% The sections closer to the root will have structural reinforcement,
% so we won't worry about them.

% Root chord normalized: c0_norm (root chord for R=1 propeller)
% Default geom formula: 0.7 * prop.x0 * 2 * sin((2*pi/prop.Nb)/2)
% You can overwrite this value manually here if needed.
%
% For prop.x0=0.18, prop.c0_norm = 0.145 is a very realistic value.
% However, as explained above, we are going to use prop.x0=0.05, and for
% that, a good value is:

prop.c0_norm = sqrt(2)*0.08;

% Mesh Setup
prop.n = 80;
prop.mesh = linspace(prop.x0, 1, prop.n);
prop.dx = gradient(prop.mesh);

% Mach tip limit:
prop.M_tip_max = 0.65;

% Propeller material characteristics:
prop.rho_mat = 1240;        % Material density
prop.sigma_p_MPa = 40;      % Yield stress [MPa]
prop.sigma_p = prop.sigma_p_MPa * 1e6;


% --------- 1.3 Brushless motor Constants (T-Motor AS2304 1500kV) ---------
motor.KV = 1500;
motor.V_batt = 7.4;    % Battery Voltage
motor.Rm = 0.251;      % Resistance [Ohms]
motor.I0 = 0.3;        % No-load current [A]
motor.I_max = 9.2;     % Max Current [A]
motor.Kt = 30 / (pi * motor.KV);


% -------------------------- 1.4 Environment ------------------------------
env.rho = 1.15;            % Air density [kg/m^3]
env.sound_speed = 340;     % Speed of sound [m/s]
env.nu = 1.5e-5;           % Kinematic viscosity


% ------------------ 1.5 Optimization Settings & Bounds -------------------
% Normalization of the objectives
% This is important if the magnitudes we use in the weighted sum are
% different. In our case, it is all efficiencies, so we don't really need
% any normalization:
opt.eta_mp_min = 0;    % Min efficiency (for objective normalization)
opt.eta_mp_max = 1;    % Max efficiency (for objective normalization)

% Optimization Bounds (Lower & Upper)
% Variable Definition (12 Variables):
% x(1:4): Twist distribution [rad] control parameters (Bezier cubic)
% x(5:7): Normalized chord distribution control points (Bezier cubic)
%         c_norm is normalized with the actual radius: c_norm = c / R
%         We use normalized chord in x because we know c0_norm, but not c0.
% x(8):    Dimensional Radius ( x(8)=R )
% x(9:12): RPM for each flight segment

% x    = [a0,  a1,  a2,   a3,  c1n, c2n,          c3n,   R,   rpm1,  rpm2,  rpm3,  rpm4]
opt.lb = [0.3,   0,   0, 0.05,   0,   0,            0, 0.04,  1000,  1000,  1000,  1000]; % Lower bound
opt.ub = [1.5, 1.5, 1.5,  1.0, 0.5, 0.5, prop.c0_norm,  0.2, 12000, 12000, 12000, 12000];  % Upper bound
% Note: Mach_tip limit is active in constraints, but we guide the search here.

% Variables normalization will be done as x_norm = (x-lb) ./ (ub-lb)
opt.lb_norm = zeros(1,12);
opt.ub_norm = ones(1,12);

% Define Initial Guess for Warm Start (can be omitted)
% Decide whether to use Initial Guess or not FOR GA:
% Initial guess will be used if it is non-empty AND use_x_init = true:
opt.use_x_norm_init_GA = false;
% All the decision logic is in the code 'runOptimization'

% This initial guess is the result of optimization for 2s lipo, Tref(3)=1.8
 opt.x_norm_init_GA = [ ...
          1, 0.214304, 0.372971, 0.221167, ... % Twist parameters norm
          0.6372084, 0.284983, 0, ...          % Chord parameters norm-norm
          0.376749, ...                        % Radius norm
          0.2497851, 0.2645628, 0.513772, 0.41865 ... % RPM norm
          ];

% opt.x_norm_init_GA = (opt.x_init_GA - opt.lb) ./ (opt.ub - opt.lb);


% ------------- 1.6 Optimization Constraints settings & switches ----------
% Stress constraints safety Factor: 
% sf = 2 (Load) * 1.5 (Material) = 3.75
opt.constr.stress_sf = 3;

% Thrust constraint tolerance for GA:
% The function evaluateObjectivesAndConstraints for GA rewards satifying
% the constraints exactly, but a constraint tolerance is needed for the GA
% in order to find initial feasible solutions, and in order to be able to
% perform more mutations.
opt.constr.GA_thrust_tol = 0.3; % 30% tolerance bands around T_ref for GA

% Constraint Control (1 = On, 0 = Off)

% Vector Order: [Slow, Cruise, Climb, Fast]
opt.constr_flags.thrust  = [1, 1, 1, 1]; % Check Thrust targets
opt.constr_flags.M_tip   = [0, 0, 1, 1]; % Check Mach limit (Usually only needed for Climb and Fast)
opt.constr_flags.stress  = [0, 0, 1, 0]; % Check Stress (Usually max load is Climb)
opt.constr_flags.current = [1, 1, 1, 1]; % Check Motor Current
opt.constr_flags.voltage = [1, 1, 1, 1]; % Check Battery Voltage



%% ========================================================================
%  SECTION 2: RUN OPTIMIZATION
%  ========================================================================
%  Runs GA (Global) and SQP (Local) using the parameters defined above.
%  Saves results to 'results/SQP_final_result.mat'.

run_optimization; 


%% ========================================================================
%  SECTION 3: SQP CONSTRAINTS AUDIT
%  ========================================================================
fprintf('\n================== SQP CONSTRAINT AUDIT ===================\n');
% If variables don't exist in workspace, load them from the file
if ~exist('x_norm_final_SQP', 'var') || ~exist('mission', 'var')
    resultsFile = 'results/SQP_final_result.mat';
    if exist(resultsFile, 'file')
        load(resultsFile, 'x_norm_final_SQP', 'fopt_final_SQP', 'signature');
    else
        error('  > Results file not found (%s). Please run main_runner.m previous sections first.', resultsFile);
    end
end

% Re-evaluate final design with 'SQP' mode
[~, c, ceq] = evaluateObjectivesAndConstraints(x_norm_final_SQP, mission, prop, motor, env, opt, 'SQP');
tol = 1e-2; % 1% tolerance

% --- 1. Thrust Targets (Equality) ---
fprintf('THRUST TARGETS (Error):\n');
segs = {'Slow', 'Cruise', 'Climb', 'Fast'};
k_eq = 1; % Counter for the ceq vector index

for i = 1:4
    if opt.constr_flags.thrust(i)
        % Constraint is active, so we read the next value from ceq
        err = ceq(k_eq) * 100;
        if abs(ceq(k_eq)) <= tol
            stat = 'SATISFIED';
        else
            stat = '**VIOLATED**';
        end
        fprintf('  %-15s: %7.3f%%  [%s]\n', segs{i}, err, stat);
        k_eq = k_eq + 1; % Increment ceq index only when a constraint was present
    else
        % Constraint was disabled in settings
        fprintf('  %-15s:       -    [SKIPPED/OFF]\n', segs{i});
    end
end

% --- 2. Physical Limits (Inequality) ---
fprintf('PHYSICAL LIMITS (Usage):\n');

% Dynamically build the names list to match how 'c' is built in evaluateObjectivesAndConstraints
names = {};
segment_names = {'Slow', 'Cruise', 'Climb', 'Fast'};

% Note: The order MUST match exactly the order in evaluateObjectivesAndConstraints.m
% 1. Loop through segments for Current, Voltage, M_tip, Stress
for i = 1:4
    % A. Motor Current
    if opt.constr_flags.current(i)
        names{end+1} = sprintf('%s I_mot', segment_names{i});
    end
    % B. Battery Voltage
    if opt.constr_flags.voltage(i)
        names{end+1} = sprintf('%s V_batt', segment_names{i});
    end
    % C. Mach Number
    if opt.constr_flags.M_tip(i)
        names{end+1} = sprintf('%s M_tip', segment_names{i});
    end
    % D. Structural Stress
    if opt.constr_flags.stress(i)
        names{end+1} = sprintf('%s Stress', segment_names{i});
    end
end

% Check for consistency
if numel(c) ~= numel(names)
    warning('Mismatch between constraint vector length (%d) and generated names (%d). Check logic consistency.', numel(c), numel(names));
end

% Print Results
for i = 1:numel(c)
    % Usage %: 100% means exactly on the limit (0). <100% means safe.
    usage = (c(i) + 1) * 100; 
    if c(i) > tol
        stat = '**VIOLATED**';
    elseif c(i) > -1e-3
        stat = 'ACTIVE (ON LIMIT)';
    else
        stat = 'SATISFIED';
    end
    % Use fixed width for cleaner alignment
    fprintf('  %-15s: %6.1f%%   [%s]\n', names{i}, usage, stat);
end



%% ========================================================================
%  SECTION 4: ANALYSIS & PLOTTING
%  ========================================================================
fprintf('\n==================== NUMERICAL RESULTS ====================\n');
if ~exist('x_norm_final_SQP', 'var') || ~exist('mission', 'var')
    resultsFile = 'results/SQP_final_result.mat';
    if exist(resultsFile, 'file')
        load(resultsFile, 'x_norm_final_SQP', 'fopt_final_SQP', 'signature');
    else
        error('  > Results file not found (%s). Please run main_runner.m previous sections first.', resultsFile);
    end
end

% --- Write the numerical results in the command window ---
writeResults(x_norm_final_SQP, mission, prop, motor, env, opt);

% STOP LOGGING
diary off;
fprintf('Full command window output saved to %s\n', logFile);


% --- Visualization ---
% The function plotGeometryAndPerformances will save the plots as .jpg
fprintf('\n===================== PLOT GENERATION =====================\n');
reply = input('\n  > Generate plots for this solution? (y/n): ','s');
if strcmpi(reply,'y')
    fprintf('  > Generating plots...\n');
    plotGeometryAndPerformances(x_norm_final_SQP, mission, prop, motor, env, opt);
else
    fprintf(' > Plot generation skipped\n');
end






