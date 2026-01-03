%% ------------------- Plot Geometry and Performances ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: plotGeometryAndPerformances
%   Project: propeller-optimization-mission-motor
%   Date: December 2025
%   Author: Miguel Frade
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Comprehensive post-optimization visualization suite. This function
%       reconstructs the optimal propeller geometry from the design vector,
%       re-evaluates it against the mission segments using BEMT, and
%       generates detailed structural and aerodynamic plots.
%
%       It automatically creates a 'Results/' directory and saves the
%       following high-resolution figures as JPGs:
%
%       1.  3D Propeller Rendering (Patch-based visualization with Hub)
%       2.  Geometry Dashboard (Chord/Twist dist., 2D Planform, 3D Mesh)
%       3.  Airfoil Aerodynamic Polars (Cl/Cd vs Alpha across Re range)
%       4.  Blade Element Velocity Triangles (Inflow geometry at 0.75R)
%       5.  Aerodynamic Distributions (Lambda, Alpha, and Stall Margins)
%       6.  Thrust Loading Distribution (dT/dx)
%       7.  Tangential Force Distribution (dFt/dx)
%       8.  Centrifugal Stress Distribution
%       9.  Bending Stress Distribution (Principal axes)
%       10. Total Stress Distribution (Von Mises/Max Principal)
%       11. Critical Airfoil Section Stress Analysis (Detailed structural view)
%       12. Performance Clouds (Eta, Ct, Cp vs Advance Ratio J)
%
%   Inputs:
%       x_norm_final : [1xN double] Normalized optimized design vector.
%       mission      : [struct]     Flight mission parameters (velocities, names).
%       prop         : [struct]     Propeller fixed parameters (airfoil, materials).
%       motor        : [struct]     Motor characteristics.
%       env          : [struct]     Environmental constants (density, viscosity).
%       opt          : [struct]     Optimization bounds (lb, ub) for denormalization.
%
%   Dependencies:
%       - getPropGeometry.m
%       - runSABEMMT.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plotGeometryAndPerformances(x_norm_final, mission, prop, motor, env, opt)

    x_final = x_norm_final .* (opt.ub - opt.lb) + opt.lb;
    
    num_segments = numel(mission.v_ref);
    geometry_indices = numel(x_final)-num_segments;

    % Reconstruct Geometry
    [R_opt, c_dist_opt, theta_deg_opt] = getPropGeometry(x_final, prop);

    % Re-Evaluate all segments to get full details
    outputs = cell(1,num_segments);
    for i = 1:num_segments
        v_curr = mission.v_ref(i);

        % rpm_indices = [9, 10, 11, 12];
        rpm_curr = x_final(geometry_indices + i);
        
        % Use the user-defined EvaluateModels function
        outputs{i} = runSABEMMT(v_curr, rpm_curr, prop.Nb, R_opt, c_dist_opt, theta_deg_opt, prop, env);
    end

    % For extracting the geometry, we want to use a non-empty output:
    idx = find(~cellfun(@isempty, outputs), 1, 'first');
    if ~isempty(idx)
        nonEmptyOutput = outputs{idx};
    else
        nonEmptyOutput = [];  % or handle the "all empty" case explicitly
    end
    
    % Define titles for each column
    titles = cell(1, 4);
    
    for i = 1:4
        % Construct the 3 lines (3 lines for each flight segment title)
        mission_name = mission.names{i};
        line1 = sprintf('\\underline{\\textbf{%s}}', mission_name);
        line2 = sprintf('$v=%.0f$, $T=%.2f$', mission.v_ref(i), outputs{i}.T);
        line3 = sprintf('$RPM=%.0f$, $\\eta_p=%.2f$', outputs{i}.rpm, outputs{i}.eta_p);
        
        % Store combined lines
        titles{i} = {line1; line2; line3};
    end
 
    
    % ================= [TOP BLOCK: START] =================
    % Don't show the figures until they are all saved.
    % Tell MATLAB to always restore visibility, even if an error occurs.

    % Store current default settings to restore them later
    defaultVis = get(0, 'DefaultFigureVisible');
    
    % Force all subsequent figures to be invisible initially
    set(0, 'DefaultFigureVisible', 'off');
    
    % Safety mechanism: If code crashes, restore visibility automatically
    cleanupVis = onCleanup(@() set(0, 'DefaultFigureVisible', defaultVis));
    % ================== [TOP BLOCK: END] ==================




    % =======================================================================
    % 1. Entire propeller Geometry in 3D (Hub matches Blade Color)
    % =======================================================================
    figure('Name', '1. Propeller geometry in 3D', 'Color', 'w', ...
           'Renderer', 'opengl'); 
    hold on; grid on; axis equal;
    
    % 1. Setup Geometry Inputs
    mesh = nonEmptyOutput.mesh;
    x0 = nonEmptyOutput.mesh(1);
    R = nonEmptyOutput.R;
    r_phys_blade = R .* mesh; 
    c_blade = nonEmptyOutput.c;
    theta_rad_blade = nonEmptyOutput.theta;

    x_hub = 0.07; % I want the hub to be always drawed up to 7% of the radius
    % The aerodynamic sections can start at greater x (x0 approx. 0.2), but
    % a smooth transition will be drawn from the hub to the aerodynamic blade.

    x_u = prop.airfoil_data.x_u;
    x_l = prop.airfoil_data.x_l;
    y_u = prop.airfoil_data.y_u;
    y_l = prop.airfoil_data.y_l;
    x_poly = [x_u, fliplr(x_l)]; 
    y_poly = [y_u, fliplr(y_l)];
    x_axis_norm = 0.25; % Stacking Axis
    
    % Pre-allocate Grids
    num_points_foil = length(x_poly);
    
    % 2. Generate Aerodynamic Blade Mesh (The original blade)  
    num_blade_sections = length(r_phys_blade);
    X_blade = zeros(num_points_foil, num_blade_sections);
    Y_blade = zeros(num_points_foil, num_blade_sections);
    Z_blade = zeros(num_points_foil, num_blade_sections);
    for i = 1:num_blade_sections
        radius = r_phys_blade(i);
        chord  = c_blade(i);
        twist  = -theta_rad_blade(i); 
        
        % Scale and Center Airfoil
        xx_local = (x_poly - x_axis_norm) * chord;
        yy_local = y_poly * chord;
        
        % Apply Twist Rotation
        x_rot = xx_local * cos(twist) - yy_local * sin(twist);
        z_rot = xx_local * sin(twist) + yy_local * cos(twist);
        
        % Map to Global Coordinates
        X_blade(:, i) = x_rot;
        Y_blade(:, i) = radius;
        Z_blade(:, i) = z_rot;
    end

    % ---------------------------------------------------------
    % 2B. Generate Transition Shank (Hub to Blade Root)
    % ---------------------------------------------------------
    X_trans = []; Y_trans = []; Z_trans = [];
    
    if x0 > x_hub
        % Settings for the transition
        num_trans_sections = 15; 
        r_trans = linspace(0.91*x_hub*R, x0*R, num_trans_sections);
        % I use 0.91*x_hub*R to make sure there is no gap between hub and
        % transition shank. Don't use less than 0.91, because the blade
        % "enters" the hub and it is not aesthetic.
        
        % Get geometry of the blade root (Target Shape)
        c_root = c_blade(1);
        twist_root = -theta_rad_blade(1);
        
        % 1. Define Root Airfoil Local Coordinates (Target)
        xx_root_local = (x_poly - x_axis_norm) * c_root;
        yy_root_local = y_poly * c_root;
        
        % 2. Define Hub Circle Local Coordinates (Source)
        % We project the airfoil points onto a circle to ensure 
        % point-to-point correspondence (topology matching).
        
        % Calculate thickness of root to sizing the shank cylinder
        root_thickness = max(y_poly) - min(y_poly);
        shank_radius = (root_thickness * c_root) * 1.8; % <--- Radius of intersection (can be modified)
        
        % Calculate angles of the airfoil points relative to centroid
        angles = atan2(yy_root_local, xx_root_local);
        
        % Create circle points based on those angles
        xx_circ_local = shank_radius .* cos(angles);
        yy_circ_local = shank_radius .* sin(angles);
        
        % Initialize grids
        X_trans = zeros(num_points_foil, num_trans_sections);
        Y_trans = zeros(num_points_foil, num_trans_sections);
        Z_trans = zeros(num_points_foil, num_trans_sections);
        
        for i = 1:num_trans_sections
            % Interpolation factor (0 = Hub, 1 = Blade Root)
            % Using a sinusoidal easing for visual smoothness
            s_lin = (i-1) / (num_trans_sections-1);
            s = sin(s_lin * pi/2); % Erase-out
            
            radius = r_trans(i);
            
            % Interpolate Shapes (Morph Circle -> Airfoil)
            xx_current = xx_circ_local * (1-s) + xx_root_local * s;
            yy_current = yy_circ_local * (1-s) + yy_root_local * s;
            
            % Apply Twist
            % We keep the twist constant (equal to root twist) for the 
            % shank so it enters the hub cleanly.
            twist = twist_root; 
            
            x_rot = xx_current * cos(twist) - yy_current * sin(twist);
            z_rot = xx_current * sin(twist) + yy_current * cos(twist);
            
            X_trans(:, i) = x_rot;
            Y_trans(:, i) = radius;
            Z_trans(:, i) = z_rot;
        end
        
        % Remove the last section of transition to avoid duplicate points 
        % at the join with the main blade
        X_trans = X_trans(:, 1:end-1);
        Y_trans = Y_trans(:, 1:end-1);
        Z_trans = Z_trans(:, 1:end-1);
    end
    
    % ---------------------------------------------------------
    % 2C. Combine Grids
    % ---------------------------------------------------------
    X_grid = [X_trans, X_blade];
    Y_grid = [Y_trans, Y_blade];
    Z_grid = [Z_trans, Z_blade];
    
    num_sections = size(X_grid, 2);

    % --- Aesthetic Settings (UPDATED) ---
    upper_color = [0.15 0.35 0.75];   % Blue (upper surface)
    lower_color = [0.80 0.55 0.20];   % Bronze (lower surface)
    
    % Set Hub Color to match Upper Blade Color
    prop_color  = upper_color;         
    
    spec_strength = 0.6;            
    diff_strength = 0.8;            
    
    % --- 3. Draw Blades using PATCH ---
    faces = [];
    vertices = [];
    colors = [];
    
    nv = 0;
    
    for i = 1:num_sections-1
        for j = 1:num_points_foil-1
            % Quad vertices (current blade)
            v1 = [X_grid(j,i),   Y_grid(j,i),   Z_grid(j,i)];
            v2 = [X_grid(j+1,i), Y_grid(j+1,i), Z_grid(j+1,i)];
            v3 = [X_grid(j+1,i+1), Y_grid(j+1,i+1), Z_grid(j+1,i+1)];
            v4 = [X_grid(j,i+1), Y_grid(j,i+1), Z_grid(j,i+1)];
    
            vertices = [vertices; v1; v2; v3; v4];
            faces = [faces; nv+1 nv+2 nv+3 nv+4];
    
            % Upper vs lower by airfoil index
            if j <= length(x_u)
                colors = [colors; upper_color];
            else
                colors = [colors; lower_color];
            end
    
            nv = nv + 4;
        end
    end
    
    patch('Faces',faces,'Vertices',vertices, ...
          'FaceVertexCData',colors, ...
          'FaceColor','flat', ...
          'CDataMapping','direct', ...
          'EdgeColor','none', ...
          'FaceLighting','gouraud', ...
          'BackFaceLighting','reverselit');
    
    % --- Second blade (mirrored) ---
    vertices2 = vertices;
    vertices2(:,1:2) = -vertices2(:,1:2);
    
    patch('Faces',faces,'Vertices',vertices2, ...
          'FaceVertexCData',colors, ...
          'FaceColor','flat', ...
          'CDataMapping','direct', ...
          'EdgeColor','none', ...
          'FaceLighting','gouraud', ...
          'BackFaceLighting','reverselit');
    
    
    % 4. Precision Hub Generation
    % Step A: Find the EXACT vertical extent of the blade root
    % NOTE: We use the transition root (first col of Grid), not the aero root
    root_z_points = Z_grid(:, 1); 
    z_max = max(root_z_points);
    z_min = min(root_z_points);
    
    % Step B: Define Hub Dimensions based on these limits
    hub_height = (z_max - z_min) * 1.05; 
    hub_center_z = (z_max + z_min) / 2;
    
    hub_radius_outer = (R * x_hub); 
    shaft_radius_inner = (R * 0.03); 
    
    % Step C: Draw Outer Hub (Uses prop_color, now Blue)
    [Hx, Hy, Hz] = cylinder(hub_radius_outer, 100); 
    Hz = (Hz - 0.5) * hub_height + hub_center_z; 
    surf(Hx, Hy, Hz, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
    
    % Step D: Draw Inner Hole
    [Sx, Sy, Sz] = cylinder(shaft_radius_inner, 50);
    Sz = (Sz - 0.5) * hub_height * 1.02 + hub_center_z; 
    surf(Sx, Sy, Sz, 'FaceColor', 'k', 'EdgeColor', 'none', ...
         'DiffuseStrength', 0, 'SpecularStrength', 0);

    % --- Step E: Draw Hub Covers (Annulus Caps) ---
    % Create a polar grid for the top/bottom caps
    r_cap = linspace(shaft_radius_inner, hub_radius_outer, 2); 
    th_cap = linspace(0, 2*pi, 80);
    [R_cap, TH_cap] = meshgrid(r_cap, th_cap);
    
    X_cap = R_cap .* cos(TH_cap);
    Y_cap = R_cap .* sin(TH_cap);
    
    % Top Cap (Blue)
    Z_cap_top = ones(size(X_cap)) * (hub_center_z + hub_height/2);
    surf(X_cap, Y_cap, Z_cap_top, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
         
    % Bottom Cap (Blue)
    Z_cap_bot = ones(size(X_cap)) * (hub_center_z - hub_height/2);
    surf(X_cap, Y_cap, Z_cap_bot, 'FaceColor', prop_color, 'EdgeColor', 'none', ...
         'SpecularStrength', spec_strength, 'DiffuseStrength', diff_strength, ...
         'FaceLighting', 'gouraud');
    
    % 5. Illumination
    % shading flat; 
    lighting gouraud;
    delete(findall(gcf,'Type','light'));
    
    light('Position', [R, R, 2*R], 'Style', 'local', 'Color', [1 1 1]);
    light('Position', [-R, -R, -R], 'Style', 'local', 'Color', [0.6 0.6 0.6]);
    light('Position', [R, -R, 0], 'Style', 'local', 'Color', [0.7 0.7 0.7]);
    
    % 6. Plot Formatting (for LaTeX)
    xlabel('Tangential coordinate, $x$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Spanwise coordinate, $y$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Thrust axis, $z$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    
    mainTitle = sprintf('\\textbf{1. Propeller geometry in 3D} ($R=%.2f$ cm)', R*100);
    t = title(mainTitle, 'Interpreter', 'latex', 'FontSize', 30);
    t.Units = 'normalized';
    t.Position(2) = t.Position(2) + 0.05;   % increase vertical offset

    
    box on; 
    grid on;
    axis on; 
    set(gca, 'GridAlpha', 0.3); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    view(-45, 20); 
    
    limit_span = R * 1.1;
    limit_axial = max(hub_height, max(abs(Z_grid(:)))) * 1.5;
    xlim([-limit_span/2 limit_span/2]); 
    ylim([-limit_span limit_span]);     
    zlim([-limit_axial limit_axial]); 
    
    camlight('headlight');






    % =======================================================================
    % 2. COMPLETE GEOMETRY ANALYSIS (3 Subplots in 1 Figure)
    % =======================================================================
    figure('Name', '2. Propeller Geometry Analysis', 'Position', [100 50 1200 900], 'Color', 'w');
    
    % --- Prepare Data Common to All Plots ---
    % 1. Extract data
    x_u = prop.airfoil_data.x_u;  % u --> upper (extrados)
    y_u = prop.airfoil_data.y_u;
    x_l = prop.airfoil_data.x_l;
    y_l = prop.airfoil_data.y_l;

    R = nonEmptyOutput.R;
    mesh_vec = nonEmptyOutput.mesh; % Non-dimensional (0 to 1)
    r_vec = mesh_vec * R;       % Dimensional radius [m]
    c = nonEmptyOutput.c;
    theta_deg = (180/pi) * nonEmptyOutput.theta;
    n = numel(mesh_vec);
    
    % ---------- SUBPLOT 1 (Top): Chord and Pitch Distribution ------------
    subplot(2, 2, [1, 2]); % Spans the entire top row
    
    % --- Limits Logic (Preserved) ---
    c_ylim = [0, max(c) * 1.1];
    theta_min = min(theta_deg);
    theta_max = max(theta_deg);
    buffer = 0.05 * (theta_max - theta_min);
    theta_ylim = [theta_min - buffer, theta_max + buffer];
    
    % --- Left Axis: Chord ---
    yyaxis left
    plot(mesh_vec, c, 'b-', 'LineWidth', 2, 'DisplayName', 'Chord $c$');
    ylabel('Chord, $c$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(c_ylim);
    set(gca, 'YColor', 'b'); % Ensure axis text is blue
    
    % --- Right Axis: Pitch ---
    yyaxis right
    plot(mesh_vec, theta_deg, 'r--', 'LineWidth', 2, 'DisplayName', 'Pitch $\theta$');
    ylabel('Pitch, $\theta$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
    ylim(theta_ylim);
    set(gca, 'YColor', 'r'); % Ensure axis text is red
    
    % --- Common Axis Settings ---
    xlim([mesh_vec(1), 1]);
    xlabel('Radial non-dimensional coordinate, $x = r/R$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. Radial Distribution of Chord and Pitch}', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('show', 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % --------- SUBPLOT 2 (Bottom Left): 2D Planform (Upper View) ---------
    subplot(2, 2, 3);
    hold on; grid on; axis equal; box on;
    
    % Calculate Leading/Trailing Edges based on 1/4 chord alignment
    X_LE = -0.25 * c;
    X_TE =  0.75 * c;
    
    % Plot Lines
    plot(r_vec, X_LE, 'b-', 'LineWidth', 2);
    plot(r_vec, X_TE, 'r-', 'LineWidth', 2);
    plot(r_vec, zeros(size(r_vec)), 'k--', 'LineWidth', 1); % Quarter chord line
    
    % Fill the blade shape
    fill([r_vec, fliplr(r_vec)], [X_LE, fliplr(X_TE)], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    
    % Formatting
    xlabel('Radius, $r$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Chordwise Position [m]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{B. 2D Planform View}', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Custom Legend (Manual construction for clarity)
    % We create dummy entries for the legend to avoid cluttering the plot
    h1 = plot(nan, nan, 'b-', 'LineWidth', 2);
    h2 = plot(nan, nan, 'r-', 'LineWidth', 2);
    h3 = plot(nan, nan, 'k--', 'LineWidth', 1);
    legend([h1 h2 h3], {'Leading Edge', 'Trailing Edge', '1/4 Chord'}, ...
        'Interpreter', 'latex', 'Location', 'SouthWest', 'FontSize', 10);
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % ----------- SUBPLOT 3 (Bottom Right): 3D Blade Geometry -------------
    subplot(2, 2, 4);
    
    % --- Data Preparation for 3D Surface ---
    % Make airfoil data local and normalized
    x_profile = [x_u, fliplr(x_l)];
    y_profile = [y_u, fliplr(y_l)];
    x_profile = x_profile - min(x_profile);
    x_profile = x_profile ./ max(x_profile);
    Np = length(x_profile);
    
    % Initialize Matrices
    X_surf = zeros(Np, n);
    Y_surf = zeros(Np, n);
    Z_surf = zeros(Np, n);
    
    % Build 3D coordinates
    for j = 1:n
        c_j = c(j);
        th = deg2rad(theta_deg(j));
        
        % Scale Airfoil
        xs = c_j * (x_profile - 0.25);
        ys = c_j * y_profile;
        
        % Rotate for Pitch (Twist)
        % Note: X is Radial, Y is Chordwise (Tangential), Z is Thickness/Vertical
        % Adjusting rotation to match standard "propeller lying flat" or "vertical" view
        
        % Standard rotation logic:
        xr = xs*cos(-th) - ys*sin(-th);
        yr = xs*sin(-th) + ys*cos(-th);
        
        X_surf(:, j) = r_vec(j); % Radial Axis
        Y_surf(:, j) = xr;       % Chord/Tangential Axis
        Z_surf(:, j) = yr;       % Thickness/Twist Axis
    end
    
    % Plot Surface
    surf(X_surf, Y_surf, Z_surf, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 1]);
    
    % Formatting
    axis equal; grid on;
    xlabel('Radius, $r$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Chordwise, $y$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Thickness, $z$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{C. 3D Blade Geometry}', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Lighting for 3D effect
    camlight; lighting gouraud;
    view(3); % Standard isometric view
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

    % Figure Super Title
    mainTitle = sprintf('\\textbf{2. Propeller Geometry Analysis} ($R=%.2f$ cm)', R*100);
    sgtitle(mainTitle, 'Interpreter', 'latex', 'FontSize', 15);


    




    % =======================================================================
    % 3. Plot the airfoil's aerodynamic curves
    % =======================================================================
    
    % ---- Plot cl, cd vs Angle of Attack (Multi-Re) ----
    figure('Name', '3. Aerodynamic Aerodynamic curves', 'Color', 'w'); clf;
    
    airfoil_data = prop.airfoil_data;

    colors = jet(length(airfoil_data.Re));
    legend_list = {};
    
    % --- Subplot 1: Lift Coefficient ---
    subplot(1,2,1); hold on; grid on;
    for i = 1:length(airfoil_data.Re)
        plot(airfoil_data.alpha_deg, airfoil_data.cl_matrix(:,i), '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        % Format Re using scientific notation or plain integers in LaTeX
        legend_list{end+1} = sprintf('$Re = %d$', airfoil_data.Re(i));
    end
    
    xlabel('Angle of Attack, $\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Lift Coefficient, $c_l$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Lift Curve} ($c_l$ vs $\alpha$)', 'Interpreter', 'latex', 'FontSize', 14);
    
    % Legend only needs to be called once if shared, but usually placed on the first plot
    legend(legend_list, 'Interpreter', 'latex', 'Location', 'Best', 'FontSize', 10);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % --- Subplot 2: Drag Coefficient ---
    subplot(1,2,2); hold on; grid on;
    for i = 1:length(airfoil_data.Re)
        plot(airfoil_data.alpha_deg, airfoil_data.cd_matrix(:,i), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    
    xlabel('Angle of Attack, $\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Drag Coefficient, $c_d$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Drag Polar} ($c_d$ vs $\alpha$)', 'Interpreter', 'latex', 'FontSize', 14);
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % Add a super title for Figure 3
    sgtitle('\textbf{3. Airfoil Aerodynamic curves}', 'Interpreter', 'latex');
    





    % =====================================================================
    % 4. Blade-element inflow geometry @x=0.75 for each mission segment
    % =====================================================================
    figure('Name', '4. Blade-Element inflow geometry comparison', 'Position', [100 100 1500 900]); 
    
    for i = 1:num_segments
        % -----------------------------------------------------------------
        % SUBPLOT: 2x2 Configuration
        % -----------------------------------------------------------------
        subplot(2, 2, i);
        hold on; axis equal; grid on; box on;
        
        % -- 1. Select Radial Station (75% Span) --
        [~, idx] = min(abs(outputs{i}.r - 0.75 * outputs{i}.R));
        
        % Extract local variables
        c_local     = outputs{i}.c(idx);
        theta_local = outputs{i}.theta(idx); 
        phi_local   = outputs{i}.phi(idx);
        alpha_local = outputs{i}.alpha(idx);
        
        % -- 2. Process Airfoil Geometry --
        x_foil = [prop.airfoil_data.x_u(:); flipud(prop.airfoil_data.x_l(:))];
        y_foil = [prop.airfoil_data.y_u(:); flipud(prop.airfoil_data.y_l(:))];
        
        % Rotate Airfoil by -theta_local 
        theta_rot = -theta_local; 
        R_mat = [cos(theta_rot), -sin(theta_rot); sin(theta_rot), cos(theta_rot)];
        coords_rot = R_mat * [x_foil'*c_local; y_foil'*c_local]; 
        x_rot = coords_rot(1, :);
        y_rot = coords_rot(2, :);
        
        % -- 3. Calculate Zero Lift Line (ZLL) --
        Fc = griddedInterpolant({prop.airfoil_data.alpha_deg*(pi/180), prop.airfoil_data.Re}, ...
                                prop.airfoil_data.cl_matrix, 'linear', 'linear');
        Re_loc = outputs{i}.Re_local(idx); 
        find_alpha0 = @(a) Fc(a, Re_loc);
        alpha_0 = fzero(find_alpha0, 0); 
        theta_zll = theta_local - alpha_0; % Note that alpha_0 < 0, thats why we subtract it (to actually add it)
                                           % theta_zll is a clockwise angle.
                                           % Later, we will rotate by -theta_zll
        
        % -- 4. Draw Reference Lines --
        
        % A. Rotor Plane
        line([-c_local*0.6, c_local*1.1], [0, 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.8);
        text(c_local*1.1, 0, 'Rotor Plane', 'FontSize', 8, 'VerticalAlignment', 'bottom', ...
             'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
        
        % B. Chord Line
        len_chord_fwd = c_local * 1.1; 
        line([0, len_chord_fwd*cos(theta_rot)], [0, len_chord_fwd*sin(theta_rot)], ...
            'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.0);
            
        len_chord_back = c_local * 0.6; 
        line([0, -0.9*len_chord_back*cos(theta_rot)], [0, -0.9*len_chord_back*sin(theta_rot)], ...
            'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.0);
        
        % C. Zero Lift Line (ZLL)
        len_zll_fwd = c_local * 0.6; 
        len_zll_aft = c_local * 0.9; 
        
        x_zll = [-len_zll_fwd, len_zll_aft] * cos(-theta_zll);
        y_zll = [-len_zll_fwd, len_zll_aft] * sin(-theta_zll);
        
        plot(x_zll, y_zll, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.2);
        
        text(x_zll(1), y_zll(1), 'ZLL ', 'FontSize', 7, 'Color', [0.5 0.5 0.5], 'Interpreter', 'latex', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'); 
        
        % D. Resultant Velocity Direction Line (W)
        phi_draw = -phi_local; 
        len_wind_fwd = c_local * 1.2;
        len_wind_back = c_local * 0.6;
        
        line([-0.9*len_wind_back*cos(phi_draw), 0.95*len_wind_fwd*cos(phi_draw)], ...
             [-0.9*len_wind_back*sin(phi_draw), 0.95*len_wind_fwd*sin(phi_draw)], ...
             'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.2);
        
        x_lbl = 0.6*len_wind_fwd * cos(phi_draw);
        y_lbl = 1.05*len_wind_fwd * sin(phi_draw);
        text(x_lbl, y_lbl, ' Resultant velocity', 'Color', [0, 0.4470, 0.7410], ...
             'FontSize', 9, 'Interpreter', 'latex', ...
             'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
        
        % -- 5. Draw Airfoil --
        fill(x_rot, y_rot, [0.2 0.2 0.2], 'FaceAlpha', 0.1, 'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.5);
        
        % -- 6. Velocity Triangle --
        
        % Visual Scale
        vis_len_Omegar = c_local * 0.7; 
        vis_len_V = vis_len_Omegar * tan(abs(phi_local));
        
        % === CHANGE 1: MOVED TRIANGLE LEFT ===
        % Changed x from -0.1 to -0.35 to shift it left
        vec_origin = [-c_local*0.35, -c_local*0.35]; 
        
        % Omega*r (Horizontal)
        quiver(vec_origin(1), vec_origin(2), vis_len_Omegar, 0, ...
            'Color', 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
            
        % V + vi (Vertical Down)
        quiver(vec_origin(1)+vis_len_Omegar, vec_origin(2), 0, -vis_len_V, ...
            'Color', [0.85, 0.32, 0.09], 'LineWidth', 1.5, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
            
        % u_R (Resultant)
        quiver(vec_origin(1), vec_origin(2), vis_len_Omegar, -vis_len_V, ...
            'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, 'MaxHeadSize', 0.3, 'AutoScale', 'off'); 
        
        % Labels
        text(vec_origin(1) + vis_len_Omegar/2, vec_origin(2), '$\Omega r$', ...
            'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 9);
        
        text(vec_origin(1) + vis_len_Omegar, vec_origin(2) - vis_len_V/2, ' $(v+v_i)$', ...
            'Color', [0.85, 0.32, 0.09], 'HorizontalAlignment', 'left', 'Interpreter', 'latex', 'FontSize', 9);
        
        text(vec_origin(1) + vis_len_Omegar*0.4, vec_origin(2) - vis_len_V*0.6, '$u_R$', ...
            'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'FontSize', 10, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'right');
        
        % Phi Angle in Triangle
        r_tri_arc = vis_len_Omegar * 0.4;
        t_tri_phi = linspace(0, phi_draw, 20);
        plot(vec_origin(1) + r_tri_arc*cos(t_tri_phi), vec_origin(2) + r_tri_arc*sin(t_tri_phi), ...
             'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.8);
        text(vec_origin(1) + r_tri_arc*1.15, vec_origin(2) + r_tri_arc*0.5*sin(phi_draw), '$\phi$', ...
             'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'FontSize', 8, 'HorizontalAlignment', 'left');
        
        % -- 7. Detailed Geometry Angles --
        r_theta = c_local * 0.50;
        r_phi   = c_local * 0.90;
        
        % A. Theta (Pitch)
        t_theta = linspace(0, theta_rot, 30);
        plot(r_theta * cos(t_theta), r_theta * sin(t_theta), 'k-', 'LineWidth', 0.8);
        text(r_theta * cos(theta_rot/2), r_theta * sin(theta_rot/2), ' $\theta$', ...
            'Color', 'k', 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 9);
        
        % B. Phi (Inflow)
        t_phi = linspace(0, phi_draw, 30);
        plot(r_phi * cos(t_phi), r_phi * sin(t_phi), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.8);
        text(r_phi * cos(phi_draw/2), r_phi * sin(phi_draw/2), ' $\phi$', ...
            'Color', [0, 0.4470, 0.7410], 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 9);
        
        % === CHANGE 2: ALPHA ON THE LEFT ===
        % Draw alpha at the Leading Edge, but on the LEFT side (upstream)
        r_alpha = c_local * 0.30; 
        
        % Extensions to the left are shifted by pi (180 deg)
        angle_chord_left = theta_rot + pi;
        angle_wind_left  = phi_draw + pi;
        
        % Draw Arc
        t_alpha = linspace(angle_wind_left, angle_chord_left, 20);
        plot(r_alpha * cos(t_alpha), r_alpha * sin(t_alpha), 'r-', 'LineWidth', 1.2);
        
        % Draw Text
        mid_alpha = (angle_wind_left + angle_chord_left) / 2;
        % Using 1.4 multiplier to push text slightly away from arc
        tx_alpha = 1.4*r_alpha * cos(mid_alpha);
        ty_alpha = 1.4*r_alpha * sin(mid_alpha);
        
        text(tx_alpha, ty_alpha, '$\alpha$', ...
            'Color', 'r', 'Interpreter', 'latex', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
        
        % -- 8. Formatting Subplot --
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 12);
        
        xlabel('Rotor Plane $x$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
        ylabel('Axial Direction $y$ [m]', 'Interpreter', 'latex', 'FontSize', 10);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 9);
        
        % -- ADJUST LIMITS --
        % Added tx_alpha to ensure alpha text is not cut off on the left
        all_x = [x_rot, vec_origin(1), x_lbl, x_zll(1), tx_alpha, c_local*1.1];
        all_y = [y_rot, vec_origin(2)-vis_len_V, y_lbl, y_zll(1)];
        
        margin = c_local * 0.15;
        xlim([min(all_x)-margin, max(all_x)+margin]);
        ylim([min(all_y)-margin, max(all_y)+margin]);
    end

    % Global Title (sgtitle)
    % Uses idx from the last loop iteration (assuming r/R station is same for all)
    sgtitle(sprintf('\\textbf{4. Blade-Element inflow geometry comparison at } $r/R = %.2f$', ...
            outputs{end}.r(idx)/outputs{end}.R), 'Interpreter', 'latex');
    


    



    % =======================================================================
    % 5. Aerodynamics of the blade: Lambda and Alpha distributions
    %    (with Re-dependent stall alpha limits from Fc surface + parabolic fit)
    % =======================================================================
    
    % --- Lambda and Alpha vs x ---
    figure('Name', '5. Lambda and Alpha comparison between flight segments', ...
           'Position', [100 100 1500 600]);
    
    % -----------------------------------------------------------------------
    % 1) Precompute stall limits for each segment (from Fc over alpha x Re),
    %    and collect global y-limits for consistent scaling across subplots
    % -----------------------------------------------------------------------
    Fc = prop.airfoil.Fc;
    
    vals_all = [];
    alpha_stall_max_all = [];
    alpha_stall_min_all = [];
    
    for i = 1:num_segments
    
        % Segment data
        lambda = outputs{i}.lambda;
        alpha  = outputs{i}.alpha;
        Re_local = outputs{i}.Re_local;
    
        % 1. Extract alpha grid points from Fc
        grid_alphas = Fc.GridVectors{1};
        num_alphas  = length(grid_alphas);
    
        % 2. Evaluate Fc for every alpha against every Re in this segment
        % cl_surface size: (num_alphas x numel(Re_local))
        cl_surface = Fc({grid_alphas, Re_local});
    
        % 3. Discrete indices for max/min Cl for each Re
        [~, idx_clmax] = max(cl_surface, [], 1);
        [~, idx_clmin] = min(cl_surface, [], 1);
    
        % 4. Refine using Parabolic Fit around the extrema (sub-grid)
        alpha_stall_clmax = zeros(size(Re_local));
        alpha_stall_clmin = zeros(size(Re_local));
    
        for k = 1:numel(Re_local)
    
            % ---- CLMAX refinement ----
            idx = idx_clmax(k);
            if idx > 1 && idx < num_alphas
                x_sub = grid_alphas(idx-1:idx+1);
                y_sub = cl_surface(idx-1:idx+1, k);
    
                p = polyfit(x_sub, y_sub, 2);      % y = a x^2 + b x + c
                a = p(1); b = p(2);
    
                if abs(a) > eps
                    alpha_refined = -b/(2*a);
    
                    % sanity: keep within one grid step of discrete peak
                    if abs(alpha_refined - grid_alphas(idx)) > (grid_alphas(idx) - grid_alphas(idx-1))
                        alpha_stall_clmax(k) = grid_alphas(idx);
                    else
                        alpha_stall_clmax(k) = alpha_refined;
                    end
                else
                    alpha_stall_clmax(k) = grid_alphas(idx);
                end
            else
                alpha_stall_clmax(k) = grid_alphas(idx);
            end
    
            % ---- CLMIN refinement ----
            idx = idx_clmin(k);
            if idx > 1 && idx < num_alphas
                x_sub = grid_alphas(idx-1:idx+1);
                y_sub = cl_surface(idx-1:idx+1, k);
    
                p = polyfit(x_sub, y_sub, 2);
                a = p(1); b = p(2);
    
                if abs(a) > eps
                    alpha_refined = -b/(2*a);
    
                    % sanity: keep within one grid step of discrete minimum
                    if abs(alpha_refined - grid_alphas(idx)) > (grid_alphas(idx) - grid_alphas(idx-1))
                        alpha_stall_clmin(k) = grid_alphas(idx);
                    else
                        alpha_stall_clmin(k) = alpha_refined;
                    end
                else
                    alpha_stall_clmin(k) = grid_alphas(idx);
                end
            else
                alpha_stall_clmin(k) = grid_alphas(idx);
            end
    
        end
    
        % Store into outputs for later use/inspection
        outputs{i}.alpha_stall_clmax_local = alpha_stall_clmax;
        outputs{i}.alpha_stall_clmin_local = alpha_stall_clmin;
    
        % Collect for global y-limits
        vals_all = [vals_all, lambda(:).', alpha(:).'];
        alpha_stall_max_all = [alpha_stall_max_all, alpha_stall_clmax(:).'];
        alpha_stall_min_all = [alpha_stall_min_all, alpha_stall_clmin(:).'];
    
    end
    
    % Global y-limits include lambda/alpha AND stall limits AND (optional) fixed airfoil limits
    global_min = min([vals_all, alpha_stall_min_all, alpha_stall_max_all, prop.airfoil.alpha_min, prop.airfoil.alpha_max]);
    global_max = max([vals_all, alpha_stall_min_all, alpha_stall_max_all, prop.airfoil.alpha_min, prop.airfoil.alpha_max]);
    
    padding = 0.1 * (global_max - global_min);
    if padding < 1e-6, padding = 0.1; end
    Aero_ylim = [global_min - padding, global_max + padding];
    
    % -----------------------------------------------------------------------
    % 2) Plotting loop: each subplot uses the “first code” content
    % -----------------------------------------------------------------------
    for i = 1:num_segments
        subplot(1, num_segments, i);
        hold on;
    
        % --- LEFT AXIS (Radians) ---
        yyaxis left
        f1 = plot(outputs{i}.mesh, outputs{i}.lambda, 'b', ...
            'LineWidth', 2, 'DisplayName', '$\lambda$');
        f2 = plot(outputs{i}.mesh, outputs{i}.alpha, 'r-', ...
            'LineWidth', 2, 'DisplayName', '$\alpha$');
    
        % Stall limits from Re-effect (same x mesh as the segment)
        f3 = plot(outputs{i}.mesh, outputs{i}.alpha_stall_clmax_local, 'r:', ...
            'LineWidth', 2, 'DisplayName', '$\alpha_{stall,max}$');
        f4 = plot(outputs{i}.mesh, outputs{i}.alpha_stall_clmin_local, 'r:', ...
            'LineWidth', 2, 'DisplayName', '$\alpha_{stall,min}$');
    
        ylim(Aero_ylim);
        ylabel('Value [rad]', 'Interpreter', 'latex', 'FontSize', 12);
        set(gca, 'YColor', 'k'); % keep left axis text black
    
        % --- RIGHT AXIS (Degrees & Red Color) ---
        yyaxis right
        ylim(rad2deg(Aero_ylim));
        ylabel('Value [deg]', 'Interpreter', 'latex', 'FontSize', 12);
        ax = gca;
        ax.YAxis(2).Color = 'r';
    
        % --- Styling & Titles ---
        xlabel('Radial coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        grid on;
    
        % Legend: show one entry for the stall-limit band (collapse f3/f4)
        legend([f1, f2, f3], { ...
            '$\lambda$', ...
            '$\alpha$', ...
            sprintf([' $\\alpha$ stall limits considering\n', ...
                     ' $Re$ number effect\n']) ...
            }, ...
            'Location', 'best', ...
            'Interpreter', 'latex', ...
            'FontSize', 6);
    
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
        hold off;
    end
    
    % Super-title
    sgtitle('\textbf{5. }$\mathbf{\lambda}$\textbf{ and }$\mathbf{\alpha}$\textbf{ vs }$\mathbf{x}$', ...
            'Interpreter', 'latex');


    





    % =======================================================================
    % 6. Force Distributions throughout the blade
    % =======================================================================
    
    % --- Thrust Distribution (dT/dx) ---
    figure('Name', '6. Thrust Distribution Comparison', 'Position', [100 100 1500 600]);
    
    % Find global Y-limits for consistent scaling
    dTdx_all = [];
    for i = 1:num_segments
        dTdx_all = [dTdx_all, outputs{i}.dTdx];
    end
    dTdx_max = max(dTdx_all);
    dTdx_min = min(dTdx_all);
    dTdx_ylim = [dTdx_min * 1.1, dTdx_max * 1.1];
    
    for i = 1:num_segments
        subplot(1, num_segments, i);
        plot(outputs{i}.mesh, outputs{i}.dTdx, 'LineWidth', 2);
        xlabel('Radial non-dimensional coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\frac{dT}{dx}$ [N/m]', 'Interpreter', 'latex', 'FontSize', 12);
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylim(dTdx_ylim);
        grid on;
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    end
    
    sgtitle('\textbf{6. Thrust Distribution} ($dT/dx$)', 'Interpreter', 'latex');
    
    
    % --- Tangential Force Distribution (dFt/dx) ---
    figure('Name', '7. Tangential Force Distribution Comparison', 'Position', [100 100 1500 600]);

    % Find global Y-limits for consistent scaling
    dFtdx_all = [];
    for i = 1:num_segments
        dFtdx_all = [dFtdx_all, outputs{i}.dFtdx];
    end
    dFtdx_max = max(dFtdx_all);
    dFtdx_min = min(dFtdx_all);
    dFtdx_ylim = [dFtdx_min * 1.1, dFtdx_max * 1.1];

    for i = 1:num_segments
        subplot(1, num_segments, i);
        plot(outputs{i}.mesh, outputs{i}.dFtdx, 'LineWidth', 2);
        xlabel('Radial non-dimensional coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\frac{dF_t}{dx}$ [N/m]', 'Interpreter', 'latex', 'FontSize', 12);
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylim(dFtdx_ylim);
        grid on;
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    end

    sgtitle('\textbf{7. Tangential Force Distribution} ($\frac{dF_t}{dx}$)', 'Interpreter', 'latex');
    
    



    % =======================================================================
    % 7. Stresses distributions throughout the blade
    % =======================================================================
    
    % --- Shared Definitions for the Reinforcement Line ---
    x_reinf_end = prop.x_finish_reinforcement;
    % Use a Cell Array for multiline LaTeX support
    reinfLabelStr = {'\textbf{End of structural reinforcement}', 'Focus only on stresses to the right'};
    % Style: Red dashed line
    reinfLineStyle = { 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5 };

    % --- Centrifugal stresses graph (\sigma_{cf}) ---
    figure('Name', '8. Centrifugal Stress Comparison', 'Position', [100 100 1500 600]);
    
    % Convert all stress to MPa and find global Y-limits for consistent scaling
    sigma_cf_all = [];
    for i = 1:num_segments
        sigma_cf_all = [sigma_cf_all, outputs{i}.sigma_cf/ 1e6]; % in MPa
    end
    sigma_cf_max = max(sigma_cf_all);
    sigma_cf_ylim = [0, sigma_cf_max * 1.1];
    
    for i = 1:num_segments
        subplot(1, num_segments, i);
        sigma_cf_MPa = outputs{i}.sigma_cf / 1e6;
        plot(outputs{i}.mesh, sigma_cf_MPa, 'b-', 'LineWidth', 2);
        hold on
        xline(x_reinf_end, ...
              'Label', reinfLabelStr, ...           % Explicit Label argument
              reinfLineStyle{:}, ...                % Expand style properties
              'Interpreter', 'latex', ...           % LaTeX support
              'LabelVerticalAlignment', 'bottom', ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 10, ...
              'HandleVisibility', 'off');
        hold off
        xlabel('Radial non-dimensional coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\sigma_{\mathrm{centrifugal}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylim(sigma_cf_ylim);
        grid on;
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    end
    sgtitle('\textbf{8. Centrifugal Stress Distribution}', 'Interpreter', 'latex');
    
    
    % --- Bending stresses graph (\sigma_{bend}) ---
    figure('Name', '9. Maximum Bending Stress Comparison', 'Position', [100 100 1500 600]);
    
    % Convert all stress to MPa and find global Y-limits for consistent scaling
    sigma_bend_Ix_all = [];
    sigma_bend_Iy_all = [];
    for i = 1:num_segments
        sigma_bend_Ix_all = [sigma_bend_Ix_all, outputs{i}.sigma_bend_max_Ix / 1e6];
        sigma_bend_Iy_all = [sigma_bend_Iy_all, outputs{i}.sigma_bend_max_Iy / 1e6];
    end
    sigma_bend_max = max([sigma_bend_Ix_all, sigma_bend_Iy_all]);
    sigma_bend_min = min([sigma_bend_Ix_all, sigma_bend_Iy_all]);
    sigma_bend_ylim = [sigma_bend_min * 1.1, sigma_bend_max * 1.1];
    
    for i = 1:num_segments
        subplot(1, num_segments, i);
        plot(outputs{i}.mesh, outputs{i}.sigma_bend_max_Ix / 1e6, 'b', 'LineWidth', 2, ...
             'DisplayName', '$\sigma_{\mathrm{bend}, Ix, \mathrm{max}}$')
        hold on
        plot(outputs{i}.mesh, outputs{i}.sigma_bend_max_Iy / 1e6, 'r--', 'LineWidth', 2, ...
             'DisplayName', '$\sigma_{\mathrm{bend}, Iy, \mathrm{max}}$')
        xline(x_reinf_end, ...
              'Label', reinfLabelStr, ...
              reinfLineStyle{:}, ...
              'Interpreter', 'latex', ...
              'LabelVerticalAlignment', 'bottom', ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 10, ...
              'HandleVisibility', 'off');
        hold off
        xlabel('Radial non-dimensional coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\sigma_{\mathrm{bend, max}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
        legend('Location', 'best', 'Interpreter', 'latex');
        grid on
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylim(sigma_bend_ylim);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    end
    sgtitle('\textbf{9. Maximum Bending Stress Distribution}', 'Interpreter', 'latex');
    
    
    % --- Total stresses graph (\sigma_{total max}) ---
    figure('Name', '10. Maximum Total Stress Comparison', 'Position', [100 100 1500 600]);
    
    % Convert all stress to MPa and find global Y-limits for consistent scaling
    sigma_total_all = [];
    for i = 1:num_segments
        sigma_total_all = [sigma_total_all, outputs{i}.sigma_total_max / 1e6]; % in MPa
    end
    sigma_total_max = max(sigma_total_all);
    sigma_total_ylim = [0, sigma_total_max * 1.1];
    
    for i = 1:num_segments
        subplot(1, num_segments, i);
        sigma_total_MPa = outputs{i}.sigma_total_max / 1e6;
        plot(outputs{i}.mesh, sigma_total_MPa, 'b-', 'LineWidth', 2);
        hold on
        xline(x_reinf_end, ...
              'Label', reinfLabelStr, ...
              reinfLineStyle{:}, ...
              'Interpreter', 'latex', ...
              'LabelVerticalAlignment', 'bottom', ...
              'LabelHorizontalAlignment', 'left', ...
              'FontSize', 10, ...
              'HandleVisibility', 'off');
        hold off
        xlabel('Radial non-dimensional coordinate, $x$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$\sigma_{\mathrm{total, max}}$ [MPa]', 'Interpreter', 'latex', 'FontSize', 12);
        title(titles{i}, 'Interpreter', 'latex', 'FontSize', 14);
        ylim(sigma_total_ylim);
        grid on;
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    end
    
    sgtitle('\textbf{10. Maximum Total Stress Distribution}', 'Interpreter', 'latex');





    % =====================================================================
    % 8. AIRFOIL Stress distributions (Critical section per flight segment)
    % =====================================================================
    
    % Create the figure
    figure('Name', '11. Airfoil Stress Distribution (Critical Sections)', ...
           'Position', [50 100 1800 500], 'Color', 'w', 'Visible', 'off');
    
    tlo = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    colormap(turbo);
    
    % Pre-allocate for common color scaling
    all_max_stresses = [];
    all_min_stresses = [];
    
    % Store data to avoid re-calculating inside the plotting loop if we want global limits
    plot_data = cell(1, num_segments);
    
    for k = 1:num_segments
        % --- A. Extract Variables for this segment ---
        out = outputs{k};
        
        dTdx  = out.dTdx;
        dFtdx = out.dFtdx;
        r     = out.r;
        dx    = out.dx;
        theta = out.theta;
        c     = out.c;
        rpm   = out.rpm;
        omega = rpm * 2 * pi / 60;
        
        % Airfoil structural properties (normalized)
        x_pts = prop.structural.x_pts;
        y_pts = prop.structural.y_pts;
        Area_sec_norm = prop.structural.A_norm;
        I_x_sec_norm  = prop.structural.I_x_norm;
        I_y_sec_norm  = prop.structural.I_y_norm;
        I_xy_sec_norm = prop.structural.I_xy_norm;
        
        % Scale properties to local chord c(r)
        Area_sec = Area_sec_norm .* (c.^2);
        I_x_sec  = I_x_sec_norm  .* (c.^4);
        I_y_sec  = I_y_sec_norm  .* (c.^4);
        I_xy_sec = I_xy_sec_norm .* (c.^4);
        
        % --- B. Structural Analysis (Loads & Moments) ---
        dT = dTdx .* dx; 
        dFt = dFtdx .* dx;
        dT_vec = dT(:)'; dFt_vec = dFt(:)'; r_vec = r(:)';
        
        % Integrate Loads (Tip to Root)
        S1_dT  = fliplr(cumsum(fliplr(dT_vec)));
        S1_dFt = fliplr(cumsum(fliplr(dFt_vec)));
        Sr_dT  = fliplr(cumsum(fliplr(dT_vec .* r_vec)));
        Sr_dFt = fliplr(cumsum(fliplr(dFt_vec .* r_vec)));
        
        % Bending Moments in Rotor Axes
        % M_T points to Leading Edge (negative X_rotor direction?) -> Check convention
        % Convention: x_rotor points to TE. Thrust creates moment about Tangential axis?
        % Standard Propeller convention:
        % M_thrust bending is about the chordwise axis (flapping).
        % M_torque bending is about the thrust axis (lagging) -> No, about thickness.
        
        M_T  = (Sr_dT - r_vec .* S1_dT);     % Moment due to Thrust
        M_Ft = (Sr_dFt - r_vec .* S1_dFt);   % Moment due to Drag/Tangential
        
        % Rotor Axes: x points to TE, y points Up (Thrust) - based on geometry plots
        M_x_rotor = -M_T;  
        M_y_rotor = M_Ft;  
        
        % Transform to Section Principal Axes
        % Rotation from Rotor to Airfoil is -theta (Clockwise)
        cos_nt = cos(-theta); sin_nt = sin(-theta);
        M_x =  M_x_rotor .* cos_nt + M_y_rotor .* sin_nt;
        M_y = -M_x_rotor .* sin_nt + M_y_rotor .* cos_nt;
        
        % --- C. Centrifugal Stress ---
        rho_mat_omega_sq = prop.rho_mat * omega^2;
        integrand = rho_mat_omega_sq .* Area_sec .* r;
        cum_from_root = cumtrapz(r, integrand); 
        F_cf = cum_from_root(end) - cum_from_root;
        
        safe_A = Area_sec; safe_A(safe_A < 1e-8) = 1e-8;
        sigma_cf = F_cf ./ safe_A; 
        
        % --- D. Calculate Total Stress at All Points ---
        sigma_total_pts = zeros(length(x_pts), length(r));
        sigma_max_per_sec = zeros(size(r));
        
        denom = I_x_sec .* I_y_sec - I_xy_sec.^2;
        denom(abs(denom) < 1e-30) = 1e-30;
        
        for j = 1:length(r)
             K_x = - (M_y(j) * I_x_sec(j) - M_x(j) * I_xy_sec(j)) / denom(j);
             K_y =   (M_x(j) * I_y_sec(j) - M_y(j) * I_xy_sec(j)) / denom(j);
             
             % Scaled geometry
             px = x_pts * c(j);
             py = y_pts * c(j);
             
             s_bend = K_x .* px + K_y .* py;
             sigma_total_pts(:, j) = sigma_cf(j) + s_bend;
             sigma_max_per_sec(j) = max(sigma_total_pts(:, j));
        end
        
        % --- E. Identify Critical Section ---
        [max_stress_val, idx_crit] = max(sigma_max_per_sec);
        
        % Store for plotting
        plot_data{k}.idx_crit = idx_crit;
        plot_data{k}.sigma_pts = sigma_total_pts(:, idx_crit);
        plot_data{k}.M_x = M_x(idx_crit);
        plot_data{k}.M_y = M_y(idx_crit);
        plot_data{k}.c   = c(idx_crit);
        plot_data{k}.theta = theta(idx_crit);
        plot_data{k}.r_R   = out.mesh(idx_crit);
        plot_data{k}.theta_p = prop.structural.theta_p; % Assuming constant airfoil
        
        all_max_stresses = [all_max_stresses, max(plot_data{k}.sigma_pts)];
        all_min_stresses = [all_min_stresses, min(plot_data{k}.sigma_pts)];
    end
    
    % Global Color Limits (MPa)
    cax = [min(all_min_stresses), max(all_max_stresses)] / 1e6;
    
    % --- Plotting Loop ---
    ax_handles = gobjects(1,4);
    for k = 1:num_segments
        ax = nexttile(tlo);
        ax_handles(k) = ax;
        
        % Prepare Title
        seg_name = mission.names{k};
        r_loc = plot_data{k}.r_R;
        t_str = sprintf('%s (crit. section at r/R=%.2f)', seg_name, r_loc);
        
        % Call Helper
        plotAirfoilStressSection(ax, plot_data{k}, t_str, cax, prop.structural);
    end
    
    % Shared Colorbar
    cb = colorbar(ax_handles(end));
    cb.Location = 'eastoutside';
    cb.Label.String = 'Normal Stress [MPa]';
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    
    % --- Create title and move it upwards ---
    % We use 'annotation' instead of 'sgtitle' for manual positioning.
    % Position format: [x_start, y_start, width, height] (Normalized 0 to 1)
    % y_start = 0.94 places it very close to the top edge.

    % Define the string (LaTeX syntax)
    titleText = '\textbf{11. Airfoil Stress Distribution at critical section for each mission segment}';
    
    % Create the annotation to act as the title
    annotation('textbox', [0, 0.90, 1, 0.05], ...  % [x, y, w, h] -> y=0.90 puts it at the top
        'String', titleText, ...
        'Interpreter', 'latex', ...      % <--- Renders it as LaTeX (matching your request)
        'FontSize', 16, ...              % <--- Standard title size (usually 14-16)
        'EdgeColor', 'none', ...         % No border
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
        











    % =====================================================================
    % 9. Performance plots of the optimal propeller: eta_p(J), ct(J), cp(J)
    % =====================================================================

    % When the performance plots of a propeller are obtained experimentally
    % in a wind tunel, the Measurement Procedure (The "Sweep") is:
    % Since the coefficients depend on the Advance Ratio J = V/nD, you must
    % vary J from 0 (static) to a value where thrust typically drops to 
    % zero. This is usually done by holding the RPM constant and changing 
    % the wind speed.
    % 1) Set Fixed RPM: The propeller is spun up to a constant rotational 
    %    speed (n). This is often chosen to match a specific Reynolds 
    %    number or the operational RPM of the real aircraft. 
    % 2) Static Test (J=0): The wind tunnel fan is off (V=0). Thrust and 
    %    Torque are measured. This gives the "static thrust" point on the 
    %    far left of the graph.
    % 3) Velocity Sweep: The wind tunnel speed (V) is gradually increased 
    %    in steps while the propeller controller maintains the fixed RPM.
    % 4) Data Acquisition: At each velocity step, the computer records:
    %    Freestream Velocity (V), Rotational Speed (n), Thrust (T), 
    %    Torque (Q), Air Density (rho) (calculated from P and Temp)
    % 5) Data Processing (Calculating the Coefficients): Once the raw data 
    %    (T, Q, V, n) is collected for the entire sweep, the 
    %    non-dimensional coefficients are calculated for each data point 
    %    using the standard formulas
    %
    % In this code, since the "experiments" are being run by a computer, we
    % don't care about reducing the number of experiments. We don't need to
    % fix RPM and only vary v.
    % We will vary both RPM and v, and that will provide a cloud of points
    % that really represent the behaviour of the propeller at different
    % Reynolds numbers (not just a pre-defined mission profile of the
    % propeller with fixed RPM and varying v).

    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Running BEMT for many v, rpm to generate PERFORMANCES PLOTS...\n'); 
    fprintf('(The code will try some unreal v, rpm. It is normal to see BEMT warning).\n');
    fprintf('-------------------------------------------------------------------------\n');

    pause(4); % Wait 4 seconds to let the user see that Warning is normal

    R = nonEmptyOutput.R;
    c = nonEmptyOutput.c;
    theta_deg = (180/pi)*nonEmptyOutput.theta;

    D = 2 * R;

    % --- 1. Define the velocity and rpm sweep ---
    max_rpm_tested = max(x_final(geometry_indices+1:end));

    v_vec = linspace(0, max(mission.v_ref), 30);
    rpm_vec = linspace(0.25*max_rpm_tested, max_rpm_tested, 30);

    n_rot_vec = rpm_vec./60;  % rev/s

    J_min_tested = 0;
    J_max_tested = max(v_vec) / ( min(n_rot_vec)*D );

    num_points = length(v_vec);
    num_rpms = length(rpm_vec);
    
    % Pre-allocate result vectors
    J_vec     = zeros(1, num_points*num_rpms);
    eta_p_vec = zeros(1, num_points*num_rpms);
    Ct_vec    = zeros(1, num_points*num_rpms);
    Cp_vec    = zeros(1, num_points*num_rpms);
    
    fprintf('Running Sweep (J = %.2f to J = %.2f) with %d points...\n', ...
            J_min_tested, J_max_tested, num_points*num_rpms);
    
    % --- 2. Execution Loop ---
    D = 2 * R;
    i_end = 0;
    for i = 1:num_points
        v_current = v_vec(i);
        
        for j = 1:num_rpms
            rpm_current = rpm_vec(j);
            n_rot_current = n_rot_vec(j);

            J_vec(i_end*num_rpms + j) = v_current / (n_rot_current*D);
            try
                % Execute BEMT
                res = runSABEMMT(v_current, rpm_current, prop.Nb, R, c, theta_deg, prop, env);
                
                % Store results
                eta_p_vec(i_end*num_rpms + j) = res.eta_p;
                Ct_vec(i_end*num_rpms + j)    = res.Ct_std; 
                Cp_vec(i_end*num_rpms + j)    = res.Cp_std;
                
            catch
                eta_p_vec(i_end*num_rpms + j) = NaN;
                Ct_vec(i_end*num_rpms + j)    = NaN;
                Cp_vec(i_end*num_rpms + j)    = NaN;
            end
        end
        i_end = i;
    end
    
    % --- 3. Process, Sort, and Scatter Plot (Separate Figures) ---
    
    % A. Filtering: Identify valid indices (remove NaN and Inf)
    valid_mask = ~isnan(J_vec) & ~isinf(J_vec) & ...
                 ~isnan(eta_p_vec) & ~isinf(eta_p_vec) & ...
                 ~isnan(Ct_vec) & ~isnan(Cp_vec) & ...
                 (J_vec < J_max_tested); 
    
    % Extract clean vectors
    J_final     = J_vec(valid_mask);
    eta_p_final = eta_p_vec(valid_mask);
    Ct_final    = Ct_vec(valid_mask);
    Cp_final    = Cp_vec(valid_mask);
    
    % B. Sorting
    [J_sorted, sort_order] = sort(J_final);
    eta_p_sorted = eta_p_final(sort_order);
    Ct_sorted    = Ct_final(sort_order);
    Cp_sorted    = Cp_final(sort_order);

    % Calculate J targets for the 4 conditions
    D = 2*R;
    J_targets = zeros(1,4);
    for k=1:4
        n_k = outputs{k}.rpm/60;
        J_targets(k) = mission.v_ref(k) / (n_k * D);
    end
    labels = {'Slow', 'Cruise', 'Climb', 'Fast'};

    
    % =======================================================================
    % COMBINED FIGURE OF PERFORMANCE PLOTS
    % =======================================================================
    figure('Name', '12. Performance Clouds', 'Position', [50 100 1800 600], 'Color', 'w');
    
    % --- SUBPLOT 1: Efficiency ---
    subplot(1, 3, 1);
    % Adjust Position (Make shorter and move down)
    pos = get(gca, 'Position');
    pos(4) = 0.53;  % Reduced from 0.60
    pos(2) = 0.10;  % Lowered from 0.15
    set(gca, 'Position', pos);
    
    hold on; grid on;
    scatter(J_sorted, eta_p_sorted, 10, [0 0.4470 0.7410], 'filled'); % Blue
    
    % Plot 4 vertical lines
    for k = 1:4
        xline(J_targets(k), '--r', 'LineWidth', 1.5, ...
            'Label', sprintf('%s ($J=%.2f$)', labels{k}, J_targets(k)), ...
            'LabelVerticalAlignment', 'bottom', ...
            'LabelOrientation', 'aligned', ...
            'Interpreter', 'latex', 'FontSize', 9);
    end
    
    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Efficiency, $\eta_p$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{A. Efficiency Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]); ylim([0, 1]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % --- SUBPLOT 2: Thrust Coefficient ---
    subplot(1, 3, 2);
    % Adjust Position (Make shorter and move down)
    pos = get(gca, 'Position');
    pos(4) = 0.53;  % Reduced from 0.60
    pos(2) = 0.10;  % Lowered from 0.15
    set(gca, 'Position', pos);
    
    hold on; grid on;
    scatter(J_sorted, Ct_sorted, 10, [0.8500 0.3250 0.0980], 'filled'); % Orange
    
    % Plot 4 vertical lines
    for k = 1:4
        xline(J_targets(k), '--r', 'LineWidth', 1.5, ...
            'Label', sprintf('%s', labels{k}), ... 
            'LabelVerticalAlignment', 'top', ...
            'LabelOrientation', 'aligned', ...
            'Interpreter', 'latex', 'FontSize', 9);
    end
    
    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Thrust Coeff., $C_t$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{B. Thrust Coefficient Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % --- SUBPLOT 3: Power Coefficient ---
    subplot(1, 3, 3);
    % Adjust Position (Make shorter and move down)
    pos = get(gca, 'Position');
    pos(4) = 0.53;  % Reduced from 0.60
    pos(2) = 0.10;  % Lowered from 0.15
    set(gca, 'Position', pos);
    
    hold on; grid on;
    scatter(J_sorted, Cp_sorted, 10, [0.9290 0.6940 0.1250], 'filled'); % Yellow
    
    % Plot 4 vertical lines
    for k = 1:4
        xline(J_targets(k), '--r', 'LineWidth', 1.5, ...
            'Label', sprintf('%s', labels{k}), ...
            'LabelVerticalAlignment', 'top', ...
            'LabelOrientation', 'aligned', ...
            'Interpreter', 'latex', 'FontSize', 9);
    end
    
    xlabel('Advance Ratio, $J$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Power Coeff., $C_p$', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{C. Power Coefficient Cloud}', 'Interpreter', 'latex', 'FontSize', 14);
    xlim([0, floor(J_max_tested * 10) / 10]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    
    % =======================================================================
    % ADD MAIN TITLE (ABOVE BOXES)
    % =======================================================================
    mainTitle = '\textbf{12. Performance plots of the optimal propeller: $\eta_p(J)$, $C_t(J)$, $C_p(J)$}';
    
    annotation('textbox', [0, 0.90, 1, 0.08], ... % Placed at very top
        'String', mainTitle, ...
        'Interpreter', 'latex', ...
        'FontSize', 16, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');

    % =======================================================================
    % DRAW 4 HEADER BOXES
    % =======================================================================
    
    % Adjusted box_y to be lower (was 0.82) to make room for the title above
    box_y = 0.74; 
    box_h = 0.12;
    box_w = 0.21;
    start_x = 0.04;
    gap = 0.03;
    
    for k = 1:4
        % Calculate X position
        x_pos = start_x + (k-1)*(box_w + gap);
        
        annotation('textbox', [x_pos, box_y, box_w, box_h], ...
            'String', titles{k}, ... % Ensure 'titles' variable is defined in your workspace
            'Interpreter', 'latex', ...
            'FontSize', 11, ... 
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'EdgeColor', 'black', ...
            'LineWidth', 1, ...
            'BackgroundColor', 'white', ...
            'FitBoxToText', 'off'); 
    end
    
    fprintf('Total points simulated: %d\n', length(J_vec));
    fprintf('Valid points plotted:   %d\n', length(J_sorted));





    % ====================================================================
    %  10. SAVE ALL PLOTS AS JPG TO 'RESULTS' FOLDER
    %  ====================================================================
    % ================= [BOTTOM BLOCK: START] =================
    fprintf('--- Finalizing and Saving Plots ---\n');
    
    % 1. Create results directory if it doesn't exist
    if ~exist('results', 'dir'), mkdir('results'); end
    
    % 2. Find all figures created by the code above
    allFigs = findobj('Type', 'figure');
    
    % 3. Sort them by Figure Number (1, 2, 3...) so they save in order
    if ~isempty(allFigs)
        [~, idx] = sort([allFigs.Number]);
        allFigs = allFigs(idx);
    end

    % 4. SAVE LOOP (Figures remain INVISIBLE here)
    % We use 'exportgraphics' or 'print' because 'saveas' requires visibility.
    fprintf('Saving images (this may take a moment)...\n');
    
    for i = 1:length(allFigs)
        hFig = allFigs(i);
        
        % A. Resize the invisible figure to be large (Full HD)
        % This ensures the saved image has good proportions/layout
        set(hFig, 'Units', 'pixels');
        set(hFig, 'Position', [100 100 1920 1080]); 
        
        % B. Construct Filename
        fName = hFig.Name;
        if isempty(fName), fName = sprintf('Figure_%d', hFig.Number); end
        safeName = regexprep(strrep(fName, ' ', '_'), '[^a-zA-Z0-9_]', '');
        fullPath = fullfile('results', [safeName, '.jpg']);
        
        % C. Save using specific renderers that support invisible figures
        try
            % OPTION 1: Modern MATLAB (R2020a+) - Best Quality
            % exportgraphics handles invisible OpenGL figures perfectly.
            exportgraphics(hFig, fullPath, 'Resolution', 300); 
        catch
            % OPTION 2: Legacy Fallback (R2019b and older)
            % 'print' forces a render even if invisible. 
            % -r150 defines 150 DPI resolution.
            print(hFig, fullPath, '-djpeg', '-r150');
        end
        
        fprintf('Saved: %s\n', fullPath);
    end

    % 5. SHOW LOOP (Restore visibility now that saving is done)
    fprintf('Displaying figures...\n');
    
    % Restore global default
    set(0, 'DefaultFigureVisible', 'on');
    
    for i = length(allFigs):-1:1
        hFig = allFigs(i);
        set(hFig, 'Visible', 'on');
        
        % Maximize the window on screen for the user
        try
            hFig.WindowState = 'maximized';
        catch
            % Fallback for older versions
            set(hFig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        end
        drawnow;
    end
    
    fprintf('Done. Check the ''results'' folder.\n');
    % ================= [BOTTOM BLOCK: END] =================
    
end







% ======================================================================
% Local helpers (filled arrowheads for elegant plotting)
% ======================================================================
function hArrow(ax, P0, P1, lw, col, headL, headW)
% Draw a straight arrow from P0 to P1 with a closed, filled head at P1.
    P0 = P0(:); P1 = P1(:);
    % Added 'ax' as the first argument to line()
    line(ax, [P0(1) P1(1)], [P0(2) P1(2)], 'Color', col, 'LineWidth', lw);
    hHeadPatch(ax, P0, P1, headL, headW, col);
end

function hDoubleArrow(ax, P0, P1, lw, col, headL, headW)
% Draw a straight line from P0 to P1 with closed, filled heads at both ends.
    P0 = P0(:); P1 = P1(:);
    line(ax, [P0(1) P1(1)], [P0(2) P1(2)], 'Color', col, 'LineWidth', lw);
    hHeadPatch(ax, P0, P1, headL, headW, col);
    hHeadPatch(ax, P1, P0, headL, headW, col);
end

function hHeadPatch(ax, P0, P1, headL, headW, col)
% Draw a closed filled triangular arrowhead at P1 pointing along P0->P1.
    P0 = P0(:); P1 = P1(:);
    v = P1 - P0;
    nv = norm(v);
    if nv < 1e-12
        return;
    end
    u = v / nv;
    p = [-u(2); u(1)];

    tip  = P1;
    base = P1 - headL*u;
    p1 = base + (headW/2)*p;
    p2 = base - (headW/2)*p;

    % Added 'ax' as the first argument to patch()
    patch(ax, [tip(1) p1(1) p2(1)], [tip(2) p1(2) p2(2)], col, ...
          'EdgeColor', col, 'FaceColor', col);
end

% ======================================================================
% Helper: Plot Single Stress Section (Fully Detailed Version)
% ======================================================================
function plotAirfoilStressSection(ax, data, titleStr, cax_MPa, struc)
    % axes(ax); % REMOVED to prevent pop-ups
    hold(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
    
    % --- 1. Extract Data ---
    x_pts = struc.x_pts;
    y_pts = struc.y_pts;
    c = data.c;
    theta = data.theta;
    theta_p = data.theta_p;
    sigma_MPa = data.sigma_pts / 1e6;
    Mx = data.M_x;
    My = data.M_y;
    
    % --- 2. Geometry & Rotation ---
    % Rotation matrix (Clockwise by theta to move from Rotor -> Airfoil)
    Rcw = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];
    
    % Airfoil coordinates in Rotor Frame
    x_sec = x_pts * c;
    y_sec = y_pts * c;
    XY = Rcw * [x_sec(:)'; y_sec(:)'];
    xR = XY(1, :)';
    yR = XY(2, :)';
    
    % Unit vectors for Airfoil Axes in Rotor Frame
    e1 = Rcw * [1; 0]; % x_airfoil direction
    e2 = Rcw * [0; 1]; % y_airfoil direction
    
    % Dynamic Limits (Scale view to fit airfoil)
    r_max = max(hypot(xR, yR));
    lim_span = 1.35 * r_max;
    xlim(ax, [-lim_span, lim_span]);
    ylim(ax, [-lim_span, lim_span]);
    
    % Arrow parameters
    headL = 0.06 * lim_span;
    headW = 0.035 * lim_span;
    
    % --- 3. Draw Rotor Axes (Reference) ---
    line(ax, [-lim_span, lim_span], [0 0], 'Color', 'k', 'LineWidth', 1);
    line(ax, [0 0], [-lim_span, lim_span], 'Color', 'k', 'LineWidth', 1);
    text(ax, 0.9*lim_span, 0.03*lim_span, '$x_{rot}$', 'Interpreter', 'latex', 'FontSize', 10);
    text(ax, 0.03*lim_span, 0.9*lim_span, '$y_{rot}$', 'Interpreter', 'latex', 'FontSize', 10);
    
    % --- 4. Draw Airfoil Axes & Labels ---
    axis_len = 0.85 * lim_span;
    colLoc = [0.4 0.4 0.4]; % Gray
    
    % Plot lines/arrows for local axes
    hArrow(ax, [0;0], axis_len*e1, 1.2, colLoc, headL, headW);
    hArrow(ax, [0;0], axis_len*e2, 1.2, colLoc, headL, headW);
    
    text(ax, 1.2*axis_len*e1(1), axis_len*e1(2), '$x_{airfoil}$', ...
         'Interpreter', 'latex', 'Color', colLoc, 'FontSize', 10, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(ax, axis_len*e2(1), axis_len*e2(2), '$y_{airfoil}$', ...
         'Interpreter', 'latex', 'Color', colLoc, 'FontSize', 10, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
         
    % --- 5. Draw Principal Axes ---
    % Rotate e1, e2 by theta_p to get principal directions
    e_xp =  cos(theta_p)*e1 + sin(theta_p)*e2;
    e_yp = -sin(theta_p)*e1 + cos(theta_p)*e2;
    
    colPr = [0 0 1]; % Blue
    L_xp = 0.6 * axis_len;
    L_yp = 0.45 * axis_len;
    
    hArrow(ax, [0;0], L_xp*e_xp, 1.2, colPr, headL, headW);
    hArrow(ax, [0;0], L_yp*e_yp, 1.2, colPr, headL, headW);
    
    % Offset text slightly to avoid overlap
    text(ax, L_xp*e_xp(1)*1.05, -L_xp*e_xp(2)*1.05, '$x_{princ}$', ...
         'Interpreter', 'latex', 'Color', colPr, 'FontSize', 9);
    
    % --- 6. Principal Angle Arc (theta_p) ---
    % Arc from -theta (Airfoil X) to (-theta + theta_p)
    ang_start = -theta;
    ang_end   = ang_start + theta_p;
    r_arc_p = 0.25 * lim_span; 
    
    t_p = linspace(ang_start, ang_end, 30);
    plot(ax, r_arc_p*cos(t_p), r_arc_p*sin(t_p), 'b', 'LineWidth', 1);
    
    text(ax, r_arc_p*1.4*cos((ang_start+ang_end)/2), ...
             r_arc_p*1.4*sin((ang_start+ang_end)/2), ...
             sprintf('$\\theta_p=%.1f^\\circ$', rad2deg(theta_p)), ...
             'Color', 'b', 'Interpreter', 'latex', 'FontSize', 9, ...
             'HorizontalAlignment', 'center');

    % --- 7. Geometric Pitch Angle Arc (theta) ---
    % Draw on the LEFT side (Tail-to-Tail) for clarity
    % Vector pointing to Rotor -X
    v_rot_neg = [-1; 0];
    % Vector pointing to Airfoil -X (Start of the angle measurement)
    v_sec_neg = -e1;
    
    arc_R = 1.15 * r_max; % Radius outside the airfoil
    
    % A. Draw the Negative Axis Extension (NEW)
    % Draw a dotted line along -x_airfoil so the arc has a visual target
    plot(ax, [0, 1.2*arc_R*v_sec_neg(1)], [0, 1.2*arc_R*v_sec_neg(2)], ...
         ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);

    % B. Calculate angles for arc
    a_start = atan2(v_rot_neg(2), v_rot_neg(1)); % pi
    a_end   = atan2(v_sec_neg(2), v_sec_neg(1));
    
    % Ensure shortest path logic for the arc drawing
    da = a_end - a_start;
    while da > pi,  da = da - 2*pi; end
    while da < -pi, da = da + 2*pi; end
    
    t_theta = linspace(a_start, a_start+da, 30);
    plot(ax, arc_R*cos(t_theta), arc_R*sin(t_theta), 'k-', 'LineWidth', 1);
    
    % Label
    theta_deg = rad2deg(theta);
    mid_t = a_start + da/2;
    text(ax, 1.1*arc_R*cos(mid_t), 1.1*arc_R*sin(mid_t), ...
         sprintf('$\\theta=%.1f^\\circ$', theta_deg), ...
         'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'right');

    % --- 8. Moment Vector Plotting ---
    Msec = [Mx; My];
    vMx = Mx * e1; 
    vMy = My * e2;
    vM  = vMx + vMy;
    Mmag = norm(vM);
    
    if Mmag > 1e-9
        % Scale arrow to be ~70% of viewport
        m_scale = 0.70 * lim_span / Mmag;
        VM  = m_scale * vM;
        VMx = m_scale * vMx;
        VMy = m_scale * vMy;
        
        colMom  = [0 0 0];       % Black Resultant
        colComp = [0.6 0.6 0.6]; % Grey Components
        
        % Component Arrows (Projections)
        hArrow(ax, [0;0], VMx, 1.5, colComp, headL, headW);
        hArrow(ax, [0;0], VMy, 1.5, colComp, headL, headW);
        
        % Construction lines (Dotted)
        plot(ax, [VMx(1), VM(1)], [VMx(2), VM(2)], ':', 'Color', colComp, 'LineWidth', 1);
        plot(ax, [VMy(1), VM(1)], [VMy(2), VM(2)], ':', 'Color', colComp, 'LineWidth', 1);
        
        % Resultant Arrow
        hArrow(ax, [0;0], VM, 2.0, colMom, headL*1.2, headW*1.2);
        
        % Label Resultant
        text(ax, VM(1), VM(2), '$\mathbf{M}_{res}$', ...
             'Interpreter', 'latex', 'FontSize', 11, ...
             'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 1);
    end

    % --- 9. Stress Scatter & Outline ---
    scatter(ax, xR, yR, 20, sigma_MPa, 'filled', 'MarkerEdgeColor', 'none');
    caxis(ax, cax_MPa);
    
    % Airfoil outline
    plot(ax, [xR; xR(1)], [yR; yR(1)], 'k-', 'LineWidth', 0.5);
    
    % Max Stress Marker
    [max_sig, idx_max] = max(sigma_MPa);
    plot(ax, xR(idx_max), yR(idx_max), 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    % --- 10. Stats Box (Updated decimals) ---
    % Matlab uses scientific notation only when needed:
    stats = {
        sprintf('$\\sigma_{max}=%.1f\\,\\mathrm{MPa}$', max_sig)
        sprintf('$M_x=%.4g\\,\\mathrm{Nm}$', Mx)
        sprintf('$M_y=%.4g\\,\\mathrm{Nm}$', My)
    };

    text(ax, -0.95*lim_span, -0.80*lim_span, stats, ...
         'Interpreter', 'latex', 'FontSize', 9, ...
         'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);
         
    title(ax, titleStr, 'Interpreter', 'latex', 'FontSize', 11);
    grid(ax, 'on');
    ax.TickLabelInterpreter = 'latex';
end

