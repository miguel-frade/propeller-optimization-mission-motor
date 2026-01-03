%% ------------------ Generate propeller geometry -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: getPropGeometry
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Reconstructs the physical propeller geometry (Radius, Chord 
%       distribution, and Twist distribution) from the normalized 
%       optimization variables (4 parameters of Bezier Cubics).
%       Scales normalized chord variables with R and uses shape functions 
%       (Bezier Cubics) to get the chord distribution and theta_deg 
%       distribution of the propeller.
%
%   Inputs:
%       x           - Design vector
%       prop_const  - Includes R_ref, c0_norm constraints, and mesh
%
%   Outputs:
%       c_dist      - Vector of chord lengths along the blade [m]
%       theta_deg   - Vector of twist angles along the blade [deg]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, c_dist, theta_deg] = getPropGeometry (x, prop_const)
    R = x(8); % vector x has the dimensional radius (NOT normalized)
    ksi = (prop_const.mesh - prop_const.x0) / (1 - prop_const.x0);
    omksi = 1 - ksi;
    theta = x(1).*omksi.^3 + 3*x(2).*ksi.*omksi.^2 + 3*x(3).*ksi.^2.*omksi + x(4).*ksi.^3;
    theta_deg = theta .* (180/pi);
    c0 = R * prop_const.c0_norm;
    c1 = R * x(5); c2 = R * x(6); c3 = R * x(7);
    c_dist = c0.*omksi.^3 + 3*c1.*ksi.*omksi.^2 + 3*c2.*ksi.^2.*omksi + c3.*ksi.^3;
end