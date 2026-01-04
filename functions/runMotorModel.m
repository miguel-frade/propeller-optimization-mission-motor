%% --------------- Motor model for T-Motor AS2304 1500kv ------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: runMotorModel
%   Project: propeller-optimization-mission-motor
%   Author: Miguel Frade
%   Affiliation: Universidad Politecnica de Madrid (At time of publication)
%   Date: December 2025
%   License: Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)
%
%   Description:
%       Analytical DC Brushless motor model. Calculates the required
%       Current and Voltage to drive a specific mechanical Load (Torque)
%       at a specific Speed (RPM).
%
%   Model Parameters (based on T-Motor AS2304):
%       - KV, Resistance (Rm), No-load current (I0)
%
%   Inputs:
%       Q_prop      - Torque required by propeller [Nm]
%       rpm         - Rotational speed [RPM]
%       motor_const - Motor constants struct
%
%   Outputs:
%       out.eta_m   - Motor efficiency
%       out.I_mot   - Required Current [A]
%       out.V_req   - Required Voltage [V]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this model, we just compute ideal outputs. The constraints of the
% battery not being able to provide sufficient voltage, the max current,
% etc. are specified in evaluateObjectivesAndConstraints.m

function out = runMotorModel (Q_prop, rpm, motor_const)
    
    % ESC efficiency is not considered in this code (assume eta_ESC = 1)

    % 1. Calculate Current Required to produce this Torque
    % Q = Kt * (I - I0)  => I = Q/Kt + I0
    I_mot = (Q_prop / motor_const.Kt) + motor_const.I0;
    
    % 2. Calculate Voltage Required to spin at this RPM with this Load
    % V_req = (RPM / KV) + I_mot * Rm
    V_req = (rpm / motor_const.KV) + (I_mot * motor_const.Rm);
    
    % 3. Calculate Motor Efficiency
    % P_mech_out = Q * omega
    omega = rpm * pi/30;
    P_out = Q_prop * omega;
    
    % P_elec_in = V_req * I_mot (Power consumed from battery/ESC)
    P_in = V_req * I_mot;
    
    if P_in > 0
        eta_m = P_out / P_in;
    else
        eta_m = 0;
    end

    out.eta_m = eta_m;
    out.I_mot = I_mot;
    out.V_req = V_req;
end