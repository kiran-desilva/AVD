% Wing Lift Analysis
%ieo18
%10th Nov 21 -
%AVD
clear
clc
%% Drag Analysis
%Parasite drag (zero-lift) & lift-induced

%% Inputs from other scripts
S_ratio=6; %=S_wet/S_ref

%% Re numbers
% Re number of each component @ cruise, TO, Landing, alternate cruise,
% loiter
%Re(1): cruise
%Re(2): TO
%Re(3): Landing
%Re(4): Alternate cruise
%Re(5): Loiter
rho=[1,1,1,1,1];
u=[1,1,1,1,1];
visc=[1,1,1,1,1];
l=[1,1,1,1,1]; %characteristic component lengths
aero_analysis.drag.Re=rho.*u.*l./visc; %array with all the Re numbers


%% M_dd
%Drag divergence Mach number

%% Parasite Drag

% C_D_0 initial estimate (not very accurate) - light aircaraft, twin engine
aero_analysis.drag.C_d0_estimate=0.0045*S_ratio;

%% Parasite Drag: Component buildup

%% Save outputs
