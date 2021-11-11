% Wing Lift Analysis
%ieo18
%10th Nov 21 -
%AVD

%% Drag Analysis
%Parasite drag (zero-lift) & lift-induced

%% Inputs from other scripts
S_ratio=6; %=S_wet/S_ref


%% M_dd
%Drag divergence Mach number

%% Parasite Drag

% C_D_0 initial estimate (not very accurate) - light aircaraft, twin engine
aero_analysis.drag.C_d0_estimate=0.0045*S_ratio;

%% Parasite Drag: Component buildup