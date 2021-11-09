clear 
clc

load("sizing.mat");
load("parameters.mat");

%% Takeoff 

Vinit = 0;
V_TO = 1.1 * sizing.takeoff.v_stall; %FAR-25

perf.Ka = (1.225 / (2 * W / Sref)) * (mu * Cl - Cd0 - (Cl^2) / (pi * AR *e));
perf.Kt = T / W - mu;
perf.Sg = (1 / (2 * g * perf.Ka)) * log(abs((perf.Kt + perf.Ka * Vtakeoff^2) / (perf.Kt + perf.Ka * Vinit^2));

%%
