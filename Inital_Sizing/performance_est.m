
load("sizing.mat");
load("parameters.mat");

%% Takeoff 

Vinit = 0;
V_TO = 1.1 * sizing.takeoff.v_stall; %FAR-25
mu = 0.03; % raymer dry asphalt runway
rolltime = 3; %raymer

perf.Ka = (1.225 / (2 * sizing.W0 / design_S_ref)) * (mu * sizing.takeoff.cl_max - sizing.takeoff.cd0 - (sizing.takeoff.cl_max^2) / (pi * sizing.AR *sizing.takeoff.e));
perf.Kt = design_t_w - mu;
perf.Sg = (1 / (2 * 9.81 * perf.Ka)) * log(abs((perf.Kt + perf.Ka * V_TO^2) / (perf.Kt + perf.Ka * Vinit^2)));

perf.Sr = 3 * V_TO; 

H_obs = 35 / 3.2808;
R = (1.15 * sizing.takeoff.v_stall)^2 / (0.2 * 9.81); %V_TR = 1.15*stall speed
perf.Str = sqrt(R^2 - (R - H_obs)^2); 

sizing.takeoff.cd = sizing.takeoff.cd0 + (sizing.takeoff.cl_max^2) /  (pi * sizing.AR * sizing.takeoff.e); 
LoverD_TR = sizing.takeoff.cl_max / sizing.takeoff.cd; 
gamma_climb = asin(design_t_w - 1 / LoverD_TR);
H_TR = R * (1 - cos(gamma_climb)); 

perf.Sto = 1.15 * (perf.Sg + perf.Sr + perf.Str); %FAR25 15% safety factor
%%
