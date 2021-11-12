
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
H_TR = R * (1 - cos(gamma_climb)); %no climb segment needed 

perf.Sto = 1.15 * (perf.Sg + perf.Sr + perf.Str); %FAR25 15% safety factor all engines operative

gamma_climb_1eng =  asin(0.5 * design_t_w - 1 / LoverD_TR); %one engine operative
gamma_min = 0.024; %FAR25 minimum climb for 2-engine a/c
G = gamma_climb_1eng - gamma_min;
%CL_climb = 1.2 * 
T_av = 0.75 * 0.5 * design_t_w * sizing.W0 * ((5 + parameters.BPR) / (4 + parameters.BPR)); %one engine operative
U = 0.01 * sizing.cl_max + 0.02; 
%perf.BFL = (0.863 / (1 + 2.3 * G)) * ((design_w_s / (1.225 * 9.81 * CL_climb)) + H_obs) * (1 / (T_av / sizing.W0 - U) + 2.7) + 655;

%% Landing 

