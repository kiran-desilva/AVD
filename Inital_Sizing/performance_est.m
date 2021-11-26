
load("sizing.mat");
load("parameters.mat");

%% Takeoff 

Vinit = 0;
v_stall_takeoff = sqrt(2 * wandb.W0 / (1.225 * Sref * high_lift_design.Clmax_takeoff)); %need clmax takeoff
V_TO = 1.1 * v_stall_takeoff; %FAR-25  
mu = 0.03; % raymer dry asphalt runway 
rolltime = 3; %raymer

performance.Ka = (1.225 / (2 * wandb.W0 / wing_design.S_ref)) * (mu * high_lift_design.Clmax_takeoff - aero_analysis.drag.Cd0_takeoff - (high_lift_design.Clmax_takeoff^2) / (pi * wing_design.AR * aero_analyis.drag.e_takeoff));
performance.Kt = (powerplant.Thrust_max_takeoff / wandb.W0) - mu;
performance.Sg = (1 / (2 * 9.81 * performance.Ka)) * log(abs((performance.Kt + performance.Ka * V_TO^2) / (performance.Kt + performance.Ka * Vinit^2)));

performance.Sr = 3 * V_TO; 

H_obs = 35 / 3.2808; %meters
R = (1.15 * v_stall_takeoff)^2 / (0.2 * 9.81); %V_TR = 1.15*stall speed
performance.Str = sqrt(R^2 - (R - H_obs)^2); 

takeoff_CD = aero_analysis.drag.Cd0_takeoff + (aero_analysis.wing.takeoff_CLmax^2) /  (pi * sizing.AR *  aero_analysis.drag.e_takeoff); 
LoverD_TR = aero_analysis.wing.takeoff_CLmax / takeoff_CD; 
gamma_climb = asin(design_t_w - 1 / LoverD_TR);
H_TR = R * (1 - cos(gamma_climb)); %no climb segment needed 

performance.Sto = 1.15 * (performance.Sg + performance.Sr + performance.Str); %FAR25 15% safety factor all engines operative

gamma_climb_1eng =  asin(0.5 * design_t_w - 1 / LoverD_TR); %one engine operative
gamma_min = 0.024; %FAR25 minimum climb for 2-engine a/c
G = gamma_climb_1eng - gamma_min;
CL_climb = 0.694 * sizing.takeoff.cl_max; %Gudmundsson
T_av = 0.75 * 0.5 * design_t_w * sizing.W0 * ((5 + parameters.BPR) / (4 + parameters.BPR)); %one engine operative
U = 0.01 * sizing.cl_max + 0.02;  
perf.BFL = (0.863 / (1 + 2.3 * G)) * ((design_w_s / (1.225 * 9.81 * CL_climb)) + H_obs) * (1 / (T_av / sizing.W0 - U) + 2.7) + 655;

%% Landing 
theta_apprch = 3; %deg from gudm.
V_a = 1.3 * sizing.landing.v_stall;
V_f = 1.23 * sizing.landing.v_stall;
V_td = 1.15 * sizing.landing.v_stall;
tfr = 2; % mid sized ac?
n_land = 1.2; 
R_land = V_f^2 / ((n_land - 1) * 9.81);
H_f = R_land * (1 - cosd(theta_apprch)); %errikos slides
mul = 0.3; %errikos

perf.Sa = (H_obs - H_f) / tand(theta_apprch);
perf.Sf = 0.1512 * sizing.landing.v_stall^2 * sind(theta_apprch); 
perf.Sfr = tfr * V_td; 
perf.Ktl = design_t_w - mul; %t/w need to change?
perf.Ka_l = (1.225 / (2 * sizing.W0 / design_S_ref)) * (mul * sizing.landing.cl_max - sizing.landing.cd0 - (sizing.landing.cl_max^2) / (pi * sizing.AR *sizing.landing.e));
perf.Sb = (1 / 2 * 9.81 * perf.Ka_l) * log((perf.Ktl + perf.Ka_l * 0) / (perf.Ktl + perf.Ka_l * Vtd^2)); 

perf.SL  = 1.666 * (perf.Sa + perf.Sf + perf.Sfr + perf.Sb); %FAR25

%% Range and Endurance
%calculate new weight fractions for cruise and loiter segments
%perf.frac_cruise1 = exp(-(parameters.cruise_range_km * 1000 * cruise1_c / 3600) / (parameters.cruise_mach * 295.07 * cruise_LoverD));  %cruise 1 using breguet range
perf.frac_cruise2 = exp(-(parameters.alternate_range_km * 1000 * cruise2_c / 3600) / (parameters.cruise_mach * 295.07 * cruise_LoverD));  %cruise 2
perf.frac_loiter = exp(-(parameters.loiter_duration * 60 * loiter_c / 3600) / (loiter_LoverD)); % loiter using endurance eqn
new_Wx_W0 = 1; %initialise 

for i = 1:9
    if i ~= 3 || 6 || 7
        new_Wx_W0 = new_Wx_W0 * sizing.roskam.fuelfrac(i);
    end
end

new_Wx_W0_ext = new_Wx_W0;

Wx_W0 = (sizing.W0 - sizing.Wf) / sizing.W0; %use original W0 and Wf to calculate Wx/W0
new_Wx_W0 = new_Wx_W0 * perf.frac_cruise2 * perf.frac_loiter; % not including cruise 1 segment
perf.frac_cruise1 = Wx_W0 / (new_Wx_W0); % calculate new cruise fuel frac using original wx/w0 and the new loiter and crusie 2 fractions
perf.cruise1_Range = (parameters.cruise_mach * 295.07 / cruise1_c) * cruise_LoverD * log(perf.frac_cruise1);

%extended range
Wx_W0_ext = (sizing.W0 - (sizing.Wf + W_pld * 0.453592 * 9.81)) / sizing.W0; %extended range replacing all pld weight with fuel
new_Wx_W0_ext = new_Wx_W0_ext * perf.frac_cruise2 * perf.frac_loiter; % not including cruise 1 segment
perf.frac_cruise1_ext = Wx_W0_ext / (new_Wx_W0_ext); % calculate new cruise fuel frac using original wx/w0 and the new loiter and crusie 2 fractions
perf.cruise1_Range = (parameters.cruise_mach * 295.07 / cruise1_c) * cruise_LoverD * log(perf.frac_cruise1_ext);

%% Point performance

perf.sigma_maxalt = (2 * sizing.fraction.before_cruise / (beta * T0)) * sqrt(

Trust_lapse = @(V,h) ; %Thrust lapse model
Drag = @(V,h,W) ; 

