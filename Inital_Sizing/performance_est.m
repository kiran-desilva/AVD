
load("sizing.mat");
load("parameters.mat");

%% Takeoff 

Vinit = 0;
v_stall_takeoff = sqrt(2 * sizing.W0 / (1.225 * Sref * aero_analysis.wing.takeoff_CLmax)); %need clmax takeoff
V_TO = aero_analysis.wing.v_takeoff_ms; %FAR-25  
mu = 0.03; % raymer dry asphalt runway 
rolltime = 3; %raymer

performance.Ka = (1.225 / (2 * wandb.W0 / wing.S_ref)) * (mu * high_lift_design.Clmax_takeoff - aero_analysis.drag.Cd0_takeoff - (high_lift_design.Clmax_takeoff^2) / (pi * wing_design.AR * aero_analyis.drag.e_takeoff));
performance.Kt = (powerplant.Thrust_max_takeoff / wandb.W0) - mu;
performance.Sg = (1 / (2 * 9.81 * performance.Ka)) * log(abs((performance.Kt + performance.Ka * V_TO^2) / (performance.Kt + performance.Ka * Vinit^2)));

performance.Sr = 3 * V_TO; 

H_obs = 35 / 3.2808; %meters
R = (1.15 * v_stall_takeoff)^2 / (0.2 * 9.81); %V_TR = 1.15*stall speed
performance.Str = sqrt(R^2 - (R - H_obs)^2); %transition distance

takeoff_CD = aero_analysis.drag.Cd0_takeoff + (aero_analysis.wing.takeoff_CLmax^2) /  (pi * sizing.AR *  aero_analysis.drag.e_takeoff); 
LoverD_TR = aero_analysis.wing.takeoff_CLmax / takeoff_CD; 
gamma_climb = asin((powerplant.Thrust_max_takeoff / wandb.W0) - 1 / LoverD_TR);
H_TR = R * (1 - cos(gamma_climb)); %no climb segment needed 

performance.Sto = 1.15 * (performance.Sg + performance.Sr + performance.Str); %FAR25 15% safety factor all engines operative

gamma_climb_1eng =  asin(0.5 * design_t_w - 1 / LoverD_TR); %one engine inoperative
gamma_min = 0.024; %FAR25 minimum climb for 2-engine a/c
G = gamma_climb_1eng - gamma_min;
CL_climb = 0.694 * sizing.takeoff.cl_max; %Gudmundsson
T_av = 0.75 * powerplant.Thrust_static_takeoff * ((5 + parameters.BPR) / (4 + parameters.BPR));
U = 0.01 * aero_analysis.wing.takeoff_CLmax + 0.02;  
performance.BFL = (0.863 / (1 + 2.3 * G)) * (((powerplant.Thrust_max_takeoff / wandb.W0) / (1.225 * 9.81 * CL_climb)) + H_obs) * (1 / (T_av / wandb.W0 - U) + 2.7) + 655;

%% Landing 
theta_apprch = 3; %deg from gudm.
v_stall_landing = sqrt(2 * wandb.Wland / (1.225 * Sref * aero_analysis.wing.landing_CLmax));
V_a = 1.3 * v_stall_landing;
V_f = 1.23 * v_stall_landing;
V_td = 1.15 * v_stall_landing; %errikos slides
tfr = 2; % mid sized ac?
n_land = 1.2; 
R_land = V_f^2 / ((n_land - 1) * 9.81);
H_f = R_land * (1 - cosd(theta_apprch)); %errikos slides
mul = 0.3; %errikos

performance.Sa = (H_obs - H_f) / tand(theta_apprch);
performance.Sf = 0.1512 * v_stall_landing^2 * sind(theta_apprch); 
performance.Sfr = tfr * V_td; 
performance.Ktl = ((2 * 0.5 * powerplant.installed_thrust_lbf * 4.44822) / wandb.wlanding - mul; %0.5 max takeoff thrust setting 2 engines -  note mjst be neative
performance.Ka_l = (1.225 / (2 *e wandb.Wland / wing.Sref)) * (mul * aero_analysis.wing.landing_CLmax - aero_analysis.drag.Cd0_landing - (aero_analysis.wing.landing_CLmax^2) / (pi * wing.AR *aero_analysis.drag.e_landing));
performance.Sb = (1 / 2 * 9.81 * performance.Ka_l) * log((performance.Ktl + performance.Ka_l * 0) / (performance.Ktl + performance.Ka_l * Vtd^2)); 

performance.SL  = 1.666 * (performance.Sa + performance.Sf + performance.Sfr + performance.Sb); %FAR25

%% Range and Endurance
%calculate new weight fractions for cruise and loiter segments
%performance.frac_cruise1 = exp(-(parameters.cruise_range_km * 1000 * cruise1_c / 3600) / (parameters.cruise_mach * 295.07 * cruise_LoverD));  %cruise 1 using breguet range
performance.frac_cruise2 = exp(-(parameters.alternate_range_km * 1000 * cruise2_c / 3600) / (parameters.cruise_mach * 295.07 * aero_analysis.wing.cruise_LoverD));  %cruise 2
performance.frac_loiter = exp(-(parameters.loiter_duration * 60 * loiter_c / 3600) / (aero_analysis.wing.loiter_LoverD)); % loiter using endurance eqn
new_Wx_W0 = 1; %initialise 

for i = 1:9
    if i ~= 3 || 6 || 7
        new_Wx_W0 = new_Wx_W0 * sizing.roskam.fuelfrac(i);
    end
end

new_Wx_W0_ext = new_Wx_W0;

Wx_W0 = (wandb.W0 - wandb.Wf) / wandb.W0; 
new_Wx_W0 = new_Wx_W0 * performance.frac_cruise2 * performance.frac_loiter; % not including cruise 1 segment
performance.frac_cruise1 = Wx_W0 / (new_Wx_W0); 
performance.cruise1_Range = (parameters.cruise_mach * 295.07 / cruise1_c) * aero_analysis.wing.cruise_LoverD * log(performance.frac_cruise1);

%extended range
Wx_W0_ext = (wandb.W0 - (wandb.Wf + W_pld * 0.453592 * 9.81)) / wandb.W0; %extended range replacing all pld weight with fuel
new_Wx_W0_ext = new_Wx_W0_ext * performance.frac_cruise2 * performance.frac_loiter; % not including cruise 1 segment
performance.frac_cruise1_ext = Wx_W0_ext / (new_Wx_W0_ext); % calculate new cruise fuel frac using original wx/w0 and the new loiter and crusie 2 fractions
performance.cruise1_Range = (parameters.cruise_mach * 295.07 / cruise1_c) * aero_analysis.wing.cruise_LoverD * log(performance.frac_cruise1_ext);

%% Point performance

performance.sigma_maxalt = (2 * sizing.fraction.before_cruise / (beta * T0)) * sqrt(

Trust_lapse = @(V,h) ; %Thrust lapse model
Drag = @(V,h,W);

