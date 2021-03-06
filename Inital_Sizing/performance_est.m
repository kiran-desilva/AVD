clear 
clc

load("sizing.mat");
load("parameters.mat");
load("aero_analysis.mat");
load("wing.mat");
load("powerplant.mat");
powerplant.sfc = 0.6545;
%% Takeoff 

Vinit = 0;
%v_stall_takeoff = %sqrt(2 * sizing.W0 / (1.225 * Sref * aero_analysis.wing.takeoff_CLmax)); %need clmax takeoff
V_TO = aero_analysis.wing.v_takeoff_ms; %FAR-25  
mu = 0.03; % raymer dry asphalt runway 
rolltime = 3; %raymer

Cl_gr = aero_analysis.summary.cl_alpha_TO * (wing.i_w*(pi/180) - aero_analysis.summary.zero_AoA.TO_rad);
performance.Ka = (1.225 / (2 * sizing.W0 / wing.Sref)) * (mu * Cl_gr - aero_analysis.drag.cd0(3) - ((Cl_gr^2) / (pi * wing.Ar * aero_analysis.summary.e_wing)));
performance.Kt = (2 * powerplant.installed_thrust_lbf * 4.44822 / sizing.W0) - mu;
performance.Sg = (1 / (2 * 9.81 * performance.Ka)) * log(abs((performance.Kt + performance.Ka * V_TO^2) / (performance.Kt + performance.Ka * Vinit^2)));

performance.Sr = 3 * V_TO; 

H_obs = 35 / 3.2808; %meters
R = (1.15 * aero_analysis.wing.stall_TO_no_safety)^2 / (0.2 * 9.81); %V_TR = 1.15*stall speed
performance.Str = sqrt(R^2 - (R - H_obs)^2); %transition distance

transition_CD = aero_analysis.drag.cd0(3) + (aero_analysis.summary.cl_transition^2) /  (pi * wing.Ar *  aero_analysis.summary.e_wing); 
LoverD_TR = aero_analysis.summary.cl_transition / transition_CD; 
gamma_climb = asin((2 * powerplant.installed_thrust_lbf * 4.44822 / sizing.W0) - 1 / LoverD_TR);
H_TR = R * (1 - cos(gamma_climb)); %no climb segment needed CHECK THIS

performance.Sto = 1.15 * (performance.Sg + performance.Sr + performance.Str); %FAR25 15% safety factor all engines operative

gamma_climb_1eng =  asin((powerplant.installed_thrust_lbf* 4.44822 - Drag_model_TO(V_TO,0,sizing.W0) )/ sizing.W0); %one engine inoperative
gamma_min = 0.024; %FAR25 minimum climb for 2-engine a/c
G = gamma_climb_1eng - gamma_min;
CL_climb = 0.694 * aero_analysis.summary.cl_max_TO; %Gudmundsson
T_av = 0.75 * powerplant.installed_thrust_lbf * 4.44822 * ((5 + parameters.BPR) / (4 + parameters.BPR));
U = 0.01 * aero_analysis.summary.cl_max_TO + 0.02;  
performance.BFL = (0.863 / (1 + 2.3 * G)) * (((powerplant.installed_thrust_lbf * 4.44822 / sizing.W0) / (1.225 * 9.81 * CL_climb)) + H_obs) * (1 / (T_av / sizing.W0 - U) + 2.7) + 655;

%% Landing 

theta_apprch = 3; %deg from gudm.
%Wland_noalternatecruise = sizing.fraction.end_cruise_1 * 0.99 * sizing.W0;
v_stall_landing = aero_analysis.wing.stall_landing_no_safety;
V_a = 1.3 * v_stall_landing;
V_f = 1.23 * v_stall_landing;
V_td = 1.15 * v_stall_landing; %errikos slides
tfr = 2; % mid sized ac?
n_land = 1.2; 
R_land = V_f^2 / ((n_land - 1) * 9.81);
H_f = R_land * (1 - cosd(theta_apprch)); %errikos slides
mul = 0.5; %raymer
H_obs_land = 50 / 3.2808; %meters

Cl_brakedist = aero_analysis.summary.cl_alpha_approach * (wing.i_w*(pi/180) - aero_analysis.summary.zero_AoA.Land_rad);
performance.Sa = (H_obs_land - H_f) / tand(theta_apprch);
performance.Sf = R_land * sind(theta_apprch); 
performance.Sfr = tfr * V_td; 
performance.Ktl = ((2 * 0.15 * powerplant.installed_thrust_lbf * 4.44822) / (0.8*sizing.W0)) - mul; %0.15 max takeoff thrust setting 2 engines -  note mjst be neative
performance.Ka_l = (1.225 / (2 * 0.8*sizing.W0 / wing.Sref)) * (mul * Cl_brakedist - aero_analysis.drag.cd0(4) - (Cl_brakedist^2) / (pi * wing.Ar *aero_analysis.summary.e_wing));
performance.Sb = (1 / (2 * 9.81 * performance.Ka_l)) * log((performance.Ktl + performance.Ka_l * 0) / (performance.Ktl + performance.Ka_l * V_td^2)); 

performance.SL  = 1.666 * (performance.Sa + performance.Sf + performance.Sfr + performance.Sb); %FAR25

%% Range and Endurance
%calculate new weight fractions for cruise and loiter segments
%performance.frac_cruise1 = exp(-(parameters.cruise_range_km * 1000 * cruise1_c / 3600) / (parameters.cruise_mach * 295.07 * cruise_LoverD));  %cruise 1 using breguet range
performance.frac_cruise2 = exp(-(parameters.alternate_range_km * 1000 * powerplant.sfc / 3600) / (parameters.cruise_mach * 295.07 * aero_analysis.summary.l_d_cruise));  %cruise 2
performance.frac_loiter = exp(-(parameters.loiter_duration * 60 * powerplant.sfc / 3600) / (aero_analysis.summary.l_d_loiter)); % loiter using endurance eqn
new_Wx_W0 = 1; %initialise 

for i = 1:9
    if i ~= 3 && i ~= 6 && i ~= 7
        disp(i)
        new_Wx_W0 = new_Wx_W0 * sizing.roskam.fuelfrac(i);
    end
end

new_Wx_W0_ext = new_Wx_W0;

Wx_W0 = (sizing.W0 - sizing.Wf) / sizing.W0; 
new_Wx_W0 = new_Wx_W0 * performance.frac_cruise2 * performance.frac_loiter; % not including cruise 1 segment
performance.frac_cruise1 = Wx_W0 / (new_Wx_W0); 
performance.cruise1_Range = (parameters.cruise_mach * 295.07 / (powerplant.sfc/3600)) * aero_analysis.summary.l_d_cruise * log(1/performance.frac_cruise1);

%extended range
Wx_W0_ext = (sizing.W0 - (sizing.Wf + 15 * 9.81)) / sizing.W0; %extended range replacing all pld weight (15kg) with fuel
new_Wx_W0_ext = new_Wx_W0_ext * performance.frac_cruise2 * performance.frac_loiter; % not including cruise 1 segment
performance.frac_cruise1_ext = Wx_W0_ext / (new_Wx_W0_ext); % calculate new cruise fuel frac using original wx/w0 and the new loiter and crusie 2 fractions
performance.cruise1_Range_exted = (parameters.cruise_mach * 295.07 / (powerplant.sfc/3600)) * aero_analysis.summary.l_d_cruise * log(1/performance.frac_cruise1_ext);
extrange = [0,performance.cruise1_Range,performance.cruise1_Range_exted];
pld = [90,90,0];

figure
plot(extrange,pld)

%% Point performance

%performance.sigma_maxalt = (2 * (sizing.fraction.before_cruise * sizing.W0) / (1.439 * powerplant.installed_thrust_lbf .* 4.482 .* 2)) * sqrt(

Ps_fts = 0; %initialise

figure;
hold on

M__ = linspace(0, 1, 100);
h__ = linspace(0, 13800, 100);
Ps_wanted = [0:1:10];

Drag_model(0.75, 12300, sizing.fraction.before_cruise * sizing.W0)
[M_mat, h_mat] = meshgrid(M__, h__);
f = @(M, h) M .* atmos(h,2) .* ((thrust_lapse(h*3.28084,M,powerplant.BPR) .* powerplant.installed_thrust_lbf .* 4.482 .* 2 - Drag_model(M,h,sizing.fraction.before_cruise * sizing.W0)) ./ (sizing.fraction.before_cruise * sizing.W0));

for i = 1:length(M__)
    for j = 1:length(h__)
        %f = @(M, h) M .* atmos(h,2) .* ((((atmos(h,4)./1.225).^0.7).* 1.439 .* powerplant.installed_thrust_lbf .* 4.482 .* 2 - Drag_model(M,h,sizing.fraction.before_cruise * sizing.W0)) ./ (sizing.fraction.before_cruise * sizing.W0));

        %Ps = f(M_mat, h_mat);
        %contour(M_mat, h_mat, Ps, Ps_wanted)
        
        Ps(j,i) = f(M__(i),h__(j));
    end
    Ms(i) = sqrt((2 * sizing.fraction.before_cruise * sizing.W0)/(atmos(h__(i),4) * wing.Sref * aero_analysis.summary.cl_max_wing_clean)) / atmos(h__(i),2);%stall
end

Vx_TAS = sqrt(1.225/atmos(13716,4)) * sqrt((2 * 1 * sizing.fraction.before_cruise * sizing.W0)/(1.225 * wing.Sref)) * (1 / ( pi * aero_analysis.summary.e_wing * wing.Ar * aero_analysis.drag.cd0(1)))^0.25;
Mx = Vx_TAS / atmos(13716,2);

Ps_serv = [2.54 10000];

contour(M_mat,h_mat,Ps,Ps_wanted,'Showtext', 'on')
hold on
plot(Ms,h__,'r')
hold on
plot(0.75, 12192,'*m', 'MarkerSize', 12)
hold on
plot(0.78,12192,'xm', 'MarkerSize', 12)
hold on
plot(Mx,13716,'or', 'MarkerSize', 12)
hold on
contour(M_mat,h_mat,Ps,Ps_serv,'--')
hold off
xlabel("Mach number",'interpreter','Latex', "Fontsize",25)
ylabel("Altitude (m)",'interpreter','Latex', "Fontsize",25)
grid on

legend("Specific Excess Power","Stall line","Cruise", "M = 0.78 at cruise alt", "$V_x$ at 45000ft requirement","Service ceiling", "interpreter","Latex", "Fontsize",15,'location','southeast')
%{
for i = 1:10
    Ps_fts = Ps_fts + 500; 
    Ps = Ps_fts * 0.3048; %ft/s to m/s
    %f = @(M, h) ((M .* atmos(h,2))./ (sizing.fraction.before_cruise * sizing.W0)) .* ((thrust_lapse(h*3.28084,M,powerplant.BPR) .* powerplant.installed_thrust_lbf .* 4.482 .* 2 - Drag_model(M,h,sizing.fraction.before_cruise * sizing.W0))) - Ps;
    %fimplicit(@(M,h) f(M,h), [0, 1, 0, 13800]);
    %fimplicit(f);
    %fimplicit(f, [0 1 0 13800])
end
%}
hold off

performance
