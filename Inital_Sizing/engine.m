clear
clc

load('design.mat')
load('sizing.mat')

newtons_to_lbf = 0.224808943;
installed_thrust_lbf = design.t_w * sizing.W0 * newtons_to_lbf / 2


%% Estimating installed thrust
% source: http://www.ae.metu.edu.tr/~ae452/lecture1_propulsion.pdf

% Pressure recovery
F = 0.96 % For podded engine
C_ram = 1.35 % For subsonic engine
pressure_recovery_percent_loss = C_ram*(1 - F) 

% Bleed air
C_bleed = 2;
bleed_fraction = 0.05; % 1-5% range
bleed_air_percent_loss = C_bleed*(bleed_fraction);

total_percent_loss = bleed_air_percent_loss + pressure_recovery_percent_loss
% total_percent_loss = 0.1;
uninstalled_thrust_lbf = installed_thrust_lbf/(1 - total_percent_loss)

%% Chosen engine PW615F
% Thrust -> 1350lbf
% BPR -> 2.8 [https://en.wikipedia.org/wiki/Pratt_%26_Whitney_Canada_PW600#cite_note-Engine_yearbook-1]
PW615F.uninstalled_thrust_lbf = 1350;
PW615F.BPR = 2.8;
PW615F.length_m = 1.252;
PW615F.basic_diam_m = 0.4064;
PW615F.dry_weight_kg = 140;

%% Rubber sizing
SF = uninstalled_thrust_lbf/PW615F.uninstalled_thrust_lbf % Scale factor
powerplant.required_uninstalled_thrust_lbf = uninstalled_thrust_lbf;
powerplant.BPR = PW615F.BPR;
powerplant.length_m = PW615F.length_m*(SF)^0.4;
powerplant.basic_diam_m = PW615F.basic_diam_m*(SF)^0.5;
powerplant.engine_weight_lb = PW615F.dry_weight_kg*(SF)^1.1*2.20462262;
aDvANCed_MaTeRiaLs_FuDgE_fACtoR = 0.08;
powerplant.sfc = 0.88*(1 - aDvANCed_MaTeRiaLs_FuDgE_fACtoR)*exp(-0.05*powerplant.BPR);
powerplant

save('powerplant', 'powerplant')


