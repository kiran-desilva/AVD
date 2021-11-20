clear
clc

load('design.mat')
load('sizing.mat')

newtons_to_lbf = 0.224808943;
installed_thrust_lbf = design.t_w * sizing.W0 * newtons_to_lbf / 2


%% Estimating installed thrust
% source: http://www.ae.metu.edu.tr/~ae452/lecture1_propulsion.pdf

% Pressure recovery
F = 0.98 % For podded engine
C_ram = 1.35 % For subsonic engine
pressure_recovery_percent_loss = C_ram*(1 - F) 

% Bleed air
C_bleed = 2;
bleed_fraction = 0.03; % 1-5% range
bleed_air_percent_loss = C_bleed*bleed_fraction;



