clear
clc

% Engine -  FJ442C

g = 9.81;
newtons_to_lbs = 0.2248089431;

load('sizing.mat')
try
	load('wandb.mat')
	load('fuselage_design')
catch e
	% TODO: CHECK ALL THE VARS HERE
	wandb.z_cg = 1.5;
	wandb.x_cg_front = 5.5;
	wandb.x_cg_aft = 6;

	lowest_spanwise_point.y = 6;
	lowest_spanwise_point.z = 0;

	tailstrike_point.x = 10;
	tailstrike_point.z = 1.3;
end

uc.nose_wheel.x = 0.8; % TODO:
uc.main_wheel.y = 3;
uc.main_wheel.x = 6.6; % TODO:
uc.main_wheel.z = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				THOUGHTS										%%
%%  	Single main wheel per strut								%%
%% 		Two nose wheels?										%%
%% 		Nose wheel should carry 8-15% of the weight				%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% VARIABLES I NEED
%	Furthest, lowest point spanwise on the aircraft (to calculate roll clearance)
%	Position of the tailstrike point, a.k.a furthest lowest point towards the back on the fuselage
%	(to calculate max tipback angle)


%% Overturn angle section
theta = atan(uc.main_wheel.y/(uc.main_wheel.x - uc.nose_wheel.x));
projected_side = uc.main_wheel.y * cos(theta);
uc.overturn_angle_deg = atand(wandb.z_cg/projected_side)
non_blocking_assert(uc.overturn_angle_deg < 63, 'overturn angle too large');

%% Roll angle section
roll_clearance = uc.main_wheel.z - (lowest_spanwise_point.y * tand(5) + lowest_spanwise_point.z)
non_blocking_assert(roll_clearance > 6*0.0254, 'not enough roll clearance')

%% Tipback angle section
tipback_angle = atand(tailstrike_point.z/(tailstrike_point.x - uc.main_wheel.x))
cg_main_wheel_angle_min = atand((uc.main_wheel.x - wandb.x_cg_aft)/wandb.z_cg)
cg_main_wheel_angle_max = atand((uc.main_wheel.x - wandb.x_cg_front)/wandb.z_cg)

non_blocking_assert(min(cg_main_wheel_angle_max, cg_main_wheel_angle_min) > max(15, tipback_angle), 'angle between cg and main wheel not sufficiently large')

uc.tyres.psi = 90; % TODO:
uc.percent_weight_nose = 0.1;

B = uc.main_wheel.x - uc.nose_wheel.x;

M_f = uc.main_wheel.x - wandb.x_cg_front;
M_a = uc.main_wheel.x - wandb.x_cg_aft;

N_f = wandb.x_cg_front - uc.nose_wheel.x;
N_a = wandb.x_cg_aft - uc.nose_wheel.x;

% Checks
non_blocking_assert(M_a/B > 0.05, 'M_a/B not in valid range');
non_blocking_assert(M_f/B < 0.2, 'M_f/B not in valid range');

%main_wheel_load = (1 - uc.percent_weight_nose) * sizing.W0 * newtons_to_lbs / 2;
%nose_wheel_load = uc.percent_weight_nose*sizing.W0 * newtons_to_lbs / 2;

far25_safety_margin = 1.07;
max_static_load = @(W) far25_safety_margin*W*N_a/B;
max_static_load_nose = @(W) far25_safety_margin*W*M_f/B;
min_static_load_nose = @(W) far25_safety_margin*W*M_a/B;
dynamic_breaking_load_nose = @(W) 10*far25_safety_margin*wandb.z_cg*W/(g*B);

W0_lbs = sizing.W0*newtons_to_lbs;
max_main_wheel_load_lbs = max_static_load(W0_lbs) / 2;
max_nose_wheel_load_lbs = max_static_load_nose(W0_lbs) / 2;
min_nose_wheel_load_lbs = min_static_load_nose(W0_lbs) / 2;
dynamic_nose_load_lbs = dynamic_breaking_load_nose(W0_lbs) / 2;
%{
max_main_wheel_load_lbs = max_static_load(main_wheel_load)
max_nose_wheel_load_lbs = max_static_load_nose(nose_wheel_load);
min_nose_wheel_load_lbs = min_static_load_nose(nose_wheel_load)
dynamic_nose_load_lbs = dynamic_breaking_load_nose(nose_wheel_load);
%}

max_nose_wheel_load_lbs = max_nose_wheel_load_lbs + dynamic_nose_load_lbs

uc.main_wheel.tyres.name = "Type III - 8.50-10"
uc.nose_wheel.tyres.name = "Type III - 5.00-4"

% Initial wheel size from Raymer pg. 344
%% Calculate main wheel diamater based on Raymer regression
%{
uc.main_wheel.diameter_cm = 8.3*main_wheel_load^0.251;
uc.main_wheel.width_cm = 3.5*main_wheel_load^0.216;

% in inches
wheel_diam_in = uc.main_wheel.diameter_cm/2.54
wheel_width_in = uc.main_wheel.width_cm/2.54
%}

uc.main_gear_shock_struts = 1;
uc.main_wheel.diameter_cm = 26.3*2.54;
uc.main_wheel.width_cm = 5; % TODO:
uc.main_wheel.rolling_radius_cm = 10.4*2.54;
uc.nose_wheel.diameter_cm = 13.25*2.54;
uc.nose_wheel.width_cm = 5; % TODO:
uc.nose_wheel.rolling_radius_cm = 5.2*2.54;

% Oleo sizing
landing_speed_ms = 10*0.3048;
shock_absorber_efficiency = 0.85; % TODO:
tire_efficiency = 0.47 % TODO:
gear_load_factor = 3; % TODO:

uc.main_wheel.oleo.stroke_m = calc_oleo_stroke(landing_speed_ms, shock_absorber_efficiency, gear_load_factor, tire_efficiency, uc.main_wheel.diameter_cm, uc.main_wheel.rolling_radius_cm)
uc.nose_wheel.oleo.stroke_m = calc_oleo_stroke(landing_speed_ms, shock_absorber_efficiency, gear_load_factor, tire_efficiency, uc.nose_wheel.diameter_cm, uc.nose_wheel.rolling_radius_cm)

uc.main_wheel.oleo.static_deflection = 2/3*uc.main_wheel.oleo.stroke_m;
uc.nose_wheel.oleo.static_deflection = 2/3*uc.nose_wheel.oleo.stroke_m;

uc.main_wheel.oleo.min_overall_length_m = 2.5*uc.main_wheel.oleo.stroke_m;
uc.nose_wheel.oleo.min_overall_length_m = 2.5*uc.nose_wheel.oleo.stroke_m;

uc.nose_wheel.oleo.diameter_in = 0.04*sqrt((max_nose_wheel_load_lbs + dynamic_nose_load_lbs));
uc.main_wheel.oleo.diameter_in = 0.04*sqrt(max_main_wheel_load_lbs);

uc.main_wheel.oleo
uc.nose_wheel.oleo

uc.main_gear_length_in = uc.main_wheel.diameter_cm*0.5/2.54 + uc.main_wheel.oleo.min_overall_length_m*39.37;
uc.nose_gear_length_in = uc.nose_wheel.diameter_cm*0.5/2.54 + uc.nose_wheel.oleo.min_overall_length_m*39.37;

%% Calculate frontal undercarriage area for Isobell
% Assume a box
main_wheel_box_width_m = max(uc.main_wheel.width_cm, uc.main_wheel.oleo.diameter_in*2.54)/100;
nose_wheel_box_width_m = max(uc.nose_wheel.width_cm, uc.nose_wheel.oleo.diameter_in*2.54)/100;

main_wheel_box_area_m_sq = uc.main_gear_length_in*2.54/100*main_wheel_box_width_m;
nose_wheel_box_area_m_sq = uc.nose_gear_length_in*2.54/100*nose_wheel_box_width_m;

uc.uc_frontal_area_m_sq = nose_wheel_box_area_m_sq + 2*main_wheel_box_area_m_sq;
uc.ratio_ngear = gear_load_factor;

save('uc', 'uc')

function oleo_stroke_m = calc_oleo_stroke(landing_speed_ms, shock_absorber_efficiency, gear_load_factor, tire_efficiency, wheel_diam_cm, wheel_rolling_rad_cm)
	tire_stroke_m = (0.5*(wheel_diam_cm) - wheel_rolling_rad_cm)/100;
	oleo_stroke_m = (landing_speed_ms^2/(2*9.81*shock_absorber_efficiency*gear_load_factor)) - tire_efficiency/ shock_absorber_efficiency * tire_stroke_m + 0.0254;
	non_blocking_assert(oleo_stroke_m > 0.2, 'Oleo stroke should be at least 20cm, will proceed with 20cm');
	oleo_stroke_m = max(oleo_stroke_m, 0.2);
end

function non_blocking_assert(cond, msg)
	if ~cond
		disp("WARNING!!!!")
		disp(msg);
		response = input("Do you still wish to continue? (y/n)" + newline, 's');
		if response == 'n'
			disp(); % crash the program
		end
	end
end