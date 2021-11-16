clear
clc

g = 9.81;
newtons_to_lbs = 0.2248089431;

load('sizing.mat')
try
	load('wandb.mat')
catch e
	wandb.x_cg_front = 2.4;
	wandb.x_cg_aft = 2.6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%				THOUGHTS										%%
%%  	Single main wheel per strut								%%
%% 		Two nose wheels?										%%
%% 		Nose wheel should carry 8-15% of the weight				%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uc.tyres.psi = 60; % TODO:
uc.percent_weight_nose = 0.1;

uc.x_back_wheel = 2.8; % TODO:
uc.x_front_wheel = 0.8; % TODO:


B = uc.x_back_wheel - uc.x_front_wheel;
H = 1; % TODO:

M_f = uc.x_back_wheel - wandb.x_cg_front;
M_a = uc.x_back_wheel - wandb.x_cg_aft;

N_f = wandb.x_cg_front - uc.x_front_wheel;
N_a = wandb.x_cg_aft - uc.x_front_wheel;

% Checks
assert(M_a/B > 0.05, 'M_a/B not in valid range');
assert(M_f/B < 0.2, 'M_f/B not in valid range');

main_wheel_load = (1 - uc.percent_weight_nose) * sizing.W0 * newtons_to_lbs / 2;
nose_wheel_load = uc.percent_weight_nose*sizing.W0 * newtons_to_lbs / 2;

far25_safety_margin = 1.07;
max_static_load = @(W) far25_safety_margin*W*N_a/B;
max_static_load_nose = @(W) far25_safety_margin*W*M_f/B;
min_static_load_nose = @(W) far25_safety_margin*W*M_a/B;
dynamic_breaking_load_nose = @(W) 10*far25_safety_margin*H*W/(g*B);

max_main_wheel_load = max_static_load(main_wheel_load)
max_nose_wheel_load = max_static_load_nose(nose_wheel_load)
min_nose_wheel_load = min_static_load_nose(nose_wheel_load)
dynamic_nose_load = dynamic_breaking_load_nose(nose_wheel_load)

% Initial wheel size from Raymer pg. 344
%% Calculate main wheel diamater based on Raymer regression
uc.main_wheel.diameter_cm = 8.3*main_wheel_load^0.251;
uc.main_wheel.width_cm = 3.5*main_wheel_load^0.216;

% in inches
wheel_diam_in = uc.main_wheel.diameter_cm/2.54
wheel_width_in = uc.main_wheel.width_cm/2.54