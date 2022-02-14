clear
clc
close all
addpath(fullfile('.', 'helper_funcs')); 

airfoilcoords.top = [1.00000     0.00000;
0.95033     0.00986;
0.90066     0.01979;
0.85090     0.02974;
0.80103     0.03935;
0.75107     0.04847;
0.70101     0.05686;
0.65086     0.06440;
0.60064     0.07085;
0.55035     0.07602;
0.50000     0.07963;
0.44962     0.08139;
0.39923     0.08139;
0.34884     0.07971;
0.29846     0.07658;
0.24881     0.07193;
0.19781     0.06562;
0.14757     0.05741;
0.09746     0.04672;
0.07247     0.04010;
0.04757     0.03227;
0.02283     0.02234;
0.01059     0.01588;
0.00580     0.01236;
0.00347     0.01010;
0.00000     0.00000];

airfoilcoords.bottom = [ 0 0;
0.00653     -0.00810;
0.00920     -0.00956;
0.01441     -0.01160;
0.02717     -0.01490;
0.05243     -0.01963;
0.07753     -0.02314;
0.10254     -0.02604;
0.15243     -0.03049;
0.20219     -0.03378;
0.25189     -0.03613;
0.30154     -0.03770;
0.35116     -0.03851;
0.40077     -0.03855;
0.45038     -0.03759;
0.50000     -0.03551;
0.54965     -0.03222;
0.59936     -0.02801;
0.64914     -0.02320;
0.69899     -0.01798;
0.74893     -0.01267;
0.79897     -0.00751;
0.84910     -0.00282;
0.89934     0.00089;
0.94967     0.00278;
1.00000     0.00000];

load('limiting_loadcase_distributions');
c_root = 1.767;
c_tip = 0.533;

top_af_fit = fit(airfoilcoords.top(:, 1), airfoilcoords.top(:, 2), 'smoothingspline');
bottom_af_fit = fit(airfoilcoords.bottom(:, 1), airfoilcoords.bottom(:, 2), 'smoothingspline');

x = linspace(0, 1, 100);
figure
hold on;
plot(x, top_af_fit(x), 'k');
plot(x, bottom_af_fit(x), 'k');
axis equal;
grid on;
hold off;

geometry.semispan = 8.94/2;
geometry.sweep_deg = 18.17;
geometry.spar.front_x_c = 0.2;
geometry.spar.rear_x_c = 0.77;
geometry.airfoil = NaN; %% TODO:
geometry.A0 = 0.0655;

geometry.c = @(y) -(c_root - c_tip)/geometry.semispan*y + c_root;

geometry.box_width_func = @(y) abs(geometry.spar.front_x_c - geometry.spar.rear_x_c)*geometry.c(y);
front_web_height = abs(top_af_fit(geometry.spar.front_x_c) - bottom_af_fit(geometry.spar.front_x_c));
rear_web_height = abs(top_af_fit(geometry.spar.rear_x_c) - bottom_af_fit(geometry.spar.rear_x_c));
geometry.web_height_func = @(y) max(front_web_height, rear_web_height)*geometry.c(y); 	%% TODO: Play with the max function, try min or avg instead
																						%% 	  max should however give largest stress

material.E = 70E9;
material.poisson_r = 1/3;

%% TODO:
design_params.stringer_pitch = 500E-3;
design_params.stringer_thickness = 1E-3;
design_params.stringer_web_height = 36E-3;
design_params.flange_to_web_ratio = 0.3;

bending_moment_dist = fit(limiting_loadcase_distributions.points', limiting_loadcase_distributions.bm', 'smoothingspline');

%design_params.stringer_pitch = 1e-3;
rib_stringer_func(geometry, material, design_params, bending_moment_dist, true);
return;

stringer_pitch_param_space = linspace(2E-3, 50E-3, 50);
% total_area = rib_stringer_func(geometry, material, design_params, bending_moment_dist, false)

area_arr = [];
for p = stringer_pitch_param_space
	test_space = design_params;
	test_space.stringer_pitch = p;
	area_arr = [area_arr, rib_stringer_func(geometry, material, test_space, bending_moment_dist, false)];
end

[min_area, min_idx] = min(area_arr);

figure;
hold on;
plot(stringer_pitch_param_space, area_arr);
scatter(stringer_pitch_param_space(min_idx), min_area, 20, 'k', 'x');
grid on;
title("Stringer Optimisation")
xlabel("Stringer Pitch [m]")
ylabel("Total Volume [m^3]")
hold off;

design_params.stringer_pitch = stringer_pitch_param_space(min_idx);
rib_stringer_func(geometry, material, design_params, bending_moment_dist, true);