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

%Material library has different units need to chage it
load materialLib
material = materialLib{1};
% material.E = 70E9;
% material.v = 1/3;


%% TODO:
design_params.stringer_pitch = 0.0364;
design_params.stringer_thickness = 3E-4;
design_params.stringer_web_height = 36E-3;
design_params.flange_to_web_ratio = 0.3;

bending_moment_dist = fit(limiting_loadcase_distributions.points', limiting_loadcase_distributions.bm', 'smoothingspline');

design_params.stringer_pitch = 100e-3;
%rib_stringer_func(geometry, material, design_params, bending_moment_dist, true)
%return;
mode = 1;

if mode == 0
    stringer_pitch_param_space = linspace(10E-3, 300E-3, 50);
    stringer_thickness_param_space = linspace(1E-4, 5E-3, 20);
    stringer_web_height_param_space = linspace(10E-3, 50E-3, 20);
    % total_area = rib_stringer_func(geometry, material, design_params, bending_moment_dist, false)

    optimisation.point_arr = [];
    optimisation.area_arr = [];
    point_arr = [];
    area_arr = [];
    parfor (p_idx = 1:numel(stringer_pitch_param_space), 16)
        test_space = design_params;
        for t = stringer_thickness_param_space
            for w = stringer_web_height_param_space
                test_space.stringer_pitch = stringer_pitch_param_space(p_idx);
                test_space.stringer_thickness = t;
                test_space.stringer_web_height = w;
                try
                    out = rib_stringer_func(geometry, material, test_space, bending_moment_dist, false);
                catch ME
                    %disp(['b = ', num2str(test_space.stringer_pitch), ' t_s = ', num2str(test_space.stringer_thickness), ' w = ', num2str(test_space.stringer_web_height)]);
                    %warning(ME.message);
                    out.total_volume = NaN;
                end

                point_arr = [point_arr, [stringer_pitch_param_space(p_idx); t; w]];
                area_arr = [area_arr, out.total_volume];
            end
        end
    end
    optimisation.point_arr = point_arr;
    optimisation.area_arr = area_arr;
elseif mode == 1
    f = @(x) optimiser_func(x, geometry, material, design_params, bending_moment_dist);
    %options = optimoptions(@fmincon,'Display','iter')
    options = optimset('TolCon',1e-18,'TolX',1e-19,'PlotFcns',@optimplotfval, 'UseParallel', true);
    x = fmincon(f, [0.0869; 0.0011; 0.0226], [], [], [], [], [5e-4; 5e-4; 5e-4], [1; 100e-3; 30e-3], [], options)
    
    design_params.stringer_pitch = x(1);
    design_params.stringer_thickness = x(2);
    design_params.stringer_web_height = x(3);
    disp(design_params);
    out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, true);
    improvePlot(gcf)
    return;
else 
    load('optimisation')
end
[optimisation.min_area, optimisation.min_idx] = min(area_arr);
save('optimisation')

figure;
hold on;
% plot(stringer_pitch_param_space, area_arr);
scatter3(optimisation.point_arr(2, :), optimisation.point_arr(3, :), optimisation.area_arr, 50, optimisation.point_arr(1, :), 'filled');
scatter3(optimisation.point_arr(2, optimisation.min_idx), optimisation.point_arr(3, optimisation.min_idx), optimisation.min_area, 100, 'k', 'x');
c = colorbar;
grid on;
title("Stringer Optimisation")
xlabel("Stringer Thickness [m]")
ylabel("Stringer Web Height [m]")
zlabel("Total Volume [m^3]")
hold off;

design_params.stringer_pitch = optimisation.point_arr(1, optimisation.min_idx);
design_params.stringer_thickness = optimisation.point_arr(2, optimisation.min_idx);
design_params.stringer_web_height = optimisation.point_arr(3, optimisation.min_idx);
disp(design_params);
out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, true);
improvePlot(gcf)

figure;
plot(out.rib_array, out.F_array);
xlabel("Spanwise Station [m]");
ylabel("F Factor");
grid on;


function output_vol = optimiser_func(x, geometry, material, design_params, bending_moment_dist)
    design_params.stringer_pitch = x(1);
    design_params.stringer_thickness = x(2);
    design_params.stringer_web_height = x(3);
    %disp(design_params)
    try
        out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, false);
        output_vol = out.total_volume;
    catch ME
        output_vol = NaN;
    end
    
end