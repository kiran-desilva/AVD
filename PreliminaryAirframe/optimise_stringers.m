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

                                                                                        
save('af_geometry', 'geometry')
%Material library has different units need to chage it
load materialLib
material = materialLib{1};
% material.E = 70E9;
% material.v = 1/3;


%% TODO:
design_params.stringer_pitch = 0.0869;
design_params.stringer_thickness = 0.0011;
design_params.stringer_web_height = 0.0226;
design_params.flange_to_web_ratio = 0.3;

bending_moment_dist = fit(limiting_loadcase_distributions.points', limiting_loadcase_distributions.bm', 'smoothingspline');

design_params.stringer_pitch = 100e-3;
%rib_stringer_func(geometry, material, design_params, bending_moment_dist, true)
%return;
mode = 1;

if mode == 0
    stringer_pitch_param_space = linspace(10E-3, 300E-3, 60);
    stringer_thickness_param_space = linspace(1E-4, 5E-3, 60);
    stringer_web_height_param_space = linspace(10E-3, 50E-3, 60);
    % total_area = rib_stringer_func(geometry, material, design_params, bending_moment_dist, false)

    optimisation.point_arr = [];
    optimisation.area_arr = [];
    point_arr = [];
    area_mesh_arr = zeros(numel(stringer_pitch_param_space), numel(stringer_thickness_param_space), numel(stringer_web_height_param_space));
    area_arr = [];
    num_t = numel(stringer_thickness_param_space);
    num_w = numel(stringer_web_height_param_space);
    parfor (p_idx = 1:numel(stringer_pitch_param_space), 16)
        test_space = design_params;
        for t_idx = 1:num_t
            for w_idx = 1:num_w
                test_space.stringer_pitch = stringer_pitch_param_space(p_idx);
                test_space.stringer_thickness = stringer_thickness_param_space(t_idx);
                test_space.stringer_web_height = stringer_web_height_param_space(w_idx);
                try
                    out = rib_stringer_func(geometry, material, test_space, bending_moment_dist, false);
                catch ME
                    %disp(['b = ', num2str(test_space.stringer_pitch), ' t_s = ', num2str(test_space.stringer_thickness), ' w = ', num2str(test_space.stringer_web_height)]);
                    %warning(ME.message);
                    out.total_volume = NaN;
                end

                point_arr = [point_arr, [stringer_pitch_param_space(p_idx); stringer_thickness_param_space(t_idx); stringer_web_height_param_space(w_idx)]];
                area_mesh_arr(p_idx, t_idx, w_idx) = out.total_volume;
                area_arr = [area_arr, out.total_volume];
            end
        end
    end
    optimisation.pitch_space = stringer_pitch_param_space;
    optimisation.thickness_space = stringer_thickness_param_space;
    optimisation.height_space = stringer_web_height_param_space;
    
    optimisation.point_arr = point_arr;
    optimisation.area_arr = area_arr;
    optimisation.area_mesh_arr = area_mesh_arr;
elseif mode == 1
    f = @(x) optimiser_func(x, geometry, material, design_params, bending_moment_dist, @fudger);
    %options = optimoptions(@fmincon,'Display','iter')
    options = optimset('TolCon',1e-18,'TolX',1e-19,'PlotFcns',@optimplotfval, 'UseParallel', false);
    x = fmincon(f, [100e-3; 0.0011; 0.0226], [], [], [], [], [5e-4; 5e-4; 5e-4], [1; 100e-3; 30e-3], [], options)
    
    design_params.stringer_pitch = x(1);
    design_params.stringer_thickness = x(2);
    design_params.stringer_web_height = x(3);
    disp(design_params);
    out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, true, @fudger);
    improvePlot(gcf);
    wing_layout = out;
    wing_layout.stringer_pitch = design_params.stringer_pitch;
    wing_layout.stringer_thickness = design_params.stringer_thickness;
    wing_layout.stringer_web_height = design_params.stringer_web_height;
    disp(['Max stress ', num2str(max(out.sigma_0)/1e6) ' MPa (Compressive)'])
    disp(['Material tensile yield stress ', num2str(material.tensile_yield/1e6) ' MPa'])
    disp(['Difference (Yield - max stress): ', num2str((material.tensile_yield - max(out.sigma_0))/1e6) ' MPa'])
    
    save('wing_layout', 'wing_layout')
  
    ctip = 0.4083;
    croot = 0.8166;

    geometry.semispan = 2.3887/2;
    
    geometry.sweep_deg = 26.65;
    geometry.spar.front_x_c = 0.15;
    geometry.spar.rear_x_c = 0.68;
    geometry.airfoil = NaN; %% TODO:
    geometry.A0 = 0.0757;
    
    geometry.c = @(y) ((ctip - croot)/geometry.semispan) * (abs(y) - geometry.semispan) + ctip; 
    
    airfoilcoords.top = [0.000000  0.000000;
      0.005000  0.009780;
      0.007500  0.011790;
      0.012500  0.014900;
      0.025000  0.020350;
      0.050000  0.028100;
      0.075000  0.033940;
      0.100000  0.038710;
      0.150000  0.046200;
      0.200000  0.051730;
      0.250000  0.055760;
      0.300000  0.058440;
      0.350000  0.059780;
      0.400000  0.059810;
      0.450000  0.057980;
      0.500000  0.054800;
      0.550000  0.050560;
      0.600000  0.045480;
      0.650000  0.039740;
      0.700000  0.033500;
      0.750000  0.026950;
      0.800000  0.020290;
      0.850000  0.013820;
      0.900000  0.007860;
      0.950000  0.002880;
      1.000000  0.000000;
      0.000000  0.000000;];
  
   airfoilcoords.bottom = [0.005000 -0.009780;
      0.007500 -0.011790;
      0.012500 -0.014900;
      0.025000 -0.020350;
      0.050000 -0.028100;
      0.075000 -0.033940;
      0.100000 -0.038710;
      0.150000 -0.046200;
      0.200000 -0.051730;
      0.250000 -0.055760;
      0.300000 -0.058440;
      0.350000 -0.059780;
      0.400000 -0.059810;
      0.450000 -0.057980;
      0.500000 -0.054800;
      0.550000 -0.050560;
      0.600000 -0.045480;
      0.650000 -0.039740;
      0.700000 -0.033500;
      0.750000 -0.026950;
      0.800000 -0.020290;
      0.850000 -0.013820;
      0.900000 -0.007860;
      0.950000 -0.002880;
      1.000000  0.000000];
  
    top_af_fit = fit(airfoilcoords.top(:, 1), airfoilcoords.top(:, 2), 'smoothingspline');
    bottom_af_fit = fit(airfoilcoords.bottom(:, 1), airfoilcoords.bottom(:, 2), 'smoothingspline'); 

    geometry.box_width_func = @(y) abs(geometry.spar.front_x_c - geometry.spar.rear_x_c)*geometry.c(y);
    front_web_height = abs(top_af_fit(geometry.spar.front_x_c) - bottom_af_fit(geometry.spar.front_x_c));
    rear_web_height = abs(top_af_fit(geometry.spar.rear_x_c) - bottom_af_fit(geometry.spar.rear_x_c));
    geometry.web_height_func = @(y) max(front_web_height, rear_web_height)*geometry.c(y);
    
    material = materialLib{1};
    
    load('HorizontalTail.mat');
    
    design_params.stringer_pitch = 0.0869;
    design_params.stringer_thickness = 0.0011;
    design_params.stringer_web_height = 0.01;
    design_params.flange_to_web_ratio = 0.3;
    
    tp_bending_moment_dist = fit(HorizontalTail.y', HorizontalTail.BendingmomentFO', 'smoothingspline');
    
    f_tail = @(x) optimiser_func(x, geometry, material, design_params, tp_bending_moment_dist, @fudger2);
    %options = optimoptions(@fmincon,'Display','iter')
    options = optimset('TolCon',1e-7,'TolX',1e-7,'PlotFcns',@optimplotfval, 'UseParallel', false);
    x = fmincon(f_tail, [120e-3; 0.001; 0.012], [], [], [], [], [9e-4; 9e-4; 9e-4], [1; 100e-3; 30e-3], [], options);

    design_params.stringer_pitch = x(1);
    design_params.stringer_thickness = x(2);
    design_params.stringer_web_height = x(3);
    disp(design_params);
    
    out = rib_stringer_func(geometry, material, design_params, tp_bending_moment_dist, true, @fudger2);
    improvePlot(gcf);
    
    tail_layout = out;
    tail_layout.stringer_pitch = design_params.stringer_pitch;
    tail_layout.stringer_thickness = design_params.stringer_thickness;
    tail_layout.stringer_web_height = design_params.stringer_web_height;
    
    disp(['Max stress ', num2str(max(out.sigma_0)/1e6) ' MPa (Compressive)'])
    disp(['Material tensile yield stress ', num2str(material.tensile_yield/1e6) ' MPa'])
    disp(['Difference (Yield - max stress): ', num2str((material.tensile_yield - max(out.sigma_0))/1e6) ' MPa'])
    save('tail_layout', 'tail_layout')
    
    return;
else 
    load('optimisation')
    
    [X, Y, Z] = meshgrid(optimisation.thickness_space, optimisation.pitch_space, optimisation.height_space);
    [min_ar, idx] = min(optimisation.area_mesh_arr, [], [3],'linear');

    figure;
    contourf(X(idx), Y(idx), min_ar, 5);
end

[optimisation.min_area, optimisation.min_idx] = min(area_arr);
save('optimisation', 'optimisation')

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
out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, true, @fudger);
improvePlot(gcf)

figure;
plot(out.rib_array, out.F_array);
xlabel("Spanwise Station [m]");
ylabel("F Factor");
grid on;

function output_vol = optimiser_func(x, geometry, material, design_params, bending_moment_dist, fudger)
    design_params.stringer_pitch = x(1);
    design_params.stringer_thickness = x(2);
    design_params.stringer_web_height = x(3);
    try
        out = rib_stringer_func(geometry, material, design_params, bending_moment_dist, false, fudger);
        output_vol = out.total_volume*material.rho;
    catch ME
        %disp(ME.message)
        output_vol = NaN;
    end
    
end

function fudge = fudger2(y)
    if y < 0.6
        fudge = 1.1;
    else
        fudge = 0.87;
    end
end

function fudge = fudger(y)
    if y < 3.5
        fudge = 0.95;
    else
        fudge = 0.85;
    end
end