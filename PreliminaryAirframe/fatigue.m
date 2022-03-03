clear;
clc;
close all

load('limiting_loadcase_distributions');
load('wing_layout');
load('af_geometry');

fig_path = fullfile('./Figures/wing');

bending_moment_dist = fit(limiting_loadcase_distributions.points', limiting_loadcase_distributions.bm', 'smoothingspline');

comp_load_per_length = @(y) bending_moment_dist(y)./(geometry.box_width_func(y).*geometry.web_height_func(y));
test_space = linspace(0, geometry.semispan, 1000)';

delta = 1e-10;
panel_thickness_arr = repelem(wing_layout.panel_thickness, 2);
rib_arr = repelem(wing_layout.rib_array, 2);
rib_arr(2:2:end) = rib_arr(2:2:end) + delta;
rib_arr = [0, rib_arr(1:end-1)];

panel_thickness_fit = fit(rib_arr', panel_thickness_arr', 'linearinterp');
%plot(panel_thickness_fit, rib_arr, roundpanel_thickness_arr);

stresses = comp_load_per_length(test_space)./panel_thickness_fit(test_space);
[~, idx] = max(abs(stresses));
max_stress_1g = stresses(idx)/2.8 % TODO: CHECK if it should be divided by 2.5
max_stress_psi = max_stress_1g/6895

plot(test_space, round(panel_thickness_fit(test_space), 4)/1000, 'b');
grid on;
xlabel('y [m]');
ylabel('t [mm]');
improvePlot(gcf);
saveas(gcf, fullfile(fig_path, 'panel_thickness_wing'), 'epsc')
