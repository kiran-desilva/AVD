clear
clc
close all


load('loadcase.mat')
load('cl_fit.mat')

spanwise_disc = linspace(0, 8.94/2, 1001);
delta = spanwise_disc(2) - spanwise_disc(1);
spanwise_disc = spanwise_disc(1:end-1) + delta;

fuel_weight_kg = 486.748216106;
%% NOTE: LOAD CASE 2 SEEMS MOST LIMITING
tmp = loadcase{1};
A1 = wing_load(1.5*2.5, tmp.v, fuel_weight_kg, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{1}.cm, false, true, "Load Case 1: ");

tmp = loadcase{2};
A2 = wing_load(1.5*2.5, tmp.v, fuel_weight_kg, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{2}.cm, false, true, "Load Case 2: ");

A3 = wing_load(1.5*2.5, tmp.v, fuel_weight_kg, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{1}.cm, true, true, "Landing: ");

torque_dists = [A1.torque; A2.torque; A3.torque];
[~, max_torque_idx] = max(abs(torque_dists), [], [1], 'linear');

bm_dists = [A1.bm; A2.bm; A3.bm];
[~, max_bm_idx] = max(abs(bm_dists), [], [1], 'linear');

shear_dists = [A1.shear; A2.shear; A3.shear];
[~, max_shear_idx] = max(abs(shear_dists), [], [1], 'linear');

limiting_loadcase_distributions = A1;
limiting_loadcase_distributions.torque = torque_dists(max_torque_idx);
limiting_loadcase_distributions.bm = bm_dists(max_bm_idx);
limiting_loadcase_distributions.shear = shear_dists(max_shear_idx);

subplot(3, 1, 1);
plot(limiting_loadcase_distributions.points, limiting_loadcase_distributions.shear);
grid on;
xlabel("y [m]")
ylabel("Shear Force [N]")
subplot(3, 1, 2);
plot(limiting_loadcase_distributions.points, limiting_loadcase_distributions.bm)
grid on;
xlabel("y [m]")
ylabel("Bending Moment [Nm]")
subplot(3, 1, 3);
plot(limiting_loadcase_distributions.points, limiting_loadcase_distributions.torque)
grid on;
xlabel("y [m]")
ylabel("Torque [Nm]")
improvePlot(gcf)

save('limiting_loadcase_distributions', 'limiting_loadcase_distributions');

