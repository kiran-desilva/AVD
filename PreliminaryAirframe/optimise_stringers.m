clear
clc
close all

load('limiting_loadcase_distributions');
c_root = 1.767;
c_tip = 0.533;

syms y;
geometry.semispan = 8.94/2;
geometry.c(y) = -(c_root - c_tip)/geometry.semispan*y + c_root;
geometry.box_width_func(y) = 0*y;
geometry.web_height_func(y) = 0*y;
geometry.sweep_deg = 18.17;
geometry.spar.front_x_c = 0.1;
geometry.spar.rear_x_c = 0.77;
geometry.airfoil = NaN; %% TODO:

material.E = 70E9;
material.poisson_r = 1/3;

%% TODO:
design_params.stringer_pitch = 20E-3;
design_params.stringer_thickness = 1E-3;
design_params.stringer_web_height = 36E-3;
design_params.flange_to_web_ratio = 0.3;

bending_moment_dist = spline(limiting_loadcase_distributions.points, limiting_loadcase_distributions.bm);

rib_stringer_func(geometry, material, design_params, bending_moment_dist, true);