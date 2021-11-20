clear
clc

load('design.mat')
load('sizing.mat')

newtons_to_lbf = 0.224808943;
installed_thrust_lbf = design.t_w * sizing.W0 * newtons_to_lbf

