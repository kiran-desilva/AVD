clear
clc
close all


load('loadcase.mat')
load('cl_fit.mat')

spanwise_disc = linspace(0, 8.94/2, 1001);
delta = spanwise_disc(2) - spanwise_disc(1);
spanwise_disc = spanwise_disc(1:end-1) + delta;

tmp = loadcase{1};
A = wing_load(1.5, tmp.v, 100, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{1}.cm, false, true, "Load Case 1: ");

tmp = loadcase{2};
A = wing_load(1.5, tmp.v, 100, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{1}.cm, false, true, "Load Case 2: ");

A = wing_load(1.5, tmp.v, 100, polyval(cl_fit.poly, spanwise_disc), spanwise_disc, loadcase{1}.cm, true, true, "Landing: ");