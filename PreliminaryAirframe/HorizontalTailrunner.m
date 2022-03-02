clear
clc
close all

load('fuselageLoading');
load('Horizontal_tail_load.m')

%%
%Va_flight
[LoaddistVA ,ShearforceVA, BendingmomentVA, TorqueVA] = HTAILCALC(fuselageLoading.Va_flight.Lt*4.2);

[LoaddistVD ,ShearforceVD, BendingmomentVD, TorqueVD] = HTAILCALC(fuselageLoading.Vd_flight.Lt*3.75);

[HorizontalTail.LoaddistFO ,HorizontalTail.ShearforceFO, HorizontalTail.BendingmomentFO, HorizontalTail.TorqueFO] = HTAILCALC(fuselageLoading.front_off.Lt*1.5);

[LoaddistOEI ,ShearforceOEI, BendingmomentOEI, TorqueOEI] = HTAILCALC(fuselageLoading.oei.Lt*1.5);

save('HorizontalTail.mat', 'HorizontalTail');