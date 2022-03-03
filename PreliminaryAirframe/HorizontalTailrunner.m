clear
clc
close all

load('fuselageLoading.mat');

%%
%Va_flight
[LoaddistVA ,ShearforceVA, BendingmomentVA, TorqueVA] = Horizontal_tail_load((fuselageLoading.Va_flight.Lt)*4.2);

[LoaddistVD ,ShearforceVD, BendingmomentVD, TorqueVD] = Horizontal_tail_load((fuselageLoading.Vd_flight.Lt)*3.75);

[HorizontalTail.LoaddistFO ,HorizontalTail.ShearforceFO, HorizontalTail.BendingmomentFO, HorizontalTail.TorqueFO] = Horizontal_tail_load((fuselageLoading.front_off.Lt)*1.5);

%[LoaddistOEI ,ShearforceOEI, BendingmomentOEI, TorqueOEI] = Horizontal_tail_load((fuselageLoading.oei.load)*1.5);

save('HorizontalTail.mat', 'HorizontalTail');

    s_h = 2.3887; %m
    y = linspace(-s_h/2,s_h/2,100);

    figure
    plot(y, HorizontalTail.LoaddistFO)
    xlabel("y (m)", 'interpreter', 'Latex')
    ylabel("Horizontal Tail Load (N)", 'interpreter', 'Latex')
    grid on

    figure
    plot(y, HorizontalTail.ShearforceFO)
    xlabel("y (m)", 'interpreter', 'Latex')
    ylabel("Horizontal Tail Shear Force (N)", 'interpreter', 'Latex')
    grid on

    figure
    plot(y, HorizontalTail.BendingmomentFO)
    xlabel("y (m)", 'interpreter', 'Latex')
    ylabel("Horizontal Tail Bending Moments (Nm)", 'interpreter', 'Latex')
    grid on
    
    figure
    plot(y, HorizontalTail.TorqueFO)
    xlabel("y (m)", 'interpreter', 'Latex')
    ylabel("Horizontal Tail Torque (Nm)", 'interpreter', 'Latex')
    grid on