clear
clc
close all

load('fuselageLoading.mat');
fig_path = fullfile('./Figures/horizontaltail/loads');

%%
%Va_flight
[LoaddistVA ,ShearforceVA, BendingmomentVA, TorqueVA] = Horizontal_tail_load((fuselageLoading.Va_flight.Lt)*4.2);

[LoaddistVD ,ShearforceVD, BendingmomentVD, TorqueVD] = Horizontal_tail_load((fuselageLoading.Vd_flight.Lt)*3.75);

[HorizontalTail.LoaddistFO ,HorizontalTail.ShearforceFO, HorizontalTail.BendingmomentFO, HorizontalTail.TorqueFO] = Horizontal_tail_load((fuselageLoading.front_off.Lt)*1.5);

%[LoaddistOEI ,ShearforceOEI, BendingmomentOEI, TorqueOEI] = Horizontal_tail_load((fuselageLoading.oei.load)*1.5);

s_h = 2.3887; %m
y = linspace(-s_h/2,s_h/2,100);
HorizontalTail.y = y; 

save('HorizontalTail.mat', 'HorizontalTail');

    figure
    subplot(4,1,1)
    plot(y, HorizontalTail.LoaddistFO)
    %xlabel("y (m)", 'interpreter', 'Latex')
    title("Nose off: Load (N)", 'interpreter', 'Latex')
    grid on
    improvePlot(gcf)

    subplot(4,1,2)
    plot(y, HorizontalTail.ShearforceFO)
   % xlabel("y (m)", 'interpreter', 'Latex')
    title("Nose off: Shear Force (N)", 'interpreter', 'Latex')
    grid on
    improvePlot(gcf)

    subplot(4,1,3)
    plot(y, HorizontalTail.BendingmomentFO)
   % xlabel("y (m)", 'interpreter', 'Latex')
    title("Nose off: Bending Moments (Nm)", 'interpreter', 'Latex')
    grid on
    improvePlot(gcf)

    subplot(4,1,4)
    plot(y, HorizontalTail.TorqueFO)
    %label("y (m)", 'interpreter', 'Latex')
    title("Nose off: Torque (Nm)", 'interpreter', 'Latex')
    grid on
    improvePlot(gcf)

    
    saveas(gcf, fullfile(fig_path, 'Htailloads'), 'epsc')
