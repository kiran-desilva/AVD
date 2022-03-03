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
plot(y, ShearforceVA)
hold on
plot(y, ShearforceVD)
hold on
plot(y, HorizontalTail.ShearforceFO)
hold off
grid on
improvePlot(gcf)
xlabel("Span (m)")
ylabel("Shear Force (N)")
legend("Va", "Vd", "Front off")
saveas(gcf, fullfile(fig_path, 'tailshearforce'), 'epsc')

figure 
plot(y, BendingmomentVA)
hold on
plot(y, BendingmomentVD)
hold on
plot(y, HorizontalTail.BendingmomentFO)
hold off
grid on
improvePlot(gcf)
xlabel("Span (m)")
ylabel("Bending Moment (Nm)")
legend("Va Load", "Vd Load", "Front off Load")
saveas(gcf, fullfile(fig_path, 'tailbendingmoments'), 'epsc')

figure 
plot(y, TorqueVA)
hold on
plot(y, TorqueVD)
hold on
plot(y, HorizontalTail.TorqueFO)
hold off
grid on
improvePlot(gcf)
xlabel("Span (m)")
ylabel("Torque (Nm)")
legend("Va Load", "Vd Load", "Front off Load")
saveas(gcf, fullfile(fig_path, 'tailtorque'), 'epsc')



    figure
%     subplot(4,1,1)
%     plot(y, HorizontalTail.LoaddistFO)
%     %xlabel("y (m)", 'interpreter', 'Latex')
%     title("Nose off: Load (N)", 'interpreter', 'Latex')
%     grid on
%     improvePlot(gcf)

    subplot(3,1,1)
    plot(y, HorizontalTail.ShearforceFO)
    xlabel("y (m)")
    ylabel("SF (N)")
    grid on
    improvePlot(gcf)

    subplot(3,1,2)
    plot(y, HorizontalTail.BendingmomentFO)
    xlabel("y (m)")
    ylabel("BM (Nm)")
    grid on
    improvePlot(gcf)

    subplot(3,1,3)
    plot(y, HorizontalTail.TorqueFO)
    xlabel("y (m)")
    ylabel("Torque (Nm)")
    grid on
    improvePlot(gcf)

    
    saveas(gcf, fullfile(fig_path, 'FOtailload'), 'epsc')
