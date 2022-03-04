clear
clc
close all

load('fuselageLoading.mat');
fig_path = fullfile('./Figures/horizontaltail/loads');

%%
%Va_flight
[HorizontalTail.LoaddistVA ,HorizontalTail.ShearforceVA, HorizontalTail.BendingmomentVA, HorizontalTail.TorqueVA] = Horizontal_tail_load((fuselageLoading.Va_flight.Lt)*4.2);

[LoaddistVD ,ShearforceVD, BendingmomentVD, TorqueVD] = Horizontal_tail_load((fuselageLoading.Vd_flight.Lt)*3.75);

[LoaddistFO ,ShearforceFO, BendingmomentFO, TorqueFO] = Horizontal_tail_load((fuselageLoading.front_off.Lt)*1.5);

%[LoaddistOEI ,ShearforceOEI, BendingmomentOEI, TorqueOEI] = Horizontal_tail_load((fuselageLoading.oei.load)*1.5);

s_h = 2.3887; %m
y = linspace(-s_h/2,s_h/2,100);
HorizontalTail.y = y; 

save('HorizontalTail.mat', 'HorizontalTail');


figure 
plot(y, HorizontalTail.ShearforceVA, 'r')
hold on
plot(y, ShearforceVD, 'g')
hold on
plot(y, ShearforceFO, 'b')
hold off
grid on
improvePlot(gcf)
xlabel("Span [m]")
ylabel("Shear Force [N]")
legend("Load Case 1", "Load Case 2", "Load Case 3")
saveas(gcf, fullfile(fig_path, 'tailshearforce'), 'epsc')

figure 
plot(y, HorizontalTail.BendingmomentVA, 'r')
hold on
plot(y, BendingmomentVD, 'g')
hold on
plot(y, BendingmomentFO, 'b')
hold off
grid on
improvePlot(gcf)
xlabel("Span [m]")
ylabel("Bending Moment [Nm]")
legend("Load Case 1", "Load Case 2", "Load Case 3")
saveas(gcf, fullfile(fig_path, 'tailbendingmoments'), 'epsc')

figure 
plot(y, HorizontalTail.TorqueVA, 'r')
hold on
plot(y, TorqueVD, 'g')
hold on
plot(y, TorqueFO, 'b')
hold off
grid on
improvePlot(gcf)
xlabel("Span [m]")
ylabel("Torque [Nm]")
legend("Load Case 1", "Load Case 2", "Load Case 3")
saveas(gcf, fullfile(fig_path, 'tailtorque'), 'epsc')



    figure
%     subplot(4,1,1)
%     plot(y, HorizontalTail.LoaddistFO)
%     %xlabel("y (m)", 'interpreter', 'Latex')
%     title("Nose off: Load (N)", 'interpreter', 'Latex')
%     grid on
%     improvePlot(gcf)

    subplot(3,1,1)
    plot(y, HorizontalTail.ShearforceVA)
    xlabel("y [m]")
    ylabel("SF [N]")
    grid on
    improvePlot(gcf)

    subplot(3,1,2)
    plot(y, HorizontalTail.BendingmomentVA)
    xlabel("y [m]")
    ylabel("BM [Nm]")
    grid on
    improvePlot(gcf)

    subplot(3,1,3)
    plot(y, HorizontalTail.TorqueVA)
    xlabel("y [m]")
    ylabel("Torque [Nm]")
    grid on
    improvePlot(gcf)

    
    saveas(gcf, fullfile(fig_path, 'VAtailload'), 'epsc')
