%% Aerodynamic Analysis
%  Summary script (with plots)
clear
clc
close all

load 'aero_analysis.mat'
load 'tailplane.mat'
load 'wing.mat'
load 'sizing.mat'

%(1): wing
%(2): HLD
%(3): horizontal stabiliser (tailplane)


x=[-0.1:0.01:1]; %radians

%% Extract data

aero_analysis.summary.cl_max_wing_clean=aero_analysis.wing.Cl_max_wing;
aero_analysis.summary.cl_max_approach= aero_analysis.wing.Cl_max_wing+aero_analysis.wing.HLD.delta_cl_max(2);
aero_analysis.summary.cl_max_TO=aero_analysis.wing.Cl_max_wing+aero_analysis.wing.HLD.delta_cl_max(1);

aero_analysis.summary.cl_alpha_wing_cruise_clean=aero_analysis.wing.Cl_alpha(1);
aero_analysis.summary.cl_alpha_wing_max_clean=aero_analysis.wing.Cl_alpha(2);
aero_analysis.summary.cl_alpha_wing_TO=aero_analysis.wing.HLD.cl_alpha_flaps(1);
aero_analysis.summary.cl_alpha_wing_approach=aero_analysis.wing.HLD.cl_alpha_flaps(2);
aero_analysis.summary.cl_alpha_wing_loiter_clean=aero_analysis.wing.Cl_alpha(5);
aero_analysis.summary.cl_alpha_tail_cruise_clean=aero_analysis.tail.Cl_alpha(1);
aero_analysis.summary.cl_alpha_tail_max_clean=aero_analysis.tail.Cl_alpha(2);
aero_analysis.summary.cl_alpha_tail_TO=aero_analysis.tail.Cl_alpha(3);
aero_analysis.summary.cl_alpha_tail_approach=aero_analysis.tail.Cl_alpha(4);
aero_analysis.summary.cl_alpha_tail_loiter_clean=aero_analysis.tail.Cl_alpha(5);

aero_analysis.summary.cl_alpha_cruise_clean=aero_analysis.summary.cl_alpha_wing_cruise_clean+aero_analysis.summary.cl_alpha_tail_cruise_clean*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_max_clean=aero_analysis.summary.cl_alpha_wing_max_clean+aero_analysis.summary.cl_alpha_tail_max_clean*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_TO=aero_analysis.summary.cl_alpha_wing_TO+aero_analysis.summary.cl_alpha_tail_TO*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_approach=aero_analysis.summary.cl_alpha_wing_approach+aero_analysis.summary.cl_alpha_tail_approach*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_loiter_clean=aero_analysis.summary.cl_alpha_wing_loiter_clean+aero_analysis.summary.cl_alpha_tail_loiter_clean*tailplane.horizontal.s/wing.Sref;

zero_AoA=aero_analysis.wing.zero_aoa;

for i=1:length(x)
    cl_alpha_wing_cruise_clean(i)=aero_analysis.summary.cl_alpha_wing_cruise_clean*x(i); %rad
    cl_alpha_wing_max_clean(i)=aero_analysis.summary.cl_alpha_wing_max_clean*x(i); %rad
    cl_alpha_wing_TO(i)=aero_analysis.summary.cl_alpha_wing_TO*x(i); %rad
    cl_alpha_wing_approach(i)=aero_analysis.summary.cl_alpha_wing_approach*x(i); %rad
    cl_alpha_wing_loiter_clean(i)=aero_analysis.summary.cl_alpha_wing_loiter_clean*x(i); %rad
end


%% L/D
%aero_analysis.summary.cl_cruise=
%aero_analysis.summary.cl_loiter=
%aero_analysis.summary.L_D_cruise=aero_analysis.summary.cl_cruise/aero_analysis.drag.cd_total(1);
%aero_analysis.summary.L_D_loiter=aero_analysis.summary.cl_loiter/aero_analysis.drag.cd_total(5);
aero_analysis.summary.e_wing=aero_analysis.induced_drag.e_theoretical(1);
%aero_analysis.summary.e_aircraft=aero_analysis.induced_drag.e_theoretical(1)+aero_analysis.induced_drag.e_theoretical(2)*tailplane.horizontal.s/wing.Sref;
%aero_analysis.summary.e_v2=aero_analysis.induced_drag.wing.e(1)+aero_analysis.induced_drag.tail.e(1)*tailplane.horizontal.s/wing.Sref;
%sizing.fraction.before_cruise
aero_analysis.summary.cl_cruise=(2*sizing.fraction.before_cruise*sizing.W0)/(aero_analysis.wing.rho(1)*(aero_analysis.wing.Mach(1)*aero_analysis.wing.air_velc(1))^2*aero_analysis.wing.HLD.s_ref);
aero_analysis.summary.cl_loiter=(2*sizing.fraction.before_loiter*sizing.W0)/(aero_analysis.wing.rho(2)*(aero_analysis.wing.Mach(2)*aero_analysis.wing.air_velc(2))^2*aero_analysis.wing.HLD.s_ref);
aero_analysis.summary.l_d_cruise=aero_analysis.summary.cl_cruise/aero_analysis.drag.cd_total(1);
aero_analysis.summary.l_d_loiter=aero_analysis.summary.cl_loiter/aero_analysis.drag.cd_total(2);
%% Drag polar

induced_AR=[wing.Ar, tailplane.horizontal.Ar];
cl_wing=[0:0.001:2];
e_polar=aero_analysis.induced_drag.e_theoretical;

cd_wing_polar=(cl_wing).^2/(pi*induced_AR(1)*e_polar(1));
cd_total_polar=aero_analysis.drag.cd0(1)+cd_wing_polar+aero_analysis.drag.wave(1);
l_d_max_wing_cruise=12.85.*cl_wing;
figure
plot(cd_total_polar,cl_wing)
hold on
plot(cl_wing,l_d_max_wing_cruise)
hold off
grid on
grid minor
xlim([0 0.15])

%l_d_cruise=aero_analysis.summary.L_D_cruise.*cl_wing_polar;
%l_d_cruise=13.2*cl_wing_polar;
% 
% figure
% plot(cd_total_polar,cl_total_polar)
% hold on
% plot(cl_wing_polar,l_d_cruise)
% %plot(aero_analysis.drag.cd_total(1),aero_analysis.summary.cl_cruise, 'o') 
% hold off
% grid on
% grid minor
% xlim([0 0.15])

%% Plots
figure
x=x*180/pi;
yline(aero_analysis.summary.cl_max_wing_clean, '--k')
hold on
yline(aero_analysis.summary.cl_max_approach, '--b')
yline(aero_analysis.summary.cl_max_TO, '--r')
%xline(0, '--m')
plot(x+zero_AoA,cl_alpha_wing_cruise_clean)
plot(x+zero_AoA, cl_alpha_wing_max_clean)
plot(x+zero_AoA+aero_analysis.wing.HLD.delta_alpha(1),cl_alpha_wing_TO, 'm')
plot(x+zero_AoA+aero_analysis.wing.HLD.delta_alpha(2),cl_alpha_wing_approach, 'k')
plot(x+zero_AoA,cl_alpha_wing_loiter_clean)
legend('Cl_max clean', 'Cl_max_approach', 'Cl_max_TO', 'Clean Cruise', 'Max Cruise', 'Takeoff', 'Approach', 'Loiter')

aero_analysis.summary.zero_aoa.TO=zero_AoA+aero_analysis.wing.HLD.delta_alpha(1);
aero_analysis.summary.zero_aoa.landing=zero_AoA+aero_analysis.wing.HLD.delta_alpha(2);

aero_analysis.summary.y_intercept_approach=0.3538;
aero_analysis.summary.y_intercept_TO=0.3488;

aero_analysis.summary.cl_transition=0.9*aero_analysis.summary.cl_max_TO;
ylim ([0 3])
xlim([-15 25])
title 'Wing'
xlabel '\alpha [degrees]'
hold off
grid on
grid minor



save('aero_analysis.mat', 'aero_analysis')
%Need to add annotations
%Fix legend
%Thicker lines
%Colours
%Terminate line after cross line

save('aero_analysis.summary','aero_analysis')

