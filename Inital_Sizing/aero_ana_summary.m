%% Aerodynamic Analysis
%  Summary script (with plots)
clear
clc

load 'aero_analysis.mat'
load 'tailplane.mat'
load 'wing.mat'

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

aero_analysis.summary.cl_alpha_loiter_clean=aero_analysis.summary.cl_alpha_wing_cruise_clean+aero_analysis.summary.cl_alpha_tail_cruise_clean*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_max_clean=aero_analysis.summary.cl_alpha_wing_max_clean+aero_analysis.summary.cl_alpha_tail_max_clean*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_TO=aero_analysis.summary.cl_alpha_wing_TO+aero_analysis.summary.cl_alpha_tail_TO*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_approach=aero_analysis.summary.cl_alpha_wing_approach+aero_analysis.summary.cl_alpha_tail_approach*tailplane.horizontal.s/wing.Sref;
aero_analysis.summary.cl_alpha_tail_loiter_clean=aero_analysis.summary.cl_alpha_wing_loiter_clean+aero_analysis.summary.cl_alpha_tail_loiter_clean*tailplane.horizontal.s/wing.Sref;

zero_AoA=aero_analysis.wing.zero_aoa;

for i=1:length(x)
    cl_alpha_wing_cruise_clean(i)=aero_analysis.summary.cl_alpha_wing_cruise_clean*x(i); %rad
    cl_alpha_wing_max_clean(i)=aero_analysis.summary.cl_alpha_wing_max_clean*x(i); %rad
    cl_alpha_wing_TO(i)=aero_analysis.summary.cl_alpha_wing_TO*x(i); %rad
    cl_alpha_wing_approach(i)=aero_analysis.summary.cl_alpha_wing_approach*x(i); %rad
    cl_alpha_wing_loiter_clean(i)=aero_analysis.summary.cl_alpha_wing_loiter_clean*x(i); %rad
end

figure
x=x*180/pi;
yline(aero_analysis.summary.cl_max_wing_clean, '--k')
hold on
yline(aero_analysis.summary.cl_max_approach, '--b')
yline(aero_analysis.summary.cl_max_TO, '--r')
plot(x+zero_AoA,cl_alpha_wing_cruise_clean)
plot(x+zero_AoA, cl_alpha_wing_max_clean)
plot(x+zero_AoA+aero_analysis.wing.HLD.delta_alpha(1),cl_alpha_wing_TO)
plot(x+zero_AoA+aero_analysis.wing.HLD.delta_alpha(2),cl_alpha_wing_approach)
plot(x+zero_AoA,cl_alpha_wing_loiter_clean)
legend('Cl_max clean', 'Cl_max_approach', 'Cl_max_TO', 'Clean Cruise', 'Max Cruise', 'Takeoff', 'Approach', 'Loiter')

ylim ([0 3])
xlim([-15 25])
title 'Wing'
xlabel '\alpha [degrees]'
hold off

save('aero_analysis.mat', 'aero_analysis')
%Need to add annotations
%Fix legend
%Thicker lines
%Colours
%Terminate line after cross line

save('aero_analysis.summary','aero_analysis')