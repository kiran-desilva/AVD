clear 
clc

%Asymetric load case with one engine inoperative assuming engines at max
%thrust upon failure for worst case scenario

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_ac_v = 9.7224; %m vertical stab xac
Y_engine = 1.1; %m
T_1eng_max = 4950; %installed sea level thrust of engine (assuming Princess didnt fuck this)

s_v = 1.1818; %vertical tail span
AR_v = 1.1;

%Directional stability in OEI
F_v = (T_1eng_max * Y_engine) / (X_ac_v - X_cg); %force on vertical tail plane

F_vperm = F_v / s_v; %rectangular lift distribution

ShearForce_vt = @(z) F_vperm * (z - s_v); %shear force and bending moment for vertical tail
BendingMoment_vt = @(z) -F_vperm * (s_v - z).^2 / 2;

figure
z = linspace(0,s_v,100);
plot(z, ShearForce_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Shear Force (N)", 'interpreter', 'Latex')
grid on

figure
plot(z, BendingMoment_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Bending Moment (Nm)", 'interpreter', 'Latex')
grid on


