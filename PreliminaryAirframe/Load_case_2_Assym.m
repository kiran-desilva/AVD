clear 
clc
close all

%Asymetric load case with one engine inoperative assuming engines at max
%thrust upon failure for worst case scenario

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_ac_W = 4.3620;
X_ac_v = 9.7224; %m vertical stab xac
Y_engine = 1.1; %m
T_1eng_max = 4950; %installed sea level thrust of engine

s_v = 1.1818; %vertical tail span
AR_v = 1.1;
z = linspace(0,s_v,1000);

%Directional stability in OEI
F_v = (T_1eng_max * Y_engine) / (X_ac_v - X_ac_W); %force on vertical tail plane

F_vperm = F_v / s_v; %rectangular lift distribution

VerticalTail.ShearForce_vt = @(z) F_vperm * (z - s_v); %shear force and bending moment for vertical tail
VerticalTail.BendingMoment_vt = @(z) -F_vperm * (s_v - z).^2 / 2;

Verticaltail.Shearforce = VerticalTail.ShearForce_vt(z);

%Vertical tail torsion

ctip = 0.8166;
croot = 1.0015;
ztip = 2.3887/2; 
c_V = @(z) ((ctip - croot)/ztip) * (z - ztip) + ztip; 

VerticalTail.VTailLoad = F_v / s_v; %Force per m
VerticalTail.VTailLoad_1000 = F_v / 1000; %force at each of 1000 stations
VerticalTail.VTorsiondist = Torsion(z, c_V, VerticalTail.VTailLoad_1000); 

VerticalTail.Torque = sum(VerticalTail.VTorsiondist) - cumsum(VerticalTail.VTorsiondist) + VerticalTail.VTorsiondist; 

figure; 
subplot(5,1,1)
plot(z, VerticalTail.VTorsiondist)
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("OEI: Vertical Tail Combined Torque (Nm)", 'interpreter', 'Latex')
grid on
subplot(5,1,2)
plot(z, VerticalTail.Torque)
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("OEI: Vertical Tail Torque Distribution (Nm)", 'interpreter', 'Latex')
grid on
subplot(5,1,3)
tailloadplot(1:1000) = VerticalTail.VTailLoad;
plot(z, tailloadplot)
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("OEI: Vertical Tail Combined Load (N)", 'interpreter', 'Latex')
grid on
subplot(5,1,4)
plot(z, VerticalTail.ShearForce_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("OEI: Shear Force (N)", 'interpreter', 'Latex')
grid on
subplot(5,1,5)
plot(z, VerticalTail.BendingMoment_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("OEI: Bending Moment (Nm)", 'interpreter', 'Latex')
grid on





save("VerticalTail.mat",'VerticalTail');

function [Torsion] = Torsion(z, c_V, VTailLoad)
    CoF = (0.7-0.15) .* c_V(z) / 2 + 0.15 .* c_V(z); %center of flexure midway between spars
    Torsion = VTailLoad .* (CoF - c_V(z)/4); %no M0, ac at c/4 for each section
end
