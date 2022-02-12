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

%Directional stability in OEI
F_v = (T_1eng_max * Y_engine) / (X_ac_v - X_ac_W); %force on vertical tail plane

F_vperm = F_v / s_v; %rectangular lift distribution

ShearForce_vt = @(z) F_vperm * (z - s_v); %shear force and bending moment for vertical tail
BendingMoment_vt = @(z) -F_vperm * (s_v - z).^2 / 2;

figure
z = linspace(0,s_v,1000);
plot(z, ShearForce_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Shear Force (N)", 'interpreter', 'Latex')
grid on

figure
plot(z, BendingMoment_vt(z))
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Bending Moment (Nm)", 'interpreter', 'Latex')
grid on

%Vertical tail torsion

coordsouter = [ 0.000000  0.000000;
  0.005000  0.012080;
  0.007500  0.014560;
  0.012500  0.018420;
  0.025000  0.025280;
  0.050000  0.035040;
  0.075000  0.042400;
  0.100000  0.048420;
  0.150000  0.057850;
  0.200000  0.064800;
  0.250000  0.069850;
  0.300000  0.073190;
  0.350000  0.074820;
  0.400000  0.074730;
  0.450000  0.072240;
  0.500000  0.068100;
  0.550000  0.062660;
  0.600000  0.056200;
  0.650000  0.048950;
  0.700000  0.041130;
  0.750000  0.032960;
  0.800000  0.024720;
  0.850000  0.016770;
  0.900000  0.009500;
  0.950000  0.003460;
  1.000000  0.000000;
  0.000000  0.000000;
  0.005000  -.012080;
  0.007500  -.014560;
  0.012500  -.018420;
  0.025000  -.025280;
  0.050000  -.035040;
  0.075000  -.042400;
  0.100000  -.048420;
  0.150000  -.057850;
  0.200000  -.064800;
  0.250000  -.069850;
  0.300000  -.073190;
  0.350000  -.074820;
  0.400000  -.074730;
  0.450000  -.072240;
  0.500000  -.068100;
  0.550000  -.062660;
  0.600000  -.056200;
  0.650000  -.048950;
  0.700000  -.041130;
  0.750000  -.032960;
  0.800000  -.024720;
  0.850000  -.016770;
  0.900000  -.009500;
  0.950000  -.003460;
  1.000000  0.000000];

%plot(coords(:,1), coords(:,2))
%axis equal

section = polyshape(coordsouter(:,1),coordsouter(:,2));
plot(section)
axis equal

[xbar,ybar] = centroid(section);

ctip = 0.8166;
croot = 1.0015;
ztip = 2.3887/2; 
c_V = @(z) ((ctip - croot)/ztip) * (z - ztip) + ztip; 

VTailLoad = F_v / 1000; %Force at each of 1000 stations
VTorsiondist = Torsion(z, c_V, VTailLoad, xbar); 
figure
plot(z, VTorsiondist)
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Vertical Tail Torque distribution (Nm)", 'interpreter', 'Latex')
grid on

Torque = sum(VTorsiondist) - cumsum(VTorsiondist) + VTorsiondist; 

figure
plot(z, Torque)
xlabel("z (m)", 'interpreter', 'Latex')
ylabel("Vertical Tail Torque (Nm)", 'interpreter', 'Latex')
grid on

function [Torsion] = Torsion(z, c_V, VTailLoad, xbar)
    CoF = (0.7-0.15) .* c_V(z) / 2 + 0.15 .* c_V(z); %center of flexure midway between spars
    Torsion = VTailLoad .* (CoF - c_V(z)/4); %no M0, ac at c/4 for each section
end
