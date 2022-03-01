clear
clc
close all

%loading on horizontal tail as a factor of unit lift force

s_h = 2.3887; %m
U = 3.75 * 137.2501; %FAR 25 3.75 * VD
AR_h = 3.9;
rho = 1.225; 
Sref_h = 1.9283; 
Sref_w = 10.32;

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_acHtail = 10.3474; 
X_acWing = 4.3620; 
WH = 17.5 * 0.453592 * 9.81; %horizontal tail weight
ULF = 3.75; % ultimate load factor
W0 = 3128.2 * 9.81; %takeoff weight
Lwing = ULF * W0; 
Vdive = 137.25;
Cm0 = -0.0266;
M0w = 0.5 * rho * (Vdive*ULF)^2 * Sref_w * Cm0;

L_Htail = (Lwing * (X_cg - X_acWing) + M0w) / (X_acHtail - X_cg); 

% y = linspace(-s_h/2,s_h/2,100);
% gamma0 = (8 * s_h * 1) / (pi * AR_h * rho * U * Sref_h); %L = 1
% gamma = @(y) gamma0 * sqrt(1 - (y ./ (s_h/2)).^2); 
% TailLoadinit = rho * U .* gamma(y); 
% TailLoadintegral = trapz(y,TailLoadinit); 
% TailLoad = @(y) L_Htail * rho * U * gamma(y) / TailLoadintegral; %unit overall tail lift
% HorizontalTail.TailLoad = TailLoad(y);

y = linspace(-s_h/2,s_h/2,100);
gamma0 = (8 * s_h * L_Htail) / (pi * AR_h * rho * U * Sref_h); %L = 1
gamma = @(y) gamma0 * sqrt(1 - (y ./ (s_h/2)).^2); 
%TailLoadinit = rho * U .* gamma(y); 
%TailLoadintegral = trapz(y,TailLoadinit); 
TailLoad = @(y) L_Htail * rho * U * gamma(y) / TailLoadintegral; %unit overall tail lift
HorizontalTail.TailLoad = TailLoad(y);


figure
plot(y, TailLoad(y))
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Load (N)", 'interpreter', 'Latex')
grid on

R = WH - L_Htail; %vertical reaction from vertical stab

[shearforce] = ShearBending(y,s_h,R,TailLoad,WH);
HorizontalTail.Shearforce = shearforce;
bendingmoment = cumtrapz(y,shearforce);
HorizontalTail.bendingmoment = bendingmoment; 

figure
plot(y, shearforce)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Shear Force (N)", 'interpreter', 'Latex')
grid on

figure
plot(y, bendingmoment)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Bending Moments (Nm)", 'interpreter', 'Latex')
grid on

coordsouter = [0.000000  0.000000;
  0.005000  0.009780;
  0.007500  0.011790;
  0.012500  0.014900;
  0.025000  0.020350;
  0.050000  0.028100;
  0.075000  0.033940;
  0.100000  0.038710;
  0.150000  0.046200;
  0.200000  0.051730;
  0.250000  0.055760;
  0.300000  0.058440;
  0.350000  0.059780;
  0.400000  0.059810;
  0.450000  0.057980;
  0.500000  0.054800;
  0.550000  0.050560;
  0.600000  0.045480;
  0.650000  0.039740;
  0.700000  0.033500;
  0.750000  0.026950;
  0.800000  0.020290;
  0.850000  0.013820;
  0.900000  0.007860;
  0.950000  0.002880;
  1.000000  0.000000;
  0.000000  0.000000;
  0.005000 -0.009780;
  0.007500 -0.011790;
  0.012500 -0.014900;
  0.025000 -0.020350;
  0.050000 -0.028100;
  0.075000 -0.033940;
  0.100000 -0.038710;
  0.150000 -0.046200;
  0.200000 -0.051730;
  0.250000 -0.055760;
  0.300000 -0.058440;
  0.350000 -0.059780;
  0.400000 -0.059810;
  0.450000 -0.057980;
  0.500000 -0.054800;
  0.550000 -0.050560;
  0.600000 -0.045480;
  0.650000 -0.039740;
  0.700000 -0.033500;
  0.750000 -0.026950;
  0.800000 -0.020290;
  0.850000 -0.013820;
  0.900000 -0.007860;
  0.950000 -0.002880;
  1.000000  0.000000];

%plot(coords(:,1), coords(:,2))
%axis equal

section = polyshape(coordsouter(:,1),coordsouter(:,2));
figure
plot(section)
axis equal

[xbar,ybar] = centroid(section);

ctip = 0.4083;
croot = 0.8166;
ytip = 2.3887/2; 
c_H = @(y) ((ctip - croot)/ytip) * (abs(y) - ytip) + ctip; 

HorizontalTail.HTorsiondist = Torsion(y, c_H, TailLoad, xbar, WH);
figure
plot(y, HorizontalTail.HTorsiondist)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Torque distribution (Nm)", 'interpreter', 'Latex')
grid on

HorizontalTail.Torque = sum(HorizontalTail.HTorsiondist) - cumsum(HorizontalTail.HTorsiondist) + HorizontalTail.HTorsiondist; 

figure
plot(y, HorizontalTail.Torque)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Torque (Nm)", 'interpreter', 'Latex')
grid on

save('HorizontalTail.mat');

function [Torsion] = Torsion(y, c_H, TailLoad, xbar, WH)
    CoF = (0.68-0.15) .* c_H(y) / 2 + 0.15 .* c_H(y); %center of flexure midway between spars
    Torsion = TailLoad(y) .* (CoF - c_H(y)/4) + WH * (xbar - CoF); %no M0, ac at c/4 for each section
end

function [ShearforceH] = ShearBending(y,s_h,R,TailLoad,WH)
    for i = 1:length(y)
        Y = linspace(-s_h/2, y(i), 1000);
        %Moments = TailMoments(Y, (y(i) - y), TailLoad);
        if y(i) > 0 
             ShearforceH(i) = -(R - WH +  trapz(Y, TailLoad(Y))); 
             %BendingmomentH(i) = -((-R + WH)*y(i) +  trapz(Y, Moments(Y)));
        else
            ShearforceH(i) = -(trapz(Y, TailLoad(Y))); 
            %BendingmomentH(i) = -(trapz(Y, Moments(Y)));
        end
    end
     
end

