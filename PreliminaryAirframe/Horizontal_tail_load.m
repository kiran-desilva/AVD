clear
clc

%loading on horizontal tail as a factor of unit lift force

s_h = 2.3887; %m
U = 100;
AR_h = 3.9;
rho = 1.225; 
Sref_h = 1.9283; 
Sref_w = 10.32;

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_acHtail = 10.3474; 
X_acWing = 4.3620; 
ULF = 3.75; % ultimate load factor
W0 = 3128.2 * 9.81; %takeoff weight
Lwing = ULF * W0;
Vdive = 137.25;
Cm0 = -0.0266;
M0w = 0.5 * rho * Vdive^2 * Sref_w * Cm0;

L_Htail = (Lwing * (X_cg - X_acWing) + M0w) / (X_acHtail - X_cg); 

y = linspace(-s_h/2,s_h/2,100);
gamma0 = (8 * s_h * 1) / (pi * AR_h * rho * U * Sref_h); %L = 1
gamma = @(y) gamma0 * sqrt(1 - (y ./ (s_h/2)).^2); 
TailLoadinit = rho * U .* gamma(y); 
TailLoadintegral = trapz(y,TailLoadinit); 
TailLoad = @(y) L_Htail * rho * U * gamma(y) / TailLoadintegral; %unit overall tail lift 

figure
plot(y, TailLoad(y))
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Load (N)", 'interpreter', 'Latex')
grid on


