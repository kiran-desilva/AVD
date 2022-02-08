clear
clc

%loading on horizontal tail for unit lift force

s_h = 2.3887; %m
U = 100;
L = 1;
AR_h = 3.9;
rho = 1.225; 
Sref_h = 1.9283; 

y = linspace(-s_h/2,s_h/2,100);
gamma0 = (8 * s_h * 1) / (pi * AR_h * rho * U * Sref_h); %L = 1
gamma = @(y) gamma0 * sqrt(1 - (y ./ (s_h/2)).^2); 
TailLoadinit = rho * U .* gamma(y); 
TailLoadintegral = trapz(y,TailLoadinit); 
TailLoad = @(y) L * rho * U * gamma(y) / TailLoadintegral; %unit overall tail lift 

figure
plot(y, TailLoad(y))
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Load (N)", 'interpreter', 'Latex')
grid on


