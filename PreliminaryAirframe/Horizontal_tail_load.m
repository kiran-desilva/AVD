clear
clc

%loading on horizontal tail for unit lift force

s_h = 2.3887; %m
U = ??
L = 1;
AR_h = 3.9;
rho
Sref_h = 1.9283; 

gamma0 = (8 * s_h * L) / (pi * AR_h * rho * U * Sref_h); 
gamma = @(y) gamma0 * sqrt(1 - (y / s)^2)

y = linspace(0,s_h,100);
