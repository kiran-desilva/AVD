clear
clc

load('sizing.mat')
load('parameters.mat')
load('tailplane.mat')

%% initial variables

eta_h = 1; % t-tail efficency
Sh = tailplane.initial.horizontal.s;
Clh = -0.5;
Zh = 1;
lh = tailplane.initial.horizontal.l;

Sw = 19.51; %cessna mustag CHANGE
xacw = 1;
Clw = .33;
cmac = 1.48;
wing_sweep_25 = 22;
wing_twist = 2;
wing_AR = 7.8;
wing_b = 13.16;


xcg = 1;

Zt = .6;



Cm0_airfoil = 1;

Cm0w = ((Cm0_airfoil* (wing_AR * (cosd(wing_sweep_25)^2)) / (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 * wing_twist)) * (wing.Cl_alpha/wing.Cl_alpha_m0); 


ka = (1/wing_AR) - (1/(1 + (wing_AR^1.7)));
klamb = (10-(3*wing.lambda))/7;
kH = (1- abs(Zh/wing_b))/(((2*lh)/(wing_b))^(1/3));

deda = 4.44 * ((ka * klamb * kH * sqrt(phi_25_w))^1.19)*(Cl_alpha_m/Cl_alpha_m0);