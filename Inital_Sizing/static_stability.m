clear
clc

load('sizing.mat')
load('parameters.mat')
load('tailplane.mat')

%% initial variables

wing.Cmac = 1.4;

Sw = 19.51; %cessna mustag CHANGE
xacw = 1;
Clw = .33;
Cmac = wing.Cmac;
wing_sweep_25 = 22;
wing_twist = 2;
wing_AR = 7.8;
wing_b = 13.16;
wing_Cl_alpha = ;
wing_xac_bar = ;
wing_Croot = 1.6;

eta_h = 1; % t-tail efficency
Sh = tailplane.initial.horizontal.s;
Clh = -0.5;
Zh = 1;
lh = tailplane.initial.horizontal.l;
tail_h_Cl_alpha = ; %64012 cl alpha
tail_h_xac_bar = (tailplane.initial.horizontal.Xac_from_tip + lh + xacw) / wing.Cmac;

Lf = 12; %fueselage length
Wf = 1.52; % fueslage max width ( diameter)

Kf_gen; % interpolate kf_gen data

wing_root_quater_chord = (0.25*wing_Croot)

fuselage_Cm_alpha = Kf*(Lf*(Wf^2))/(Cmac * Sw);

xcg = 1;

Zt = .6;

Cm0_airfoil = 1;


required_static_margin = .05*Cmac;
%temporaries
sh_sw = sh/sw;



ka = (1/wing_AR) - (1/(1 + (wing_AR^1.7)));
klamb = (10-(3*wing.lambda))/7;
kH = (1- abs(Zh/wing_b))/(((2*lh)/(wing_b))^(1/3));

deda = 4.44 * ((ka * klamb * kH * sqrt(phi_25_w))^1.19)*(Cl_alpha_m/Cl_alpha_m0);

syms xcg_xac_bar
%tailplane stability analysis
static_stab_constraint_sh_sw(xcg_xac_bar) = (wing_Cl_alpha * xcg_xac_bar)/(tail_h_Cl_alpha*eta_h*(1-deda)*((lh/Cmac)-xcg_xac_bar));
control_stab_constraint_sh_sw(xcg_xac_bar) = 


Cm0w = ((Cm0_airfoil* (wing_AR * (cosd(wing_sweep_25)^2)) / (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 * wing_twist)) * (wing.Cl_alpha/wing.Cl_alpha_m0); 




Xnp_bar = ((wing_Cl_alpha * wing_xac_bar) - fuselage_Cm_alpha + (eta_h * sh_sw * tail_h_Cl_alpha * (1-deda)*tail_h_xac_bar))/(wing_Cl_alpha + (eta_h*tail_h_Cl_alpha*(1-deda)*sh_sw));


%% TRIMMED CG ANALYSIS

syms aoa


