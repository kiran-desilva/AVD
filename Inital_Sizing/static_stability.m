clear
clc

load('sizing.mat')
load('parameters.mat')
load('tailplane.mat')
load('wing.mat')
load('locations.mat')
load('fuse.mat')

%% initial variables



wing.Cmac = 1.4;

Sw = 19.51; %cessna mustag CHANGE
xacw = locations.x_ac_w;
Clw = .33;
Cmac = wing.Cmac;
wing_sweep_25 = wing.sweep_25;
wing_twist = 2;
wing_AR = wing.Ar;
wing_b = wing.b;
wing_Cl_alpha = ;%%%%%
wing_xac_bar = xacw/Cmac;
wing_Croot = wing.Croot;

eta_h = 1.0 ; % t-tail efficency -> errikos slides
Sh = tailplane.horizontal.s; % this prob should be have set in postions
Clh = -0.5;
zh = z_above_fuselage_bar * Cmac;

lh = tailplane.horizontal.l;

tail_h_Cl_alpha = ; %64012 cl alpha
tail_h_xac_bar = locations.x_ac_h/Cmac;

Lf = convlen(fuse.total_fuselage_length,'in','m'); %fueselage length
Wf = convlen(fuse.d_f,'in','m'); % fueslage max width ( diameter)

%%%%REGIONS%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kf_gen; % interpolate kf_gen data
wing_root_quater_chord_percent_length = ((0.25*wing_Croot) + locations.x_wing)/Lf;
Kf = Kf_fit(wing_root_quater_chord_percent_length);

fuselage_Cm_alpha = Kf*(Lf*(Wf^2))/(Cmac * Sw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcg = 1;
Zt = .6;

Cm0_airfoil = 1;


required_static_margin = .05*Cmac;
%temporaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms lh_syms zh_syms kh deda
%set to constant for now
zh_syms = zh;
% lh_syms = lh;

ka = (1/wing_AR) - (1/(1 + (wing_AR^1.7)));
klamb = (10-(3*wing.lambda))/7;
kh(lh_syms) = (1- abs(zh_syms/wing_b))/(((2*lh_syms)/(wing_b))^(1/3));

deda(lh_syms) = 4.44 * ((ka * klamb * kh * sqrt(wing_sweep_25))^1.19)*(Cl_alpha_m/Cl_alpha_m0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% syms xcg_xac_bar 
%tailplane stability analysis
% static_stab_constraint_sh_sw(xcg_xac_bar) = (wing_Cl_alpha * xcg_xac_bar)/(tail_h_Cl_alpha*eta_h*(1-deda)*((lh/Cmac)-xcg_xac_bar));
% control_stab_constraint_sh_sw(xcg_xac_bar) = 

syms xnp_bar wing_xac_bar_syms sh_sw_syms tail_h_xac_bar_syms lh_eq

lh_eq(tail_h_xac_bar_syms,wing_xac_bar_syms) = tail_h_xac_bar_syms - wing_xac_bar_syms;

xnp_bar(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms) = ((wing_Cl_alpha * wing_xac_bar_syms) - fuselage_Cm_alpha + (eta_h * sh_sw_syms * tail_h_Cl_alpha * (1-deda(lh_eq))*tail_h_xac_bar_syms))/(wing_Cl_alpha + (eta_h*tail_h_Cl_alpha*(1-deda(lh_eq))*sh_sw_syms));


%% TRIMMED CG ANALYSIS

syms aoa


Cm0w = ((Cm0_airfoil* (wing_AR * (cosd(wing_sweep_25)^2)) / (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 * wing_twist)) * (wing.Cl_alpha/wing.Cl_alpha_m0); 


