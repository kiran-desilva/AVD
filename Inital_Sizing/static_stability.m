clear
clc

load('sizing.mat')
load('parameters.mat')
load('tailplane.mat')
load('wing.mat')
load('locations.mat')
load('fuse.mat')
load('aero_analysis')

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kf_gen; % interpolate kf_gen data
wing_root_quater_chord_percent_length = ((0.25*wing_Croot) + locations.x_wing)/Lf;
Kf = Kf_fit(wing_root_quater_chord_percent_length);

fuselage_Cm_alpha = Kf*(Lf*(Wf^2))/(Cmac * Sw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xcg_wrt_ac = @(xacw,xach) x;

xcg = @(xacw,xach,wf_tail,wf_wing,Pax) 5.6;


Zt = .6;



required_static_margin = .05*Cmac;
%temporaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms lh_syms zh_syms kh deda cla_eta_syms
%set to constant for now
zh_syms = zh;
% lh_syms = lh;

ka = (1/wing_AR) - (1/(1 + (wing_AR^1.7)));
klamb = (10-(3*wing.lambda))/7;
kh(lh_syms) = (1- abs(zh_syms/wing_b))/(((2*lh_syms)/(wing_b))^(1/3));

deda(lh_syms,eta_mach_syms) = 4.44 * ((ka * klamb * kh * sqrt(wing_sweep_25))^1.19)*(cla_eta_syms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CM0W%%%%%%%
syms Cm0w 

Cm0_airfoil = -0.075;
Cm0w(cla_eta_syms) = ((Cm0_airfoil* (wing_AR * (cosd(wing_sweep_25)^2)) / (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 * wing.twsit)) * (cla_eta_syms); 

%%%%%%%%%%%%%%%%

%%%%REGIONS%%%%%%
region.regions = ["Takeoff","Landing","Cruise Start","Cruise End"];
region.cla_eta = [aero_analysis.wing.eta([3 4 1 1])];
region.deda = subs(deda,eta_mach_syms,region.cla_eta);
regions.fuel_fraction = [sizing.fraction.before_take_off,sizing.fraction.before_alternate_cruise,sizing.fraction.before_cruise,sizing.fraction.end_cruise_1];
region.xcg = []; % -> cams function
region.cl_alpha_w = [ aero_analysis.wing.HLD.cl_alpha_flaps([1 2]) aero_analysis.wing.Cl_alpha([1 1]) ];
region.Cm0w = double(subs(Cm0w,cla_eta_syms,region.cla_eta));
regions.C_thrust = [~,~,CD/Cmac,CD/Cmac]; % T/(q*sw*cbar) -> takeoff = 100% thrust, landing = 25% thrust
regions.Cl = []; % design lift coefficents at specific regions
regions.rho = [1.225,1.225,0.3016,0.3016];
regions.a = [340.3,340.3,295.0696,295.0696];
regions.m = [.1,.1,.75,.75];
regions.v = [aero_analysis.wing.v_takeoff_ms,aero_analysis.wing.v_landing_ms,221.3022,221.3022];
regions.q = 0.5 * (regions.rho .* (regions.v.^2));
regions.alpha_0_w = [];
regions.alpha_0_h = [];



%%%%%%%%%%%%%%%%%

% syms xcg_xac_bar 
%tailplane stability analysis
% static_stab_constraint_sh_sw(xcg_xac_bar) = (wing_Cl_alpha * xcg_xac_bar)/(tail_h_Cl_alpha*eta_h*(1-deda)*((lh/Cmac)-xcg_xac_bar));
% control_stab_constraint_sh_sw(xcg_xac_bar) = 


syms xnp_bar wing_xac_bar_syms sh_sw_syms tail_h_xac_bar_syms lh_eq sm wing_Cl_alpha_syms, wf_tail

lh_eq(tail_h_xac_bar_syms,wing_xac_bar_syms) = tail_h_xac_bar_syms - wing_xac_bar_syms;

xnp_bar(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms) = ((wing_Cl_alpha_syms * wing_xac_bar_syms) - fuselage_Cm_alpha + (eta_h * sh_sw_syms * tail_h_Cl_alpha * (1-deda(lh_eq))*tail_h_xac_bar_syms))/(wing_Cl_alpha_syms + (eta_h*tail_h_Cl_alpha*(1-deda(lh_eq))*sh_sw_syms));

sm(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms) = xnp_bar - (xcg(wing_xac_bar_syms,tail_h_xac_bar_syms,));

figure
hold on
%lb sm > 0%
fimplicit(@(x,y) sm(y,x,locations.x_ac_h,regions.cl_alpha_w(3)),'--') 
% ub sm < 20%
fimplicit(@(x,y) sm(y,x,locations.x_ac_h,regions.cl_alpha_w(3)) - .2,'--') 
% ideal lb sm > 4%
fimplicit(@(x,y) sm(y,x,locations.x_ac_h,regions.cl_alpha_w(3)) - .04,'--') 
%ideal ub sm < 7%
fimplicit(@(x,y) sm(y,x,locations.x_ac_h.regions.cl_alpha_w(3)) - .07,'--') 

legend('Sm > 0%','Sm < 20%','Sm > 4%','Sm < 7%')
ylabel('{S_h}/{S_w}')
xlabel('{x_ac_w}/{Cmac}')
grid on


%% TRIMMED CG ANALYSIS

syms CL_w CL_h C

CL_w = 


function [] = plot_sm_constraints()

