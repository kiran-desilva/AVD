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
alpha_0_w = degtorad(-3);

eta_h = 1.0 ; % t-tail efficency -> errikos slides
Sh = tailplane.horizontal.s; % this prob should be have set in postions
Clh = -0.5;
zh = z_above_fuselage_bar * Cmac;

lh = tailplane.horizontal.l;

tail_h_Cl_alpha = aero_analysis.tail.Cl_alpha; %64012 cl alpha
tail_h_xac_bar = locations.x_ac_h/Cmac;
alpha_0_h = degtorad(0);

Lf = convlen(fuse.total_fuselage_length,'in','m'); %fueselage length
Wf = convlen(fuse.d_f,'in','m'); % fueslage max width ( diameter)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms wing_xac_bar_syms kfx

Kf_gen; % interpolate kf_gen data

wing_root_quater_chord_percent_length = @(xac) (((0.25*wing_Croot) - wing.Xac_from_tip) + xac)/Lf;
Kf_fit = poly2sym(Kf_fit_coef,kfx);
Kf(wing_xac_bar_syms) =  Kf_fit(wing_root_quater_chord_percent_length(wing_xac_bar_syms.*Cmac));

fuselage_Cm_alpha(wing_xac_bar_syms) = Kf*(Lf*(Wf^2))/(Cmac * Sw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xcg_wrt_ac = @(xacw,xach) x;

xcg = @(xacw,xach,wf_wing,wf_tail,payload_factor) 5.6;


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

deda(lh_syms,cla_eta_syms) = 4.44 * ((ka * klamb * kh * sqrt(wing_sweep_25))^1.19)*(cla_eta_syms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CM0W%%%%%%%
syms Cm0w 

Cm0_airfoil = -0.075;
Cm0w(cla_eta_syms) = ((Cm0_airfoil* (wing_AR * (cosd(wing_sweep_25)^2)) / (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 * wing.twsit)) * (cla_eta_syms); 

%%%%%%%%%%%%%%%%

%%%%REGIONS%%%%%%
region.regions = ["Takeoff","Landing","Cruise Start","Cruise End"];
region.cla_eta = [aero_analysis.wing.eta([3 4 1 1])];
% region.deda = subs(deda,cla_eta_syms,region.cla_eta);

region.fuel_fraction = [sizing.fraction.before_take_off,sizing.fraction.before_alternate_cruise,sizing.fraction.before_cruise,sizing.fraction.end_cruise_1];
region.xcg = []; % -> cams function
region.cl_alpha_w = [ aero_analysis.wing.HLD.cl_alpha_flaps([1 2]) aero_analysis.wing.Cl_alpha([1 1]) ];
region.Cm0w = double(subs(Cm0w,cla_eta_syms,region.cla_eta));
region.C_thrust = [~,~,CD/Cmac,CD/Cmac]; % T/(q*sw*cbar) -> takeoff = 100% thrust, landing = 25% thrust
region.CL = []; % design lift coefficents at specific regions
region.rho = [1.225,1.225,0.3016,0.3016];
region.a = [340.3,340.3,295.0696,295.0696];
region.m = [.1,.1,.75,.75];
region.v = [aero_analysis.wing.v_takeoff_ms,aero_analysis.wing.v_landing_ms,221.3022,221.3022];
region.q = 0.5 * (regions.rho .* (regions.v.^2));
% region.alpha_0_w = degtorad([-3 -3 -3 -3]);
% region.alpha_0_h = degtorad([0 0 0 0]);



%%%%%%%%%%%%%%%%%

% syms xcg_xac_bar 
%tailplane stability analysis
% static_stab_constraint_sh_sw(xcg_xac_bar) = (wing_Cl_alpha * xcg_xac_bar)/(tail_h_Cl_alpha*eta_h*(1-deda)*((lh/Cmac)-xcg_xac_bar));
% control_stab_constraint_sh_sw(xcg_xac_bar) = 


syms xnp_bar wing_xac_bar_syms sh_sw_syms tail_h_xac_bar_syms lh_eq sm wing_Cl_alpha_syms wf_tail_syms wf_wing_syms payload_factor_syms

xcg_bar_eq(wing_xac_bar_syms,tail_h_xac_bar_syms,wf_wing_syms,wf_tail_syms,payload_factor_syms) = (xcg(wing_xac_bar_syms,tail_h_xac_bar_syms,wf_wing_syms,wf_tail_syms,payload_factor_syms)./Cmac);

lh_eq(tail_h_xac_bar_syms,wing_xac_bar_syms) = tail_h_xac_bar_syms - wing_xac_bar_syms;

xnp_bar(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms) = ((wing_Cl_alpha_syms * wing_xac_bar_syms) - fuselage_Cm_alpha(wing_xac_bar_syms) + (eta_h * sh_sw_syms * tail_h_Cl_alpha * (1-deda(lh_eq))*tail_h_xac_bar_syms))/(wing_Cl_alpha_syms + (eta_h*tail_h_Cl_alpha*(1-deda(lh_eq))*sh_sw_syms));
%payload conf
sm(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms,wf_wing_syms,wf_tail_syms,payload_factor_syms) = xnp_bar - xcg_bar_eq;

figure
hold on

% takeoff max fuel and passengers
%lb sm > 0%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(1),1,1,1),'--','color','red') 
% ub sm < 20%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(1),1,1,1) - .2,'--','color','red') 
% ideal lb sm > 4%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(1),1,1,1) - .04,'--','color','red') 
%ideal ub sm < 7%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(1),1,1,1) - .07,'--','color','red') 

%landing empty fuel and passengers?
%lb sm > 0%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(2),0,0,1),'--','color','red') 
% ub sm < 20%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(2),0,0,1) - .2,'--','color','red') 
% ideal lb sm > 4%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(2),0,0,1) - .04,'--','color','red') 
%ideal ub sm < 7%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(2),0,0,1) - .07,'--','color','red') 

%cruise moderate fuel and full passengers
%lb sm > 0%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(3),0,0,1),'--','color','red') 
% ub sm < 20%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(3),0,0,1) - .2,'--','color','red') 
% ideal lb sm > 4%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(3),0,0,1) - .04,'--','color','red') 
%ideal ub sm < 7%
fimplicit(@(x,y) sm(y,x,tail_h_xac_bar,region.cl_alpha_w(3),0,0,1) - .07,'--','color','red') 


legend('Sm > 0%','Sm < 20%','Sm > 4%','Sm < 7%')
ylabel('{S_h}/{S_w}')
xlabel('{x_ac_w}/{Cmac}')
grid on


%% TRIMMED CG ANALYSIS

syms CL_w CL_h alpha iw_syms ih_syms C_thrust_syms



Cmalphaf = fuselage_Cm_alpha(locations.x_ac_w);
sh_sw = Sh/Sw;

CL_w(alpha,iw_syms,wing_Cl_alpha_syms) = wing_Cl_alpha_syms*(alpha + iw_syms - alpha_0_w);
CL_h(alpha,ih_syms,iw_syms,cla_eta_syms) = tail_h_Cl_alpha*( ((alpha + iw_syms - alpha_0_w)*(1-deda(lh,cla_eta_syms))) + (ih_syms-iw_syms) - (alpha_0_h - alpha_0_w));


Cmcg(alpha,ih_syms,iw_syms) = -CL_w*(wing_xac_bar_syms - xcg_bar_eq) + Cm0w + (Cmalphaf*alpha) - (eta_h*CL_h*sh_sw*(tail_h_xac_bar - xcg_bar_eq)) + (Zt*C_thrust_syms);
CL(alpha,ih_syms,iw_syms) = CL_w + (eta_h*sh_sw*CL_h);


for idx = 1:length(region.regions)

   Cmcg_function = subs(Cmcg,[cla_eta_syms,wing_Cl_alpha_syms,wing_xac_bar_syms,tail_h_xac_bar_syms],[region.cla_eta(idx),region.cl_alpha_w(idx),wing_xac_bar,tail_h_xac_bar])
   % sub in weights
   Cmcg_function = subs(Cmcg_function,[wf_wing_syms,wf_tail_syms,payload_factor_syms],[])


    %x(1) = alpha, x(2) = ih
    F(1) = @(x) double(subs(Cmcg,[wing_C]))