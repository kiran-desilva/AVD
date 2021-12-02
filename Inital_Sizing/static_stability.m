clear
clc

load('sizing.mat')
load('parameters.mat')
load('tailplane.mat')
load('wing.mat')
load('locations.mat')
load('fuse.mat')
load('aero_analysis')
load('design')
load('wing')
load('wandb')

%% initial variables





Sw = wing.Sref; %cessna mustag CHANGE
xacw = locations.x_ac_w;
Cmac = wing.Cmac;
wing_sweep_25 = wing.sweep_25;
wing_twist = wing.twist;
wing_AR = wing.Ar;
wing_b = wing.b;
wing_xac_bar = xacw/Cmac;
wing_Croot = wing.Croot;
alpha_0_w = degtorad(-3);

eta_h = 1.0 ; % t-tail efficency -> errikos slides
Sh = tailplane.horizontal.s; % this prob should be have set in postions
zh = tailplane.horizontal.z_above_fuselage_bar .* Cmac;

lh = tailplane.horizontal.l;

tail_h_Cl_alpha = aero_analysis.tail.Cl_alpha(1); 
tail_h_xac_bar = locations.x_ac_h/Cmac;
alpha_0_h = degtorad(0);

Lf = convlength(fuse.total_fuselage_length,'in','m'); %fueselage length
Wf = convlength(fuse.d_f,'in','m'); % fueslage max width ( diameter)

% +ve if thrustline under cg
Zt = -0.716;

required_static_margin = .05*Cmac;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms wing_xac_bar_syms kfx

Kf_gen; % interpolate kf_gen data

wing_root_quater_chord_percent_length = @(xac) (((0.25*wing_Croot) - wing.Xac_from_tip) + xac)./Lf;
Kf_fit = poly2sym(Kf_fit_coef,kfx);
Kf(wing_xac_bar_syms) =  subs(Kf_fit,kfx,wing_root_quater_chord_percent_length(wing_xac_bar_syms.*Cmac));

fuselage_Cm_alpha(wing_xac_bar_syms) = Kf.*(Lf.*(Wf.^2))./(Cmac .* Sw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xacw - m
% wf_fuel - fuel fraction
% payload factor - 0 - 1
xcg = @(xacw,wf_fuel,payload_factor) wandb.x_cg_function(xacw,wf_fuel,payload_factor)*0.3048;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms lh_syms zh_syms kh deda cla_eta_syms
%set to constant for now
zh_syms = zh;
% lh_syms = lh;

ka = (1./wing_AR) - (1./(1 + (wing_AR.^1.7)));
klamb = (10-(3.*wing.lambda))./7;
kh(lh_syms) = (1- abs(zh_syms./wing_b))./(((2.*lh_syms)./(wing_b)).^(1/3));

deda(lh_syms,cla_eta_syms) = 4.44 .* ((ka .* klamb .* kh .* sqrt(cosd(wing_sweep_25))).^1.19).*(cla_eta_syms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CM0W%%%%%%%
syms Cm0w 

Cm0_airfoil = -0.075;
Cm0w(cla_eta_syms) = ((Cm0_airfoil.* (wing_AR .* (cosd(wing_sweep_25)^2)) ./ (wing_AR+(2*cosd(wing_sweep_25)))) - (0.01 .* wing.twist)) .* (cla_eta_syms); 

%%%%%%%%%%%%%%%%

%%%%REGIONS%%%%%%
region.regions = ["Takeoff","Landing","Cruise Start","Cruise End"];
region.cla_eta = [aero_analysis.wing.eta([3 4 1 1])];
% region.deda = subs(deda,cla_eta_syms,region.cla_eta);

region.fuel_fraction = [sizing.fraction.before_take_off,sizing.fraction.before_alternate_cruise,sizing.fraction.before_cruise,sizing.fraction.end_cruise_1];

% region.xcg = []; % -> cams function
region.cl_alpha_w = [ aero_analysis.wing.HLD.cl_alpha_flaps([1 2]) aero_analysis.wing.Cl_alpha([1 1]) ];
region.cl_alpha_h = [aero_analysis.tail.Cl_alpha([3 4 1 1])];
region.Cm0w = double(subs(Cm0w,cla_eta_syms,region.cla_eta));


region.rho = [1.225,1.225,0.3016,0.3016];
% region.a = [340.3,340.3,295.0696,295.0696];
% region.m = [.1,.1,.75,.75];
region.v = [aero_analysis.wing.v_takeoff_ms,aero_analysis.wing.v_landing_ms,221.3022,221.3022];
region.q = 0.5 * (region.rho .* (region.v.^2));

region.CL = (region.fuel_fraction*sizing.W0)./(region.q * Sw); % design lift coefficents at specific regions

CD_cruise = aero_analysis.drag.cd_total(1);

region.C_thrust = [design.t0/(region.q(1)*Sw*Cmac),(design.t0/(region.q(1)*Sw*Cmac)) * 0.5,CD_cruise/Cmac,CD_cruise/Cmac]; % T/(q*sw*cbar) -> takeoff = 100% thrust, landing = 50% thrust
% region.alpha_0_w = degtorad([-3 -3 -3 -3]);
% region.alpha_0_h = degtorad([0 0 0 0]);



%%%%%%%%%%%%%%%%%

% syms xcg_xac_bar 
%tailplane stability analysis
% static_stab_constraint_sh_sw(xcg_xac_bar) = (wing_Cl_alpha * xcg_xac_bar)/(tail_h_Cl_alpha*eta_h*(1-deda)*((lh/Cmac)-xcg_xac_bar));
% control_stab_constraint_sh_sw(xcg_xac_bar) = 


syms xnp_bar wing_xac_bar_syms sh_sw_syms tail_h_xac_bar_syms lh_eq sm wing_Cl_alpha_syms tail_Cl_alpha_syms wf_fuel_syms payload_factor_syms

xcg_bar_eq(wing_xac_bar_syms,wf_fuel_syms,payload_factor_syms) = (sym(xcg(wing_xac_bar_syms.*Cmac,wf_fuel_syms,payload_factor_syms)./Cmac));

lh_eq(tail_h_xac_bar_syms,wing_xac_bar_syms) = tail_h_xac_bar_syms - wing_xac_bar_syms;

xnp_bar(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms) = ((wing_Cl_alpha_syms .* wing_xac_bar_syms) - fuselage_Cm_alpha(wing_xac_bar_syms) + (eta_h .* sh_sw_syms .* tail_Cl_alpha_syms .* (1-deda(lh_eq,cla_eta_syms)).*tail_h_xac_bar_syms)) ./ (wing_Cl_alpha_syms + (eta_h.*tail_Cl_alpha_syms.*(1-deda(lh_eq,cla_eta_syms)).*sh_sw_syms));
%payload conf
sm(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms) = xnp_bar(sh_sw_syms,wing_xac_bar_syms,tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms) - xcg_bar_eq(wing_xac_bar_syms,wf_fuel_syms,payload_factor_syms);

disp('static margin')

for i = 1:length(region.regions)
    region.Cm0w(i) = double(Cm0w(region.cla_eta(i)))
    region.de_da(i) = double(deda(lh,region.cla_eta(i)));
    region.Kn_off(i) = double(sm(Sh/Sw,locations.x_ac_w/Cmac,locations.x_ac_h/Cmac,region.cl_alpha_w(i),region.cl_alpha_h(i),region.cla_eta(i),region.fuel_fraction(i),1));
    region.Kn_on(i) = region.Kn_off(i) - 0.02;
    region.xnp_bar(i) = double(xnp_bar(Sh/Sw,locations.x_ac_w/Cmac,locations.x_ac_h/Cmac,region.cl_alpha_w(i),region.cl_alpha_h(i),region.cla_eta(i)));
end

region

figure
hold on

xmin = 3;
xmax = 4.5;
ymin = 0;
ymax = 0.2;


% takeoff max fuel and passengers
f1 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(1),region.cl_alpha_h(1),region.cla_eta(1),sizing.fraction.before_take_off,1]));
plot_sm_constraints(f1,'red')

%landing empty fuel and passengers
f2 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(2),region.cl_alpha_h(2),region.cla_eta(2),sizing.fraction.end,1]));
plot_sm_constraints(f2,'blue')

%cruise start fuel and full passengers
f3 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(3),region.cl_alpha_h(3),region.cla_eta(3),sizing.fraction.before_cruise,1]));
plot_sm_constraints(f3,'green')

%cruise end fuel and full passengers
f4 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(4),region.cl_alpha_h(4),region.cla_eta(4),sizing.fraction.end_cruise_1,1]));
plot_sm_constraints(f4,'magenta')

yline(Sh/Sw,'--','linewidth',1.5)
xline(xacw/Cmac,'--','linewidth',1.5)

% legend('Sm > 0%','Sm < 20%','Sm > 4%','Sm < 7%')
legend('Takeoff - Passengers','','','','Landing - Passengers','','','','Cruise Start - Passengers','','','','Cruise End - Passengers','','','')

ylabel('{S_h}/{S_w}')
xlabel('{x_ac_w}/{Cmac}')
grid on
grid minor
xlim([xmin,xmax])
ylim([ymin,ymax])

improvePlot(gcf)

figure
hold on

% takeoff max fuel and passengers
f1 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(1),region.cl_alpha_h(1),region.cla_eta(1),sizing.fraction.before_take_off,0]));
plot_sm_constraints(f1,'red')

%landing empty fuel and passengers
f2 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(2),region.cl_alpha_h(2),region.cla_eta(2),sizing.fraction.end,0]));
plot_sm_constraints(f2,'blue')

%cruise start fuel and full passengers
f3 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(3),region.cl_alpha_h(3),region.cla_eta(3),sizing.fraction.before_cruise,0]));
plot_sm_constraints(f3,'green')

%cruise end fuel and full passengers
f4 = matlabFunction(subs(sm,[tail_h_xac_bar_syms,wing_Cl_alpha_syms,tail_Cl_alpha_syms,cla_eta_syms,wf_fuel_syms,payload_factor_syms],[tail_h_xac_bar,region.cl_alpha_w(4),region.cl_alpha_h(4),region.cla_eta(4),sizing.fraction.end_cruise_1,0]));
plot_sm_constraints(f4,'magenta')


yline(Sh/Sw,'--','linewidth',1.5)
xline(xacw/Cmac,'--','linewidth',1.5)

% legend('Sm > 0%','Sm < 20%','Sm > 4%','Sm < 7%')
legend('Takeoff - No Passengers','','','','Landing - No Passengers','','','','Cruise Start - No Passengers','','','','Cruise End - No Passengers','','','')

ylabel('{S_h}/{S_w}')
xlabel('{x_ac_w}/{Cmac}')
grid on
grid minor
xlim([xmin,xmax])
ylim([ymin,ymax])

improvePlot(gcf)



%% TRIMMED CG ANALYSIS

syms CL_w CL_h alpha iw_syms ih_syms C_thrust_syms



% Cmalphaf = fuselage_Cm_alpha(wing_xac_bar);
sh_sw = Sh/Sw;

CL_w(alpha,iw_syms,wing_Cl_alpha_syms) = wing_Cl_alpha_syms.*(alpha + iw_syms - alpha_0_w);
CL_h(alpha,ih_syms,iw_syms,cla_eta_syms) = tail_h_Cl_alpha.*( ((alpha + iw_syms - alpha_0_w).*(1-deda(lh,cla_eta_syms))) + (ih_syms-iw_syms) - (alpha_0_h - alpha_0_w));


Cmcg(alpha,ih_syms,iw_syms) = -CL_w(alpha,iw_syms,wing_Cl_alpha_syms).*(wing_xac_bar_syms - xcg_bar_eq(wing_xac_bar_syms,wf_fuel_syms,payload_factor_syms)) + Cm0w(cla_eta_syms) + (fuselage_Cm_alpha(wing_xac_bar_syms).*alpha) - (eta_h.*CL_h(alpha,ih_syms,iw_syms,cla_eta_syms).*sh_sw.*(tail_h_xac_bar - xcg_bar_eq(wing_xac_bar_syms,wf_fuel_syms,payload_factor_syms))) + (Zt.*C_thrust_syms);
CL(alpha,ih_syms,iw_syms) = CL_w(alpha,iw_syms,wing_Cl_alpha_syms) + (eta_h.*sh_sw.*CL_h(alpha,ih_syms,iw_syms,cla_eta_syms));

load('optimized_incidence')
wing.i_w = iw_required;
for idx = 1:length(region.regions)

    %sub in aerodymaic constants
    Cmcg_function = subs(Cmcg,[cla_eta_syms,wing_Cl_alpha_syms,wing_xac_bar_syms,tail_h_xac_bar_syms],[region.cla_eta(idx),region.cl_alpha_w(idx),wing_xac_bar,tail_h_xac_bar]);
    % sub in for weight configuration -> assuming always fully loaded with payload
    Cmcg_function = subs(Cmcg_function,[wf_fuel_syms,payload_factor_syms],[region.fuel_fraction(idx),1]);
    % sub in thrust settings
    Cmcg_function = subs(Cmcg_function,[C_thrust_syms],[region.C_thrust(idx)]);
    % expose the variables we want
    Cmcg_function(alpha,ih_syms,iw_syms) = matlabFunction(Cmcg_function,'vars',[alpha,ih_syms,iw_syms]);

    CL_function = subs(CL,[cla_eta_syms,wing_Cl_alpha_syms],[region.cla_eta(idx),region.cl_alpha_w(idx)]);
    CL_function(alpha,ih_syms,iw_syms) = matlabFunction(CL_function,'vars',[alpha,ih_syms,iw_syms]);

    %x(1) = alpha, x(2) = ih
    F = @(x) [double(Cmcg_function(x(1),x(2),wing.i_w * (pi/180))); double(CL_function(x(1),x(2),wing.i_w * (pi/180))) - region.CL(idx)];
 

    x0 = [0,0];
    res = fsolve(F,x0);

    region.aoa(idx) = res(1) * (180/pi);
    region.ih(idx) = res(2) * (180/pi);

end

region

%sub in aerodymaic constants
Cmcg_function = subs(Cmcg,[cla_eta_syms,wing_Cl_alpha_syms,wing_xac_bar_syms,tail_h_xac_bar_syms],[region.cla_eta(3),region.cl_alpha_w(3),wing_xac_bar,tail_h_xac_bar]);
% sub in for weight configuration -> assuming always fully loaded with payload
Cmcg_function = subs(Cmcg_function,[wf_fuel_syms,payload_factor_syms],[region.fuel_fraction(idx),1]);
% sub in thrust settings
Cmcg_function = subs(Cmcg_function,[C_thrust_syms],[region.C_thrust(3)]);
% expose the variables we want
Cmcg_function(alpha,ih_syms,iw_syms) = matlabFunction(Cmcg_function,'vars',[alpha,ih_syms,iw_syms]);

CL_function = subs(CL,[cla_eta_syms,wing_Cl_alpha_syms],[region.cla_eta(3),region.cl_alpha_w(3)]);
CL_function(alpha,ih_syms,iw_syms) = matlabFunction(CL_function,'vars',[alpha,ih_syms,iw_syms]);


F = @(x) [double(Cmcg_function(0,x(1),x(2))); double(CL_function(0,x(1),x(2))) - region.CL(3)];
x0 = [0,0];
res = fsolve(F,x0);
ih_required = res(1) * (180/pi)
iw_required = res(2) * (180/pi)

save('optimized_incidence','ih_required','iw_required')



function [] = plot_sm_constraints(func,color)
    %lb sm > 0%
    fimplicit(@(x,y) func(y,x),'-','color',color,'linewidth',1.2) 
    % ub sm < 20%
    fimplicit(@(x,y) func(y,x) - .2,'-','color',color,'linewidth',1.2) 
    % ideal lb sm > 4%
    fimplicit(@(x,y) func(y,x) - .04,'--','color',color) 
    %ideal ub sm < 7%
    fimplicit(@(x,y) func(y,x) - .07,'--','color',color) 
end