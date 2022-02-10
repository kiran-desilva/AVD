clear
clc
load('voldist.mat')

%Constants
lbs_to_n = 4.4482216153;
ft_to_m = 0.3048;
ulf = 1.5*2.5;

%%Weights and CG from Report UNITS IN LBS AND FT

i = 0;

% i = i + 1;
% component_weight(i) = 189.6; %wing
% component_cg(i,1) = 15.32;
% component_cg(i,2) = -2.75;

i = i + 1;
component.weight(i) = 17.1; % vertical tail
component.cg(i,1) = 32.57;
component.cg(i,2) = 4.53;

i = i + 1;
component.weight(i) = 17.5; % horizontal tail
component.cg(i,1) = 35.72;
component.cg(i,2) = 6.99;

% i = i + 1;
% component.weight(i) = 723.9; % fuselage
% component.cg(i,1) = 18.1;
% component.cg(i,2) = 0;

% i = i + 1;
% component.weight(i) = 160.0; % main landing gear
% component.cg(i,1) = 18.24;
% component.cg(i,2) = -2.75;

i = i + 1;
component.weight(i) = 71.5; %nose landing gear
component.cg(i,1) = 4.92;
component.cg(i,2) = -2.75;

i = i + 1;
component.weight(i) = 128.6; % nacelle
component.cg(i,1) = 30.17;
component.cg(i,2) = 0;

i = i + 1;
component.weight(i) = 710.5; %engines
component.cg(i,1) = 30.17;
component.cg(i,2) = 0;

i = i + 1;
component.weight(i) = 85.2; %engine controls
component.cg(i,1) = 30.17;
component.cg(i,2) = 0;

i = i + 1;
component.weight(i) = 35.8; %engine starter
component.cg(i,1) = 30.17;
component.cg(i,2) = 0;

i = i + 1;
component.weight(i) = 60.2; % fuel system
component.cg(i,1) = 15.14;
component.cg(i,2) = -2.75;

i = i + 1;
component.weight(i) = 359.2; % flight controls
component.cg(i,1) = 4.5;
component.cg(i,2) = -0.66;

i = i + 1;
component.weight(i)= 108.1; %Instruments
component.cg(i,1) = 1.67;
component.cg(i,2) = -1;

i = i + 1;
component.weight(i) = 303.1; %Avionics
component.cg(i,1) = 1.67;
component.cg(i,2) = -1;

i = i + 1;
component.weight(i) = 55.7; % hydraulics
component.cg(i,1) = 1.83;
component.cg(i,2) = -1.67;

i = i + 1;
component.weight(i) = 530; % electrical system
component.cg(i,1) = 5;
component.cg(i,2) = -1.67;

i = i + 1;
component.weight(i) = 109.8; %air-con
component.cg(i,1) = 22;
component.cg(i,2) = -1.67;

i = i + 1;
component.weight(i) = 13.8; %anti-ice
component.cg(i,1) = 22;
component.cg(i,2) = -1.67;

i = i + 1;
component.weight(i) = 171.8; %furnishings
component.cg(i,1) = 20;
component.cg(i,2) = -2.5;

% i = i + 1;
% component.weight(i) = 1601.1;% full fuel weight
% component.cg(i,1) = 15.14;
% component.cg(i,2) = -2.75;

%%PAYLOAD AND CREW%% -> maybe model these as distributed loads?
payload_factor = 1;

i = i + 1;
component.weight(i) = 1243.4*payload_factor; % full passengers and crew??
component.cg(i,1) = 15.68;
component.cg(i,2) = -1.25;

i = i + 1;
component.weight(i) = 198.4*payload_factor; % full luggage
component.cg(i,1) = 8.92;
component.cg(i,2) = -1.25;
%%


%%CONVERSION TO SI%%
component.weight = component.weight*lbs_to_n;
component.cg = component.cg*ft_to_m;

%%FUEL%%

%aux fuel tank cg in line with the main wing fuel tanks
fuel_weight_percent = 1;

aux_fuel_tank_vol = 0.3087; %m^3
aux_full_fuel_mass = aux_fuel_tank_vol*775;%kg
aux_full_fuel_weight = aux_full_fuel_mass*9.81;%N
tot_full_fuel_weight = 1601.1*lbs_to_n;
aux_fuel_weight_fraction = aux_full_fuel_weight/tot_full_fuel_weight;
wing_fuel_weight_fraction = 1 - aux_fuel_weight_fraction;

current_fuel_weight = fuel_weight_percent * tot_full_fuel_weight;
aux_fuel_weight = current_fuel_weight*aux_fuel_weight_fraction;
wing_fuel_weight = current_fuel_weight*wing_fuel_weight_fraction;

fuel_cg(1) = 15.14*ft_to_m;
fuel_cg(2) = -2.75*ft_to_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%Fuselage Weight Distribution%%
fuselage_total_weight = 723.9*lbs_to_n;

%distribution is calculated from stl which has end of cabin as datum 
%so we need to offset the volume distribution

offset = min(fsolve(@(x) voldist.fit(x), min(voldist.range)));

fuse_unit_dist = @(x) (voldist.fit(x+offset)./max(voldist.totalvol));
fuse_weight_dist = @(x) fuse_unit_dist(x) * fuselage_total_weight;

%%MAIN WING%% -> assuming chord at root - this needs to be changed later
wing.weight = 189.6*lbs_to_n;
wing.cg = 15.32*ft_to_m;
wing.x = 3.3156;
wing.croot = 1.7674;
x_ac_w = 4.3620;

front_spar_c = 0.1;
rear_spar_c = 0.77;

front_spar_x = wing.x + front_spar_c*wing.croot;
rear_spar_x = wing.x + rear_spar_c*wing.croot;


%%Horizontal Tailplane CONSTANTS%%
x_ac_h = 10.3474;

%%MAIN LANDING GEAR%%
mlg.weight = 160*lbs_to_n;
mlg.cg= 5.16305;

%LOAD DISTRIBUTION%

%%Calculating sectional loading of fuselage 
load_dist.x =  linspace(0,11.73,1000)';
load_dist.spacing = abs(load_dist.x(2)-load_dist.x(1));
load_dist.load = zeros(size(load_dist.x));

%%load distribution for Fuselage Only%%

%add fuselage distributed load modelled as a distributuion of point loads i.e multiplied by discretization
load_dist.load = load_dist.load + (fuse_weight_dist(load_dist.x) * load_dist.spacing); 
%add aux fuel tank load
load_dist.load = add_point_load(load_dist.load,load_dist.x,aux_fuel_weight,fuel_cg(1));
%add component weights
for i=1:length(component.weight)
    load_dist.load = add_point_load(load_dist.load,load_dist.x,component.weight(i),component.cg(i,1)); 
end

%intertial weight acts downards so make load negative
load_dist.load = -load_dist.load;

[load_dist.Fi,load_dist.Fi_cg] = point_load_from_dist(load_dist.x,load_dist.load); % calculate the effective total load and cg of the distributed load
%%%%%%%%%%%

%functors for force and moment distributions 
sectional_load = @(Ff,Fr,Ft) load_dist.load + dist_from_point(load_dist.x,Ff,front_spar_x) + dist_from_point(load_dist.x,Fr,rear_spar_x) + dist_from_point(load_dist.x,Ft,x_ac_h);
%shear force is simply addtion of the point load distribution
shear_force = @(s_load) cumsum(s_load);
bending_moment = @(x_dist,s_force) cumtrapz(x_dist,s_force);

%find required tail force for trim by taking global equilbrium of aircraft
%add wing weight, fuel in wing,main landing gear note the negative signs as wieght acts downwards
load_dist.global.load = load_dist.load + dist_from_point(load_dist.x,-wing.weight,wing.cg) + dist_from_point(load_dist.x,-wing_fuel_weight,fuel_cg(1)) + dist_from_point(load_dist.x,-mlg.weight,mlg.cg);
[load_dist.global.Fi,load_dist.global.Fi_cg] = point_load_from_dist(load_dist.x,load_dist.global.load);

%force eq
% Lw + Lt + gFi = 0
%moment eq
% x_ac_w*Lw + x_ac_h*Lt + gFi_cg * gFi = 0
res = [1 1;x_ac_w x_ac_h]\[-load_dist.global.Fi;-load_dist.global.Fi*load_dist.global.Fi_cg];
sym_flight.Lw = res(1); % wing lift force
sym_flight.Lt = res(2); %tail lift force

%calculating local equilbirum considering only the fuselage loads now
%force eq
%Ff + Fr + Lt + Fi = 0
%moment eq
%f_x * Ff + r_x * Fr + x_ac_h*Lt + Fi_cg * cg = 0
res = [1 1;front_spar_x rear_spar_x]\[-load_dist.Fi-sym_flight.Lt;(-load_dist.Fi*load_dist.Fi_cg) + (-sym_flight.Lt * x_ac_h)];
sym_flight.Ff = res(1);
sym_flight.Fr = res(2);


% %calculating local equilbirum considering only the fuselage loads now
% %force eq
% %Ff + Fr + Lt + Fi = 0
% %moment eq
% %f_x * Ff + r_x * Fr + x_ac_h*Lt + Fi_cg * cg = 0
% res = [1 1;front_spar_x rear_spar_x]\[-load_dist.Fi-sym_flight.Lt;(-load_dist.Fi*load_dist.Fi_cg) + (-sym_flight.Lt * x_ac_h)];
% sym_flight.Ff = res(1);
% sym_flight.Fr = res(2);


sym_flight.x = load_dist.x;
sym_flight.load = sectional_load(sym_flight.Ff,sym_flight.Fr,sym_flight.Lt);
sym_flight.Fi = load_dist.Fi;
sym_flight.Fi_cg = load_dist.Fi_cg;
sym_flight.shear = shear_force(sym_flight.load);
sym_flight.bending_moment = bending_moment(sym_flight.x,sym_flight.shear);
sym_flight.shear_fit = fit(sym_flight.x,sym_flight.shear,'linearinterp');
sym_flight.bending_moment_fit = fit(sym_flight.x,sym_flight.bending_moment,'linearinterp');


fuselage_load_plots(sym_flight);


fuselageLoading.sym_flight = sym_flight;

%SOLVE FOR FRONT OFF CASE -> NO AERO LOADS 

%find required tail force for trim by taking global equilbrium of aircraft
res = [1 1;mlg.cg x_ac_h]\[-load_dist.global.Fi;-load_dist.global.Fi*load_dist.global.Fi_cg];
front_off.Fuc = res(1); % main landing gear force
front_off.Lt = res(2); %tail lift force

%calculating local equilbirum considering only the fuselage loads now
res = [1 1;front_spar_x rear_spar_x]\[-load_dist.Fi-front_off.Lt;(-load_dist.Fi*load_dist.Fi_cg) + (-front_off.Lt * x_ac_h)];
front_off.Ff = res(1);
front_off.Fr = res(2);

front_off.x = load_dist.x;
front_off.load = sectional_load(front_off.Ff,front_off.Fr,front_off.Lt);
front_off.Fi = load_dist.Fi;
front_off.Fi_cg = load_dist.Fi_cg;
front_off.shear = shear_force(front_off.load);
front_off.bending_moment = bending_moment(front_off.x,front_off.shear);
front_off.shear_fit = fit(front_off.x,front_off.shear,'linearinterp');
front_off.bending_moment_fit = fit(front_off.x,front_off.bending_moment,'linearinterp');


fuselage_load_plots(front_off);


fuselageLoading.front_off = front_off;


save('fuselageLoading.mat','fuselageLoading');