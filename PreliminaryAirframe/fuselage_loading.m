clear
clc
load('voldist.mat')

lbs_to_n = 4.4482216153;
ft_to_m = 0.3048;

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

%%WING SPAR CONSTANTS%% -> assuming chord at root - this needs to be changed later
x_wing = 3.3156;
croot = 1.7674;

front_spar_c = 0.1;
rear_spar_c = .77;

front_spar_x = x_wing + front_spar_c*croot;
rear_spar_x = x_wing + rear_spar_c*croot;


%%Horizontal Tailplane CONSTANTS%%
x_ac_h = 10.3474;

%LOAD DISTRIBUTION%


%!!!! POSITIVE DIRECTION DOWNWARDS!!!!

%%Calculating sectional loading of fuselage 
load_dist.x =  linspace(0,11.73,10000)';
load_dist.load = zeros(size(load_dist.x));
%add fuselage distributed load
load_dist.load = load_dist.load + (fuse_weight_dist(load_dist.x)); 
%add aux fuel tank load
load_dist.load = add_point_load(load_dist.load,load_dist.x,aux_fuel_weight,fuel_cg(1));
%add component weights
for i=1:length(component.weight)
    load_dist.load = add_point_load(load_dist.load,load_dist.x,component.weight(i),component.cg(i,1)); % gravity is down
end

[load_dist.Fi,load_dist.Fi_cg] = point_load_from_dist(load_dist.x,load_dist.load);

%CALCULATE TAIL LOAD

%force equilbrium
force_eq = @(Ff,Fr,Ft) sum(load_dist.load) - Ff - Fr - Ft;
%define ccw as positive
%moments taken from ac datum i.e the nose
moment_eq = @(Ff,Fr,Ft) (-sum((load_dist.load).*(load_dist.x))) + (Ff*front_spar_x) + (Fr*rear_spar_x) + (Ft*x_ac_h);

shear_force = @(Ff,Fr,Ft) cumsum(load_dist.load + dist_from_point(load_dist.x,-Ff,front_spar_x) + dist_from_point(load_dist.x,-Fr,rear_spar_x) + dist_from_point(load_dist.x,-Ft,x_ac_h));

bending_moment = @(Ff,Fr,Ft) cumtrapz(load_dist.x,shear_force(Ff,Fr,Ft));

%enforce BC at end of bending moment distribution
index_at = @(expr,idx) expr(idx);

F = @(x) [force_eq(x(1),x(2),x(3)); moment_eq(x(1),x(2),x(3)); index_at(bending_moment(x(1),x(2),x(3)),end)]

res = fsolve(@(x) F(x),[0 0 0]);
Ff = res(1);
Fr = res(2);
Ft = res(3);


sym_flight.x = load_dist.x;
sym_flight.load = load_dist.load + dist_from_point(load_dist.x,Ff,front_spar_x) + dist_from_point(load_dist.x,Fr,rear_spar_x) + dist_from_point(load_dist.x,Ft,x_ac_h);
sym_flight.Fi = load_dist.Fi;
sym_flight.Fi_cg = load_dist.Fi_cg;
sym_flight.shear = shear_force(Ff,Fr,Ft);
sym_flight.bending_moment = bending_moment(Ff,Fr,Ft);
sym_flight.shear_fit = fit(sym_flight.x,sym_flight.shear,'linearinterp');
sym_flight.bending_moment_fit = fit(sym_flight.x,sym_flight.bending_moment,'linearinterp');

fuselage_load_plots(sym_flight);


fuselageLoading.sym_flight = sym_flight;

%%SOLVE FOR FRONT OFF CASE -> NO AERO LOADS 
% constants for offset
x_gear = 4.6941;
a = x_gear - front_spar_x;
b = x_gear - rear_spar_x;

Ff_frontoff = @(Fuc) Fuc*(1-(a/(b-a)));
Fr_frontoff = @(Fuc) Fuc*(a/(b-a));

F = @(x) [force_eq(Ff_frontoff(x(1)),Fr_frontoff(x(1)),x(2));moment_eq(Ff_frontoff(x(1)),Fr_frontoff(x(1)),x(2))];

res = fsolve(@(x) F(x),[0 0]);
Fuc = res(1);
Ft = res(2);
Ff = Ff_frontoff(Fuc);
Fr = Fr_frontoff(Fuc);

front_off.x = load_dist.x;
front_off.load = load_dist.load + dist_from_point(load_dist.x,Ff,front_spar_x) + dist_from_point(load_dist.x,Fr,rear_spar_x) + dist_from_point(load_dist.x,Ft,x_ac_h);
front_off.Fi = load_dist.Fi;
front_off.Fi_cg = load_dist.Fi_cg;
front_off.shear = shear_force(Ff,Fr,Ft);
front_off.bending_moment = bending_moment(Ff,Fr,Ft);
front_off.shear_fit = fit(front_off.x,front_off.shear,'linearinterp');
front_off.bending_moment_fit = fit(front_off.x,front_off.bending_moment,'linearinterp');

fuselage_load_plots(front_off);

fuselageLoading.front_off = front_off;


save('fuselageLoading.mat','fuselageLoading');