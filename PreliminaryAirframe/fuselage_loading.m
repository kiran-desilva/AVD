clear
clc
load('voldist.mat')

lbs_to_n = 4.4482216153;
ft_to_m = 0.3048;

%%Weights and CG from Report UNITS IN LBS AND FT

i = 0;

% i = i + 1;
% weight(i) = 189.6; %wing
% cg(i,1) = 15.32;
% cg(i,2) = -2.75;

i = i + 1;
weight(i) = 17.1; % vertical tail
cg(i,1) = 32.57;
cg(i,2) = 4.53;

i = i + 1;
weight(i) = 17.5; % horizontal tail
cg(i,1) = 35.72;
cg(i,2) = 6.99;

% i = i + 1;
% weight(i) = 723.9; % fuselage
% cg(i,1) = 18.1;
% cg(i,2) = 0;

i = i + 1;
weight(i) = 160.0; % main landing gear
cg(i,1) = 18.24;
cg(i,2) = -2.75;

i = i + 1;
weight(i) = 71.5; %nose langind gear
cg(i,1) = 4.92;
cg(i,2) = -2.75;

i = i + 1;
weight(i) = 128.6; % nacelle
cg(i,1) = 30.17;
cg(i,2) = 0;

i = i + 1;
weight(i) = 710.5; %engines
cg(i,1) = 30.17;
cg(i,2) = 0;

i = i + 1;
weight(i) = 85.2; %engine controls
cg(i,1) = 30.17;
cg(i,2) = 0;

i = i + 1;
weight(i) = 35.8; %engine starter
cg(i,1) = 30.17;
cg(i,2) = 0;

i = i + 1;
weight(i) = 60.2; % fuel system
cg(i,1) = 15.14;
cg(i,2) = -2.75;

i = i + 1;
weight(i) = 359.2; % flight controls
cg(i,1) = 4.5;
cg(i,2) = -0.66;

i = i + 1;
weight(i)= 108.1; %Instruments
cg(i,1) = 1.67;
cg(i,2) = -1;

i = i + 1;
weight(i) = 303.1; %Avionics
cg(i,1) = 1.67;
cg(i,2) = -1;

i = i + 1;
weight(i) = 55.7; % hydraulics
cg(i,1) = 1.83;
cg(i,2) = -1.67;

i = i + 1;
weight(i) = 530; % electrical system
cg(i,1) = 5;
cg(i,2) = -1.67;

i = i + 1;
weight(i) = 109.8; %air-con
cg(i,1) = 22;
cg(i,2) = -1.67;

i = i + 1;
weight(i) = 13.8; %anti-ice
cg(i,1) = 22;
cg(i,2) = -1.67;

i = i + 1;
weight(i) = 171.8; %furnishings
cg(i,1) = 20;
cg(i,2) = -2.5;

% i = i + 1;
% weight(i) = 1601.1;% full fuel weight
% cg(i,1) = 15.14;
% cg(i,2) = -2.75;

% i = i + 1;
% weight(i) = 1243.4; % full passengers and crew??
% cg(i,1) = 15.68;
% cg(i,2) = -1.25;

% i = i + 1;
% weight(i) = 198.4; % full luggage
% cg(i,1) = 8.92;
% cg(i,2) = -1.25;

%%CONVERSION TO SI%%
weight = weight*lbs_to_n;
cg = cg*ft_to_m;


%aux fuel tank cg in line with the main wing fuel tanks
aux_fuel_tank_vol = 0.3087; %m^3
aux_fuel_mass = aux_fuel_tank_vol*775;%kg
aux_fuel_weight = aux_fuel_mass*9.81;%N
tot_fuel_weight = 1601.1*lbs_to_n;
wing_fuel_weight = tot_fuel_weight - aux_fuel_weight;
fuel_cg(1) = 15.14;
fuel_cg(2) = -2.75;

%%Fuselage Weight Distribution
fuselage_total_weight = 723.9*lbs_to_n;

%distribution is calculated from stl which has end of cabin as datum 
%so we need to offset the volume distribution

offset = min(fsolve(@(x) voldist.fit(x), min(voldist.range)));

fuse_unit_dist = @(x) (voldist.fit(x+offset)./max(voldist.totalvol));
fuse_mass_dist = @(x) fuse_unit_dist(x) * fuselage_total_weight;


%%Calculating sectional loading of fuselage

discretized_x = linspace(0,11.73,10000);
discretized_loading = fuse_mass_dist(discretized_x); %first add fuselage distributed load

%loop thru all components
for i=1:length(weight)
    %find closest station
    [~,idx] = min(abs(discretized_x - cg(i,1)));
    discretized_loading(idx) = discretized_loading(idx) + weight(i); % add point load to current value
end


%%Sectional Loading

figure
hold on

plot(discretized_x,discretized_loading)

xlabel('X coordinate from Tip (m)')
ylabel('Sectional Loading (N/m)')
grid on
grid minor

%%Shear force
sectional_loading_fit = fit(discretized_x',discretized_loading,'linearinterp');
discretized_shear_force = integrate(sectional_loading_fit,discretized_x,0);

figure
hold on

plot(discretized_x,discretized_shear_force)

xlabel('X coordinate from Tip (m)')
ylabel('Shear Force (N)')
grid on
grid minor

%%Bending Moment

shear_force_fit = fit(discretized_x',discretized_shear_force','linearinterp');
discretized_bending_moment = integrate(shear_force_fit,discretized_x,0);

figure
hold on

plot(discretized_x,discretized_bending_moment)

xlabel('X coordinate from Tip (m)')
ylabel('Bending Moment (Nm)')
grid on
grid minor
