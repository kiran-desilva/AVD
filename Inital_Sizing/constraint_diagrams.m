clear
clc
clf

u = symunit; % used for units
subindex = @(A, idx) A(idx); %% function for anonymous indexing

open("parameters.mat");
open("sizing.mat");

%% shit we want
%% l/d max - l_dmax
%% specific fuel consumption - sfc
%% gross weight - w0
%% operating empty weight - we
%% fuel weight - wf

%{

We don't seem to be using any of these??

sizing.ld_max = 10;
sizing.ld_cruise = 0.866*sizing.ld_max;
sizing.sfc = 10;

sizing.service_ceiling = 40000; % ft

sizing.w_landing_w_total_max = 0.8;
%}

sizing.absolute_ceiling = 45000; % ft
sizing.AR = 5;
sizing.maxTakeoffWeight = 5000*9.81; % TODO: -> from chorley are we sure this is in Newtons
sizing.cl_max = 1.5; % TODO:        
sizing.e = 0.8; %% e will change with the addition of hld
sizing.k = 1/(sizing.AR*pi*sizing.e); % This eq. is used in Gudmundsson p. 64 -> idk if this is valid 
sizing.cd_min = 0.05; % TODO: 
sizing.cd_0 = 0.05; % TODO:
sizing.sref = 73.73; % TODO:
sizing.n = 1/cosd(40); % max load factor as defined by far
sizing.runway_length = 1200; %% in meters
sizing.engine_number = 2;
sizing.rho_0 = 1.225; % reference density

q = @(V_inf, rho) 0.5*rho*V_inf^2;

%% takeoff 
sizing.takeoff.rho = 1.225;
sizing.takeoff.cl_max = 1.6; % TODO: -> raymer?
sizing.takeoff.v_stall = sqrt((sizing.maxTakeoffWeight/sizing.sref) * 2/(sizing.takeoff.rho * sizing.takeoff.cl_max));
sizing.takeoff.v_inf = 1.3*sizing.takeoff.v_stall; % we should check if this is reasonable... Errikos did it in his excel but idk -> from FAR25
sizing.takeoff.q = q(sizing.takeoff.v_inf, sizing.takeoff.rho);

%% climb
sizing.climb.rho = sizing.takeoff.rho;
sizing.climb.v_inf= sizing.takeoff.v_inf; % we should check if this is reasonable... Errikos did it in his excel but idk
sizing.climb.q = q(sizing.climb.v_inf, sizing.climb.rho);
sizing.climb.dh_dt = 12200/(30*60); % TODO



%% cruise
[~, sizing.cruise.a, ~, sizing.cruise.rho] = atmosisa(distdim(40000, 'ft', 'm'));
sizing.cruise.v_inf= 0.75*sizing.cruise.a;
sizing.cruise.q = q(sizing.cruise.v_inf, sizing.cruise.rho);


%% loiter
[~, ~, ~, sizing.loiter.rho] = atmosisa(distdim(5000, 'ft', 'm'));
sizing.loiter.v_inf = sizing.cruise.v_inf*(1/3)^(0.25); %TODO: this would be the exact equation for VminD, we need the loiter weight for that tho sqrt(1.225/sizing.loiter.rho)*pow((4*sizing.k*sizing.loiter.w^2)/(pi*sizing.cd_0*sizing.AR), 0.25);
sizing.loiter.q = q(sizing.loiter.v_inf, sizing.loiter.rho);

%% landing
sizing.landing.rho = sizing.takeoff.rho;
sizing.landing.cl_max = sizing.takeoff.cl_max; % TODO:
sizing.landing.obstacle_height = 183; %m
sizing.landing.kr = 0.66; % with thrust reversers










%% What's next?
%% Make sure aircraft can complete the stipulated design mission profile
%% Make sure aircraft capable of achieving performance targets like
%% Ceilings (absolute, service, combat) -> T/W at height, v has to be bigger than 1/(L/D)_max
%% Maximum speed
%% Time to climb / Rates of climb
%% Sustained turn rates / radii
%% Level (axial) acceleration
%% Takeoff and Landing distances (TODA & LDA) -> T/W, stall speed 
%% Stall speed -> function of Cl_max (Cl_max changes in takeoff and clean config) and wing loading (W/ S)
%% Make sure aircraft meets airworthiness requirements

%% FAR25 - https://www.engineerstoolkit.com/Airworthiness%20Standards%20%20FAA%20FAR%20Part%2025.pdf
%% take off speed = 1.1* stall speed


syms wing_loading thrust_to_weight V_inf rho turn_height_m % use syms for constraints
syms alpha beta % alpha -> scales weight to initial weight, beta -> scales thrust to sea level thrust

%% Computing climb constraint

climb_constraint(wing_loading) = sizing.climb.dh_dt/sizing.climb.v_inf + sizing.climb.q/wing_loading*sizing.cd_min + sizing.k/sizing.climb.q*wing_loading;

%% Cruise constraint

cruise_constraint(wing_loading, alpha, beta) = alpha/beta*(sizing.cruise.q*sizing.cd_min/(alpha*wing_loading) + alpha*sizing.k*wing_loading/sizing.cruise.q);

%% Take-off distance constraint

takeoff_bfl_constraint(wing_loading) = ( wing_loading*(0.297-0.019*sizing.engine_number) ) / ( sizing.runway_length * (sizing.takeoff.rho / sizing.rho_0) * sizing.takeoff.cl_max) ; 


take_off_constraint(wing_loading) = double(separateUnits(unitConvert(37.5*(u.ft^3/u.lbf), u.m^3/u.N)))*wing_loading/(sizing.takeoff.cl_max*sizing.runway_length);
% take_off_constraint(wing_loading) = 37.5*wing_loading/(sizing.takeoff.cl_max*sizing.runway_length);

%% service ceiling constraint
sizing.service_ceiling_rho = sizing.cruise.rho;
sizing.service_climb_velocity_at_ceiling = 100; % TODO: what is the rate of climb at the ceiling alt we want? Look at far mby...

thing = sqrt(sizing.k/(3*sizing.cd_min));
service_ceiling_constraint(wing_loading) = sizing.service_climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.service_ceiling_rho) + 4*sqrt(sizing.k*sizing.cd_min/3);

%% absolute ceiling constraint

[~, ~, ~, sizing.absolute_ceiling_rho] = atmosisa(distdim(sizing.absolute_ceiling,'ft','km'));

sizing.absolute_climb_velocity_at_ceiling = 0; % at absolute ceiling, climb rate should be zero

thing = sqrt(sizing.k/(3*sizing.cd_min));
absolute_ceiling_constraint(wing_loading) = sizing.absolute_climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.service_ceiling_rho) + 4*sqrt(sizing.k*sizing.cd_min/3);

%% turn constraint

turn_constraint = @(wing_loading, turn_height_m, V_inf) q(V_inf, subindex(atmosisa(turn_height_m), 4))*(sizing.cd_min/wing_loading + sizing.k*wing_loading*(sizing.n/q(V_inf, subindex(atmosisa(turn_height_m), 4)))^2);

%% landing constraint
%% page 111 in roskam
%% 
%% landing_distance_thrust_reversal_constraint = sizing.runway_length - 
% TODO: Add constraint with thrust reversals?
landing_constraint_wing_loading = sizing.runway_length*sizing.landing.rho*sizing.landing.cl_max/(0.6*1.3*double(separateUnits(unitConvert(u.ft/u.kts^2, u.m/(u.m/u.s)^2))))^2;

%% stall constraint

%TODO: stall_constraint
stall_constraint = @(cl_max,q_stall) q_stall*cl_max;


hold on

weight_loading_interval = [0, 4000];
fplot(climb_constraint, weight_loading_interval);
fplot(@(wing_loading) cruise_constraint(wing_loading, 1, 1), weight_loading_interval);
fplot(take_off_constraint, weight_loading_interval);
fplot(service_ceiling_constraint, weight_loading_interval);
fplot(absolute_ceiling_constraint, weight_loading_interval);
fplot(@(wing_loading) turn_constraint(wing_loading, distdim(40000, 'ft', 'm'), sizing.cruise.v_inf), weight_loading_interval);
xline(landing_constraint_wing_loading)

legend('Climb', 'Cruise', 'Take-Off', 'Service Ceiling', 'Absolute Ceiling', 'Turn', 'Landing');

hold off;



