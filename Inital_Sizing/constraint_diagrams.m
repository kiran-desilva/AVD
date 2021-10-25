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
sizing.AR = 7.8;
sizing.maxTakeoffWeight = 5000*9.81; % TODO: -> from chorley are we sure this is in Newtons

%% clean
sizing.cl_max = 1.6; % TODO: 
sizing.e = 0.6364; %% e will change with the addition of hld 
sizing.cd_0 = 0.02; % TODO:

%% cd 0 taken from errikos' slides
sizing.undercarriage_e = -0.05;
sizing.undercarriage_cd0 = 0.02;

sizing.t_o_flaps_e = -0.05;
sizing.t_o_flaps_cd0 = 0.02;

sizing.landing_flaps_e = -0.1;
sizing.landing_flaps_cd0 = 0.07;

sizing.k = 1/(sizing.AR*pi*sizing.e); % This eq. is used in Gudmundsson p. 64 -> idk if this is valid 


sizing.sref = 73.73; % TODO:
sizing.n = 1/cosd(40); % max load factor as defined by far
sizing.runway_length = 1200; %% in meters
sizing.engine_number = 2;
sizing.rho_0 = 1.225; % reference density

q = @(V_inf, rho) 0.5*rho*V_inf^2;


%% takeoff 
sizing.takeoff.rho = 1.225;


sizing.takeoff.cl_max = 1.9; %  roskam
sizing.takeoff.e = sizing.e + sizing.undercarriage_e + sizing.t_o_flaps_e;
sizing.takeoff.cd_0 = sizing.cd_0 + sizing.undercarriage_cd0 + sizing.t_o_flaps_cd0;

%% landing
sizing.landing.rho = sizing.takeoff.rho;
sizing.landing.cl_max = 2.2; % roskam:
sizing.landing.obstacle_height = 183; %m
sizing.landing.kr = 0.66; % with thrust reversers

%% landing constraint
%% page 111 in roskam
%% 
%% landing_distance_thrust_reversal_constraint = sizing.runway_length - 
% TODO: Add constraint with thrust reversals?
landing_constraint_wing_loading_roskam = sizing.runway_length*sizing.landing.rho*sizing.landing.cl_max/((0.6*1.3)^2*double(separateUnits(unitConvert(u.ft/u.kts^2, u.m/(u.m/u.s)^2))));
landing_constraint_wing_loading = (sizing.runway_length - sizing.landing.obstacle_height)*sizing.landing.cl_max/0.51;
landing_constraint_wing_loading_trev = (sizing.runway_length - sizing.landing.obstacle_height)*sizing.landing.cl_max/(0.51*0.66);

takeoff_constraint_wing_loading = landing_constraint_wing_loading / 0.72;

% sizing.takeoff.v_stall = sqrt((sizing.maxTakeoffWeight/sizing.sref) * 2/(sizing.takeoff.rho * sizing.takeoff.cl_max));
sizing.landing.v_stall = sqrt((2*landing_constraint_wing_loading)/(sizing.landing.rho*sizing.landing.cl_max));
sizing.takeoff.v_stall = sqrt((2*takeoff_constraint_wing_loading)/(sizing.takeoff.rho*sizing.takeoff.cl_max));
sizing.v_stall = sizing.takeoff.v_stall*sqrt(sizing.takeoff.cl_max/sizing.cl_max);
sizing.takeoff.v_inf = 1.3*sizing.takeoff.v_stall; % we should check if this is reasonable... Errikos did it in his excel but idk -> from FAR25
sizing.takeoff.q = q(sizing.takeoff.v_inf, sizing.takeoff.rho);


sizing.takeoff.initial_climb.q = q(1.2*sizing.takeoff.v_stall, sizing.takeoff.rho);
sizing.takeoff.initial_climb.e = sizing.e + sizing.t_o_flaps_e;
sizing.takeoff.initial_climb.cd_0 = sizing.cd_0 + sizing.t_o_flaps_cd0;

sizing.takeoff.transition.q = sizing.takeoff.initial_climb.q;
sizing.takeoff.transition.e = sizing.e + sizing.t_o_flaps_e + sizing.undercarriage_e;
sizing.takeoff.transition.cd_0 = sizing.cd_0 + sizing.t_o_flaps_cd0 + sizing.undercarriage_cd0;


sizing.takeoff.second_segment.q = sizing.takeoff.initial_climb.q;
sizing.takeoff.second_segment.e = sizing.e + sizing.t_o_flaps_e;
sizing.takeoff.second_segment.cd_0 = sizing.cd_0 + sizing.t_o_flaps_cd0;


sizing.takeoff.en_route_climb.q = q(1.25*sizing.v_stall, sizing.takeoff.rho);
sizing.takeoff.en_route_climb.e = sizing.e;
sizing.takeoff.en_route_climb.cd_0 = sizing.cd_0;

sizing.landing.first.q = q(1.3*sizing.landing.v_stall, sizing.landing.rho);
sizing.landing.first.e = sizing.e + sizing.landing_flaps_e + sizing.undercarriage_e;
sizing.landing.first.cd_0 = sizing.cd_0 + sizing.landing_flaps_cd0 + sizing.undercarriage_cd0;


sizing.landing.second.q = q(1.5*1.1*sizing.landing.v_stall, sizing.landing.rho);
sizing.landing.second.e = sizing.e + sizing.t_o_flaps_e + sizing.undercarriage_e;
sizing.landing.second.cd_0 = sizing.cd_0 + sizing.t_o_flaps_cd0 + sizing.undercarriage_cd0;


%% climb
sizing.climb.rho = sizing.takeoff.rho;
sizing.climb.v_inf= sizing.takeoff.v_inf*1.5;
sizing.climb.q = q(sizing.climb.v_inf, sizing.climb.rho);
sizing.climb.dh_dt = double(separateUnits(unitConvert(2500*u.ft/u.min, u.m/u.s))); % TODO



%% cruise
[~, sizing.cruise.a, ~, sizing.cruise.rho] = atmosisa(distdim(40000, 'ft', 'm'));
sizing.cruise.v_inf= 0.75*sizing.cruise.a;
sizing.cruise.q = q(sizing.cruise.v_inf, sizing.cruise.rho);


%% loiter
[~, ~, ~, sizing.loiter.rho] = atmosisa(distdim(5000, 'ft', 'm'));
sizing.loiter.v_inf = sizing.cruise.v_inf*(1/3)^(0.25); %TODO: this would be the exact equation for VminD, we need the loiter weight for that tho sqrt(1.225/sizing.loiter.rho)*pow((4*sizing.k*sizing.loiter.w^2)/(pi*sizing.cd_0*sizing.AR), 0.25);
sizing.loiter.q = q(sizing.loiter.v_inf, sizing.loiter.rho);





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


syms wing_loading thrust_to_weight V_inf rho turn_height_m  climb_grad k cd_0 q_var % use syms for constraints
syms alpha beta % alpha -> scales weight to initial weight, beta -> scales thrust to sea level thrust

%% Computing climb constraint

% climb_constraint(wing_loading) = sizing.climb.dh_dt/sizing.climb.v_inf + sizing.climb.q/wing_loading*sizing.cd_0 + sizing.k/sizing.climb.q*wing_loading;
climb_constraint(wing_loading, climb_grad) = climb_grad + sizing.climb.q/wing_loading*sizing.cd_0 + sizing.k/sizing.climb.q*wing_loading;

%% Cruise constraint

cruise_constraint(wing_loading, alpha, beta) = alpha/beta*(sizing.cruise.q*sizing.cd_0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/sizing.cruise.q);

%% Take-off distance constraint

takeoff_bfl_constraint(wing_loading) = ( wing_loading*(0.297-0.019*sizing.engine_number) ) / ( sizing.runway_length * (sizing.takeoff.rho / sizing.rho_0) * sizing.takeoff.cl_max) ; 


take_off_constraint(wing_loading) = double(separateUnits(unitConvert(37.5*(u.ft^3/u.lbf), u.m^3/u.N)))*wing_loading/(sizing.takeoff.cl_max*sizing.runway_length);
% take_off_constraint(wing_loading) = 37.5*wing_loading/(sizing.takeoff.cl_max*sizing.runway_length);

% max velocity constraint
max_velocity_q = q(sizing.cruise.a*0.78, sizing.cruise.rho);
max_velocity_constraint(wing_loading, alpha, beta) =alpha/beta*(max_velocity_q*sizing.cd_0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/max_velocity_q);
%  max_velocity_q/wing_loading*sizing.cd_0 + sizing.k/max_velocity_q*wing_loading;

%% service ceiling constraint
sizing.service_ceiling_rho = sizing.cruise.rho;
sizing.service_climb_velocity_at_ceiling = 6; % TODO: what is the rate of climb at the ceiling alt we want? Look at far mby...

thing = sqrt(sizing.k/(3*sizing.cd_0));
service_ceiling_constraint(wing_loading) = sizing.service_climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.service_ceiling_rho) + 4*sqrt(sizing.k*sizing.cd_0/3);

%% absolute ceiling constraint

[~, ~, ~, sizing.absolute_ceiling_rho] = atmosisa(distdim(sizing.absolute_ceiling,'ft','km'));

% sizing.absolute_climb_velocity_at_ceiling = 0; % at absolute ceiling, climb rate should be zero

% thing = sqrt(sizing.k/(3*sizing.cd_0));
% absolute_ceiling_constraint(wing_loading,alpha,beta) = sizing.absolute_climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.service_ceiling_rho) + 4*sqrt(sizing.k*sizing.cd_0/3);

Vx(wing_loading) = sqrt( (2/sizing.absolute_ceiling_rho) * (wing_loading) * sqrt(sizing.k/sizing.cd_0) * 1);

absolute_ceiling_constraint(wing_loading,alpha,beta) = simplify((alpha/beta) * ( ( (0.5*sizing.absolute_ceiling_rho*(Vx^2)*sizing.cd_0)/(wing_loading) ) + ( ((1^2)*(wing_loading))/(0.5*sizing.absolute_ceiling_rho*(Vx^2)*pi*sizing.AR*sizing.e) ) ));

%% turn constraint

turn_constraint = @(wing_loading, turn_height_m, V_inf) q(V_inf, atmos(turn_height_m, 4))*(sizing.cd_0/wing_loading + sizing.k*wing_loading*(sizing.n/q(V_inf, atmos(turn_height_m, 4)))^2);


%% stall constraint

%TODO: stall_constraint
stall_constraint = @(cl_max,q_stall) q_stall*cl_max;

%% Climb OEI
climb_constraint_oei(wing_loading, climb_grad, q_var, k, cd_0) = 2*(climb_grad + q_var/wing_loading*cd_0 + k/q_var*wing_loading);

hold on

weight_loading_interval = [1, 12000];
% fplot(@(wing_loading) climb_constraint(wing_loading, 3.2/100), weight_loading_interval);
fplot(@(wing_loading) cruise_constraint(wing_loading, 0.8296, 0.24), weight_loading_interval);
fplot(take_off_constraint, weight_loading_interval);
% fplot(service_ceiling_constraint, weight_loading_interval);
%fplot(absolute_ceiling_constraint, weight_loading_interval);
fplot(@(wing_loading) absolute_ceiling_constraint(wing_loading,0.8296, 0.25))
fplot(@(wing_loading) turn_constraint(wing_loading, distdim(40000, 'ft', 'm'), sizing.cruise.v_inf), weight_loading_interval);
xline(landing_constraint_wing_loading,'color','red');
xline(landing_constraint_wing_loading_trev,'color','cyan');
% xline(sizing.maxTakeoffWeight/sizing.sref);
xline(landing_constraint_wing_loading_roskam, 'color', 'magenta');

fplot(@(wing_loading) max_velocity_constraint(wing_loading, 0.8296, 0.25), weight_loading_interval);

k_func = @(e) 1/(pi*sizing.AR*e);
% fplot(@(wing_loading) climb_constraint_oei(wing_loading, 1.2/100, sizing.takeoff.initial_climb.q, k_func(sizing.takeoff.initial_climb.e), sizing.takeoff.initial_climb.cd_0))
% fplot(@(wing_loading) climb_constraint_oei(wing_loading, 0, sizing.takeoff.transition.q, k_func(sizing.takeoff.transition.e), sizing.takeoff.transition.cd_0))
% fplot(@(wing_loading) climb_constraint_oei(wing_loading, 2.4/100, sizing.takeoff.second_segment.q, k_func(sizing.takeoff.second_segment.e), sizing.takeoff.second_segment.cd_0))
% fplot(@(wing_loading) climb_constraint_oei(wing_loading, 1.2/100, sizing.takeoff.en_route_climb.q, k_func(sizing.takeoff.en_route_climb.e), sizing.takeoff.en_route_climb.cd_0))
% fplot(@(wing_loading) 0.5*climb_constraint_oei(wing_loading, 3.2/100, sizing.landing.first.q, k_func(sizing.landing.first.e), sizing.landing.first.cd_0))
% fplot(@(wing_loading) climb_constraint_oei(wing_loading, 2.1/100, sizing.landing.second.q, k_func(sizing.landing.second.e), sizing.landing.second.cd_0))

% pls tell me i did the dumb dumb
optimum_w_over_s = 0.5*sizing.cruise.rho*((sizing.cruise.v_inf/3^0.25)^2)*sqrt(pi*sizing.AR*sizing.cd_0/sizing.k)
% i think u did?
% that aint gonna change anything :(
% who knows lol
% nah it's so high cus v_inf^2 is like 40k...
% I dont think we can cruise at Vmd...
% hmmmmm yeah maybe we are curising too quickly
% lemme read a bit more of the notes...
%%  can we relate t/w to vinf?
% ive asked tanvi to jsut check we arent being dumb
% yeah cool, thanks

%ok so.. if we want max range, we don't cruise at Vmd, but at ~1.3*Vmd, so that should be our cruise speed

yline(double(climb_constraint_oei(takeoff_constraint_wing_loading, 1.2/100, sizing.takeoff.initial_climb.q, k_func(sizing.takeoff.initial_climb.e), sizing.takeoff.initial_climb.cd_0)), 'r--')
yline(double(climb_constraint_oei(takeoff_constraint_wing_loading, 0, sizing.takeoff.transition.q, k_func(sizing.takeoff.transition.e), sizing.takeoff.transition.cd_0)), 'b--')
yline(double(climb_constraint_oei(takeoff_constraint_wing_loading, 2.4/100, sizing.takeoff.second_segment.q, k_func(sizing.takeoff.second_segment.e), sizing.takeoff.second_segment.cd_0)), 'g--')
yline(double(climb_constraint_oei(takeoff_constraint_wing_loading, 1.2/100, sizing.takeoff.en_route_climb.q, k_func(sizing.takeoff.en_route_climb.e), sizing.takeoff.en_route_climb.cd_0)), 'c--')
yline(0.5*double(climb_constraint_oei(landing_constraint_wing_loading, 3.2/100, sizing.landing.first.q, k_func(sizing.landing.first.e), sizing.landing.first.cd_0)), 'm--')
yline(double(climb_constraint_oei(landing_constraint_wing_loading, 2.1/100, sizing.landing.second.q, k_func(sizing.landing.second.e), sizing.landing.second.cd_0)), 'k--')

ylim([0 1]);
xlim([0,8000]);
xlabel('Wing Loading [Nm^{-2}]');
ylabel('Thrust to Weight ratio (errikos had smt dif here, CHECK!!!!)', 'color', 'red');
grid on;
legend(...%'Climb',...
       'Cruise',... 
       'Take-Off',... %    'Service Ceiling',...
       'Absolute Ceiling',... 
       'Turn',...
       'Landing Raymer',...
       'Landing Raymer Trev',... 
       'Landing Roskam',... 
       'Max Velocity',...
       'OEI initial climb',...
       'OEI transition',...
       'OEI second segement',...
       'OEI en-route climb',...
       'OEI first landing',...
	   'OEI second landing');

xline(optimum_w_over_s);

hold off;

improveplot(gcf)

