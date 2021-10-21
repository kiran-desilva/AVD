
open("parameters.mat");

%% shit we want
%% l/d max - l_dmax
%% specific fuel consumption - sfc
%% gross weight - w0
%% operating empty weight - we
%% fuel weight - wf

sizing.ld_max = 10;
sizing.ld_cruise = 0.866*sizing.ld_max;
sizing.sfc = 10;
sizing.maxTakeoffWeight = 1000;
sizing.cd_min = nan;

%% shit from errikos' excel
sizing.cl_max = nan;
sizing.AR = nan;
sizing.cd_0 = nan;
sizing.e = nan;
sizing.k = nan;  % idk what k is, Gudmundsson p. 59 says lift-induced drag constant
sizing.field_length = nan;

sizing.w_landing_w_total_max = nan;

sizing.takeoff.rho = nan;
sizing.cruise.rho = nan;
sizing.loiter.rho = nan;
sizing.climb.rho = nan;


sizing.takeoff.v_inf = nan;
sizing.cruise.v_inf= nan;
sizing.loiter.v_inf= nan;
sizing.climb.v_inf= nan;

sizing.climb.dh_dt = nan;

sizing.v_stall_max = nan;

sizing.w_total_S_max_landing = nan;
sizing.w_total_S_max_stall = nan;


%% shit from errikos' excel
sizing.cl_max = nan;
sizing.AR = nan;
sizing.cd_0 = nan;
sizing.e = nan;
sizing.field_length = nan;

sizing.w_landing_w_total_max = nan;
sizing.rho_landing = nan;
sizing.v_stall_max = nan;

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


syms wing_loading thrust_to_weight V_inf rho % use syms for constraints
q(V_inf, rho) = 0.5*rho*V_inf^2;

%% Computing climb constraint
climb = sizing.climb;

q_climb = q(climb.v_inf, climb.rho);
thrust_to_weight(wing_loading) = climb.dh_dt/climb.v_inf + q_climb/wing_loading*sizing.cd_min + sizing.k/q*wing_loading;

%% Cruise constraint
cruise = sizing.cruise;

q_cruise = q(cruise.v_inf, cruise.rho);
thrust_to_weight(wing_loading) = q_cruise*sizing.cd_min/wing_loading + sizing.k*wing_loading/q_cruise;

%% Take-off distance constraint
takeoff = sizing.takeoff;

q_takeoff = q(takeoff.v_inf, takeoff.rho);
thrust_to_weight(wing_loading) = nan;

%% Operating ceiling constraint
sizing.ceiling_rho = nan;
sizing.climb_velocity_at_ceiling = nan; % what is the rate of climb at the ceiling alt we want? Look at far mby...

thing = sqrt(sizing.k/(3*sizing.cd_min));
thrust_to_weight(wing_loading) = sizing.climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.ceiling_rho) + 4*sqrt(sizing.k*sizing.cd_min/3);

