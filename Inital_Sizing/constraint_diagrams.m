clear
clc
clf

u = symunit; % used for units
subindex = @(A, idx) A(idx); %% function for anonymous indexing

load("parameters.mat");
load("sizing.mat");


parameters.absolute_ceiling = 45000; % ft


%% clean
sizing.cl_max = 1.6;

sizing.k = 1/(sizing.AR*pi*sizing.e); % This eq. is used in Gudmundsson p. 64 -> idk if this is valid 


%% cd 0 taken from errikos' slides
sizing.undercarriage_e = -0.05;
sizing.undercarriage_cd0 = 0.02;

sizing.t_o_flaps_e = -0.05;
sizing.t_o_flaps_cd0 = 0.02;

sizing.landing_flaps_e = -0.1;
sizing.landing_flaps_cd0 = 0.07;

sizing.n = 1/cosd(40); % max load factor as defined by far


q = @(V_inf, rho) 0.5*rho*V_inf^2;


%% takeoff 
sizing.takeoff.rho = parameters.rho_0;

sizing.takeoff.cl_max = 1.9; %  roskam
sizing.takeoff.e = sizing.e + sizing.undercarriage_e + sizing.t_o_flaps_e;
sizing.takeoff.cd0 = sizing.cd0 + sizing.undercarriage_cd0 + sizing.t_o_flaps_cd0;

%% landing
sizing.landing.rho = sizing.takeoff.rho;
sizing.landing.cl_max = 2.4; % roskam:
sizing.landing.obstacle_height = 183; %m
sizing.landing.kr = 0.66; % with thrust reversers

%% landing constraint
%% page 111 in roskam
%% 
% TODO: Add constraint with thrust reversals?


landing_constraint_wing_loading_roskam = (1/parameters.target_landing_weight)*(parameters.runway_length*0.6)*sizing.landing.rho*sizing.landing.cl_max/((0.6*1.3)^2*double(separateUnits(unitConvert(u.ft/u.kts^2, u.m/(u.m/u.s)^2))));
landing_constraint_wing_loading = (1/parameters.target_landing_weight)*((parameters.runway_length*0.6) - sizing.landing.obstacle_height)*sizing.landing.cl_max/0.51;
landing_constraint_wing_loading_trev = (1/parameters.target_landing_weight)*((parameters.runway_length*0.6) - sizing.landing.obstacle_height)*sizing.landing.cl_max/(0.51*0.66);

takeoff_constraint_wing_loading = landing_constraint_wing_loading;

% sizing.takeoff.v_stall = sqrt((sizing.maxTakeoffWeight/sizing.sref) * 2/(sizing.takeoff.rho * sizing.takeoff.cl_max));
sizing.landing.v_stall = sqrt((2*landing_constraint_wing_loading)/(sizing.landing.rho*sizing.landing.cl_max));
sizing.takeoff.v_stall = sqrt((2*takeoff_constraint_wing_loading)/(sizing.takeoff.rho*sizing.takeoff.cl_max));
sizing.v_stall = sizing.takeoff.v_stall*sqrt(sizing.takeoff.cl_max/sizing.cl_max);
sizing.takeoff.v_inf = 1.3*sizing.takeoff.v_stall; % we should check if this is reasonable... Errikos did it in his excel but idk -> from FAR25
sizing.takeoff.q = q(sizing.takeoff.v_inf, sizing.takeoff.rho);


sizing.takeoff.initial_climb.q = q(1.2*sizing.takeoff.v_stall, sizing.takeoff.rho);
sizing.takeoff.initial_climb.e = sizing.e + sizing.t_o_flaps_e;
sizing.takeoff.initial_climb.cd0 = sizing.cd0 + sizing.t_o_flaps_cd0;

sizing.takeoff.transition.q = sizing.takeoff.initial_climb.q;
sizing.takeoff.transition.e = sizing.e + sizing.t_o_flaps_e + sizing.undercarriage_e;
sizing.takeoff.transition.cd0 = sizing.cd0 + sizing.t_o_flaps_cd0 + sizing.undercarriage_cd0;


sizing.takeoff.second_segment.q = sizing.takeoff.initial_climb.q;
sizing.takeoff.second_segment.e = sizing.e + sizing.t_o_flaps_e;
sizing.takeoff.second_segment.cd0 = sizing.cd0 + sizing.t_o_flaps_cd0;


sizing.takeoff.en_route_climb.q = q(1.25*sizing.v_stall, sizing.takeoff.rho);
sizing.takeoff.en_route_climb.e = sizing.e;
sizing.takeoff.en_route_climb.cd0 = sizing.cd0;

sizing.landing.first.q = q(1.3*sizing.landing.v_stall, sizing.landing.rho);
sizing.landing.first.e = sizing.e + sizing.landing_flaps_e + sizing.undercarriage_e;
sizing.landing.first.cd0 = sizing.cd0 + sizing.landing_flaps_cd0 + sizing.undercarriage_cd0;


sizing.landing.second.q = q(1.5*1.1*sizing.landing.v_stall, sizing.landing.rho);
sizing.landing.second.e = sizing.e + sizing.t_o_flaps_e + sizing.undercarriage_e;
sizing.landing.second.cd0 = sizing.cd0 + sizing.t_o_flaps_cd0 + sizing.undercarriage_cd0;


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
[~, sizing.loiter.a, ~, sizing.loiter.rho] = atmosisa(distdim(5000, 'ft', 'm'));
%sizing.loiter.v_inf = sizing.cruise.v_inf*(1/3)^(0.25); %TODO: this would be the exact equation for VminD, we need the loiter weight for that tho sqrt(1.225/sizing.loiter.rho)*pow((4*sizing.k*sizing.loiter.w^2)/(pi*sizing.cd0*sizing.AR), 0.25);
sizing.loiter.v_inf = 0.28*sizing.loiter.a;
sizing.loiter.q = q(sizing.loiter.v_inf, sizing.loiter.rho);

%%diversion

[~, ~, ~, sizing.diversion.rho] = atmosisa(distdim(40000, 'ft', 'm'));
sizing.diversion.v_inf = sizing.cruise.v_inf; %TODO: this would be the exact equation for VminD, we need the loiter weight for that tho sqrt(1.225/sizing.loiter.rho)*pow((4*sizing.k*sizing.loiter.w^2)/(pi*sizing.cd0*sizing.AR), 0.25);
sizing.diversion.q = q(sizing.diversion.v_inf, sizing.diversion.rho);


%% FAR25 - https://www.engineerstoolkit.com/Airworthiness%20Standards%20%20FAA%20FAR%20Part%2025.pdf
%% take off speed = 1.1* stall speed


syms wing_loading thrust_to_weight V_inf rho turn_height_m  climb_grad k cd0 q_var % use syms for constraints
syms alpha beta % alpha -> scales weight to initial weight, beta -> scales thrust to sea level thrust

%% Computing climb constraint

% climb_constraint(wing_loading) = sizing.climb.dh_dt/sizing.climb.v_inf + sizing.climb.q/wing_loading*sizing.cd0 + sizing.k/sizing.climb.q*wing_loading;
climb_constraint(wing_loading, climb_grad) = climb_grad + sizing.climb.q/wing_loading*sizing.cd0 + sizing.k/sizing.climb.q*wing_loading;

%% Cruise constraint

cruise_constraint(wing_loading, alpha, beta) = alpha/beta*(sizing.cruise.q*sizing.cd0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/sizing.cruise.q);

loiter_constraint(wing_loading, alpha, beta) = alpha/beta*(sizing.loiter.q*sizing.cd0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/sizing.loiter.q);

diversion_constraint(wing_loading, alpha, beta) = alpha/beta*(sizing.diversion.q*sizing.cd0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/sizing.diversion.q);

%% Take-off distance constraint

takeoff_bfl_constraint(wing_loading) = ( wing_loading*(0.297-0.019*parameters.engine_number) ) / ( parameters.runway_length * (sizing.takeoff.rho / parameters.rho_0) * sizing.takeoff.cl_max) ; 


take_off_constraint(wing_loading) = double(separateUnits(unitConvert(37.5*(u.ft^3/u.lbf), u.m^3/u.N)))*wing_loading/(sizing.takeoff.cl_max*parameters.runway_length);
% take_off_constraint(wing_loading) = 37.5*wing_loading/(sizing.takeoff.cl_max*parameters.runway_length);

% max velocity constraint
max_velocity_q = q(sizing.cruise.a*parameters.max_mach, sizing.cruise.rho);
max_velocity_constraint(wing_loading, alpha, beta) =alpha/beta*(max_velocity_q*sizing.cd0/(alpha*wing_loading) + alpha*sizing.k*wing_loading/max_velocity_q);
%  max_velocity_q/wing_loading*sizing.cd0 + sizing.k/max_velocity_q*wing_loading;

%% absolute ceiling constraint

[~, ~, ~, parameters.absolute_ceiling_rho] = atmosisa(distdim(parameters.absolute_ceiling,'ft','km'));

% sizing.absolute_climb_velocity_at_ceiling = 0; % at absolute ceiling, climb rate should be zero

% thing = sqrt(sizing.k/(3*sizing.cd0));
% absolute_ceiling_constraint(wing_loading,alpha,beta) = sizing.absolute_climb_velocity_at_ceiling / sqrt(wing_loading*2*thing/sizing.service_ceiling_rho) + 4*sqrt(sizing.k*sizing.cd0/3);

Vx(wing_loading) = sqrt( (2/parameters.absolute_ceiling_rho) * (wing_loading) * sqrt(sizing.k/sizing.cd0) * 1);

% Vimd(wing_loading) = sqrt((2/1.225) * wing_loading) * ((sizing.k / (pi*sizing.AR * sizing.cd0))^(1/4))

% Vx = EquivalentToTrue(Vimd,parameters.absolute_ceiling)

absolute_ceiling_constraint(wing_loading,alpha,beta) = simplify((alpha/beta) * ( ( (0.5*parameters.absolute_ceiling_rho*(Vx^2)*sizing.cd0)/(alpha*wing_loading) ) + ( ((1^2)*(alpha*wing_loading))/(0.5*parameters.absolute_ceiling_rho*(Vx^2)*pi*sizing.AR*sizing.e) ) ));


%% turn constraint

turn_constraint = @(wing_loading, turn_height_m, V_inf) q(V_inf, atmos(turn_height_m, 4))*(sizing.cd0/wing_loading + sizing.k*wing_loading*(sizing.n/q(V_inf, atmos(turn_height_m, 4)))^2);


%% stall constraint

%TODO: stall_constraint
stall_constraint = @(cl_max,q_stall) q_stall*cl_max;

%% Climb OEI
climb_constraint_oei(wing_loading, climb_grad, q_var, k, cd0) = 2*(climb_grad + q_var/wing_loading*cd0 + k/q_var*wing_loading);


hold on

linewidth = 2;

weight_loading_interval = [1, 12000];

fplot(@(wing_loading) cruise_constraint(wing_loading, sizing.fraction.before_cruise, 0.25), weight_loading_interval,'LineWidth',linewidth);
fplot(@(wing_loading) loiter_constraint(wing_loading, sizing.fraction.before_loiter, 0.7), weight_loading_interval,'LineWidth',linewidth);
fplot(@(wing_loading) diversion_constraint(wing_loading, sizing.fraction.before_alternate_cruise, 0.25), weight_loading_interval,'LineWidth',linewidth);

fplot(take_off_constraint, weight_loading_interval,'LineWidth',linewidth);

fplot(@(wing_loading) absolute_ceiling_constraint(wing_loading, sizing.fraction.before_cruise, 0.20),'LineWidth',linewidth)
fplot(@(wing_loading) sizing.fraction.before_cruise*turn_constraint(wing_loading, distdim(40000, 'ft', 'm'), sizing.cruise.v_inf), weight_loading_interval,'LineWidth',linewidth);
xline(landing_constraint_wing_loading,'color','red','LineWidth',linewidth);
xline(landing_constraint_wing_loading_trev,'LineWidth',linewidth);

xline(landing_constraint_wing_loading_roskam, 'color', 'magenta','LineWidth',linewidth);

fplot(@(wing_loading) max_velocity_constraint(wing_loading, sizing.fraction.before_cruise, 0.25), weight_loading_interval,'LineWidth',linewidth);

k_func = @(e) 1/(pi*sizing.AR*e);

% pls tell me i did the dumb dumb
% optimum_w_over_s = 0.5*sizing.cruise.rho*((sizing.cruise.v_inf/(3^0.25))^2)*sqrt(pi*sizing.AR*sizing.cd0/sizing.k)
optimum_w_over_s = 0.5*sizing.cruise.rho*((109.8/(3^0.25))^2)*sqrt(pi*sizing.AR*sizing.cd0/sizing.k);
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

yline(sizing.fraction.before_take_off*double(climb_constraint_oei(takeoff_constraint_wing_loading, 1.2/100, sizing.takeoff.initial_climb.q, k_func(sizing.takeoff.initial_climb.e), sizing.takeoff.initial_climb.cd0)), 'r--','LineWidth',linewidth)
yline(sizing.fraction.before_take_off*double(climb_constraint_oei(takeoff_constraint_wing_loading, 0, sizing.takeoff.transition.q, k_func(sizing.takeoff.transition.e), sizing.takeoff.transition.cd0)), 'b--','LineWidth',linewidth)
yline(sizing.fraction.before_take_off*double(climb_constraint_oei(takeoff_constraint_wing_loading, 2.4/100, sizing.takeoff.second_segment.q, k_func(sizing.takeoff.second_segment.e), sizing.takeoff.second_segment.cd0)), 'g--','LineWidth',linewidth)
yline(sizing.fraction.before_take_off*double(climb_constraint_oei(takeoff_constraint_wing_loading, 1.2/100, sizing.takeoff.en_route_climb.q, k_func(sizing.takeoff.en_route_climb.e), sizing.takeoff.en_route_climb.cd0)), 'c--','LineWidth',linewidth)
yline(sizing.fraction.end*0.5*double(climb_constraint_oei(landing_constraint_wing_loading, 3.2/100, sizing.landing.first.q, k_func(sizing.landing.first.e), sizing.landing.first.cd0)), 'm--','LineWidth',linewidth)
yline(sizing.fraction.end*double(climb_constraint_oei(landing_constraint_wing_loading, 2.1/100, sizing.landing.second.q, k_func(sizing.landing.second.e), sizing.landing.second.cd0)), 'k--','LineWidth',linewidth)

design_w_s = 2973;
design_t_w = .3845;
design_S_ref = 36.7e3 /design_w_s;
t0 = design_t_w * 36.7e3;

ylim([0 1]);
xlim([0,11000]);
xlabel('Wing Loading [Nm^{-2}]');
ylabel('T_0/W_0');
grid on;


plot(design_w_s, design_t_w, 'x', 'MarkerSize', 20, 'color','black');
% plot(design_w_s, max_velocity_constraint(design_w_s, sizing.fraction.before_cruise, 0.25), 'o', 'MarkerSize', 20, 'color','black');

legend('Cruise',...
       'Loiter',...
       'Diversion',...
       'Take-Off',...
       'Absolute Ceiling',... 
       'Turn',...
       'Landing Raymer w/o Trev',...
       'Landing Raymer Trev',...
       'Landing Roskam w/o Trev',... 
       'Max Velocity',...
       'OEI initial climb',...
       'OEI transition',...
       'OEI second segement',...
       'OEI en-route climb',...
       'AEO go-around climb',...
	'OEI go-around climb',...
       'Design Point');

improvePlot(gcf)
hold off

design.t_w = design_t_w;
design.w_s = design_w_s;
design.sref = design_S_ref;
design.t0 = t0;

save('design','design')
