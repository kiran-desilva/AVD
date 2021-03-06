clc
clear
metres_to_ft = 3.28084;

load('uc.mat')
load('wing.mat')
load('powerplant.mat')
load('tailplane.mat')
load('sizing.mat')
load('fuse.mat')
load('locations.mat')
load('control_surface.mat')
load('cg.mat')
load('aero_analysis.mat')

%%%WING FUELTANK CG%%%%%
fueltank.x_cg_from_tip = 1.3; % m
fueltank.z_cg_from_tip = 0.027; % m up 


AR = sizing.AR;                         % Wing aspect ratio
AR_h = tailplane.horizontal.Ar;         % Horizontal tailplane aspect ratio
AR_v = tailplane.vertical.Ar;           % Vertical tailplane aspect ratio
B_h = tailplane.horizontal.b * 3.28084; % Horizontal tailplane span (ft)
B_w = wing.b * 3.28084;                 % Wing span (ft)
F_w = 0;                  % Fuselage width at horizontal tail intersection (ft)
H_t_H_v = 1               % Location of horizontal tailplane on vertical tail. 0.0 for fuselage mounted horizontal tail; 1.0 for T-tail
K_buf = 1.02;             % 1.02 for short ranges; 5.68 for very long ranges
K_door = 1;               % 1.0 if no cargo door; 1.06 for one side cargo door; 1.12 for two sidecargo doors; 1.12 for aft clamshell door; 1.25 for two side and anaft clamshell cargo doors
K_lav = 0.39;              % 1.11 for long range aircraft; 0.31 for short range aircraft; 3.9 for business jets.
K_Lg = 1;                 % 1.12  for fuselage mounted landing gear; 1.0 otherwise
K_mp = 1;                 % 1.126 for kneeling main gear; 1.0 otherwise
K_ng = 1;                 % 1.017 for pylon mounted nacelle; 1.0 otherwise
K_np = 1;                 % 1.15 for kneeling nose-gear; 1.0 otherwise
K_r = 1;                  % 1.133 for reciprocating engines; 1.0 otherwise
K_tp = 1;                 % 0.793 for turboprop; 1.0 otherwise
K_uht = 1;                % 1.143 for all-moving tail; 1.0 otherwise
L = fuse.structural_length/12;    % Fuselage structural length (ft)
L_a = 47                 % Electrical routing distance; generators to avionics to cockpit (ft)
L_ec = 94                % Engine controls routing distance; engine to cockpit - total if multiengine (ft)
L_f = fuse.total_fuselage_length/12; % Total fuselage length (ft)
L_ht = (locations.x_ac_h - locations.x_ac_w)* 3.28084 ;  % Length from wing aerodynamic centre to horizontal tailplane aerodynamic centre (ft)
L_m = uc.main_gear_length_in;        % Main landing gear length (inches)
L_n = uc.nose_gear_length_in;        % Nose landing gear length (inches)
L_vt = (locations.x_ac_v - locations.x_ac_w) * 3.28084;  % Length from wing aerodynamic centre to vertical tailplane aerodynamic centre (ft)
N_c = 2;                  % Number of crew
N_en = 2;                 % Number of engines
N_f = 4;                  % Number of functions performed by controls; typically 4 ??? 7 (pitch, yaw, roll, throttle, HLDs, undercarriage)
N_gear = 2.7%uc.ratio_ngear;  % The ratio of the average total load applied to the main gear during landing to the max landing weight; typically 2.7???3 for commercial aircraft
N_gen = 2;                % Number of generators; typically = Nen
N_l = N_gear * 1.5;       % Ultimate landing gear load factor. 1.5 ?? Ngear
N_Lt = powerplant.nacelle_length_ft;   % Nacelle length (ft)
N_m = 2;                  % Number of mechanical function performed by controls; typically 0 ??? 2 (HLDs, stability controls)
N_mss = uc.main_gear_shock_struts * 2; % Number of main gear shock struts
N_mw = 2;                 % Number of main wheels
N_nw = 2;                 % Number of nose wheels
N_p = 6;                  % Total number of persons onboard (crew + passengers)
N_seat = 4;               % Number of seats of given type
N_t = 3;                  % Total number of fuel tanks
N_w = powerplant.nacelle_width_ft; % Nacelle width (ft)
N_z = 2.7%2.75*1.5;          % Ultimate load factor; 1.5?? limit load factor (2.7-3)
R_kva = 40               % System electrical rating; typically 40 ??? 60 for transports (kVA)
S_cs = control_surface.total_area_ft2;        % Total control surface area (ft^2)
S_csw = control_surface.wingmounted_area_ft2; % Area of wing mounted control surfaces (ft^2)
S_e  = control_surface.elevator_area_ft2;     % Elevator area (ft^2)
S_f  = 50.192 * 10.7639;     % Fuselage wetted area (ft^2) (21.94m^2)
S_ht = tailplane.initial.horizontal.s * 10.7639;   % Horizontal tailplane area (ft^2)
S_n = powerplant.nacelle_wetted_area_ft_sq;    % Nacelle wetted area (ft^2)
S_vt = tailplane.initial.vertical.s * 10.7639; % Vertical tailplane area (ft^2)
S_w = wing.Sref * 10.7639;                     % Reference wing area (ft^2)
tc_root = 0.12;       % Wing root thickness to chord ratio
tc_rootv = 0.15;      % Vertical tailplane root thickness to chord ratio

V_i = 257.43;          % Integral fuel tank volume (gal) -> should be 247.63

V_p  = 0;             % Self sealing tank volume (gal)
V_pr = 509.484697     % Volume of pressurized sections (ft^3) (12.93m^3)
V_s = aero_analysis.wing.v_landing_fts;       % Landing stall speed (ft/s) (from Isobel)

V_t = 257.43;          % Total volume of fuel tanks (gal) -> changed from 259.6

W_APU = 0 ;           % Uninstalled APU weight (lb) (12x13x24)
W_c = 0;              % Maximum cargo weight (lb) [NO CARGO WEIGHT NEEDED]
W_dg = convforce(sizing.W0,'n','lbf');      % Design gross weight (lb)
W_en = powerplant.engine_weight_lb;      % Engine weight (lb)
W_enc = 2.231*(W_en^0.901);              % Weight of engine and contents (lb); ??? 2.331KpKtrW0.901en - Kp is 1.4 for engine with propeller or 1.0 otherwise, Ktr is 1.18 for jet with thrust reversers or 1.0 otherwise
weights.W_f = 1601.1  % Fuel mass (lb)
W_l = W_dg*0.8;         % Landing design gross weight (lb)
W_seat = 40;            % Weight of single seat (lb); ??? 60 for flight deck seats, 32 for passenger seats, and 11 for troop seats
W_uav = 191.6;         % Uninstalled avionics weight; typically 800 ??? 1400 (lb)
lambda = wing.lambda;                % Wing taper ratio
cap_lambda = wing.sweep_25;          % Wing quarter chord sweep
cap_lambda_ht = tailplane.horizontal.sweep_25;        % Horizontal tailplane quarter chord sweep
cap_lambda_vt = tailplane.vertical.sweep_25;          % Vertical tailplane quarter chord sweep

K_ws = 0.75*((1+2*lambda)/(1+lambda)) * B_w * tand(cap_lambda/L);    % 0.75[(1 + 2??)/(1 + ??)]Bw tan ??/L
K_y = 0.3*L_ht;          % Aircraft pitching radius of gyration; ??? 0.3Lht (ft)
K_z = L_vt;              % Aircraft yaw radius of gyration; ??? Lvt (ft)
I_y = (sizing.W0/9.81)*2.20462*K_y^2   % Pitching moment of inertia; ??? W_o*Ky^2 (lb ft^2)
L_over_D = fuse.L_D_f

% Aircraft wings
weights.W_w = 0.0051*((((W_dg*N_z)^0.557)*(S_w^0.649)*(AR^0.5)*((1+lambda)^0.1)*(S_csw^0.1)) / ((cosd(cap_lambda))*(tc_root^0.4)))*0.78;

% Horizontal tailplane
weights.W_ht = 0.75*0.0379*(K_uht*(W_dg^0.639)*(N_z^0.1)*(S_ht^0.75)*(K_y^0.704)*(AR_h^0.166)*((1+S_e/S_ht)^0.1)) / (((1+F_w/B_h)^0.25)*L_ht*cosd(cap_lambda_ht));

% Vertical tailplane
weights.W_vt = 0.75*0.0026*(((1+H_t_H_v)^0.225)*(W_dg^0.556)*(N_z^0.536)*(S_vt^0.5)*(K_z^0.875)*(AR_v^0.35)) / ((L_vt^0.5)*cosd(cap_lambda_vt)*(tc_rootv^0.5));

% Fuselage
weights.W_fus = 0.85*0.328*K_door*K_Lg*((W_dg*N_z)^0.5)*(L^0.25)*(S_f^0.302)*((1+K_ws)^0.04)*(L_over_D^0.1);

% Main landing gear 
weights.W_mlg = 0.88*0.0106*K_mp*(W_l^0.888)*(N_l^0.25)*(L_m^0.4)*(N_mw^0.321)*(V_s^0.1)/(N_mss^0.5);

% Nose landing gear
weights.W_nlg = 0.88*0.032*K_np*(W_l^0.646)*(N_l^0.2)*(L_n^0.5)*(N_nw^0.45);

% Nacelle
weights.W_inl = 0.85*0.6724*K_ng*(N_Lt^0.1)*(N_w^0.294)*(N_z^0.119)*(W_enc^0.611)*(N_en^0.984)*(S_n^0.224);

% Engine controls
weights.W_ec = 5*N_en + 0.8*L_ec;

% Engine pneumatic starter
weights.W_es = 49.19*((N_en*W_en/1000)^0.541);

% Engine weight
weights.W_en = N_en * W_enc;

% Fuel system
weights.W_fs = 2.405*(V_t^0.606)*(N_t^0.5)*((1+V_p/V_t)/(1+V_i/V_t));

% Flight controls
weights.W_fc = (145.9*(N_f^0.554)*(S_cs^0.2)*((I_y*(10^(-6)))^0.07))/(1+(N_m/N_f));

% Installed APU 
weights.W_APUinst = 0%2.2*W_APU;

% Instruments
weights.W_instr = 4.509*K_r*K_tp*(N_c^0.541)*N_en*((L_f+B_w)^0.5);

% Hydraulic system
weights.W_hydr = 0.2673*N_f*((L_f+B_w)^0.937);

% Electrical system
weights.W_el = 7.291*(R_kva^0.782)*(L_a^0.346)*(N_gen^0.1);

% Avionics
weights.W_av = 1.73*(W_uav^0.983);

% Furnishings
weights.W_furn = 0.0577*(N_c^0.1)*(W_c^0.393)*(S_f^0.75) + N_seat*W_seat + K_lav*(N_p^1.33) + K_buf*(N_p^1.12);

% Air-conditioning
weights.W_ac = 62.36*(N_p^0.25)*((V_pr*10^(-3))^0.604)*(W_uav^0.1);

% Anti-icing
weights.W_ai = 0.002*W_dg;

% Passengers + crew
weights.W_p = 207.23 * 6;

weights.W_pay = convmass((15 * 6),"kg","lbm")


% Total weight with use of fudge factors -> MTOW
weights.Total_weight = weights.W_w + (weights.W_ht + weights.W_vt) + weights.W_fus + (weights.W_mlg + weights.W_nlg) + weights.W_inl + weights.W_ec + weights.W_es + weights.W_en + weights.W_f +weights.W_fs + weights.W_fc + weights.W_APUinst + weights.W_instr + weights.W_hydr + weights.W_el + weights.W_av + weights.W_furn + weights.W_ac + weights.W_ai + weights.W_p + weights.W_pay
weights.Total_weight_no_fucking = weights.W_w + (weights.W_ht + weights.W_vt) + weights.W_fus + (weights.W_mlg + weights.W_nlg) + weights.W_inl + weights.W_ec + weights.W_es + weights.W_en + weights.W_f +weights.W_fs + weights.W_fc + weights.W_APUinst + weights.W_instr + weights.W_hydr + weights.W_el + weights.W_av + weights.W_furn + weights.W_ac + weights.W_ai + weights.W_p + weights.W_pay


fuel_fraction_to_fuel_weight = @(wf_fuel) ((wf_fuel*sizing.W0)-(sizing.W0-sizing.Wf))*0.2248089431; % conversion from newton to lbf

payload_factor_func = @(payload_factor) (2+(4*payload_factor))/6;
weights.Total_weight_func = @(wf_fuel,payload_factor) weights.W_w*0.78 + (weights.W_ht + weights.W_vt)*0.75 + weights.W_fus*0.85 + (weights.W_mlg + weights.W_nlg)*0.88 + weights.W_inl*0.85 + weights.W_ec + weights.W_es + weights.W_en + fuel_fraction_to_fuel_weight(wf_fuel) +weights.W_fs + weights.W_fc + weights.W_APUinst + weights.W_instr + weights.W_hydr + weights.W_el + weights.W_av + weights.W_furn + weights.W_ac + weights.W_ai + (payload_factor_func(payload_factor)*(weights.W_p + weights.W_pay))
weights.Total_weight_func = @(wf_fuel,payload_factor) weights.W_w + (weights.W_ht + weights.W_vt) + weights.W_fus + (weights.W_mlg + weights.W_nlg) + weights.W_inl + weights.W_ec + weights.W_es + weights.W_en + fuel_fraction_to_fuel_weight(wf_fuel) +weights.W_fs + weights.W_fc + weights.W_APUinst + weights.W_instr + weights.W_hydr + weights.W_el + weights.W_av + weights.W_furn + weights.W_ac + weights.W_ai + (payload_factor_func(payload_factor)*(weights.W_p + weights.W_pay))
weights.Total_weight = weights.Total_weight_func(1,1);

save('weights','weights')

% I_y = W*K_y^2/9.81 %not sure what W is here

% x & z cg coords
cg_w = [cg.x_wing -33/12]; % wing cg
cg_ht = [cg.x_htail cg.z_htail]; % htail cg
cg_vt = [cg.x_vtail cg.z_vtail]
cg_fus = [cg.x_fuse 0]; % fuselage cg
cg_mlg = [cg.x_mlg cg.z_mlg]; % main landing gear cg
cg_nlg = [cg.x_nlg cg.z_nlg]; % '' '' etc
cg_inl = [cg.x_inl cg.z_inl];
cg_ec = [cg.x_inl cg.z_inl];
cg_es = [cg.x_inl cg.z_inl];
cg_f = [metres_to_ft*(locations.x_ac_w - wing.Xac_from_tip + fueltank.x_cg_from_tip); -33/12]; % TODO: Might need to add underbelly in
cg_fs = [metres_to_ft*(locations.x_ac_w - wing.Xac_from_tip + fueltank.x_cg_from_tip); -33/12]; % TODO: Might need to add underbelly in
cg_fc = [4.5 -0.66];
cg_APUinst = [cg.x_APUinst 0]
cg_instr = [20 -12]/12;
cg_hydr= [22 -20]/12;
% cg_el = [22 -20]/12;
cg_el = [60 -20]/12;
cg_en = [metres_to_ft*locations.x_rear_bulkhead metres_to_ft*0];
cg_av = [20 -12]/12;
cg_furn = [20 -(30/12)];
%locating passenger cg in center of wing root for easy emergency access
passenger_cg_func = @(xac) (1.465 + cg.x_wing_tip_from_ac(xac))*metres_to_ft;
cg_pass = [passenger_cg_func(locations.x_ac_w) -15/12];
cg_pilot = [convlength(1.779,'m','ft') -15/12];
cg_ac = [22 -20/12];
cg_ai = [22 -20/12];
cg_pay = [107/12,-15/12];
% cg_pay = [333/12,-15/12];
%
%returns wing_fuel_tank_cg in feet
fuel_tank_cg = @(wing_ac) metres_to_ft*(wing_ac - wing.Xac_from_tip + fueltank.x_cg_from_tip);


wandb.x_cg = (weights.W_en*cg_en(1) + weights.W_w*cg_w(1)+weights.W_ht*cg_ht(1)+weights.W_vt*cg_vt(1)+weights.W_fus*cg_fus(1)+weights.W_mlg*cg_mlg(1)+weights.W_nlg*cg_nlg(1)+weights.W_inl*cg_inl(1)+weights.W_ec*cg_es(1)+weights.W_es*cg_es(1)+weights.W_f*cg_f(1)+weights.W_fs*cg_fs(1)+weights.W_fc*cg_fc(1)+weights.W_instr*cg_instr(1)+weights.W_hydr*cg_hydr(1)+weights.W_el*cg_el(1)+weights.W_av*cg_av(1)+weights.W_furn*cg_furn(1)+weights.W_ac*cg_ac(1)+weights.W_ai*cg_ai(1)+(weights.W_p*cg_pass(1))+(weights.W_pay*cg_pay(1)))/weights.Total_weight


%	wing_ac	-	aerodynamic centre of the wing, measured from the nose of the plane, in metres
%   x_gear - in meters from tip
%	wf_fuel -	fuel weight fraction, defined as current_weight/max_takeoff_weight 
%	payload_factor -	payload fraction, varies from 0 to 1, where 0 means no passengers and no baggage (but 2 pilots)
%						and 1 means 4 passengers, 4 bags (and 2 pilots)
%   returns [ cg location in FEET]

wandb.x_cg_function = @(wing_ac,x_gear_main,x_gear_nose,wf_fuel,payload_factor) (weights.W_en*cg_en(1) + weights.W_w*cg.x_wing_from_ac_func(wing_ac)+weights.W_ht*cg_ht(1)+weights.W_vt*cg_vt(1)+weights.W_fus*cg_fus(1)+(weights.W_mlg*x_gear_main*metres_to_ft)+(weights.W_nlg*x_gear_nose*metres_to_ft)+weights.W_inl*cg_inl(1)+weights.W_ec*cg_es(1)+weights.W_es*cg_es(1)+(fuel_fraction_to_fuel_weight(wf_fuel)*fuel_tank_cg(wing_ac))+(weights.W_fs*fuel_tank_cg(wing_ac))+weights.W_fc*cg_fc(1)+weights.W_instr*cg_instr(1)+weights.W_hydr*cg_hydr(1)+weights.W_el*cg_el(1)+weights.W_av*cg_av(1)+weights.W_furn*cg_furn(1)+weights.W_ac*cg_ac(1)+weights.W_ai*cg_ai(1)+( (payload_factor*(4/6)*weights.W_p*passenger_cg_func(wing_ac)) + ((2/6)*weights.W_p*cg_pilot(1)) )+(payload_factor_func(payload_factor)*weights.W_pay*cg_pay(1)))./weights.Total_weight_func(wf_fuel,payload_factor);

wandb.x_cg = wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,1,1)

wandb.z_cg = (weights.W_en*cg_en(2) + weights.W_w*cg_w(2)+weights.W_ht*cg_ht(2)+weights.W_vt*cg_vt(2)+weights.W_fus*cg_fus(2)+weights.W_mlg*cg_mlg(2)+weights.W_nlg*cg_nlg(2)+weights.W_inl*cg_inl(2)+weights.W_ec*cg_es(2)+weights.W_es*cg_es(2)+weights.W_f*cg_f(2)+weights.W_fs*cg_fs(2)+weights.W_fc*cg_fc(2)+weights.W_instr*cg_instr(2)+weights.W_hydr*cg_hydr(2)+weights.W_el*cg_el(2)+weights.W_av*cg_av(2)+weights.W_furn*cg_furn(2)+weights.W_ac*cg_ac(2)+weights.W_ai*cg_ai(2)+weights.W_p+cg_pass(2)+weights.W_pay*cg_pay(2))/weights.Total_weight

fuel_used_at_weight_fraction = @(total_initial_weight, weight_fraction) (total_initial_weight*(1 - weight_fraction));
x_cg_at_weight_fraction = @(total_initial_weight, weight_fraction) (wandb.x_cg*total_initial_weight - fuel_used_at_weight_fraction(total_initial_weight, weight_fraction)*cg_f(1))/(total_initial_weight*weight_fraction);  

weights.total_weight_no_fuel = weights.Total_weight - weights.W_f
wandb.x_cg_drained = (wandb.x_cg*weights.Total_weight - weights.W_f*cg_f(1))/weights.total_weight_no_fuel
wandb.z_cg_drained = (wandb.z_cg*weights.Total_weight - weights.W_f*cg_f(2))/weights.total_weight_no_fuel

% %% Cg envelope calc
% % with passengers
% with_passengers = @(w_frac) x_cg_at_weight_fraction(weights.Total_weight, w_frac);
% start_cg = with_passengers(1 - weights.W_f/weights.Total_weight);
% start_weight = weights.Total_weight - weights.W_f; 
% cg_envelope = [start_cg, wandb.x_cg, with_passengers(sizing.fraction.before_cruise), with_passengers(sizing.fraction.end_cruise_1), start_cg;...
% 			   start_weight, weights.Total_weight, sizing.fraction.before_cruise*weights.Total_weight, sizing.fraction.end_cruise_1*weights.Total_weight, start_weight];

% % without passengers and baggage
% x_cg_at_weight_fraction = @(total_initial_weight, weight_fraction) (wandb.x_cg*total_initial_weight - fuel_used_at_weight_fraction(total_initial_weight, weight_fraction)*cg_f(1) - weights.W_p*4/6*)/(total_initial_weight*weight_fraction);  
% no_passengers = @(w_frac) x_cg_at_weight_fraction(weights.Total_weight - weights.W_p*4/6, w_frac);
% start_cg = no_passengers(1 - weights.W_f/weights.Total_weight);
% start_weight = weights.Total_weight - weights.W_p*4/6 - weights.W_f; 
% cg_envelope = [start_cg, no_passengers(1), no_passengers(sizing.fraction.before_cruise), no_passengers(sizing.fraction.end_cruise_1), start_cg;...
% 			   start_weight, weights.Total_weight, sizing.fraction.before_cruise*weights.Total_weight, sizing.fraction.end_cruise_1*weights.Total_weight, start_weight];

% plot(cg_envelope(1, :), cg_envelope(2, :));

fraction = [sizing.fraction.before_take_off,sizing.fraction.before_cruise,sizing.fraction.end_cruise_1,sizing.fraction.before_alternate_cruise,sizing.fraction.before_loiter,sizing.fraction.end];
wew0 = (1-(sizing.Wf/sizing.W0)); 

empty_cg = wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,wew0,0);
empty_weight = convforce(wew0 * sizing.W0,'n','lbf') - (weights.W_pay + weights.W_p);

figure
hold on
plot(0.3048*[empty_cg;wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,fraction,1)';empty_cg],4.4482*[empty_weight;weights.Total_weight_func(fraction,1)';empty_weight],'x-')
plot(0.3048*[empty_cg;wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,fraction,0)';empty_cg],4.4482*[empty_weight;weights.Total_weight_func(fraction,0)';empty_weight],'x-')

xlabel('X_{cg} [m]')
ylabel('Weight [N]')
grid on
legend('Full Payload','No Payload')
save('wandb','wandb')

figure
hold on

scatter(cg_w(1),cg_w(2),weights.W_w,'LineWidth',1.5)
scatter(cg_ht(1),cg_ht(2),weights.W_ht,'LineWidth',1.5)
scatter(cg_vt(1),cg_vt(2),weights.W_vt,'LineWidth',1.5)
scatter(cg_fus(1),cg_fus(2),weights.W_fus,'LineWidth',1.5)
scatter(cg_mlg(1),cg_mlg(2),weights.W_mlg,'LineWidth',1.5)
scatter(cg_nlg(1),cg_nlg(2),weights.W_nlg,'LineWidth',1.5)
scatter(cg_inl(1),cg_inl(2),weights.W_inl,'LineWidth',1.5)
scatter(cg_ec(1),cg_ec(2),weights.W_ec,'LineWidth',1.5)
scatter(cg_es(1),cg_es(2),weights.W_es,'LineWidth',1.5)
scatter(cg_f(1),cg_f(2),weights.W_f,'LineWidth',1.5)
scatter(cg_fs(1),cg_fs(2),weights.W_fs,'LineWidth',1.5)
scatter(cg_fc(1),cg_fc(2),weights.W_fc,'LineWidth',1.5)
%scatter(cg_APUinst(1),weights.)
scatter(cg_instr(1),cg_instr(2),weights.W_instr,'LineWidth',1.5)
scatter(cg_hydr(1),cg_hydr(2),weights.W_hydr,'LineWidth',1.5)
scatter(cg_el(1),cg_el(2),weights.W_el,'LineWidth',1.5)
scatter(cg_en(1),cg_en(2),weights.W_en,'LineWidth',1.5)
scatter(cg_av(1),cg_av(2),weights.W_av,'LineWidth',1.5)
scatter(cg_furn(1),cg_furn(2),weights.W_furn,'LineWidth',1.5)
scatter(cg_pass(1),cg_pass(2),weights.W_p * (4/6),'LineWidth',1.5)
scatter(cg_pilot(1),cg_pilot(2),weights.W_p * (2/6),'LineWidth',1.5)
scatter(cg_ac(1),cg_ac(2),weights.W_ac,'LineWidth',1.5)
scatter(cg_ai(1),cg_ai(2),weights.W_ai,'LineWidth',1.5)
scatter(cg_pay(1),cg_pay(2),weights.W_pay,'LineWidth',1.5)
plot(wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,sizing.fraction.before_take_off,1),wandb.z_cg,'x','color','black')

% xline(wandb.x_cg_function(locations.x_ac_w,cg.x_mlg/metres_to_ft,cg.x_nlg/metres_to_ft,sizing.fraction.before_take_off,1),'--','linewidth',1.2)
% xline(metres_to_ft*locations.x_wing)
% xline(metres_to_ft*(locations.x_wing + wing.Croot))

legend('Wing','Horizontal Stabilizer','Vertical Stabilizer','Fuselage','Main Landing Gear','Nose Landing Gear','Nacelle','Engine Controls','Engine Starter','Fuel','Fuel System','Flight Controls','Instruments','Hydraulics','Electrical System','Engine','Avionics','Furnishings','Passengers','Pilots','Air Conditioning','Anti-icing','Luggage','CG');
xlabel('X coordinate from nose [ft]')
ylabel('Z coordinate from fuselage centre [ft]')

grid on


improvePlot(gcf)