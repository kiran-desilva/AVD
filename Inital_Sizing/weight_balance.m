load('uc.mat')
load('wing.mat')
load('powerplant.mat')
load('tailplane.mat')
load('sizing.mat')

AR = sizing.AR;                                 % Wing aspect ratio
AR_h = tailplane.horizontal.Ar;                 % Horizontal tailplane aspect ratio
AR_v = tailplane.vertical.Ar;                   % Vertical tailplane aspect ratio
B_h = tailplane.horizontal.b;                   % Horizontal tailplane span (ft)
B_w = wing_design.wing_span_ft;                 % Wing span (ft)
F_w = 0;                  % Fuselage width at horizontal tail intersection (ft)

K_buf = 1.02;             % 1.02 for short ranges; 5.68 for very long ranges
K_door = 1;               % 1.0 if no cargo door; 1.06 for one side cargo door; 1.12 for two sidecargo doors; 1.12 for aft clamshell door; 1.25 for two side and anaft clamshell cargo doors
K_lav = 3.9;              % 1.11 for long range aircraft; 0.31 for short range aircraft; 3.9 for business jets.
K_Lg = 1;                 % 1.12  for fuselage mounted landing gear; 1.0 otherwise
K_mp = 1;                 % 1.126 for kneeling main gear; 1.0 otherwise
K_ng = 1;                 % 1.017 for pylon mounted nacelle; 1.0 otherwise
K_np = 1;                 % 1.15 for kneeling nose-gear; 1.0 otherwise
K_r = 1;                  % 1.133 for reciprocating engines; 1.0 otherwise
K_tp = 1;                 % 0.793 for turboprop; 1.0 otherwise
K_uht = 1;                % 1.143 for all-moving tail; 1.0 otherwise
L = fuselage.structural_length_ft;    % Fuselage structural length (ft)
L_a = 47                 % Electrical routing distance; generators to avionics to cockpit (ft)
L_ec = 94                % Engine controls routing distance; engine to cockpit - total if multiengine (ft)
L_f = 46.66667;          % Total fuselage length (ft)
L_ht = distdim(tailplane.initial.horizontal.Xac,'m','ft');              % Length from wing aerodynamic centre to horizontal tailplane aerodynamic centre (ft)
L_m = uc.main_gear_length_in;        % Main landing gear length (inches)
L_n = uc.nose_gear_length_in;        % Nose landing gear length (inches)
L_vt = tailplane.ac_vt;              % Length from wing aerodynamic centre to vertical tailplane aerodynamic centre (ft)
N_c = 2;                  % Number of crew
N_en = 2;                 % Number of engines
N_f = 6;                  % Number of functions performed by controls; typically 4 − 7 (pitch, yaw, roll, throttle, HLDs, undercarriage)
N_gear = uc.ratio_ngear;  % The ratio of the average total load applied to the main gear during landing to the max landing weight; typically 2.7−3 for commercial aircraft
N_gen = 2;                % Number of generators; typically = Nen
N_l = N_gear * 1.5;       % Ultimate landing gear load factor. 1.5 × Ngear
N_Lt = powerplant.nacelle_length_ft;   % Nacelle length (ft)
N_m = 2;                  % Number of mechanical function performed by controls; typically 0 − 2 (HLDs, stability controls)
N_mss = uc.main_gear_shock_struts;     % Number of main gear shock struts
N_mw = 2;                % Number of main wheels
N_nw = 2;                % Number of nose wheels
N_p = 6;                 % Total number of persons onboard (crew + passengers)
N_seat = 4;              % Number of seats of given type
N_t = 2;                 % Total number of fuel tanks
N_w = powerplant.nacelle_width_ft;    % Nacelle width (ft)
N_z = 2.75*1.5;          % Ultimate load factor; 1.5× limit load factor (2.7-3)
R_kva = 40               % System electrical rating; typically 40 − 60 for transports (kVA)
S_cs = control_surface.total_area_ft2;      % Total control surface area (ft^2)
S_csw = control_surface.wingmounted_area_ft2;      % Area of wing mounted control surfaces (ft^2)
S_e  = control_surface.elevator_area_ft2;          % Elevator area (ft^2)
S_f  = fuselage.wetted_area_ft2;     % Fuselage wetted area (ft^2)
S_ht = tailplane.initial.horizontal.s * 10.7639;   % Horizontal tailplane area (ft^2)
S_n = powerplant.nacelle_wetted_area;    % Nacelle wetted area (ft^2)
S_vt = tailplane.initial.vertical.s * 10.7639;     % Vertical tailplane area (ft^2)
S_w = wing_design.wing_area_ft2;                  % Reference wing area (ft^2)
tc_root = wing_design.tc_ratio;                 % Wing root thickness to chord ratio
tc_rootv = 0.15;         % Vertical tailplane root thickness to chord ratio
% V_i =               % Integral fuel tank volume (gal)
% V_p  =              % Self sealing tank volume (gal)
% W_pr =              % Volume of pressurized sections (ft^3)
% V_s =               % Landing stall speed (ft/s)
V_t = 259.6;             % Total volume of fuel tanks (gal)
W_APU = 0;               % Uninstalled APU weight (lb)
W_c = 0;                 % Maximum cargo weight (lb) [NO CARGO WEIGHT NEEDED]
W_dg = 3100;             % Design gross weight (lb)
W_en = powerplant.engine_weight_lb;      % Engine weight (lb)
W_enc = 2.231*(W_en^0.901);              % Weight of engine and contents (lb); ≈ 2.331KpKtrW0.901en - Kp is 1.4 for engine with propeller or 1.0 otherwise, Ktr is 1.18 for jet with thrust reversers or 1.0 otherwise
W_l = uc.landing_design_gross_weight_lb;      % Landing design gross weight (lb)
W_seat = 40;            % Weight of single seat (lb); ≈ 60 for flight deck seats, 32 for passenger seats, and 11 for troop seats
W_uav = 191.61;         % Uninstalled avionics weight; typically 800 − 1400 (lb)
lambda = wing_design.taper_ratio;                             % Wing taper ratio
cap_lambda = wing_design.wing_quarter_chord_sweep;            % Wing quarter chord sweep
cap_lambda_ht = tailplane.horizontal.sweep_25;        % Horizontal tailplane quarter chord sweep
cap_lambda_vt = tailplane.vertical.sweep_25;          % Vertical tailplane quarter chord sweep

K_ws = 0.75*((1+2*lambda)/(1+lambda)) * B_w * tan(cap_lambda/L);    % 0.75[(1 + 2λ)/(1 + λ)]Bw tan Λ/L
K_y = 0.3*L_ht;          % Aircraft pitching radius of gyration; ≈ 0.3Lht (ft)
K_z = L_vt;              % Aircraft yaw radius of gyration; ≈ Lvt (ft)
%I_y =                   % Pitching moment of inertia; ≈ W_o*Ky^2 (lb ft^2)

% Aircraft wings
W_w = 0.0051*((((W_dg*N_z)^0.557)*(S_w^0.649)*(AR^0.5)*((1+lambda)^0.1)*(S_csw^0.1)) / ((cos(cap_lambda))*(tc_root^0.4)));

% Horizontal tailplane
W_ht = 0.0379*(K_uht*(W_dg^0.639)*(N_z^0.1)*(S_ht^0.75)*(K_y^0.704)*(AR_h^0.166)*((1+S_e/S_ht)^0.1)) / (((1+F_w/B_h)^0.25)*L_ht*cos(cap_lambda_ht));

% Vertical tailplane
W_vt = 0.0026*(((1+H_t/H_v)^0.225)*(W_dg^0.556)*(N_z^0.536)*(S_vt^0.5)*(K_z^0.875)*(AR_v^0.35)) / ((L_vt^0.5)*cos(cap_lambda_vt)*(tc_root_v^0.5));

% Fuselage
W_fus = 0.328*K_door*K_Lg*((W_dg*N_z)^0.5)*(L^0.25)*(S_f^0.302)*((1+K_ws)^0.04)*(L_over_D^0.1);

% Main landing gear 
W_mlg = 0.0106*K_mp*(W_l^0.888)*(N_l^0.25)*(L_m^0.4)*(N_mw^0.321)*(V_s^0.1)/(N_mss^0.5);

% Nose landing gear
W_nlg = 0.032*K_np*(W_l^0.646)*(N_l^0.2)*(L_n^0.5)*(N_nw^0.45);

% Nacelle
W_inl = 0.6724*K_ng*(N_Lt^0.1)*(N_w^0.294)*(N_z^0.119)*(W_enc^0.611)*(N_en^0.984)*(S_n^0.224);

% Engine controls
W_ec = 5*N_en + 0.8*L_ec;

% Engine pneumatic starter
W_es = 49.19*((N_en*W_en/1000)^0.541);

% Fuel system
W_fs = 2.405*(V_t^0.606)*(N_t^0.5)*((1+V_p/V_t)/(1+V_i/V_t));

% Flight controls
W_fc = 145.9*(N_f^0.554)*(S_cs^0.2)*((I_y*10^(-6))^0.07);

% Installed APU - NOT APPLICABLE
%W_APUinst = 2.2*W_APU;

% Instruments
W_instr = 4.509*K_r*K_tp*(N_c^0.541)*N_en*((L_f+B_w)^0.5);

% Hydraulic system
W_hydr = 0.2673*N_f*((L_f+B_w)^0.937);

% Electrical system
W_el = 7.291*(R_kva^0.782)*(L_a^0.346)*(N_gen^0.1);

% Avionics
W_av = 1.73*(W_uav^0.983);

% Furnishings
W_furn = 0.0577*(N_c^0.1)*(W_c^0.393)*(S_f^0.75) + N_seat*W_seat + K_lav*(N_p^1.33) + K_buf*(N_p^1.12);

% Air-conditioning
W_ac = 62.36*(N_p^0.25)*((V_pr*10^(-3))^0.604)*(W_uav^0.1);

% Anti-icing
W_ai = 0.002*W_dg;

% Handling gear (civilian)
W_hg = 3*10^(-4)*W_dg;

% Handling gear military transport - NOT APPLICABLE
% W_hg = 2.4*A_fc

% Total weight with use of fudge factors
Total_weight = W_w*0.78 + (W_ht + W_vt)*0.75 + W_fus*0.85 + (W_mlg + W_nlg)*0.88 + W_inl*0.85 + W_ec + W_es + W_fs + W_fc + W_APUinst + W_instr + W_hydr + W_el + W_av + W_furn + W_ac + W_ai + W_hg

I_y = W*K_y^2/9.81 %not sure what W is here

cg_w = [x z]; % wing cg
cg_t = [x z]; % tail cg
cg_fus = []; % fuselage cg
cg_mlg = []; % main landing gear cg
cg_nlg = []; % '' '' etc
cg_inl = [];
cg_ec = [];
cg_es = [];
cg_fs = [];
cg_fc = [];
cg_instr = [];
cg_hydr= [];
cg_el = [];
cg_av = [];
cg_furn = [];
cg_ac = [];
cg_ai = [];
cg_hg = [];

wandb.x_cg = (W_w*cg_w(1)+(W_ht+W_vt)*cg_t(1)+W_fus*cg_fus(1)+W_mlg*cg_mlg(1)+W_nlg*cg_nlg(1)+W_inl*cg_inl(1)+W_ec*cg_es(1)+W_es*cg_es(1)+W_fs*cg_fs(1)+W_fc*cg_fc(1)+W_w*cg_w(1)+W_instr*cg_instr(1)+W_hydr*cg_hydr(1)+W_el*cg_el(1)+W_av*cg_av(1)+W_furn*cg_furn(1)+W_ac*cg_ac(1)+W_ai*cg_ai(1)+W_hg*cg_hg(1))/Total_weight

wandb.z_cg = (W_w*cg_w(2)+(W_ht+W_vt)*cg_t(2)+W_fus*cg_fus(2)+W_mlg*cg_mlg(2)+W_nlg*cg_nlg(2)+W_inl*cg_inl(2)+W_ec*cg_es(2)+W_es*cg_es(2)+W_fs*cg_fs(2)+W_fc*cg_fc(2)+W_w*cg_w(2)+W_instr*cg_instr(2)+W_hydr*cg_hydr(2)+W_el*cg_el(2)+W_av*cg_av(2)+W_furn*cg_furn(2)+W_ac*cg_ac(2)+W_ai*cg_ai(2)+W_hg*cg_hg(2))/Total_weight

save('wandb','wandb')
