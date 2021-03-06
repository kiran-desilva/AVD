% Wing Lift Analysis
%ieo18
%10th Nov 21 - 29th Nov 21
%AVD
%clear
clc

%% Load files

load('sizing.mat')
load('design.mat')
load('wing.mat')
load('fuse.mat')
load('aero_analysis.mat')
load('tailplane.mat')
%% Wing Lift Estimation

%% Inputs from other scripts
d=(fuse.d_f)/39.3700787; %diameter of the fuselage [m] from Camille 26/11
AR=tailplane.horizontal.Ar;
S_exp=1; %[m^2]
S_ref=tailplane.horizontal.s; %[m^2]
%h=0.3; %winglet height (m)
b=tailplane.horizontal.b; %wing span (m)
max_sweep_t=tailplane.horizontal.sweepmaxTC*pi/180; %sweep of wing at the chord location where the airfoil is thickest [rad]
sweep_LE=tailplane.horizontal.sweepLE*pi/180; %[rad]
cl_max_airfoil=1.6; %maximum lift of airfoil
lambda_quarter=tailplane.horizontal.sweep_25*pi/180; %quarter-chord sweep
aero_analysis.tail.t_c=tailplane.horizontal.tc;
aero_analysis.tail.te_angle=wing.sweepTE; %degrees

%% Mach number definition, wing
%aero_analysis.wing.stall_landing=sqrt(2*(sizing.W0-sizing.Wf)/(1.225*wing.Sref*2.3));
%NOTE: need to change the weights used for the stall speed!!!!
aero_analysis.tail.v_landing_ms=1.3*sqrt(2*(sizing.W0-sizing.Wf)/(1.225*wing.Sref*2.3)); %units:m/s, safety margin of 1.3
aero_analysis.tail.v_landing_fts=3.2808*aero_analysis.tail.v_landing_ms; %units: ft/s
aero_analysis.tail.v_takeoff_ms=1.2*sqrt(2*(sizing.W0)/(1.225*wing.Sref*2.3)); %units:m/s, safety margin of 1.2
aero_analysis.tail.v_takeoff_fts=3.2808*aero_analysis.tail.v_takeoff_ms; %units: ft/s 

aero_analysis.tail.dyn_visc=[0.0000143226,0.0000143226,0.0000181206,0.0000181206,0.0000174141];
taper=tailplane.horizontal.Ctip/tailplane.horizontal.Croot;
aero_analysis.tail.MAC=tailplane.horizontal.Croot*2/3*((1+taper+taper^2)/1+taper);
aero_analysis.tail.rho=[0.302,0.237,1.225,1.225,1.0552];
aero_analysis.tail.air_velc=[294.9,294.9,340.3,340.3,334.4];
aero_analysis.tail.Mach=[0.75,0.78,aero_analysis.tail.v_takeoff_ms/aero_analysis.tail.air_velc(3),aero_analysis.tail.v_landing_ms/aero_analysis.tail.air_velc(4), 0.28]; %Mach numbers of interest
aero_analysis.tail.Re=aero_analysis.tail.rho.*aero_analysis.tail.air_velc*aero_analysis.tail.MAC./aero_analysis.tail.dyn_visc; %Reynolds number at each aspect of flight for the wing (accounts for the characteristic lengths)

aero_analysis.tail.Re(1)=aero_analysis.tail.rho(1)*aero_analysis.tail.Mach(1)*aero_analysis.tail.air_velc(1)*aero_analysis.tail.MAC/aero_analysis.tail.dyn_visc(1);
aero_analysis.tail.Re(2)=aero_analysis.tail.rho(2)*aero_analysis.tail.Mach(2)*aero_analysis.tail.air_velc(2)*aero_analysis.tail.MAC/aero_analysis.tail.dyn_visc(2);
aero_analysis.tail.Re(3)=aero_analysis.tail.rho(3)*aero_analysis.tail.v_takeoff_ms*aero_analysis.tail.MAC/aero_analysis.tail.dyn_visc(3);
aero_analysis.tail.Re(4)=aero_analysis.tail.rho(4)*aero_analysis.tail.v_landing_ms*aero_analysis.tail.MAC/aero_analysis.tail.dyn_visc(4);
aero_analysis.tail.Re(5)=aero_analysis.tail.rho(5)*aero_analysis.tail.Mach(5)*aero_analysis.tail.air_velc(5)*aero_analysis.tail.MAC/aero_analysis.tail.dyn_visc(5);

%(1): cruise
%(2): max speed at max altitude
%(3): take-off
%(4): approach
%(5): loiter (5000 ft @ M=0.28)


%% Correction factors, tail
aero_analysis.tail.beta=sqrt(1-aero_analysis.tail.Mach.^2); %compressibility effects
aero_analysis.tail.cl_alpha_ratio=[0.95,0.945,0.945,0.945, 0.95]; %read off graph - DATCOM 1978
aero_analysis.tail.cl_alpha_theory=2*pi+4.7*aero_analysis.tail.t_c*(1+0.00375*0.6); %0.6 deg
aero_analysis.tail.Cl_alpha_aerofoil=(1.05./aero_analysis.tail.beta).*aero_analysis.tail.cl_alpha_ratio*aero_analysis.tail.cl_alpha_theory; %cl of the airfoil as a function of the mach number. Varies so inputs must be changed manually based on readings off graph
aero_analysis.tail.eta=(aero_analysis.tail.beta).*(aero_analysis.tail.Cl_alpha_aerofoil)./(2*pi); %fraction between both
%aero_analysis.tail.eta=[1,1,1,1,1];
aero_analysis.tail.F=1.07*(1+d/b)^2;

%% winglet
%changes AR
%aero_analysis.wing.A_effective=AR*(1+h/b)^2;
%aero_analysis.wing.updated_AR=5; %need to establish the equation for this
%use the updated AR for induced drag calculations
%% CL
%lift curve slope per radian, accurate up to Mdd and reasonably accurate almost to M=1
%from Raymer
%aero_analysis.wing.Cl_alpha=2*pi*AR*(S_exp/S_ref)*aero_analysis.wing.F./(2+sqrt(4+(((AR^2.*aero_analysis.wing.beta.^2)/aero_analysis.wing.eta.^2)*(1+((tan(max_sweep_t)^2)./aero_analysis.wing.beta.^2))))); 
aero_analysis.tail.Cl_alpha=2*pi*AR*1./(2+sqrt(4+(((AR^2.*aero_analysis.tail.beta.^2)/aero_analysis.tail.eta.^2)*(1+((tan(max_sweep_t)^2)./aero_analysis.tail.beta.^2))))); 

aero_analysis.tail.Cl_max_tail=0.9*cl_max_airfoil*cos(lambda_quarter);

aero_analysis.tail.zero_aoa=0; %[degrees] to be defined

save('aero_analysis.mat', 'aero_analysis')