% Wing Lift Analysis
%ieo18
%10th Nov 21 - 29th Nov 21
%AVD
clear
clc

%% Load files
load('sizing.mat')
load('design.mat')
load('wing.mat')
%load('tailplane.mat')
load('fuse.mat')
%load("aero_analysis.mat")
load('tailplane.mat')
%% Wing Lift Estimation

%% Inputs from other scripts
d=(fuse.d_f)/39.3700787; %diameter of the fuselage [m] from Camille 26/11
AR=wing.Ar;
S_exp=8.1292; %[m^2]
S_ref=design.sref; %[m^2]
%h=0.3; %winglet height (m)
b=wing.b; %wing span (m)
max_sweep_t=16.024*pi/180; %sweep of wing at the chord location where the airfoil is thickest [rad]
sweep_LE=wing.sweepLE*pi/180; %[rad]
cl_max_airfoil=1.6; %maximum lift of airfoil
lambda_quarter=wing.sweep_25*pi/180; %quarter-chord sweep
aero_analysis.wing.t_c=0.12;
aero_analysis.wing.te_angle=wing.sweepTE; %degrees

%% Mach number definition, wing
aero_analysis.wing.stall_landing_no_safety=sqrt(2*(sizing.W0-sizing.Wf)/(1.225*wing.Sref*2.3));
aero_analysis.wing.stall_TO_no_safety=sqrt(2*(sizing.W0)/(1.225*wing.Sref*2.3));
%NOTE: need to change the weights used for the stall speed!!!!
aero_analysis.wing.v_landing_ms=1.3*sqrt(2*(sizing.W0-sizing.Wf)/(1.225*wing.Sref*2.3)); %units:m/s, safety margin of 1.3
aero_analysis.wing.v_landing_fts=3.2808*aero_analysis.wing.v_landing_ms; %units: ft/s
aero_analysis.wing.v_takeoff_ms=1.2*sqrt(2*(sizing.W0)/(1.225*wing.Sref*2.3)); %units:m/s, safety margin of 1.2
aero_analysis.wing.v_takeoff_fts=3.2808*aero_analysis.wing.v_takeoff_ms; %units: ft/s 

aero_analysis.wing.dyn_visc=[0.0000143226,0.0000143226,0.0000181206,0.0000181206,0.0000174141];
taper=wing.Ctip/wing.Croot;
aero_analysis.wing.MAC=wing.Croot*2/3*((1+taper+taper^2)/1+taper);
aero_analysis.wing.rho=[0.302,0.237,1.225,1.225,1.0552];
aero_analysis.wing.air_velc=[294.9,294.9,340.3,340.3,334.4];
aero_analysis.wing.Mach=[0.75,0.78,aero_analysis.wing.v_takeoff_ms/aero_analysis.wing.air_velc(3),aero_analysis.wing.v_landing_ms/aero_analysis.wing.air_velc(4),0.28]; %Mach numbers of interest
%aero_analysis.wing.Re=aero_analysis.wing.rho.*aero_analysis.wing.air_velc*aero_analysis.wing.MAC./aero_analysis.wing.dyn_visc; %Reynolds number at each aspect of flight for the wing (accounts for the characteristic lengths)

aero_analysis.wing.Re(1)=aero_analysis.wing.rho(1)*aero_analysis.wing.Mach(1)*aero_analysis.wing.air_velc(1)*aero_analysis.wing.MAC/aero_analysis.wing.dyn_visc(1);
aero_analysis.wing.Re(2)=aero_analysis.wing.rho(2)*aero_analysis.wing.Mach(2)*aero_analysis.wing.air_velc(2)*aero_analysis.wing.MAC/aero_analysis.wing.dyn_visc(2);
aero_analysis.wing.Re(3)=aero_analysis.wing.rho(3)*aero_analysis.wing.v_takeoff_ms*aero_analysis.wing.MAC/aero_analysis.wing.dyn_visc(3);
aero_analysis.wing.Re(4)=aero_analysis.wing.rho(4)*aero_analysis.wing.v_landing_ms*aero_analysis.wing.MAC/aero_analysis.wing.dyn_visc(4);
aero_analysis.wing.Re(5)=aero_analysis.wing.rho(5)*aero_analysis.wing.Mach(5)*aero_analysis.wing.air_velc(5)*aero_analysis.wing.MAC/aero_analysis.wing.dyn_visc(5);

%(1): cruise
%(2): max speed at max altitude
%(3): take-off
%(4): approach
%(5): loiter

%% Correction factors, wing
aero_analysis.wing.beta=sqrt(1-aero_analysis.wing.Mach.^2); %compressibility effects
aero_analysis.wing.cl_alpha_ratio=[0.95,0.945,0.945,0.945,0.95]; %read off graph - DATCOM 1978
aero_analysis.wing.cl_alpha_theory=2*pi+4.7*aero_analysis.wing.t_c*(1+0.00375*0.6); %0.6 deg
aero_analysis.wing.Cl_alpha_aerofoil=(1.05./aero_analysis.wing.beta).*aero_analysis.wing.cl_alpha_ratio*aero_analysis.wing.cl_alpha_theory; %cl of the airfoil as a function of the mach number. Varies so inputs must be changed manually based on readings off graph
aero_analysis.wing.eta=(aero_analysis.wing.beta).*(aero_analysis.wing.Cl_alpha_aerofoil)./(2*pi); %fraction between both
%aero_analysis.wing.eta=[1,1,1,1,1];
aero_analysis.wing.F=1.07*(1+d/b)^2;

%% winglet
%changes AR
%aero_analysis.wing.A_effective=AR*(1+h/b)^2;
%aero_analysis.wing.updated_AR=5; %need to establish the equation for this
%use the updated AR for induced drag calculations
%% CL
%lift curve slope per radian, accurate up to Mdd and reasonably accurate almost to M=1
%from Raymer
%aero_analysis.wing.Cl_alpha=2*pi*AR*(S_exp/S_ref)*aero_analysis.wing.F./(2+sqrt(4+(((AR^2.*aero_analysis.wing.beta.^2)/aero_analysis.wing.eta.^2)*(1+((tan(max_sweep_t)^2)./aero_analysis.wing.beta.^2))))); 
aero_analysis.wing.Cl_alpha=2*pi*AR*1./(2+sqrt(4+(((AR^2.*aero_analysis.wing.beta.^2)/aero_analysis.wing.eta.^2)*(1+((tan(max_sweep_t)^2)./aero_analysis.wing.beta.^2))))); 

aero_analysis.wing.Cl_max_wing=0.9*cl_max_airfoil*cos(lambda_quarter);

aero_analysis.wing.zero_aoa=-3; %[degrees]
%% Supersonic range
%Mach number above which the Mach number is supersonic
aero_analysis.wing.supersonic_Mach=1/cos(sweep_LE);

%% HLD
% inputs from other scripts
aero_analysis.wing.HLD.c_fraction=1.2415;
aero_analysis.wing.HLD.s_ref=wing.Sref;
aero_analysis.wing.HLD.s_flapped=6.6422;
aero_analysis.wing.HLD.delta_hl=wing.sweepTE; %hinge lift surface

%calculations
aero_analysis.wing.HLD.delta_cl_device=1.6*(aero_analysis.wing.HLD.c_fraction); %assumed double-slotted TE. Can change - pg 415 Raymer
%initialise vectors
aero_analysis.wing.HLD.delta_cl_max=[0,0];
aero_analysis.wing.HLD.delta_cl_max(2)=0.9*aero_analysis.wing.HLD.delta_cl_device*(aero_analysis.wing.HLD.s_flapped/aero_analysis.wing.HLD.s_ref)*cos(aero_analysis.wing.HLD.delta_hl);
aero_analysis.wing.HLD.delta_cl_max(1)=aero_analysis.wing.HLD.delta_cl_max(2)*0.8; %for take-off, lift increment is about 80% of these values
aero_analysis.wing.HLD.cl_alpha_flaps=[0,0];
aero_analysis.wing.HLD.delta_alpha=[0,0];
%(1): cruise
%(2): max

%(3): take-off
%(4): approach
%==> ignore the first 2 terms of the cl_alpha_flaps
aero_analysis.wing.HLD.alpha=[-10,-15]; %change in flap angle [degrees]
for i=1:2
    aero_analysis.wing.HLD.cl_alpha_flaps(i)=aero_analysis.wing.Cl_alpha(i+2)*(1+(aero_analysis.wing.HLD.c_fraction-1)*aero_analysis.wing.HLD.s_flapped/aero_analysis.wing.HLD.s_ref);
    aero_analysis.wing.HLD.delta_alpha(i)=aero_analysis.wing.HLD.alpha(i)*(aero_analysis.wing.HLD.s_flapped/aero_analysis.wing.HLD.s_ref)*cos(aero_analysis.wing.HLD.delta_hl);
end
save('aero_analysis.mat', 'aero_analysis')

%aero_analysis.wing.Cl_max_landing=aero_analysis.wing.HLD.delta_cl_max(2)+aero_analysis.wing.Cl_max_wing;

%aero_analysis.wing.Cl_max_approach=aero_analysis.wing.HLD.delta_cl_max(1)+aero_analysis.wing.Cl_max_wing;

%% Plots
