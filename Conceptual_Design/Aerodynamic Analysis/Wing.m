% Wing Lift Analysis
%ieo18
%10th Nov 21 -
%AVD
% clear
% clc
%% Wing Lift Estimation

%% Inputs from other scripts
d=1; %will come from another script
b=1; %will come from another script
AR=6;
S_exp=0.5;
S_ref=1;
h=0.3; %winglet height (m)
b=5; %wing span (m)
max_sweep_t=20*pi/180; %sweep of wing at the chord location where the airfoil is thickest
sweep_LE=20*pi/180;
cl_max_airfoil=2.4; %maximum lift of airfoil, at similar Reynolds number (Raymer)
lambda_quarter=20*pi/180; %quarter-chord sweep

%% Mach number definition
aero_analysis.wing.Mach=[0.75,0.78]; %Mach numbers of interest
%(1): cruise
%(2): max

%% Correction factors
aero_analysis.wing.beta=sqrt(1-aero_analysis.wing.Mach.^2); %compressibility effects
aero_analysis.wing.Cl_alpha_M=[1,1];
aero_analysis.wing.eta=(aero_analysis.wing.beta).*(aero_analysis.wing.Cl_alpha_M)./(2*pi);
aero_analysis.wing.F=1.07*(1+d/b)^2;

%% winglet
%changes AR
aero_analysis.wing.A_effective=AR*(1+h/b)^2;
aero_analysis.wing.updated_AR=5; %need to establish the equation for this
%use the updated AR for induced drag calculations
%% CL
%lift curve slope per radian, accurate up to Mdd and reasonably accurate almost to M=1
%from Raymer
aero_analysis.wing.Cl_alpha=2*pi*AR*(S_exp/S_ref)*aero_analysis.wing.F./(2+sqrt(4+(((AR^2.*aero_analysis.wing.beta.^2)/aero_analysis.wing.eta.^2)*(1+((tan(max_sweep_t)^2)./aero_analysis.wing.beta.^2))))); 
aero_analysis.wing.Cl_max_wing=0.9*cl_max_airfoil*cos(lambda_quarter);


%% Supersonic range
%Mach number above which the Mach number is supersonic
aero_analysis.wing.supersonic_Mach=1/cos(sweep_LE);

%% HLD
% inputs from other scripts
aero_analysis.wing.HLD.c=1;
aero_analysis.wing.HLD.c_dash=1.2;
aero_analysis.wing.HLD.s_ref=1;
aero_analysis.wing.HLD.s_flapped=1.2;
aero_analysis.wing.HLD.delta_hl=20*pi/180; %hinge lift surface


aero_analysis.wing.HLD.delta_cl_device=1.3*(aero_analysis.wing.HLD.c_dash/aero_analysis.wing.HLD.c); %assumed Fowler. Can change - pg 415 Raymer
aero_analysis.wing.HLD.delta_cl_max=0.9*aero_analysis.wing.HLD.delta_cl_device*(aero_analysis.wing.HLD.s_flapped/aero_analysis.wing.HLD.s_ref)*cos(aero_analysis.wing.HLD.delta_hl);


