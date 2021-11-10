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
max_sweep_quarter=20*pi/180;
sweep_LE=20*pi/180;

%% Mach number definition
aero_analysis.wing.Mach=[0.75,0.78]; %Mach numbers of interest
%(1): cruise
%(2): max

%% Correction factors
aero_analysis.wing.beta=sqrt(1-aero_analysis.wing.Mach.^2); %compressibility effects
aero_analysis.wing.Cl_alpha_M=[1,1];
aero_analysis.wing.eta=(aero_analysis.wing.beta).*(aero_analysis.wing.Cl_alpha_M)./(2*pi);
aero_analysis.wing.F=1.07*(1+d/b)^2;


%% CL_alpha
%from Raymer
aero_analysis.wing.Cl_alpha=2*pi*AR*(S_exp/S_ref)*aero_analysis.wing.F./(2+sqrt(4+(((AR^2.*aero_analysis.wing.beta.^2)/aero_analysis.wing.eta.^2)*(1+((tan(max_sweep_quarter)^2)./aero_analysis.wing.beta.^2)))));

%% Supersonic range
%Mach number above which the Mach number is supersonic
aero_analysis.wing.supersonic_Mach=1/cos(sweep_LE);

%% HLD

