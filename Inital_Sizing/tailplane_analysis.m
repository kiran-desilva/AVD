clear
clc
load('tailplane.mat')

%% Correction factors, wing
aero_analysis.tail.beta=sqrt(1-aero_analysis.wing.Mach.^2); %compressibility effects
aero_analysis.tail.cl_alpha_ratio=[0.6,0.7,0.5,1]; %read off graph
aero_analysis.tail.cl_alpha_theory=2*pi+4.7*aero_analysis.wing.t_c*(1+0.00375*aero_analysis.wing.te_angle);
aero_analysis.tail.Cl_alpha_M=(1.05./aero_analysis.wing.beta).*aero_analysis.wing.cl_alpha_ratio*aero_analysis.wing.cl_alpha_theory; %cl of the airfoil as a function of the mach number. Varies so inputs must be changed manually based on readings off graph
aero_analysis.tail.eta=(aero_analysis.wing.beta).*(aero_analysis.wing.Cl_alpha_M)./(2*pi);
aero_analysis.tail.F=1.07*(1+d/b)^2;

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