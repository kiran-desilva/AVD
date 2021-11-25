clear
clc

load('sizing')

%parameters from poster
W_S= 2973 %
C_cruise=295.1 %m/s
M_cruise=0.75
V_cruise=C_cruise*M_cruise %m/s
rho_cruise= 0.302748958
W_cruise_W_total=0.842
Max_takeoff_weight=3740 %kg
%mu_cruise=2.969*10^(-7)
mu_cruise=1.45*10^(-5)
character_c=0.6096
CL_design_cruise=( (W_cruise_W_total)*(W_S) )/ (0.5*rho_cruise*V_cruise^2)

% L_over_D=Cl/cd;
% plot( Alpha, Cl)
% ylabel('Cl')
% xlabel('alpha')
% [CL_Alpha_curve, gof] = CL_ALPHA_FIT(Alpha, Cl)
% 
% Cl_design=CL_Alpha_curve(4.25)



% Cl_max=max(Cl)



%mach


sweep_25=18.1692
M_DD=0.77
M_DD_Eff=M_DD*sqrt(cosd(sweep_25))
AR=7.8
x=cosd(22)

Cl_design_a=(CL_design_cruise/(0.9*0.95))/cosd(sweep_25)


Km=1.05

%lambda_opt=0.45*(0.85)^(-0.036*sweep_angle)
lambda_opt=0.45*exp(-0.036*sweep_25)
lambda_min=0.2*(AR)^(0.25)*(cosd(sweep_25)^2)
SP=(Max_takeoff_weight*9.8)/(W_S) 
mainwing.b = sqrt(AR*SP) 
lambda_opt=lambda_min
mainwing.Croot= (2*12.34086)/(9.811*(1+lambda_opt))
mainwing.Ctip=0.51*mainwing.Croot
C_mean=(2/3) * mainwing.Croot* ((1+lambda_opt+lambda_opt^2)/(1+lambda_opt))

y_bar=(mainwing.b/6)*((1+2*lambda_opt)/(1+lambda_opt))

t_c= 0.3*cosd(sweep_25)*((1-((5+M_DD_Eff^2)/(5+(Km-0.25*Cl_design_a)^2))^3.5)*((sqrt(1-M_DD_Eff^2)/M_DD_Eff^2)))^(2/3)

Re= (rho_cruise * V_cruise * C_mean )/ mu_cruise
twist=-3

[mainwing.sweepLE] = sweep_angle(sweep_25,0,25,AR,lambda_opt)
[fig] = horizontal_stab_plot(mainwing)

% i_w_mac=1
% i_w= (y_bar/(b/2))*twist + i_w_mac
alpha_0=-3
CL_cruise=0.4809

Cl_alpha=(12--8)/(8.5--6.2)
i_w= (CL_cruise/Cl_alpha)+alpha_0 +0.4*twist

wing.sweep_25 = sweep_25;
wing.Croot = mainwing.Croot;


save('wing','wing')