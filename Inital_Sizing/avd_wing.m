clear
clc

load('parameters')
load('sizing')
load('design')

%parameters from poster
W_S= design.w_s;%
% C_cruise=295.1 %m/s
[~,C_cruise,~,rho_cruise] = atmosisa(distdim(parameters.cruise_alt_ft,'ft','m'));
M_cruise=parameters.cruise_mach;
V_cruise=C_cruise*M_cruise %m/s
% rho_cruise= 0.302748958
W_cruise_W_total=0.842 % -> needs to be changed
Max_takeoff_weight=sizing.W0; %kg
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
t_c_cruiseonly=-0.0439*atan(3.3450*M_cruise+-3.0231)+0.0986



%mach
sweep_25=18.1692
M_DD=0.75
M_DD_Eff=M_DD*sqrt(cosd(sweep_25))
AR=7.8

x=cosd(22)

Cl_design_a=(CL_design_cruise/(0.9*0.95))/cosd(sweep_25)


Km=1.05

%lambda_opt=0.45*(0.85)^(-0.036*sweep_angle)
lambda_opt=0.45*exp(-0.036*sweep_25)
lambda_min=0.2*(AR)^(0.25)*(cosd(sweep_25)^2)
SP=(Max_takeoff_weight)/(W_S) 
mainwing.b = sqrt(AR*SP) 
if lambda_min>lambda_opt
    lambda=lambda_min
else
    lambda=lambda_op
end
mainwing.lambda=lambda

mainwing.Croot= (2*SP)/(mainwing.b*(1+lambda));
mainwing.Ctip=lambda*mainwing.Croot;
C_mean=(2/3) * mainwing.Croot* ((1+lambda+lambda^2)/(1+lambda))
mainwing.C_mean=C_mean
%c_dist
m=(0-mainwing.b/2)/(mainwing.Croot/mainwing.Ctip)

y_bar=(mainwing.b/6)*((1+2*lambda)/(1+lambda))
mainwing.y_bar=y_bar
t_c= 0.3*cosd(sweep_25)*((1-((5+M_DD_Eff^2)/(5+(Km-0.25*Cl_design_a)^2))^3.5)*((sqrt(1-M_DD_Eff^2)/M_DD_Eff^2)))^(2/3)

Re= (rho_cruise * V_cruise * C_mean )/ mu_cruise
mainwing.Re=Re
twist=-3
mainwing.twist=twist

[mainwing.sweepLE] = sweep_angle(sweep_25,0,25,AR,lambda);
[mainwing.sweepTE] = sweep_angle(sweep_25,100,25,AR,lambda);

% i_w_mac=1
% i_w= (y_bar/(b/2))*twist*(x_position) + i_w_mac
alpha_0=-3
CL_cruise=0.4809
mainwing.alpha_0=alpha_0
mainwing.CL_cruise=CL_cruise

Cl_alpha=(1.6-0.4)/((16*(pi/180)--0*(pi/180))) %issue find clalpha
i_w= (CL_cruise/Cl_alpha)+alpha_0 +0.4*twist

i_mac=1
x_span=[0:mainwing.b/2];
span_twist=((twist)/(mainwing.b/2))*x_span +i_mac+1.20939
figure
plot(x_span,span_twist)

i_w_method=((twist*y_bar)/(mainwing.b/2)) +i_mac+1.20939
CLmax=1.6
CLmax_eff=0.9*CLmax*cosd(sweep_25)
delta_y=19.3*t_c %for 65 series

CLmax_takeooff = 1.9
Clmax_landing = 2.4
if CLmax_takeooff<Clmax_landing
    q=Clmax_landing
else
    q=CLmax_takeooff
end
CL_req_inc=CLmax-q
%need 1.0669

%30% of the wing must be allocated to aileron placement for sufficient controllability




m=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b);
m2=mainwing.sweepLE*(pi/180);
mainwing.HDL_PERC=0.7689;
mainwing.HDL_Croot=(1-mainwing.HDL_PERC)*mainwing.Croot;
mainwing.HDL_Ctip=(1-mainwing.HDL_PERC)*mainwing.Ctip;
[mainwing.sweepHDL]= sweep_angle(sweep_25,mainwing.HDL_PERC*100,25,AR,lambda);


sref=((mainwing.Croot+mainwing.Ctip)/2)*mainwing.b;
mainwing.sref=sref
%aileron starts at 70 of span
mainwing.aileron_start_top=0.7*mainwing.b/2;
mainwing.aileron_ypos_top=tand(mainwing.sweepHDL)*mainwing.aileron_start_top;
mainwing.aileron_Croot=(1-mainwing.HDL_PERC)*(((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*mainwing.aileron_start_top +mainwing.Croot);
saileron=0.5*(mainwing.aileron_Croot+(1-mainwing.HDL_PERC)*mainwing.Ctip)*(mainwing.b-mainwing.aileron_start_top);
mainwing.aileron_start_TE=0.7*mainwing.b/2;
mainwing.aileron_ypos_TE=tand(mainwing.sweepTE)*mainwing.aileron_start_TE;
%1.52/2
mainwing.HDL_start=1.52/2 + 0.05*mainwing.b %therefore 0.05 not including fueselage
C_HDL_start_root=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*(mainwing.HDL_start) +mainwing.Croot
C_70=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*(0.7*mainwing.b*0.5) +mainwing.Croot
sflapped=(((C_HDL_start_root+C_70)/2)*(0.7*mainwing.b-mainwing.HDL_start))
mainwing.sflapped=sflapped
mainwing.C_HDL_root=(1-mainwing.HDL_PERC)*(((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*mainwing.HDL_start +mainwing.Croot);



x_dist=[0:mainwing.b/2]
chord_dist=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*x_dist +mainwing.Croot
c_mac=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*y_bar +mainwing.Croot
mainwing.c_mac=c_mac
delta_CL=q-CLmax_eff

syms x sweep_hdl
C_dist=((2*(mainwing.Ctip-mainwing.Croot))/mainwing.b)*x+mainwing.Croot
c_flap_ext=((CL_req_inc*sref)/(0.9*1.3*sflapped))*(1/cosd(sweep_hdl))*C_dist

%delta_Cl=1.3*((c_dash)/chord_dist)

%CLeff is the 1.0669, delta_cl = 1.3cdash/c then find the other in terms of
%cf
%delta_CL=0.9*1.3*((c_dash)/chord_dist*(sflapped/sref)*cosd(mainwing.sweepHDL)
%c_dash=delta_CL*0.9*1.3*chord_dist*(sref/sflapped)*(1/cosd(mainwing.sweepHDL))
syms sweep_hdl_find
[mainwing.sweepHDL] = sweep_angle(sweep_25,75,25,AR,lambda);
HDL_coeff=1.6
k=(delta_CL*sref) / (0.9*HDL_coeff*sflapped*cosd(sweep_hdl_find)) %c'/c
% fplot(@(sweep_hdl_find) k)
fplot(k, [mainwing.sweepTE, mainwing.sweepLE])
hold on
plot([mainwing.sweepHDL, mainwing.sweepHDL],[1,2])
k=(delta_CL*sref) / (0.9*HDL_coeff*sflapped*cosd(mainwing.sweepTE))
%k=cbar/c
cbar_root=k*mainwing.Croot
HDL_chord=cbar_root-mainwing.Croot
fraction_HDLchord_over_wingchord=HDL_chord/mainwing.Croot
perc_chord=(mainwing.Croot-HDL_chord)/mainwing.Croot
mainwing.fraction_HDLchord_over_wingchord=fraction_HDLchord_over_wingchord

[fig] = hdl(mainwing)

%s_exposed=
%s_wet=s_exposed*(1.977+0.52(t_c))


mainwing.AR=AR
mainwing.sweep_25=sweep_25
mainwing.M_DD=M_DD
mainwing.M_DD_Eff=M_DD_Eff
mainwing.Cl_design_a=Cl_design_a
mainwing.Km=Km

 save('wing','wing.mat')