% Wing Lift Analysis
%ieo18
%10th Nov 21 -
%AVD
clear
clc
%all units in SI

%% Run lift matlab file
%run Wing.m; error
load('sizing.mat')
load('wing.mat')
load('uc.mat')
load("aero_analysis.mat")
load("tailplane.mat")
load("powerplant.mat")
load("fuse.mat")
%% Drag Analysis
%Parasite drag (zero-lift) & lift-induced

%% Inputs from other scripts
S_ratio=6; %=S_wet/S_ref
S=[wing.Sref,100,8.1292,((1.3985/2)/3.2808)^2*pi*2,50];
%S(1):S_ref
%S(2):A_eff
%S(3):S_exposed_wing
%S(4): frontal area of nacelles (total)
%S(5): wetted surface area of the whole plane

t_c=[0.12,tailplane.horizontal.tc,tailplane.vertical.tc]; %thickness to chord ratio
%(1): wing
%(2): horizontal stabiliser
%(3): vertical stabiliser

sweep_angles_quarter=[wing.sweep_25,tailplane.horizontal.sweep_25,tailplane.vertical.sweep_25].*pi/180;
%sweep angles @ 1/4 chord
%sweep angles(1): wing
%sweep angles(2): horizontal tail
%sweep angles(3): vertical tail
sweep_angles_half=14.5429*pi/180; %sweep of wing at 1/2 chord


d_fuselage=fuse.d_f/39.3700787; %fuselage diameter [m]
l_fuselage=fuse.structural_length/39.3700787; %fuselage length [m]

delta=[0,0,15,45,0]*pi/180; %angle of deflection of HLD [rad]
%delta(1):Cruise
%delta(2):Max
%delta(3):TO
%delta(4):Landing
%delta(5):Loiter

%uc_frontal_area_m_sq comes from P's script
A_UC_front=[0, 0, uc.uc_frontal_area_m_sq, uc.uc_frontal_area_m_sq,0]; %Area of the UC [m^2]
%(1): cruise
%(2): max
%(3): TO
%(4): Landing

%% Wave drag, wing, M_dd
%Drag divergence Mach number
%only for wing
cl_cruise=sizing.fraction.before_cruise*sizing.W0/(0.5*0.302*(aero_analysis.wing.Mach(2)*aero_analysis.wing.air_velc(2))^2*wing.Sref); %need to define
korn=0.87;
aero_analysis.drag.M_DD=(korn/cos(sweep_angles_half))-(t_c(1)/cos(sweep_angles_half)^2)-(cl_cruise/(10*cos(sweep_angles_half)^3));
aero_analysis.drag.M_cr_wing=aero_analysis.drag.M_DD-0.08;
aero_analysis.drag.M.wing=aero_analysis.wing.Mach;

aero_analysis.drag.C_d_aerofoil=zeros(1,length(aero_analysis.drag.M.wing));
for i=1:length(aero_analysis.drag.M.wing)
    if aero_analysis.drag.M.wing(i) >= aero_analysis.drag.M_cr_wing
        aero_analysis.drag.C_d_aerofoil(i)=20*(aero_analysis.drag.M.wing(i)-aero_analysis.drag.M_cr_wing)^4;
    else
        aero_analysis.drag.C_d_aerofoil(i)=0;
    end
end

aero_analysis.drag.wave=aero_analysis.drag.C_d_aerofoil*(S(3)/S(1));

%% Parasitic drag (not accounting for wave drag): skin friction, form & intereferece, leakage, protuberances

%% U/C
for i=1:length(A_UC_front)
    aero_analysis.drag.uc(i)=2.25*(A_UC_front(i)/S(1));
end

%% Flaps
b=wing.b;
b_f=3;

for i=1:length(delta)
    aero_analysis.drag.HLD(i)=0.0023*(0.4138)*delta(i);
end

%% windmilling engines
aero_analysis.drag.WE=0.3*S(4)/wing.Sref;

%% Fuselage upsweep
upsweep=11*pi/180; %upsweep of the fuselage
aero_analysis.drag.upsweep=3.83*(pi*d_fuselage^2/(4*S(1)))*(upsweep)^(2.5);

%% Combine Misc terms
aero_analysis.drag.cd_misc=zeros(1,length(delta));
for j=1:length(delta)
    aero_analysis.drag.cd_misc(j)=aero_analysis.drag.upsweep+aero_analysis.drag.WE+aero_analysis.drag.HLD(j)+aero_analysis.drag.uc(j);
end

%% Friction drag
%(1):fuselage
%(2):wing
%(3): horizontal stabiliser
%(4): vertical stabiliser
%(5): nacelle - note, there are 2 to consider (2 engines!)
d_components=[d_fuselage,b,tailplane.horizontal.Cmac,tailplane.vertical.Cmac,sqrt(4*uc.uc_frontal_area_m_sq/pi)];

aero_analysis.drag.l_components=[l_fuselage,wing.Cmac,tailplane.horizontal.Cmac,tailplane.vertical.Cmac,powerplant.nacelle_length_ft/3.281
]; %define characteristic length of each component
aero_analysis.drag.fineness=zeros(1,length(d_components));
for i=1:length(d_components)
    aero_analysis.drag.fineness(i)=aero_analysis.drag.l_components(i)/d_components(i);
end

chordwise_max_thickness=[0.399,0.4,0.35];
%(1):wing
%(2):horizontal stabiliser
%(3):vertical stabiliser

sweep_max_thickness=[16.024,tailplane.horizontal.sweepmaxTC,tailplane.vertical.sweepmaxTC]*pi/180;
%(1):wing
%(2):horizontal stabiliser
%(3):vertical stabiliser

aero_analysis.drag.FF=zeros(1,length(d_components));
aero_analysis.drag.FF(1)=1+(60/aero_analysis.drag.fineness(1)^3)+(aero_analysis.drag.fineness(1)/400);
for i=1:length(chordwise_max_thickness)
    aero_analysis.drag.FF(i+1)=(1+(0.6/chordwise_max_thickness(i))*t_c(i)+100*(t_c(i))^4)*(1.34*1^0.18*(cos(sweep_max_thickness(i))^0.28));
end
aero_analysis.drag.FF(5)=1+(0.35/aero_analysis.drag.fineness(5));

Q=[1,1,1.055,1.055,1.5]; %assuming 5.5% for horizontal and vertical stabiliser

S_wet=[28.252,14,2.5,2.5,powerplant.nacelle_wetted_area_ft_sq*2/10.764]; %wetted area, exact values to be added
%(1):fuselage
%(2):wing
%(3): horizontal stabiliser
%(4): vertical stabiliser
%(5): nacelle - note, there are 2 to consider (2 engines!)

aero_analysis.drag.Re_cutoff=38.21.*(aero_analysis.drag.l_components*3.28084/2.08E-05).^(1.053);
%note: must convert lengths to ft for the cutoff Re

%conditions
%(1):cruise
%(2):max
%(3):TO
%(4):landing
%(5): loiter

U=aero_analysis.wing.air_velc.*aero_analysis.wing.Mach; %velocity for conditions
rho=aero_analysis.wing.rho; %density for conditions
nu=aero_analysis.wing.dyn_visc; %dynamic viscosity for conditions
M=aero_analysis.wing.Mach; %Mach number

Re=zeros(length(U), length(aero_analysis.drag.l_components));
for i=1:length(aero_analysis.drag.l_components)
    for j=1:length(U)
        Re(j,i)=rho(j)*U(j)*aero_analysis.drag.l_components(i)/nu(j);
    end
end

%must evaluate if Re is greater than Re_cutoff - not the case

aero_analysis.drag.cf=zeros(length(U), length(aero_analysis.drag.l_components));
for i=1:length(aero_analysis.drag.l_components)
    for j=1:length(U)
        aero_analysis.drag.cf(j,i)=0.455/((log10(Re(j,i)))^2.58*(1+0.144*M(j)^2)^0.65);
    end
end

aero_analysis.drag.skin_comp=zeros(length(U), length(aero_analysis.drag.l_components));
aero_analysis.drag.skin_friction=zeros(1,length(U));
for j=1:length(U)
    for i=1:length(aero_analysis.drag.l_components)
        aero_analysis.drag.skin_comp(j,i)=aero_analysis.drag.cf(j,i)*aero_analysis.drag.FF(i)*Q(i)*S_wet(i);
    end
    aero_analysis.drag.skin_friction(j)=sum(aero_analysis.drag.skin_comp(j,:))/S(5);
end

%% Delta parasitic flap drag
% s_flapped=aero_analysis.wing.HLD.s_flapped;
% aero_analysis.drag.delta_cd0=zeros(1,length(U));
% aero_analysis.drag.delta_cd0(3)=0.0074*aero_analysis.wing.HLD.c_fraction*(s_flapped/wing.Sref)*(30-10);
% aero_analysis.drag.delta_cd0(4)=0.0074*aero_analysis.wing.HLD.c_fraction*(s_flapped/wing.Sref)*(65-10);

%% Leakage
%aero_analysis.drag.leakage=1; %find a method
%Assuming it's between 3-5% (currently 4%)
%% Total parasitic drag
% Skin friction, Misc, Leakage
%accounting for 4% of leakage
for j=1:length(U)
    aero_analysis.drag.cd0(j)=1.02*(aero_analysis.drag.skin_friction(j)+aero_analysis.drag.cd_misc(j));
end

%% Induced drag
% e_theoretical
sweep_25=[wing.sweep_25,tailplane.horizontal.sweep_25]; %input, sweep angle of 1/4 chord (degrees)
df_b=0.12;
ke_f=0.971;
k_e_d0=0.864;
big_k=0.38;
induced_AR=[wing.Ar, tailplane.horizontal.Ar];
a_e=-0.001521;
b_e=10.82;
cl=[0.45,0.2]; %define later

aero_analysis.induced_drag.taper=[wing.Ctip/wing.Croot, tailplane.horizontal.Ctip/tailplane.horizontal.Croot] ; %from Nadia
aero_analysis.induced_drag.lambda_delta=-0.357+0.45*exp(0.0375.*sweep_25);
aero_analysis.induced_drag.taper_and_delta=aero_analysis.induced_drag.taper-aero_analysis.induced_drag.lambda_delta;
aero_analysis.induced_drag.fourth_order_lambda=0.0524.*aero_analysis.induced_drag.taper_and_delta.^4-0.15.*aero_analysis.induced_drag.taper_and_delta.^3+0.1659.*aero_analysis.induced_drag.taper_and_delta.^2-0.0706.*aero_analysis.induced_drag.taper_and_delta+0.0119;

aero_analysis.induced_drag.e_theoretical=1./(1+aero_analysis.induced_drag.fourth_order_lambda.*induced_AR);

for j=1:length(U)
    aero_analysis.induced_drag.k_e_m(j)=((a_e)*(M(j)/0.3-1)^(b_e))+1;
end

aero_analysis.induced_drag.Q=1./(aero_analysis.induced_drag.e_theoretical.*ke_f);

aero_analysis.induced_drag.P=big_k.*aero_analysis.drag.cd0;

%assuming efficiency of tail (eta_h) is 1 since it's a t-tail
angle_flight=[-0.0168,-0.0168,9.8,8.42,-0.0168]*pi/180;
for j=1:length(U)
    aero_analysis.induced_drag.wing.e(j)=aero_analysis.induced_drag.k_e_m(j)/(aero_analysis.induced_drag.Q(1)+aero_analysis.induced_drag.P(j)*pi*induced_AR(1));
    aero_analysis.induced_drag.wing.e_V2(j)=aero_analysis.induced_drag.e_theoretical(1)*ke_f*k_e_d0*aero_analysis.induced_drag.k_e_m(j); %without knowing Cd0
    aero_analysis.induced_drag.tail.e(j)=aero_analysis.induced_drag.k_e_m(j)/(aero_analysis.induced_drag.Q(2)+aero_analysis.induced_drag.P(j)*pi*induced_AR(2));
    aero_analysis.induced_drag.tail.e_V2(j)=aero_analysis.induced_drag.e_theoretical(2)*ke_f*k_e_d0*aero_analysis.induced_drag.k_e_m(j); %without knowing Cd0
%    aero_analysis.induced_drag.wing.cd_0(j)=cl^2/(pi*induced_AR(1)*aero_analysis.induced_drag.wing.e(j));
%    aero_analysis.induced_drag.wing.cd_0_V2(j)=cl^2/(pi*induced_AR(1)*aero_analysis.induced_drag.wing.e_V2(j));
    aero_analysis.induced_drag.wing.cd_i(j)=(aero_analysis.wing.Cl_alpha(j)*angle_flight(j))^2/(pi*induced_AR(1)*aero_analysis.induced_drag.e_theoretical(1));
    aero_analysis.induced_drag.tail.cd_i(j)=(aero_analysis.tail.Cl_alpha(j)*angle_flight(j))^2/(pi*induced_AR(2)*aero_analysis.induced_drag.e_theoretical(2))*1*tailplane.horizontal.s/wing.Sref;
    aero_analysis.induced_drag.cd_i(j)=aero_analysis.induced_drag.wing.cd_i(j)+aero_analysis.induced_drag.tail.cd_i(j);
end

aero_analysis.induced_drag.HLD=0.28^2.*aero_analysis.wing.HLD.delta_cl_max.^2*cos(wing.sweep_25);
aero_analysis.induced_drag.cd_i(3)=aero_analysis.induced_drag.cd_i(3)+aero_analysis.induced_drag.HLD(1);
aero_analysis.induced_drag.cd_i(4)=aero_analysis.induced_drag.cd_i(4)+aero_analysis.induced_drag.HLD(2);

%% Total Drag


aero_analysis.drag.cd_induced_total=aero_analysis.induced_drag.cd_i;
aero_analysis.drag.cd_parasitic_total=aero_analysis.drag.cd0;
aero_analysis.drag.cd_total=aero_analysis.drag.cd_parasitic_total+aero_analysis.drag.cd_induced_total+aero_analysis.drag.wave;

save('aero_analysis.mat', 'aero_analysis')



