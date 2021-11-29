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
%% Drag Analysis
%Parasite drag (zero-lift) & lift-induced

%% Inputs from other scripts
S_ratio=6; %=S_wet/S_ref
S=[11.4246,100,5];
%S(1):S_ref
%S(2):A_eff
%S(3):S_exposed_wing

t_c=[0.12,tailplane.horizontal.tc,tailplane.vertical.tc]; %thickness to chord ratio
%(1): wing
%(2): horizontal stabiliser
%(3): vertical stabiliser

sweep_angles=[wing.sweep_25,tailplane.horizontal.sweep_25,tailplane.vertical.sweep_25].*pi/180;
%sweep angles @ 1/4 chord
%sweep angles(1): wing
%sweep angles(2): horizontal tail
%sweep angles(3): vertical tail

d_fuselage=1.524; %fuselage diameter [m]
l_fuselage=14.224; %fuselage length [m]

delta=[0,0,15,45]*pi/180; %angle of deflection of HLD [rad]
%delta(1):Cruise
%delta(2):Max
%delta(3):TO
%delta(4):Landing

%uc_frontal_area_m_sq comes from P's script
A_UC_front=[0, 0, uc.uc_frontal_area_m_sq, uc.uc_frontal_area_m_sq]; %Area of the UC [m^2]
%(1): cruise
%(2): max
%(3): TO
%(4): Landing

%% Wave drag, wing, M_dd
%Drag divergence Mach number
%only for wing
cl=1; %need to define
korn=0.87;
aero_analysis.drag.M_DD=(korn-(t_c(1)/cos(sweep_angles(1)))-(cl/(10*cos(sweep_angles(1))^2)))/cos(sweep_angles(1));
aero_analysis.drag.M_cr_wing=aero_analysis.drag.M_DD-nthroot((0.1/80),3);
aero_analysis.drag.M.wing=1;

for i=1:length(aero_analysis.drag.M.wing)
    if aero_analysis.drag.M.wing(i) >= aero_analysis.drag.M_cr_wing
        C_d_aerofoil=20*(aero_analysis.drag.M.wing-aero_analysis.drag.M_cr_wing)^4;
    else
        C_d_aerofoil=0;
    end
end

aero_analysis.drag.wave=C_d_aerofoil*(S(3)/S(1));

%% Parasitic drag (not accounting for wave drag): skin friction, form & intereferece, leakage, protuberances

%% U/C
for i=1:length(A_UC_front)
    aero_analysis.drag.uc(i)=2.25*(A_UC_front(i)/S(1));
end

%% Flaps
b=5;
b_f=3;

for i=1:length(delta)
    aero_analysis.drag.HLD(i)=0.0023*(b_f/b)*delta(i);
end

%% windmilling engines
aero_analysis.drag.WE=0.3*S(2)/S(1);

%% Fuselage upsweep
upsweep=11*pi/180; %upsweep of the fuselage
aero_analysis.drag.upsweep=3.83*(pi*d_fuselage^2/(4*S(1)))*upsweep^(2.5);

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
d_components=[d_fuselage,b,tailplane.horizontal.Cmac,tailplane.vertical.Cmac,1];

aero_analysis.drag.l_components=[20,2,3,4,5]; %define characteristic length of each component
aero_analysis.drag.fineness=zeros(1,length(d_components));
for i=1:length(d_components)
    aero_analysis.drag.fineness(i)=aero_analysis.drag.l_components(i)/d_components(i);
end

chordwise_max_thickness=[0.399,0.4,0.5];
%(1):wing
%(2):horizontal stabiliser
%(3):vertical stabiliser

sweep_max_thickness=[10,10,10]*pi/180;
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

S_wet=[1,1,1,1,1]; %wetted area, exact values to be added

aero_analysis.drag.Re_cutoff=38.21.*(aero_analysis.drag.l_components*3.28084/2.08E-05).^(1.053);
%note: must convert lengths to ft for the cutoff Re

%conditions
%(1):cruise
%(2):max
%(3):TO
%(4):landing

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

%must evaluate if Re is greater than Re_cutoff

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
    aero_analysis.drag.skin_friction(j)=sum(aero_analysis.drag.skin_comp(j,:))/S(1);
end

%% Leakage
%aero_analysis.drag.leakage=1; %find a method
%Assuming it's between 3-5% (currently 4%)
%% Total parasitic drag
% Skin friction, Misc, Leakage
%accounting for 4% of leakage
for j=1:length(U)
    aero_analysis.drag.cd0(j)=1.04*(aero_analysis.drag.skin_friction(j)+aero_analysis.drag.cd_misc(j));
end

%% Induced drag
% e_theoretical
sweep_25=18.1692; %input, sweep angle of 1/4 chord (degrees)
df_b=0.12;
ke_f=0.971;
k_e_d0=0.864;
big_k=0.38;
wing_A=7.8;
a_e=-0.001521;
b_e=10.82;
cl=2.1; %define later

aero_analysis.induced_drag.wing.taper=wing.Ctip/wing.Croot;; %from Nadia
aero_analysis.induced_drag.wing.lambda_delta=-0.357+0.45*exp(0.0375*sweep_25);
aero_analysis.induced_drag.wing.taper_and_delta=aero_analysis.induced_drag.wing.taper-aero_analysis.induced_drag.wing.lambda_delta;
aero_analysis.induced_drag.wing.fourth_order_lambda=0.0524*aero_analysis.induced_drag.wing.taper_and_delta^4-0.15*aero_analysis.induced_drag.wing.taper_and_delta^3+0.1659*aero_analysis.induced_drag.wing.taper_and_delta^2-0.0706*aero_analysis.induced_drag.wing.taper_and_delta+0.0119;

aero_analysis.induced_drag.wing.e_theoretical=1/(1+aero_analysis.induced_drag.wing.fourth_order_lambda*wing_A);

for j=1:length(U)
    aero_analysis.induced_drag.wing.k_e_m(j)=(a_e*(M(j)/0.3-1)^b_e)+1;
end

aero_analysis.induced_drag.wing.Q=1/(aero_analysis.induced_drag.wing.e_theoretical*ke_f);

aero_analysis.induced_drag.wing.P=big_k.*aero_analysis.drag.cd0;


for j=1:length(U)
    aero_analysis.induced_drag.wing.e(j)=aero_analysis.induced_drag.wing.k_e_m(j)/(aero_analysis.induced_drag.wing.Q+aero_analysis.induced_drag.wing.P(j)*pi*wing_A);
    aero_analysis.induced_drag.wing.e_V2(j)=aero_analysis.induced_drag.wing.e_theoretical*ke_f*k_e_d0*aero_analysis.induced_drag.wing.k_e_m(j); %without knowing Cd0
    aero_analysis.induced_drag.wing.cd_0(j)=cl^2/(pi*wing_A*aero_analysis.induced_drag.wing.e(j));
    aero_analysis.induced_drag.wing.cd_0_V2(j)=cl^2/(pi*wing_A*aero_analysis.induced_drag.wing.e_V2(j));
end
save('aero_analysis.mat', 'aero_analysis')

%% Total Drag

