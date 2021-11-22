% Wing Lift Analysis
%ieo18
%10th Nov 21 -
%AVD
clear
clc
%all units in SI
%% Drag Analysis
%Parasite drag (zero-lift) & lift-induced

%% Inputs from other scripts
S_ratio=6; %=S_wet/S_ref
S=[500,100,5];
%S(1):S_ref
%S(2):A_eff
%S(3): S_exposed_wing

t_c=[0.12,0.12,0.13]; %thickness to chord ratio
%(1): wing
%(2): horizontal stabiliser
%(3): vertical stabiliser

sweep_angles=[20,40,20].*pi/180;
%sweep angles(1): wing
%sweep angles(2): vertical tail
%sweep angles(3): horizontal tail

d_fuselage=5; %fuselage diameter


%% Wave drag, wing, M_dd
%Drag divergence Mach number
%only for wing
cl=1; %need to define
korn=0.87;
aero_analysis.drag.M_DD=(korn-(t_c_wing/cos(sweep_angles(1)))-(cl/(10*cos(sweep_angles(1))^2)))/cos(sweep_angles(1));
aero_analysis.drag.M_cr_wing=aero_analysis.drag.M_DD-nthroot((0.1/80),3);
aero_analysis.drag.M.wing=1;

if aero_analysis.drag.M.wing >= aero_analysis.drag.M_cr_wing
    C_d_aerofoil=20*(aero_analysis.drag.M.wing-aero_analysis.drag.M_cr_wing)^4;
else
    C_d_aerofoil=0;
end

aero_analysis.drag.wave=C_d_aerofoil*(S(3)/S(1));

%% Parasitic drag (not accounting for wave drag): skin friction, form & intereferece, leakage, protuberances
%% U/C
A_UC_front=5; %TO BE FILLED IN LATER
aero_analysis.drag.uc=2.25*(A_UC_front/S(1));

%% Flaps
b=5;
b_f=3;
delta=[15,45]*pi/180;
%delta(1):TO
%delta(2):Landing
for i=1:length(delta)
    aero_analysis.drag.HLD(i)=0.0023*(b_f/b)*delta(i);
end

%% windmilling engines
aero_analysis.drag.WE=0.3*S(2)/S(1);

%% Fuselage upsweep
upsweep=4*pi/180;
aero_analysis.drag.upsweep=3.83*(pi*d_fuselage^2/(4*S(1)))*upsweep^(2.5);

% %% Parasite Drag
% 
% % C_D_0 initial estimate (not very accurate) - light aircaraft, twin engine
% aero_analysis.drag.C_d0_estimate=0.0045*S_ratio;

%% Friction drag
%(1):fuselage
%(2):wing
%(3): horizontal stabiliser
%(4): vertical stabiliser
%(5): nacelle - note, there are 2 to consider (2 engines!)
d_components=[d_fuselage,b,1,1,1];

aero_analysis.drag.l_components=[1,2,3,4,5]; %define characteristic length of each component
aero_analysis.drag.fineness=zeros(1,length(d_components));
for i=1:length(d_components)
    aero_analysis.drag.fineness(i)=aero_analysis.drag.l_components(i)/d_components(i);
end

aero_analysis.drag.FF=zeros(1,length(d_components));
aero_analysis.drag.FF(1)=1+(60/aero_analysis.drag.fineness(1)^3)+(aero_analysis.drag.fineness(1)/400);
max_thickness=[0.3,0.4,0.5];
aero_analysis.drag.FF(2)=1;
aero_analysis.drag.FF(3)=1;
aero_analysis.drag.FF(4)=1;
aero_analysis.drag.FF(5)=1+(0.35/aero_analysis.drag.fineness(5));

Q=[1,1.5,1.2,1.2,1.5]; %exact values to be added

S_wet=[1,1,1,1,1]; %wetted area, exact values to be added

aero_analysis.drag.Re_cutoff=38.21.*(aero_analysis.drag.l_components*3.28084/2.08E-05).^(1.053);
%note: must convert lengths to ft for the cutoff Re

%conditions
%(1):cruise
%(2):TO
%(3):landing
U=[300,100,100]; %velocity for consitions
rho=[1,2,3]; %density for conditions
nu=[1,1,1]; %dynamic viscosity for conditions
M=[1,1,1]; %Mach number

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
aero_analysis.drag.leakage=1; %find a method
%% Induced drag
