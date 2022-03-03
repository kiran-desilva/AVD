%% D section of wing
% not the most efficient script in the world but seems the best way to
% input the K_s value manually in a short period of time. Not sure how
% efficient it is to digitise the plots.
clear
clc

%in PreliminaryAirframe file
load b_a_int.mat
load a_b_int.mat
load wing_layout.mat
load ribs_check.mat
load materialLib.mat
%% Wing box parameters
%b_ref=3.63; %wingbox span [m] -NEED TO VERIFY
a=ribs_check.s; %d-cell span length
b=0.2.*ribs_check.chord; %curved panel length - proportion of x/c. [NEED TO ESTABLISH THE PROPORTION]
t1=wing_layout.panel_thickness; %d-cell thickness 
%R1= 0.01.*ribs_check.chord; %D-cell radius at position of rib (taking minimum value of radius (1% of x/c). Could take the average?  
R1=ones(1,length(a));
frac_A1=b./a; %b/a if b>a
frac_B1=a./sqrt(R1.*t1); %if b>a
frac_A2=a./b; %a/b if a>b
frac_B2=b./sqrt(R1.*t1); %if a>b
% K_S=zeros(1,length(a));
K_S=[9,10,10.5,11,12,13.5, 13.5, 13, 12,10.5,9.5,8.5,7.5,7,7];
% for i=1:length(a)
%     K_S(i)=d_cell_k_s(a(i),b(i),R1(i),t1(i));
% end
E_panel=materialLib{1, 1}.E; %Young's modulus of Panel material
tau_cr=K_S.*E_panel.*(t1./b).^2; %buckling stress of curved panel in shear
tresca=145e+06; %Tresca shear yielding stress is the same as the shear yield strength [Pa]
%t_tresca=
%save ('d_section.mat', 'tau_cr','k_s','b','t1','E_panel')


