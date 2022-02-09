%% Script in order to optimise the aircraft wing-box skin-stringer panels
% AVD
% 8th Feb 2022

% ieo18

% Note: paper uses inches and lbs.
% URL: https://www.sciencedirect.com/science/article/pii/S1270963809001072?via%3Dihub
%% Inputs, to initialise the script

t_s=0.1; % thickness of the skin [inches]
b_s=5; %panel width/ stringer pitch [inches]
A_sk=t_s*b_s; %area of skin [inches^2]
SR=0.5; %Stiffening ratio. Usually set to 0.5 for preliminary designs
A_st=SR*A_sk; %area of the stringer [inches^2]
A=A_st+A_sk; %area of the panel (=skin+stringer)
BR=0.5; %NEED TO DEFINE, width ratio (=bs/be)
b_e=b_s/BR;
equal_flang=0; %0: Z-type, false; 1: Z-type, true; 2: J-type; equal flanges
integrally_machined=0; %0: false; 1: true
L=4;

if integrally_machined==0
    if t_s<=0.3
        b_a=2.08*t_s+0.68;
    else
        b_a=1.312;
    end
    t_a=0.7*t_s;
end

if integrally_machined==1
    b_a=0;
    t_a=0;
end



if equal_flang==0
    b_w=sqrt((b_e/t_s)*((A_st-0.7*b_a*t_s)/1.327));
elseif equal_flang==1
    b_w=sqrt((b_e/t_s)*(A_st-1.4*b_a*t_s));
elseif equal_flang==2
    b_w=sqrt((b_e/t_s)*((A_st-2*b_a*t_a)/1.327));  
end
        

t_w=b_w*(t_s/b_e);
b_f=0.327*b_w;
t_f=t_w;

mass_1=(A_sk*rho_sk
