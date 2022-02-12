%% Code to solve for optimum ribs distribution

%% Housekeeping
close all
clc
clear

%% Daiqing method from Excel

% Assume rib at every 0.5m

ribs_loc_1=[0:0.5:3.63];

n_ribs=length(ribs_loc_1);
c_ribs=ones(1,n_ribs)*1.1; %chord of each rib [m]
panels=n_ribs-1; %number of panels - assuming no rib at tip of wing

for i=1:panels
    s_panel(i)=ribs_loc_1(i+1)-ribs_loc_1(i); %length of panels - distance between 2 consecutive panels
end
s_panel(panels+1)=3.63-max(ribs_loc_1); % panel after the last rib - between last rib and tip

moment_1=ones(1,n_ribs)*1000; % initalise value of moment at each rib

t_e=2e-03; % effective thickness of the panel [m]
h_c=ones(1,n_ribs).*0.3; % Wingbox depth [m]

Young_E=72e+09; %Young's modulus [GPa]

I=(c_ribs.*(t_e)^3/12 + c_ribs.*(t_e).*(h_c./2).^2); %[m^4]

t_r=ones(1,n_ribs).*5e-03; %design rib thickness

crush_F=(moment_1.^2.*s_panel.*h_c.*t_e.*c_ribs)./(2*Young_E.*I.^2);
crush_stress=crush_F./(t_r.*c_ribs); %[Pa]

buckling_stress_crit=3.62*Young_E.*(t_r./h_c).^2; % [Pa]
yield_stress=350e+06; %[Pa]

new_t_r=crush_F./(yield_stress.*c_ribs);



