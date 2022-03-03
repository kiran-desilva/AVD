%% Code to output thickness of wing ribs

%% Housekeeping
close all
clc
clear

load wing_layout.mat
load materialLib.mat
%% Daiqing method from Excel

% % Assume rib at every 0.5m
ribs_check.ribs_loc_1=wing_layout.rib_array; %input that need to get; location of rib in spanwise direction
ribs_check.ribs_n=length(ribs_check.ribs_loc_1);

for i=1:ribs_check.ribs_n %chord at each spanwise location of rib
    ribs_check.chord(i)=feval(geometry.c,ribs_check.ribs_loc_1(i));
end

for i=2:length(ribs_check.ribs_loc_1) 
    ribs_check.s(i-1)=ribs_check.ribs_loc_1(i)-ribs_check.ribs_loc_1(i-1); %displacement between 2 ribs
end

ribs_check.s(ribs_check.ribs_n)=geometry.semispan-ribs_check.ribs_loc_1(ribs_check.ribs_n);

ribs_check.E_panel=materialLib{1, 1}.E; %Young's modulus of panel (???) [72 GPa]
ribs_check.t_s=design_params.stringer_thickness; %thickness of stringer
ribs_check.t_e=3*ribs_check.t_s/2; % effective thickness of the panel [m] - assuming stiffness ratio is 0.5

%input max M of each section to be considered for each rib
for i=1:ribs_check.ribs_n
    ribs_check.M(i)=feval(bending_moment_dist,ribs_check.ribs_loc_1(i));
    ribs_check.h_c(i)=feval(geometry.web_height_func,ribs_check.ribs_loc_1(i)); % wing box height %assume wing box is constant height - need to modify
end

ribs_check.I=(ribs_check.chord.*(ribs_check.t_e)^3/12 + ribs_check.chord.*(ribs_check.t_e).*(ribs_check.h_c./2).^2); %[m^4]
ribs_check.F=(ribs_check.M.^2.*ribs_check.s.*ribs_check.h_c.*ribs_check.t_e.*ribs_check.chord)./(2*ribs_check.E_panel.*ribs_check.I.^2);
wing_layout.rib_thickness=((ribs_check.F.*ribs_check.h_c.^2)./(3.62*ribs_check.E_panel.*ribs_check.chord)).^(1/3); %required rib thickness for optimal design

save wing_layout.mat

% sigma_y=350E+06; %yield stress
% design_t_r=F./(sigma_y.*chord);
% save ('wing_ribs.mat')
% n_ribs=length(ribs_loc_1);
% c_ribs=ones(1,n_ribs)*1.1; %chord of each rib [m]
% panels=n_ribs+1; %number of panels - assuming no rib at tip of wing
% 
% for i=1:panels
%     s_panel(i)=ribs_loc_1(i+1)-ribs_loc_1(i); %length of panels - distance between 2 consecutive panels
% end
% s_panel(panels+1)=3.63-max(ribs_loc_1); % panel after the last rib - between last rib and tip
% 
% moment_1=ones(1,n_ribs)*1000; % initalise value of moment at each rib
% 
% t_e=2e-03; % effective thickness of the panel [m]
% h_c=ones(1,n_ribs).*0.3; % Wingbox depth [m]
% 
% Young_E=72e+09; %Young's modulus [GPa]
% 
% I=(c_ribs.*(t_e)^3/12 + c_ribs.*(t_e).*(h_c./2).^2); %[m^4]
% 
% t_r=ones(1,n_ribs).*5e-03; %design rib thickness
% 
% crush_F=(moment_1.^2.*s_panel.*h_c.*t_e.*c_ribs)./(2*Young_E.*I.^2);
% crush_stress=crush_F./(t_r.*c_ribs); %[Pa]
% 
% buckling_stress_crit=3.62*Young_E.*(t_r./h_c).^2; % [Pa]
% yield_stress=350e+06; %[Pa]
% 
% new_t_r=crush_F./(yield_stress.*c_ribs);

% %% Method by Isobel - optimise placement of ribs w.r.t bending moment and resulting crush force
% 
% %% Step 1: Assume a spacing distribution of 0.5m (from the root)
% % Find the thickness of the first rib (at 0.5 m from the root)
% 
% M(1)=1000000; %max(moment)
% s(1)=0.5; %distance from previous rib [m]
% h_c(1)=0.4; %average wing box height between y=0 m and y=0.5 m
% t_e=5e-03; %effective thickness of skin panel [m] - 2mm
% c(1)=1.9; %[m]
% I=(c(1)*(t_e)^3/12 + c(1)*(t_e)*(h_c(1)/2)^2); %[m^4]
% E=72e+09; %72 Gpa
% F(1)=(M(1)^2*s(1)*h_c(1)*t_e*c)/(2*E*I^2);
% 
% %Find t_r 
% t_r=((F(1)*h_c(1)^2)/(3.62*E*c(1)))^(1/3) % this will be the thickness of each rib in the wing, for ease of manufacturing and also calculations
% 
% %% Step 2: using t_r find the spanwise location of the next rib
% M(2)=500000; % moment at point just after previous rib
% % assume that the h_c is constant - aspect of over-engineering here but
% % good for first iteration
% h_c(2)=0.35; % value of wing box height at y=0.5m
% 
% % equate two expressions for F and solve for s
% % assume c is the same as the previous rib (over-engineering) but then
% % iterate
% I(2)=(c(1)*(t_e)^3/12 + c(1)*(t_e)*(h_c(2)/2)^2); %moment of inertia, first guess
% s(2)=(7.24*E^2*I(2)^2*t_r^3)/(M(2)^2*h_c(2)^2*t_e);
% F(2)=(3.62*E*t_r^3*c(1))/(h_c(2)^2);
% diff=1; %initialise difference
% tol=0.001; % tolerance of 1mm
% %use s value to calculate c to then calculate I to then re-calculate s
% %value etc
% 
% d(1)=s(2)+s(1); % initialise spanwise position
% g(1)=s(2);
% i=2;
% while abs(diff) >= tol
%     c(2)=-0.34*d(i-1)+1.7674; % formula to find chord length at a given spanwise location (s)
%     I(2)=c(2)*(t_e)^3/12 + c(2)*(t_e)*(h_c(2)/2)^2; %new moment of inertia
%     g(i)=(7.24*E^2*I(2)^2*t_r^3)/(M(2)^2*h_c(2)^2*t_e); %iterated value of s
%     diff=g(i)-g(i-1); %convergence of s value
%     d(i)=s(1)+g(i); %new spanwise position of rib
%     i=i+1;
% end
% n=length(g);
% s(2)=g(n);
% disp(1)=s(1);
% disp(2)=d(n);
% 
% 
% %% Repeat process for all remaining ribs
% M(3)=250000; %value of M at location of s(2)
% M(4)=300;
% M(5)=40;
% M(6)=4;
% h_c(3)=0.3; %assume constant for now - need a formula to show variation of h_c
% h_c(4)=0.3;
% h_c(5)=0.3;
% h_c(6)=0.3;
% j=3;
% 
% 
% while max(disp)< 3.63
%     I(j)=(c(j-1)*(t_e)^3/12 + c(j-1)*(t_e)*(h_c(j)/2)^2); %moment of inertia, first guess
%     s(j)=(7.24*E^2*I(j)^2*t_r^3)/(M(j)^2*h_c(j)^2*t_e);
%     F(j)=(3.62*E*t_r^3*c(j-1))/(h_c(j)^2);
%     diff=1; %reinitialise diff
%     i=2;
%     clear d
%     while abs(diff) >= tol
%         c(j)=-0.34*sum(s)+1.7674; % formula to find chord length at a given spanwise location (s)
%         I(j)=c(j)*(t_e)^3/12 + c(j)*(t_e)*(h_c(j)/2)^2;
%         g(i)=(7.24*E^2*I(j)^2*t_r^3)/(M(j)^2*h_c(j)^2*t_e);
%         f=sum(s)
%         d(i)=f+g(i); %update spanwise position of rib
%         diff=g(i)-g(i-1);
%         i=i+1;
%         s(j)=g(i-1)
%     end
%     c(j)=-0.34*d(length(d))+1.7674;
%     s(j)=g(length(g)) %update value of s for converged value of g
%     disp(j)=d(length(g)) %total displacement from root of latest rib
%     j=j+1
% 
% end
% 
% %     I(3)=(c(2)*(t_e)^3/12 + c(2)*(t_e)*(h_c(3)/2)^2); %moment of inertia, first guess
% %     s(3)=(7.24*E^2*I(3)^2*t_r^3)/(M(3)^2*h_c(3)^2*t_e);



    


