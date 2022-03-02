clear 
clc 
close all

%Note we assume that most of the load is taken by the skin and spar web,
%and the flange doesnt take much load

%%WING
%% Calculating web thickness

load('limiting_loadcase_distributions');
load('materialLib.mat');
load('VerticalTail.mat');
load('HorizontalTail.mat');
load('wing_layout.mat');

c_root = 1.767;
c_tip = 0.533;
span = 8.94;
semispan = span/2;

chord = @(y) -(c_root - c_tip)/semispan*y + c_root; %chord distribution

coords = [1.00000     0.00000;
0.95033     0.00986;
0.90066     0.01979;
0.85090     0.02974;
0.80103     0.03935;
0.75107     0.04847;
0.70101     0.05686;
0.65086     0.06440;
0.60064     0.07085;
0.55035     0.07602;
0.50000     0.07963;
0.44962     0.08139;
0.39923     0.08139;
0.34884     0.07971;
0.29846     0.07658;
0.24881     0.07193;
0.19781     0.06562;
0.14757     0.05741;
0.09746     0.04672;
0.07247     0.04010;
0.04757     0.03227;
0.02283     0.02234;
0.01059     0.01588;
0.00580     0.01236;
0.00347     0.01010;
0.00000     0.00000;
0.00653     -0.00810;
0.00920     -0.00956;
0.01441     -0.01160;
0.02717     -0.01490;
0.05243     -0.01963;
0.07753     -0.02314;
0.10254     -0.02604;
0.15243     -0.03049;
0.20219     -0.03378;
0.25189     -0.03613;
0.30154     -0.03770;
0.35116     -0.03851;
0.40077     -0.03855;
0.45038     -0.03759;
0.50000     -0.03551;
0.54965     -0.03222;
0.59936     -0.02801;
0.64914     -0.02320;
0.69899     -0.01798;
0.74893     -0.01267;
0.79897     -0.00751;
0.84910     -0.00282;
0.89934     0.00089;
0.94967     0.00278;
1.00000     0.00000];

aerofoiltop = fit(coords(1:25,1), coords(1:25,2), 'cubicspline');
%pointstop = polyval(aerofoiltop, coords(:,1));
aerofoilbottom = fit(coords(26:51,1), coords(26:51,2), 'cubicspline');
%pointsbottom = polyval(aerofoilbottom, coords(:,1));
figure
plot(aerofoiltop, coords(1:25,1), coords(1:25,2))
hold on
plot(aerofoilbottom, coords(26:51,1), coords(26:51,2))
hold off
axis equal

%% 
y = linspace(0,(span/2),1000);
%a = %web panel spacing design
c = chord(y) * (0.77-0.1); %Distance between spars
b2front = (aerofoiltop(0.1) - aerofoilbottom(0.1)).*chord(y);
b2rear = (aerofoiltop(0.77) - aerofoilbottom(0.77)).*chord(y);
Ks = 8.1;%lookup from graph
E_sp = materialLib{3}.E;%Youngs modulus of spar web 

Shear_Load = limiting_loadcase_distributions.shear; %N
Torque_Load = limiting_loadcase_distributions.torque; %Nm

Shearflow_Shearfront = Shear_Load ./ (2 .* b2front); %N/m
Shearflow_Shearrear = Shear_Load ./ (2 .* b2rear); %N/mm note here we assume that each spar takes half of vertical load despite having different heights

Shearflow_Torque = Torque_Load ./ (0.5 * c .* (b2front + b2rear) * 2); %T/2A

Shearflow_Fspar = Shearflow_Shearfront + Shearflow_Torque;
Shearflow_Rspar = abs(Shearflow_Shearrear - Shearflow_Torque);

wing.t_fs = (Shearflow_Fspar .* b2front ./ (Ks * E_sp)).^(1/3);
wing.t_rs = (Shearflow_Rspar .* b2rear ./ (Ks * E_sp)).^(1/3); 

figure
plot(y, wing.t_fs)
hold on
plot(y,wing.t_rs)
hold off

Shearstress_Fs = Shearflow_Fspar / wing.t_fs;
Shearstress_Rs = Shearflow_Rspar / wing.t_rs; 

ShearflowTfunc = fit(y', Shearflow_Torque', 'cubicspline'); 


%% Shear flow in the skin
load('wing_layout.mat');
positions = wing_layout.rib_array - 0.1;
Shearflow_Torque2 = ShearflowTfunc(positions);
E_sk = materialLib{1}.E; %N/m^2 (71-74.6)
t_skin = wing_layout.panel_thickness_eff';
Shearflow_skin = (Shearflow_Torque2);
Shearstress_skin = Shearflow_Torque2 ./ (t_skin);
stringer_pitch = wing_layout.stringer_pitch; 
Stresscrit = 8.1*E_sk*(t_skin./stringer_pitch).^2; %Shear buckling stress (note: double check why this is)
Trescastress = materialLib{4}.shear_yield;%Tresca shear yielding stress
sigma_0 = wing_layout.sigma_0;%stringer panel critical buckling stress
sigma_crit = wing_layout.sigma_crit; %stringer panel sigma crit

Rc = sigma_0 / sigma_crit; %compressive stress ratio
Rs = Shearstress_skin / Trescastress; 

criteria_sc = Rs.^2 + Rc;

if criteria_sc > 0.99
    disp("Bad design (Criteria not satisfied)")
end

%% HORIZONTAL TAIL
%% Calculating web thickness
c_root = 0.8166;
c_tip = 0.4083;
span = 2.3887;
semispan = span/2;

chord = @(y) -(c_root - c_tip)/semispan*y + c_root; %chord distribution

coordsH = [0.000000  0.000000;
  0.005000  0.009780;
  0.007500  0.011790;
  0.012500  0.014900;
  0.025000  0.020350;
  0.050000  0.028100;
  0.075000  0.033940;
  0.100000  0.038710;
  0.150000  0.046200;
  0.200000  0.051730;
  0.250000  0.055760;
  0.300000  0.058440;
  0.350000  0.059780;
  0.400000  0.059810;
  0.450000  0.057980;
  0.500000  0.054800;
  0.550000  0.050560;
  0.600000  0.045480;
  0.650000  0.039740;
  0.700000  0.033500;
  0.750000  0.026950;
  0.800000  0.020290;
  0.850000  0.013820;
  0.900000  0.007860;
  0.950000  0.002880;
  1.000000  0.000000;
  0.000000  0.000000;
  0.005000 -0.009780;
  0.007500 -0.011790;
  0.012500 -0.014900;
  0.025000 -0.020350;
  0.050000 -0.028100;
  0.075000 -0.033940;
  0.100000 -0.038710;
  0.150000 -0.046200;
  0.200000 -0.051730;
  0.250000 -0.055760;
  0.300000 -0.058440;
  0.350000 -0.059780;
  0.400000 -0.059810;
  0.450000 -0.057980;
  0.500000 -0.054800;
  0.550000 -0.050560;
  0.600000 -0.045480;
  0.650000 -0.039740;
  0.700000 -0.033500;
  0.750000 -0.026950;
  0.800000 -0.020290;
  0.850000 -0.013820;
  0.900000 -0.007860;
  0.950000 -0.002880;
  1.000000  0.000000];

aerofoiltop = fit(coordsH(1:26,1), coordsH(1:26,2), 'cubicspline');
%pointstopH = polyval(aerofoiltop, coordsH(1:26,1));
aerofoilbottom = fit(coordsH(27:52,1), coordsH(27:52,2), 'cubicspline');
%pointsbottomH = polyval(aerofoilbottom, coordsH(1:26,1));
figure
plot(aerofoiltop, coordsH(1:26,1), coordsH(1:26,2))
hold on
plot(aerofoilbottom, coordsH(27:52,1), coordsH(27:52,2))
hold off
axis equal

%% 
s_h = 2.3887;
y = linspace(-s_h/2,s_h/2,100);
%a = %web panel spacing design
c = chord(y) * (0.68-0.15); %Distance between spars
b2front = (aerofoiltop(0.15) - aerofoilbottom(0.15)).*chord(y);
b2rear = (aerofoiltop(0.68) - aerofoilbottom(0.68)).*chord(y);
Ks = 8.1;%lookup from graph
E_sp = materialLib{3}.E;%Youngs modulus of spar web 

Shear_Load = HorizontalTail.Shearforce; %N
Torque_Load = HorizontalTail.Torque; %Nm

Shearflow_Shearfront = Shear_Load ./ (2 .* b2front); %N/m
Shearflow_Shearrear = Shear_Load ./ (2 .* b2rear); %N/mm note here we assume that each spar takes half of vertical load despite having different heights

Shearflow_Torque = Torque_Load ./ (0.5 * c .* (b2front + b2rear) * 2); %T/2A

Shearflow_Fspar = Shearflow_Shearfront + Shearflow_Torque;
Shearflow_Rspar = abs(Shearflow_Shearrear - Shearflow_Torque);

HorizontalTail.t_fs = (Shearflow_Fspar .* b2front ./ (Ks * E_sp)).^(1/3);
HorizontalTail.t_rs = (Shearflow_Rspar .* b2rear ./ (Ks * E_sp)).^(1/3); 

figure
plot(y, HorizontalTail.t_fs)
hold on
plot(y,HorizontalTail.t_rs)
hold off

Shearstress_Fs = Shearflow_Fspar / HorizontalTail.t_fs;
Shearstress_Rs = Shearflow_Rspar / HorizontalTail.t_rs; 

%% Shear flow in the skin
% E_sk = 
% t_skin = 
% Shearflow_skin = Shearflow_Torque;
% Shearstress_skin = Shearflow_Torque ./ t_skin; 
% stringer_pitch = 
% Stresscrit = 8.1*E*(t_skin/stringer_pitch)^2; %Shear buckling stress (note: double check why this is)
% Trescastress = %Tresca shear yielding stress
% sigma_0 = %stringer panel critical buckling stress
% sigma_crit = %stringer panel sigma crit
% 
% Rc = sigma_0 / sigma_crit; %compressive stress ratio
% Rs = Shearstress_skin / Trescastress; 
% 
% criteria_sc = Rs^2 + Rc;
% 
% if criteria_sc <= 0.99
%     disp("Good design (Criteria satisfied)")
% else 
%     disp("Bad design (Criteria not satisfied)")
% end

% %% VERTICAL TAIL
% %% Calculating web thickness
% 
% load('VerticalTail.mat');
% 
% c_root = 1.0015;
% c_tip = 0.8166;
% span = 1.1818;
% 
% chord = @(y) -(c_root - c_tip)/span*y + c_root; %chord distribution
% 
% coordsV = [0.000000  0.000000;
%   0.005000  0.012080;
%   0.007500  0.014560;
%   0.012500  0.018420;
%   0.025000  0.025280;
%   0.050000  0.035040;
%   0.075000  0.042400;
%   0.100000  0.048420;
%   0.150000  0.057850;
%   0.200000  0.064800;
%   0.250000  0.069850;
%   0.300000  0.073190;
%   0.350000  0.074820;
%   0.400000  0.074730;
%   0.450000  0.072240;
%   0.500000  0.068100;
%   0.550000  0.062660;
%   0.600000  0.056200;
%   0.650000  0.048950;
%   0.700000  0.041130;
%   0.750000  0.032960;
%   0.800000  0.024720;
%   0.850000  0.016770;
%   0.900000  0.009500;
%   0.950000  0.003460;
%   1.000000  0.000000;
%   0.000000  0.000000;
%   0.005000  -.012080;
%   0.007500  -.014560;
%   0.012500  -.018420;
%   0.025000  -.025280;
%   0.050000  -.035040;
%   0.075000  -.042400;
%   0.100000  -.048420;
%   0.150000  -.057850;
%   0.200000  -.064800;
%   0.250000  -.069850;
%   0.300000  -.073190;
%   0.350000  -.074820;
%   0.400000  -.074730;
%   0.450000  -.072240;
%   0.500000  -.068100;
%   0.550000  -.062660;
%   0.600000  -.056200;
%   0.650000  -.048950;
%   0.700000  -.041130;
%   0.750000  -.032960;
%   0.800000  -.024720;
%   0.850000  -.016770;
%   0.900000  -.009500;
%   0.950000  -.003460;
%   1.000000  0.000000];
% 
% aerofoiltop = polyfit(coordsV(1:26,1), coordsV(1:26,2), 20);
% pointstop = polyval(aerofoiltop, coordsV(1:26,1));
% aerofoilbottom = polyfit(coordsV(27:52,1), coordsV(27:52,2), 20);
% pointsbottom = polyval(aerofoilbottom, coordsV(1:26,1));
% figure
% plot(coordsV(1:26,1), pointstop)
% hold on
% plot(coordsV(1:26,1), pointsbottom)
% hold off
% axis equal
% 
% %% 
% z = linspace(0,s_v,1000);
% %a = %web panel spacing design
% c = chord(z) * (0.7-0.15); %Distance between spars
% b2front = (polyval(aerofoiltop, 0.15) - polyval(aerofoilbottom, 0.15)).*chord(y);
% b2rear = (polyval(aerofoiltop, 0.7) - polyval(aerofoilbottom, 0.7)).*chord(y);
% Ks = 8.1;%lookup from graph
% E_sp = materialLib{3}.E;%Youngs modulus of spar web 
% 
% Shear_Load = VerticalTail.ShearForce_vt(z); %N
% Torque_Load = VerticalTail.Torque; %Nm
% 
% Shearflow_Shearfront = Shear_Load ./ (2 .* b2front); %N/m
% Shearflow_Shearrear = Shear_Load ./ (2 .* b2rear); %N/mm note here we assume that each spar takes half of vertical load despite having different heights
% 
% Shearflow_Torque = Torque_Load ./ (0.5 * c .* (b2front + b2rear) * 2); %T/2A
% 
% Shearflow_Fspar = Shearflow_Shearfront + Shearflow_Torque;
% Shearflow_Rspar = Shearflow_Shearrear - Shearflow_Torque;
% 
% VerticalTail.t_fs = (Shearflow_Fspar .* b2front ./ (Ks * E_sp)).^(1/3);
% VerticalTail.t_rs = (Shearflow_Rspar .* b2rear ./ (Ks * E_sp)).^(1/3); 
% 
% figure
% plot(z, VerticalTail.t_fs)
% hold on
% plot(z,VerticalTail.t_rs)
% hold off
% 
% Shearstress_Fs = Shearflow_Fspar / VerticalTail.t_fs;
% Shearstress_Rs = Shearflow_Rspar / VerticalTail.t_rs; 
% 
% %% Shear flow in the skin
% % E_sk = 
% % t_skin = 
% % Shearflow_skin = Shearflow_Torque;
% % Shearstress_skin = Shearflow_Torque ./ t_skin; 
% % stringer_pitch = 
% % Stresscrit = 8.1*E*(t_skin/stringer_pitch)^2; %Shear buckling stress (note: double check why this is)
% % Trescastress = %Tresca shear yielding stress
% % sigma_0 = %stringer panel critical buckling stress
% % sigma_crit = %stringer panel sigma crit
% % 
% % Rc = sigma_0 / sigma_crit; %compressive stress ratio
% % Rs = Shearstress_skin / Trescastress; 
% % 
% % criteria_sc = Rs^2 + Rc;
% % 
% % if criteria_sc <= 0.99
% %     disp("Good design (Criteria satisfied)")
% % else 
% %     disp("Bad design (Criteria not satisfied)")
% % end

save("wing.mat", 'wing');
save("HorizontalTail.mat", 'HorizontalTail');

