clear 
clc 
close all

%% Calculating web thickness

load('limiting_loadcase_distributions');

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

aerofoiltop = polyfit(coords(1:25,1), coords(1:25,2), 20);
pointstop = polyval(aerofoiltop, coords(:,1));
aerofoilbottom = polyfit(coords(26:51,1), coords(26:51,2), 20);
pointsbottom = polyval(aerofoilbottom, coords(:,1));
figure
plot(coords(:,1), pointstop)
hold on
plot(coords(:,1), pointsbottom)
hold off
axis equal

%% 
y = linspace(0,(span/2),1000);
%a = %web panel spacing design
c = chord(y) * (0.77-0.1); %Distance between spars
b2front = polyval(aerofoiltop, chord(y)*0.1) - polyval(aerofoilbottom, chord(y)*0.1);
b2rear = polyval(aerofoiltop, chord(y)*0.77) - polyval(aerofoilbottom, chord(y)*0.77);
Ks = 8.1;%lookup from graph
E_sp = 70e9;%Youngs modulus of spar web 

Shear_Load = limiting_loadcase_distributions.shear; %N
Torque_Load = limiting_loadcase_distributions.torque; %Nm

Shearflow_Shearfront = Shear_Load ./ (2 .* b2front); %N/m
Shearflow_Shearrear = Shear_Load ./ (2 .* b2rear); %N/mm note here we assume that each spar takes half of vertical load despite having different heights

Shearflow_Torque = Torque_Load ./ (0.5 * c * (b2front + b2rear) * 2); %T/2A

Shearflow_Fspar = Shearflow_Shearfront + Shearflow_Torque;
Shearflow_Rspar = Shearflow_Shearrear - Shearflow_Torque;

t_fs = (Shearflow_Fspar .* b2front ./ (Ks * E_sp)).^(1/3);
t_rs = (Shearflow_Rspar .* b2rear ./ (Ks * E_sp)).^(1/3); 

figure
plot(y, t_fs)
hold on
plot(y,t_rs)
hold off

Shearstress_Fs = ShearflowFspar / t_fs;
Shearstress_Rs = ShearflowRspar / t_rs; 

%% Shear flow in the skin
E_sk = 
t_skin = 
Shearflow_skin = Shearflow_Torque;
Shearstress_skin = Shearflow_Torque ./ t_skin; 
stringer_pitch = 
Stresscrit = 8.1*E*(t_skin/stringer_pitch)^2; %Shear buckling stress (note: double check why this is)
Trescastress = %Tresca shear yielding stress
sigma_0 = %stringer panel critical buckling stress
sigma_crit = %stringer panel sigma crit

Rc = sigma_0 / sigma_crit; %compressive stress ratio
Rs = Shearstress_skin / Trescastress; 

criteria_sc = Rs^2 + Rc;

if criteria_sc <= 0.99
    disp("Good design (Criteria satisfied)")
else 
    disp("Bad design (Criteria not satisfied)")
end

