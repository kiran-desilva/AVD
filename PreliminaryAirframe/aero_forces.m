%% AVD Script
% steps:
% - extract data from .txt file
% - calculate Lift distribution
% - calculate shear force distribution
% - calculate bending moment distribution

%% Outputs to save
%cruise.cl_scaled: the scaled Cl for each strip section which produces the
%overall CL required for the aircraft in cruise (0.2689)

%p1: 64th order polynomial terms for Cl distribution. used air properties
%at 40 000 ft for air in order to obtain the necessary Cl.

%% Housekeeping
clear
clc
close all
load sizing.mat

%% CRUISE

a=importfile("AVL Data Files/AVL_10th_feb.dat", [21, 84]);


% Parameters from AVL
cruise.j=table2array(a(:,1));
cruise.y_le=table2array(a(:,4));
cruise.chord=table2array(a(:,7));
cruise.area=table2array(a(:,10));
cruise.c_cl=table2array(a(:,13));
cruise.ai=table2array(a(:,16));
cruise.cl_norm=table2array(a(:,19));
cruise.cl=table2array(a(:,22));
cruise.cd=table2array(a(:,25))+table2array(a(:,24));
cruise.cdv=table2array(a(:,28));
cruise.cm_c_4=table2array(a(:,29))+table2array(a(:,30));
cruise.cm_LE=table2array(a(:,32))+table2array(a(:,31));
cruise.c_p_x_c=table2array(a(:,35))+table2array(a(:,36));


% Flight conditions
cruise.Mach=0.75;
%cruise.rho=0.302; 
%cruise.air_vel=294.9; 
cruise.rho=1.225; %sea level densty, SI
cruise.air_vel=340.3; %sea level m/s, SI
cruise.vel=cruise.air_vel*cruise.Mach; %m/s
cruise.q=0.5*cruise.rho*cruise.vel^2; %dynamic pressure, SI

%% Lift
cruise.sectional_lift=cruise.q.*cruise.area.*cruise.cl; % lift of each segment
figure(1)
plot(cruise.y_le,cruise.sectional_lift, 'k', 'Linewidth', 1.25)
xlabel 'Spanwise displacement [m]'
ylabel 'Lift distribution [N]'
grid on
grid minor
cruise.lift_no_dynamic_p=cruise.area.*cruise.cl;
cruise.total_lift=sum(cruise.sectional_lift); % sum the lift of each segment

%% Shear
t=(length(cruise.j):-1:1);
cruise.shear=zeros(length(t),1);
cruise.shear(cruise.j,1)=cruise.sectional_lift(cruise.j);
cruise.shear(length(cruise.j),1)=cruise.sectional_lift(length(cruise.j),1);
cruise.shear(length(cruise.j)-1,1)=cruise.shear(length(cruise.j),1)+cruise.sectional_lift((length(cruise.j)-1),1);

for i=1:length(t)-2
    cruise.shear(length(t)-(i+1),1)=cruise.sectional_lift(length(t)-(i+1),1)+cruise.shear(length(t)-(i),1);
end
figure(2)
plot(cruise.y_le, cruise.shear, 'k', 'Linewidth', 1.25)
xlabel 'Spanwise displacement [m]'
ylabel 'Shear Force [N]'
grid on
grid minor

%% Bending Moment
cruise.sectional_moment=cruise.cm_c_4.*cruise.area.*cruise.chord*cruise.q.*cruise.chord;
cruise.moment(length(cruise.j),1)=cruise.sectional_moment(length(cruise.j),1);
cruise.moment(length(cruise.j)-1,1)=cruise.moment(length(cruise.j),1)+cruise.sectional_moment((length(cruise.j)-1),1);
for i=1:length(t)-2
    cruise.moment(length(t)-(i+1),1)=cruise.sectional_moment(length(t)-(i+1),1)+cruise.moment(length(t)-(i),1);
end

% %% Bending moment V2
% 
% 
% for i=1:length(t)-2
%     moment(length(t)-(i+1),1)=(
% end

figure(3)
plot(cruise.y_le,cruise.moment, 'k', 'Linewidth', 1.25)
xlabel 'Spanwise displacement [m]'
ylabel 'Bending moment [Nm]'
grid on
grid minor

%%
% [fitobject1]=fit(y_LE,p1,'fourier8');
% a0 = fitobject1.a0;
% a1 = fitobject1.a1;
% b1 = fitobject1.b1;
% a2 = fitobject1.a2;
% b2 = fitobject1.b2;
% a3 = fitobject1.a3;
% b3= fitobject1.b3;
% a4 = fitobject1.a4;
% b4= fitobject1.b4;
% a5 = fitobject1.a5;
% b5= fitobject1.b5;
% a6 = fitobject1.a6;
% b6= fitobject1.b6;
% a7 = fitobject1.a7;
% b7= fitobject1.b7;
% a8 = fitobject1.a8;
% b8= fitobject1.b8;
% w  = fitobject1.w;

%[fitresult,gof]=fitobject1(a0,a1,b1)

%% normalise lift coefficient
cruise.lift_no_dynamic_p=cruise.area.*cruise.cl; %lift coefficient*area for each strip. Lift/ dynamic pressure. Doesn't account for dynamic pressure
z=sum(cruise.lift_no_dynamic_p); %L/(dynamic pressure); for whole aircraft
cruise.total_area=sum(cruise.area); %total area of wing, verification that it is 4.1745 m^2
cruise.cl_total_dim=z/cruise.total_area; %CL of wing. CL=sum(Cl(i)*S(i))/S_total
%scale cl to be equal to Cl necessary at cruise (0.2689)
cruise.cl_req=(sizing.W0*sizing.fraction.before_cruise*0.5)/(0.5*0.302*294.9^2*cruise.total_area); %calculate the required CL (=W*0.5/(q*s)); using conditions of air at cruise; 0.5 since for 1 wing
cruise.ratio_cl=cruise.cl_req/cruise.cl_total_dim; %scaling factor
%% mutliply lift coefficient distribution by necessary factor

cruise.cl_scaled=cruise.cl*cruise.ratio_cl;
cl_fit.poly=polyfit(cruise.y_le,cruise.cl_scaled,13);

x = linspace(min(cruise.y_le), max(cruise.y_le), 1000);
figure;
hold on
plot(x, polyval(cl_fit.poly, x))
plot(cruise.y_le, cruise.cl_scaled)
hold off


save('cl_fit')
% test=[0:0.001:3.63]
% %plot(cruise.y_le,cruise.lift_no_dynamic_p)
% q=polyfit(cruise.y_le,cruise.lift_no_dynamic_p,64)
%                