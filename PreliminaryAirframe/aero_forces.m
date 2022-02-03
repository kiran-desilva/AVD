%% AVD Script
% steps:
% - extract data from .txt file
% - calculate Lift distribution
% - calculate shear force distribution
% - calculate bending moment distribution


%% Housekeeping
clear
clc
close all

%% CRUISE

a=importfile("C:\Users\izzye\OneDrive\Documents\GitHub\AVD\PreliminaryAirframe\AVL Data Files\AVD_AVL_CRUISE_FINAL.dat", [21, 84]);


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
cruise.rho=1.225;
cruise.air_vel=340.3;
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

disp cruise.total_lift