clear 
clc

%Asymetric load case with one engine inoperative assuming engines at max
%thrust upon failure for worst case scenario

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_ac_v = 9.7224; %m vertical stab xac
Y_engine = 1.1; %m
T_1eng_max = 4950; %installed sea level thrust of engine (assuming Princess didnt fuck this)

%Directional stability in OEI
F_v = (T_1eng_max * Y_engine) / (X_ac_v - X_cg); %force on vertical tail plane





