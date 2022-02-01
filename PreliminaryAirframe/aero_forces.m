%% AVD Script
% steps:
% - extract data from .txt file
% - calculate Lift distribution
% - calculate shear force distribution
% - calculate bending moment distribution


%% Housekeeping
clear
clc

%% Extract data

a=importfile("C:\Users\izzye\OneDrive\Documents\GitHub\AVD\PreliminaryAirframe\AVL Data Files\AVL_1.dat", [21, 84])

j=a(:,1);
y_le=a(:,4);
c_cl=a(:,5);
chord=a(:,7);
area=a(:,10);