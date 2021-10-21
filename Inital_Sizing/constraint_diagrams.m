
open("parameters.mat");

%% shit we want
%% l/d max - l_dmax
%% specific fuel consumption - sfc
%% gross weight - w0
%% operating empty weight - we
%% fuel weight - wf

sizing.ld_max = 10;
sizing.ld_cruise = 0.866*sizing.ld_max;
sizing.sfc = 10;
sizing.maxTakeoffWeight = 1000;

%% shit from errikos' excel
sizing.cl_max = nan;
sizing.AR = nan;
sizing.cd_

%% What's next?
%% Make sure aircraft can complete the stipulated design mission profile
%% Make sure aircraft capable of achieving performance targets like
%% Ceilings (absolute, service, combat) -> T/W at height, v has to be bigger than 1/(L/D)_max
%% Maximum speed
%% Time to climb / Rates of climb
%% Sustained turn rates / radii
%% Level (axial) acceleration
%% Takeoff and Landing distances (TODA & LDA) -> T/W, stall speed 
%% Stall speed -> function of Cl_max (Cl_max changes in takeoff and clean config) and wing loading (W/ S)
%% Make sure aircraft meets airworthiness requirements

%% FAR25 - https://www.engineerstoolkit.com/Airworthiness%20Standards%20%20FAA%20FAR%20Part%2025.pdf
%% take off speed = 1.1* stall speed