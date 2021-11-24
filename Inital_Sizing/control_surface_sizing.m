clear 
clc

load("tailpane.mat");

tailplaneC_elevatoredge = tailplane.Croot + 0.9 * (tailplane.Ctip - tailplane.Croot);
elevator_area = 0.5 * 0.32 * (tailplane.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.b;

