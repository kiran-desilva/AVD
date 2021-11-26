clear 
clc

load("tailplane.mat");

tailplaneC_elevatoredge = tailplane.initial.horizontal.Croot + 0.9 * (tailplane.initial.horizontal.Ctip - tailplane.initial.horizontal.Croot);
control_surface.elevator_area = 0.32 * (tailplane.initial.horizontal.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.initial.horizontal.b; %note 2 elevators

rudderC_rudderedge = tailplane.initial.vertical.Croot + 0.9 * (tailplane.initial.vertical.Ctip - tailplane.initial.vertical.Croot);
control_surface.rudder_area =  0.5 * 3 * (tailplane.initial.vertical.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.initial.vertical.b; 

aileronC_inner = tailplane.initial.vertical.Croot + 0.9 * (tailplane.initial.vertical.Ctip - tailplane.initial.vertical.Croot);
control_surface.rudder_area =  0.5 * 3 * (tailplane.initial.vertical.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.initial.vertical.b; 
