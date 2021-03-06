clear 
clc

load("tailplane.mat");
load("wing.mat");

tailplaneC_elevatoredge = tailplane.initial.horizontal.Croot + 0.9 * (tailplane.initial.horizontal.Ctip - tailplane.initial.horizontal.Croot);
control_surface.elevator_area = 0.32 * (tailplane.initial.horizontal.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.initial.horizontal.b; %note 2 elevators

rudderC_rudderedge = tailplane.initial.vertical.Croot + 0.9 * (tailplane.initial.vertical.Ctip - tailplane.initial.vertical.Croot);
control_surface.rudder_area =  0.5 * 0.3 * (tailplane.initial.vertical.Croot  + tailplaneC_elevatoredge) * 0.9 * 0.5 * tailplane.initial.vertical.b; 

aileronC_inneredge = wing.Croot + 0.7 * (wing.Ctip - wing.Croot);
aileronC_outeredge = wing.Croot + 0.95 * (wing.Ctip - wing.Croot);
control_surface.aileron_area =  (1 - wing.HDL_PERC) * (aileronC_inneredge  + aileronC_outeredge) * 0.25 * 0.5 * wing.b; 

flapC_inneredge = wing.Croot + (wing.HDL_start / (wing.b / 2)) * (wing.Ctip - wing.Croot);
flapC_outeredge = wing.Croot + 0.7 * (wing.Ctip - wing.Croot);
control_surface.flap_area =  (1 - wing.HDL_PERC) * (flapC_inneredge  + flapC_outeredge) * (0.7 - (wing.HDL_start / (wing.b / 2))) * 0.5 * wing.b; 

control_surface.total_area_ft2 = (control_surface.elevator_area + control_surface.rudder_area + control_surface.aileron_area + control_surface.flap_area) * 3.28084^2;
control_surface.wingmounted_area_ft2 = (control_surface.aileron_area + control_surface.flap_area) * 3.28084^2;
control_surface.elevator_area_ft2 = control_surface.elevator_area * 3.28084^2;

save('control_surface','control_surface')
