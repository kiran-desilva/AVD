clear

open('sizing.mat')
open('parameters.mat')

%cessana citation mustang currently
fuselage_length = 12;
cmac = 1.48;
sref = 19.51;
b = 13.16;
ar = sizing.AR;
lambdaLE = sizing.lambdaLE;

tailplane.initial.lh = 0.45*fuselage_length;
tailplane.initial.lv = 0.45*fuselage_length;

%correction for T-tail endplate effect 
correction_factor = .95; 

%from aero tn tailsizing pdf
tailplane.initial.vh = 0.7 * correction_factor; % howe
tailplane.initial.vv = 0.065 * correction_factor; 

tailplane.initial.sh = ( tailplane.initial.vh * cmac * sref ) / (tailplane.initial.lh);

tailplane.initial.sv = ( tailplane.initial.vv * b * sref ) / (tailplane.initial.lv);



save("tailplane","tailplane")

