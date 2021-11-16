clear

open('sizing.mat')
open('parameters.mat')

%cessana citation mustang currently
fuselage_length = 12;
cmac = 1.48;
sref = 19.51;
b = 13.16;
ar = sizing.AR;
wing_lambdaLE = sizing.lambdaLE;


tailplane.initial.lh = 0.45*fuselage_length;
tailplane.initial.lv = 0.45*fuselage_length;

%correction for T-tail endplate effect 
correction_factor = .95; 

%from aero tn tailsizing pdf
tailplane.initial.vh = 0.7 * correction_factor; % howe
tailplane.initial.vv = 0.065 * correction_factor; 

tailplane.initial.sh = ( tailplane.initial.vh * cmac * sref ) / (tailplane.initial.lh);

tailplane.initial.sv = ( tailplane.initial.vv * b * sref ) / (tailplane.initial.lv);

%from Obert
tailplane.initial.arh = 3.9; %this is probably too big
tailplane.initial.arv = 0.9; 

tailplane.initial.bh = sqrt(tailplane.initial.arh * tailplane.initial.sh);
tailplane.initial.bv = sqrt(tailplane.inital.arv * tailplane.initial.sv);

tailplane.initial.taperh = 0.35; % obert -> looked good to me idk
tailplane.initial.taperv = 0.7; %obert -> cant be too small as the vertical stabilizer has to support horizotnal tailplane

tailplane.initial.sweepLEh = wing_sweepLE + 5; 

save("tailplane","tailplane")

