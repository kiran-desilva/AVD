clear

load('sizing.mat');
load('parameters.mat');



%cessana citation mustang currently
fuselage_length = 12;
cmac = 1.48;
sref = 19.51;
b = 13.16;
ar = sizing.AR;
% wing_lambdaLE = sizing.sweepLE;
wing_sweepLE = 22;


tailplane.initial.horizontal.l = 0.45*fuselage_length;
tailplane.initial.vertical.l = 0.45*fuselage_length;

%correction for T-tail endplate effect 
correction_factor = .95; 

%from aero tn tailsizing pdf
tailplane.initial.horizontal.v = 0.7 * correction_factor; % howe
tailplane.initial.vertical.v = 0.065 * correction_factor; 

tailplane.initial.horizontal.s = ( tailplane.initial.horizontal.v * cmac * sref ) / (tailplane.initial.horizontal.l);

tailplane.initial.vertical.s = ( tailplane.initial.vertical.v * b * sref ) / (tailplane.initial.vertical.l);

%from Obert
tailplane.initial.horizontal.Ar = 3.9; %this is probably too big
tailplane.initial.vertical.Ar = .9; 

% tailplane.initial.bh = sqrt(tailplane.initial.arh * tailplane.initial.sh);
% tailplane.initial.bv = sqrt(tailplane.inital.arv * tailplane.initial.sv);

%taper
tailplane.initial.horizontal.lambda = 0.35; % obert -> looked good to me idk
tailplane.initial.vertical.lambda = 0.7; %obert -> cant be too small as the vertical stabilizer has to support horizotnal tailplane

tailplane.initial.horizontal.sweepLE = wing_sweepLE + 5; 
tailplane.initial.vertical.sweepLE = 40; 

tailplane.initial.horizontal.sweep_25 = sweep_angle(tailplane.initial.horizontal.sweepLE,25,0,tailplane.initial.horizontal.Ar,tailplane.initial.horizontal.lambda)
tailplane.initial.vertical.sweep_25 = sweep_angle(tailplane.initial.vertical.sweepLE,25,0,tailplane.initial.vertical.Ar,tailplane.initial.vertical.lambda)

tailplane.initial.horizontal.tc = 0.12; %64012
tailplane.initial.vertical.tc = 0.12;

tailplane.initial.horizontal = wing_geometry_calc_struct(tailplane.initial.horizontal);
tailplane.initial.vertical = wing_geometry_calc_struct(tailplane.initial.vertical);

horizontal_stab_plot(tailplane.initial.horizontal)
vertical_stab_plot(tailplane.initial.vertical)

save("tailplane","tailplane")

