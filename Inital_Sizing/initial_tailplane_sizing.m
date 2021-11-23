clear

load('sizing.mat');
load('parameters.mat');
load('design.mat')



%cessana citation mustang currently
fuselage_length = 10;
cmac = 1.4;
sref = design.sref;
b = sqrt(sizing.AR * sref);
ar = sizing.AR;
% wing_lambdaLE = sizing.sweepLE;

wing_sweepLE = 27.4549;


tailplane.initial.horizontal.l = 0.45*fuselage_length;
tailplane.initial.vertical.l = 0.45*fuselage_length;

%correction for T-tail endplate effect  and trimmable horizontal stab
correction_factor = .95 * .9; 

%from aero tn tailsizing pdf
tailplane.initial.horizontal.v = 0.7 * correction_factor; % howe
tailplane.initial.vertical.v = 0.065 * correction_factor; 

tailplane.initial.horizontal.s = ( tailplane.initial.horizontal.v * cmac * sref ) / (tailplane.initial.horizontal.l);

tailplane.initial.vertical.s = ( tailplane.initial.vertical.v * b * sref ) / (tailplane.initial.vertical.l);

%from Obert
tailplane.initial.horizontal.Ar = 3.9; %this is probably too big
tailplane.initial.vertical.Ar = 1.2; 

%taper
tailplane.initial.horizontal.lambda = 0.35; % obert -> looked good to me idk
tailplane.initial.vertical.lambda = .7; %obert -> cant be too small as the vertical stabilizer has to support horizotnal tailplane

tailplane.initial.horizontal.sweepLE = wing_sweepLE + 5; 
tailplane.initial.vertical.sweepLE = 40; 

tailplane.initial.horizontal.sweep_25 = sweep_angle(tailplane.initial.horizontal.sweepLE,25,0,tailplane.initial.horizontal.Ar,tailplane.initial.horizontal.lambda)
tailplane.initial.vertical.sweep_25 = sweep_angle(tailplane.initial.vertical.sweepLE,25,0,tailplane.initial.vertical.Ar,tailplane.initial.vertical.lambda)

tailplane.initial.horizontal.tc = 0.12; %64012
tailplane.initial.vertical.tc = 0.15; % 64015

tailplane.initial.horizontal = wing_geometry_calc_struct(tailplane.initial.horizontal,1);
% using Vertical wing function overload to account for sref being 2*Sv 
tailplane.initial.vertical = wing_geometry_calc_struct(tailplane.initial.vertical,2); % mulitplier of 2 as vertical is only half a wing
% tailplane.initial.vertical.bv = tailplane.initial.vertical.b/2; %define bv as specifc span to vertcial wing

tailplane.initial.horizontal.z_above_fuselage_bar = tailplane.initial.vertical.b / cmac;
tailplane.initial.horizontal.l_bar = tailplane.initial.horizontal.l / cmac;


names = {'Vertical';'Horizontal'};
lever = [tailplane.initial.vertical.l;tailplane.initial.horizontal.l];
initial_volume_coefficent = [0.065; 0.7];
volume_coefficent_corrected = [tailplane.initial.vertical.v;tailplane.initial.horizontal.v];
tail_sref = [tailplane.initial.vertical.s;tailplane.initial.horizontal.s];
AR = [tailplane.initial.vertical.Ar;tailplane.initial.horizontal.Ar];
taper = [tailplane.initial.vertical.lambda;tailplane.initial.horizontal.lambda];
initial_sizing_table = table(names,lever,initial_volume_coefficent,volume_coefficent_corrected,tail_sref,AR,taper)
table2latex(initial_sizing_table,'initial_sizing_table.tex')


horizontal_stab_plot(tailplane.initial.horizontal)
vertical_stab_plot(tailplane.initial.vertical)

save("tailplane","tailplane")

