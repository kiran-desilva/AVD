clear

load('sizing.mat');
load('parameters.mat');
load('design.mat');
load('wing.mat');




fuselage_length = 11.83;
cmac = wing.Cmac;
sref = design.sref;
b = sqrt(sizing.AR * sref);
ar = sizing.AR;
wing_taper = 0.3017;
% wing_lambdaLE = sizing.sweepLE;

wing_sweepLE = sweep_angle(18.1692,0,25,ar,wing_taper)


tailplane.initial.horizontal.l = 0.45*fuselage_length;
tailplane.initial.vertical.l = 0.45*fuselage_length;



%from aero tn tailsizing pdf
tailplane.initial.horizontal.v = 0.7 * .95 * .9; % howe
tailplane.initial.vertical.v = 0.065 * .95; 

tailplane.initial.horizontal.s = ( tailplane.initial.horizontal.v * cmac * sref ) / (tailplane.initial.horizontal.l);

tailplane.initial.vertical.s = ( tailplane.initial.vertical.v * b * sref ) / (tailplane.initial.vertical.l);

%from Obert
tailplane.initial.horizontal.Ar = 3.9; %this is probably too big
tailplane.initial.vertical.Ar = 1.3; 

%taper
tailplane.initial.horizontal.lambda = .5; % obert -> looked good to me idk
tailplane.initial.vertical.lambda = .7; %obert -> cant be too small as the vertical stabilizer has to support horizotnal tailplane

tailplane.initial.horizontal.sweepLE = wing_sweepLE + 5; 
tailplane.initial.vertical.sweepLE = 35; 

tailplane.initial.horizontal.sweep_25 = sweep_angle(tailplane.initial.horizontal.sweepLE,25,0,tailplane.initial.horizontal.Ar,tailplane.initial.horizontal.lambda)
tailplane.initial.vertical.sweep_25 = sweep_angle(tailplane.initial.vertical.sweepLE,25,0,tailplane.initial.vertical.Ar,tailplane.initial.vertical.lambda)

tailplane.initial.horizontal.sweepTE = sweep_angle(tailplane.initial.horizontal.sweepLE,25,100,tailplane.initial.horizontal.Ar,tailplane.initial.horizontal.lambda)
tailplane.initial.vertical.sweepTE = sweep_angle(tailplane.initial.vertical.sweepLE,25,100,tailplane.initial.vertical.Ar,tailplane.initial.vertical.lambda)

tailplane.initial.horizontal.sweepmaxTC = sweep_angle(tailplane.initial.horizontal.sweepLE,25,40,tailplane.initial.horizontal.Ar,tailplane.initial.horizontal.lambda)
tailplane.initial.vertical.sweepmaxTC= sweep_angle(tailplane.initial.vertical.sweepLE,25,35,tailplane.initial.vertical.Ar,tailplane.initial.vertical.lambda)

tailplane.initial.horizontal.tc = 0.12; %64012
tailplane.initial.vertical.tc = 0.15; % 64015


tailplane.initial.horizontal = wing_geometry_calc_struct(tailplane.initial.horizontal,1);
% using Vertical wing function overload to account for sref being 2*Sv 
tailplane.initial.vertical = wing_geometry_calc_struct(tailplane.initial.vertical,2); % mulitplier of 2 as vertical is only half a wing
% tailplane.initial.vertical.bv = tailplane.initial.vertical.b/2; %define bv as specifc span to vertcial wing

tailplane.initial.horizontal.z_above_fuselage_bar = tailplane.initial.vertical.b / cmac;
tailplane.initial.horizontal.l_bar = tailplane.initial.horizontal.l / cmac;

tailplane.horizontal = tailplane.initial.horizontal;
tailplane.vertical = tailplane.initial.vertical;

%%%%%SIZING TO T-TAIL CONFIG%%%%%%%
%get deep stall limit
[h_tail_constraint_fit,h_tail_constraint_fig] = horizontal_placement_graph;
h_tail_constraint_fig;
% currently assuming wing is on the bottom of the fueslage
fuselage_offset = 1.52;

current_horizontal_position_bar = tailplane.horizontal.z_above_fuselage_bar + (fuselage_offset/cmac);
plot(tailplane.horizontal.l_bar,current_horizontal_position_bar,'o','color','green');

% min_vertical_span_bar = h_tail_constraint_fit(tailplane.horizontal.l_bar);
% yline(min_vertical_span_bar,'--');
% if (min_vertical_span_bar > current_horizontal_position_bar)
%     disp("Bad Vertical Position -> fixing span")
%     required_span = (min_vertical_span_bar * cmac) - fuselage_offset;
%     tailplane.vertical.Ar = ((required_span)^2)/tailplane.vertical.s;
%     tailplane.vertical = wing_geometry_calc_struct(tailplane.vertical,2);
% end

% current_horizontal_position_bar = (tailplane.vertical.b + fuselage_offset)/ cmac;

% plot(tailplane.horizontal.l_bar,current_horizontal_position_bar,'o','color','red');


improvePlot(h_tail_constraint_fig)


%ensure horizontal tailplane can fit on vertical tailplane
min_vertical_Ctip = tailplane.horizontal.Croot;
if tailplane.vertical.Ctip < min_vertical_Ctip
    % correct taper ratio to fix this issue
    disp("Initial vertical stablizer Ctip too small -> fixing taper")
    syms lambda ctip
    ctip(lambda) = (2*lambda*tailplane.vertical.s)/(sqrt(tailplane.vertical.Ar * tailplane.vertical.s) * (1+lambda));
    tailplane.vertical.lambda = fsolve(@(x) double(ctip(x)) - min_vertical_Ctip,0.9);
    %tailplane.vertical.lambda = tailplane.horizontal.Croot/tailplane.vertical.Croot;
    tailplane.vertical = wing_geometry_calc_struct(tailplane.vertical,2); % recalculate paramters
end








names = {'Vertical';'Horizontal'};
lever = [tailplane.vertical.l;tailplane.horizontal.l];
initial_volume_coefficent = [0.065; 0.7];
volume_coefficent_corrected = [tailplane.vertical.v;tailplane.horizontal.v];
tail_sref = [tailplane.vertical.s;tailplane.horizontal.s];
AR = [tailplane.vertical.Ar;tailplane.horizontal.Ar];
taper = [tailplane.vertical.lambda;tailplane.horizontal.lambda];
initial_sizing_table = table(names,lever,initial_volume_coefficent,volume_coefficent_corrected,tail_sref,AR,taper)
table2latex(initial_sizing_table,'initial_sizing_table.tex')


horizontal_stab_plot(tailplane.horizontal);
vertical_stab_plot(tailplane.vertical);

% save("tailplane","tailplane")

