clc
clear

load("wing.mat")
load("locations.mat")
load("fuse.mat")
load("tailplane.mat")

%wing xcg
x_wing = locations.x_wing
x_sweep = 0.35*0.5*wing.b/tand(wing.sweepLE)
chord_35 = (wing.Ctip-wing.Croot)*0.35 +wing.Croot
spar_distance = chord_35*0.75-chord_35*0.1 %spars at 10% and 75% chord
x_chord = 0.7*spar_distance + 0.1*chord_35

wing_xcg_metres = x_wing + x_sweep + x_chord
cg.x_wing = xcg_metres * 3.28084 %in feet


%fuselage xcg
cg.x_fuse = (0.47 * fuse.total_fuselage_length)/12

%tail xcg
offset_55_c = 0.55*tailplane.vertical.b*tand(tailplane.vertical.sweepLE);
cg.x_tail = locations.x_vertical_stabilizer+offset_55_c+((tailplane.vertical.Croot - ((tailplane.vertical.Croot-tailplane.vertical.Ctip)*0.55))*.42)



save('cg','cg')




