clc
clear

load("wing.mat")
load("locations.mat")
load("fuse.mat")
load("tailplane.mat")
load("uc.mat")

%wing xcg
x_wing = locations.x_wing;
x_sweep = 0.35*0.5*wing.b/tand(wing.sweepLE);
chord_35 = (wing.Ctip-wing.Croot)*0.35 +wing.Croot;
spar_distance = chord_35*0.75-chord_35*0.1; %spars at 10% and 75% chord
x_chord = 0.7*spar_distance + 0.1*chord_35;

wing_xcg_metres = x_wing + x_sweep + x_chord;
cg.x_wing = wing_xcg_metres * 3.28084 %in feet


%fuselage xcg
cg.x_fuse = (0.47 * fuse.total_fuselage_length)/12

%vertical tail cg
offset_55_c = 0.55*tailplane.vertical.b*tand(tailplane.vertical.sweepLE);
vtail_xcg_metres = locations.x_vertical_stabilizer+offset_55_c+((tailplane.vertical.Croot - ((tailplane.vertical.Croot-tailplane.vertical.Ctip)*0.55))*.42);
cg.x_vtail = vtail_xcg_metres * 3.28084

vtail_zcg_metres = 0.55*tailplane.vertical.b + (0.84-(locations.x_vertical_stabilizer-fuse.x_afterbody)*tand(4))
cg.z_vtail = vtail_zcg_metres * 3.28084

%horizontal tail cg
offset_vtail = tailplane.vertical.b*tand(tailplane.vertical.sweepLE);
x_tailsweep = 1.2888*tand(26.6513);
tailchord_35 = (tailplane.horizontal.Ctip-tailplane.horizontal.Croot)*0.35+tailplane.horizontal.Croot;
tailspar_distance = tailchord_35*0.68-tailchord_35*0.15;
x_tailchord = 0.7*tailspar_distance + 0.1*tailchord_35;

htail_xcg_metres = x_tailchord + x_tailsweep + offset_vtail + locations.x_vertical_stabilizer;
cg.x_htail = htail_xcg_metres * 3.28084

htail_zcg_metres = tailplane.vertical.b + (0.84-(fuse.x_afterbody-locations.x_vertical_stabilizer)*tand(4)) + 12/((tailplane.horizontal.Ctip-tailplane.horizontal.Croot)*0.535+tailplane.horizontal.Croot);
cg.z_htail = htail_zcg_metres * 3.28084

%mlg
cg.z_mlg = -33/12
cg.x_mlg = locations.x_wing

%nlg
cg.z_nlg = -33/12
cg.x_nlg = uc.nose_wheel.x * 3.28084

%nacelle & engine
cg.x_inl = locations.x_rear_bulkhead * 3.28084
cg.z_inl = 1.1 * 3.28084

%APU
cg.x_APUinst = (fuse.total_fuselage_length-12)/12 





save('cg','cg')




