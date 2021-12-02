clear
clc
close all

param_gen
init_weight_sizing
constraint_diagrams

avd_wing
initial_tailplane_sizing
fuselage

initial_placement
change_wing_ac(3.39*wing.Cmac)
% horizontal_placement_graph
% load('tailplane')
% load('wing')
% fuselage_offset = 1.52;
% current_horizontal_position_bar = tailplane.horizontal.z_above_fuselage_bar + (fuselage_offset/wing.Cmac);
% plot(tailplane.horizontal.l/wing.Cmac,current_horizontal_position_bar,'o','color','green');

cg_calcs
weight_balance
static_stability
