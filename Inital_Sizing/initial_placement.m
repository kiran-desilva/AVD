clear
clc

load('wing.mat')
load('tailplane.mat')

%% x is defined at zero at nose 

total_length = convlength(462,'in','m');

%assume vertical tail plane front spar is 0.1bv from tip

vertical_stabalizer_fs = 0.1 * tailplane.vertical.Croot;
x_rear_bulkhead = total_length - convlength(100,'in','m'); % this needs to be checked

x_engine = x_rear_bulkhead;

x_vertical_stabilizer = x_rear_bulkhead - vertical_stabalizer_fs; % align vert stablizer front spar with x_rear_bulkhead
% x_vertical_stabilizer =  11.73 - tailplane.vertical.Croot 
x_ac_v = x_vertical_stabilizer + tailplane.vertical.Xac_from_tip;

x_horizontal_stabilizer = x_vertical_stabilizer + (tailplane.vertical.b * tand(tailplane.vertical.sweepLE)) + (tailplane.vertical.Ctip/2) - (tailplane.horizontal.Croot /2 ); % align center of horizontal stab with cetner of vertical tailplane
x_ac_h = x_horizontal_stabilizer + tailplane.horizontal.Xac_from_tip;

% base wing dimentioning from lever arm of horiztonal tailplane as this is where it will be crucial

x_ac_w = x_ac_h - tailplane.horizontal.l;
% calculate tip location of wing
x_wing = x_ac_w - wing.Xac_from_tip;

% gear estimate in center of chord at edge of fuselage
fuselage_radius = convlength(66/2,'in','m');
fuselage_diameter_bar = fuselage_radius/(wing.b/2);
chord_at_point =  wing.Croot - ((wing.Croot-wing.Ctip)*fuselage_diameter_bar);
x_gear = x_wing + (chord_at_point*tand(wing.sweepLE)) + (chord_at_point/2); 



locations.total_length = total_length;
locations.vertical_stabalizer_fs = vertical_stabalizer_fs;
locations.x_rear_bulkhead = x_rear_bulkhead;
locations.x_engine = x_engine;
locations.x_vertical_stabilizer = x_vertical_stabilizer;
locations.x_ac_v = x_ac_v;
locations.x_horizontal_stabilizer = x_horizontal_stabilizer;
locations.x_ac_h = x_ac_h;
locations.x_ac_w = x_ac_w;
locations.x_wing = x_wing;
locations.x_gear = x_gear;

figure
hold on
keys  = fieldnames(locations)
for i = 1:numel(keys)
    plot(locations.(keys{i}),0,'o','linewidth',2,'markersize',5)
    
end
legend(keys)
for i = 1:numel(keys)
    xline(locations.(keys{i}),'--')
end
grid on
axis equal
xlabel('m')

save('locations','locations')