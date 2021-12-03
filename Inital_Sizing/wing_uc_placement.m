clear
clc
load('locations')
load('wing')
load('uc')

struct = wing;

offset_from_center = 1.18;


% uc_lim = [4.95];
uc_lim=[5.16];



span_percent = (offset_from_center/(0.5*wing.b));
chord_offset = wing.Croot - ((wing.Croot-wing.Ctip )* span_percent);
LE_offset = offset_from_center*tand(wing.sweepLE);
TE_offset = offset_from_center*tand(wing.sweepTE);

uc_lim_transformed = wing.Croot - (uc_lim-(locations.x_wing));

figure 
hold on

rhs_wing = [0,0;
            0,struct.Croot;
            struct.b/2,struct.Croot - (0.5*struct.b * tand(struct.sweepLE));
            struct.b/2,(struct.Croot - (0.5*struct.b * tand(struct.sweepLE)) - struct.Ctip);
            0,0;
]
% lhs_wing = [0,0;
%             0,struct.Croot;
%             -struct.b/2,struct.Croot - (0.5*struct.b * tand(struct.sweepLE));
%             -struct.b/2,(struct.Croot - (0.5*struct.b * tand(struct.sweepLE)) - struct.Ctip);
%             0,0;
% ]

plot(rhs_wing(:,1),rhs_wing(:,2),'color','black');

plot(offset_from_center,mean(uc_lim_transformed),'x','color','green','markersize',10,'linewidth',1.5)

xline(offset_from_center,'--','color','blue')
xline(0.84,'-','Fuselage','color','red','linewidth',1.2)
yline(uc_lim_transformed(1),'--','color','blue')
yline(uc_lim_transformed(2),'--','color','blue')

grid on
grid minor
axis equal
