% dimensions in inches unless stated otherwise

clear

d_f = 66
fuse.d_f = d_f;


%cabin length
seat_length = 41;
leg_room = 40;
lavatory = 40;
xtra = 49.2;

cabin_length = seat_length*2 + leg_room + lavatory + xtra
fuse.cabin_length = cabin_length;

%nose length
L_D_n = 1.3;
flightdeck_length = 70
nose_cone = L_D_n*d_f - flightdeck_length
nose_length = nose_cone + flightdeck_length
fuse.nose_length = nose_length;

% %afterbody length
% L_D_a = 5
% afterbody_length = L_D_a * d_f
% theta_afterbody = atand(d_f/afterbody_length)

%afterbody length
L_D_a = 2.5
afterbody_length = L_D_a * d_f
bottom_theta_afterbody = 11
% upper_theta_afterbody = atand((d_f-(tand(bottom_theta_afterbody)*afterbody_length))/afterbody_length)
upper_theta_afterbody = 4;

fuse.x_afterbody = (nose_length + cabin_length)*0.0254; %in metres

fuse.bottom_theta_afterbody = bottom_theta_afterbody;
fuse.upper_theta_afterbody = upper_theta_afterbody;

%total fuselage length
fuse.total_fuselage_length = afterbody_length+cabin_length+nose_length
fuse.L_D_f = fuse.total_fuselage_length/d_f


%structural length
aft_bulkhead = 100
fuse.structural_length = fuse.total_fuselage_length - nose_cone - aft_bulkhead
fuse.x_aft_bulkhead = fuse.total_fuselage_length-aft_bulkhead;


save('fuse', 'fuse')
