% dimensions in inches unless stated otherwise

clear

d_f = 66


%cabin length
seat_length = 41;
leg_room = 40;
lavatory = 40;
xtra = 49.2;

cabin_length = seat_length*2 + leg_room + lavatory + xtra

%nose length
L_D_n = 1.3;
flightdeck_length = 70
nose_cone = L_D_n*d_f - flightdeck_length
nose_length = nose_cone + flightdeck_length

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


%total fuselage length
fuse.total_fuselage_length = afterbody_length+cabin_length+nose_length
L_D_f = fuse.total_fuselage_length/d_f

% %volume of pressurised sections
% cone_vol = 
% nose_vol = ... -cone_vol
% cabin_vol = pi * (d_f/2)^2 * cabin_length
% aft_vol =
% 
% vol_press = nose_vol + cabin_vol + aft_vol

%structural length
aft_bulkhead = 100
fuse.structural_length = fuse.total_fuselage_length - nose_cone - aft_bulkhead

save('fuse', 'fuse')
