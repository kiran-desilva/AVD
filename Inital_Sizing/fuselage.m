% dimensions in inches unless stated otherwise

d_f = 60


%cabin length
seat_length = 41;
leg_room = 40;
lavatory = 40;
xtra = 20;

cabin_length = seat_length*2 + leg_room + lavatory + xtra;

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
L_D_a = 4
afterbody_length = L_D_a * d_f
bottom_theta_afterbody = 11
upper_theta_afterbody = atand((60-(tand(11)*afterbody_length))/afterbody_length)

%total fuselage length
total_fuselage_length = afterbody_length+cabin_length+nose_length
L_D_f = total_fuselage_length/d_f

% %volume of pressurised sections
% cone_vol = 
% nose_vol = ... -cone_vol
% cabin_vol = pi * (d_f/2)^2 * cabin_length
% aft_vol =
% 
% vol_press = nose_vol + cabin_vol + aft_vol

