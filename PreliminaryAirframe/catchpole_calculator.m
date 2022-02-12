function K_corrected = catchpole_calculator(h_over_b, t_s_over_t, poissons_ratio)
	% TODO:
	K = NaN; % Read off of catchpole boi
	K_corrected = K*pi^2/(12*0.9*(1 - poissons_ratio^2));
	K_corrected = 3;
	disp("WARNING!! Catchpole function not properly implemented yet!!")
end