function K_corrected = catchpole_calculator(h_over_b, t_s_over_t, poissons_ratio)
	% TODO:
	K = ...; % Read off of catchpole boi
	K_corrected = K*pi^2/(12*0.9*(1 - poissons_ratio^2));
end