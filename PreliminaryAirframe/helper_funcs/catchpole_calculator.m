function K_corrected = catchpole_calculator(h_over_b, t_s_over_t, poissons_ratio)
	if h_over_b < 0 || h_over_b > 1 || t_s_over_t < 0.2 || t_s_over_t > 2.2
        error('Tried to sample catchpole diagram out of distribution');
    end

    persistent cd

	if isempty(cd)
		load catchpole_data.mat catchpole_data;
		cd = catchpole_data;
	end

	% select the right curve
	[~, idx] = min(abs(cd.t_s_over_t - t_s_over_t));

	K = cd.fit{idx}(h_over_b); % Read off of catchpole boi
	K_corrected = K*pi^2/(12*0.9*(1 - poissons_ratio^2));
end