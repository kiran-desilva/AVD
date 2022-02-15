function K_corrected = catchpole_calculator(h_over_b, t_s_over_t, poissons_ratio, use3D)
	if h_over_b < 0 || h_over_b > 1 || t_s_over_t < 0.2 || t_s_over_t > 2.2
        error('Tried to sample catchpole diagram out of distribution');
    end

    persistent cd
	persistent cdi

	if isempty(cd)
		load catchpole_data.mat catchpole_data;
		cd = catchpole_data;
	end
	if isempty(cdi)
		load catchpole_int.mat catchpoleInt3D;
		cdi = catchpoleInt3D;
	end


	% select the right curve
	if use3D
		K = cdi(h_over_b,t_s_over_t);
	else
		[~, idx] = min(abs(cd.t_s_over_t - t_s_over_t));

		K = cd.fit{idx}(h_over_b); % Read off of catchpole boi
	end
	
	K_corrected = K*pi^2/(12*0.9*(1 - poissons_ratio^2));
end