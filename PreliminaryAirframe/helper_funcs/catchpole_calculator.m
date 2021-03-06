function K_corrected = catchpole_calculator(h_over_b, t_s_over_t, poissons_ratio, use3D)
	if h_over_b < 0 || h_over_b > 0.9
        error(['Tried to sample catchpole diagram out of distribution, got h/b = ', num2str(h_over_b), ', t_s/t = ', num2str(t_s_over_t)]);
    end
    
    if t_s_over_t > 2
        K_corrected = 6 * 4 * (pi^2)/(12*(1-(poissons_ratio^2)));
        return;
    end
    
    if t_s_over_t < 0.43
        K_corrected = 4 * (pi^2)/(12*(1-(poissons_ratio^2))); % K = 1
        return;
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
	
	%K_corrected = K*4*pi^2/(12*(1 - poissons_ratio^2));
	K_corrected = K;
end