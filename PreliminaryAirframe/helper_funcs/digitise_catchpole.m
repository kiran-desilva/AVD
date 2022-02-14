clear
clc

catchpol_data_path = fullfile('..', 'digitised_plots', 'catchpole_data', 'catchpole_z_d_h_03.csv');

curve_arr = [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2.0];
fit_arr = digitise_single_catchpole_plot(catchpol_data_path, curve_arr);

catchpole_data.t_s_over_t = curve_arr;
catchpole_data.fit = fit_arr;

save catchpole_data

if false
    figure;
    hold on;
    % for i = 1:numel()
	% 	scatter( 'r');
    % end
    hold off;
	grid on;
end

function catchpole_fit_arr = digitise_single_catchpole_plot(fpath, curve_arr)
	data = readtable(fpath);

	for i = 1:numel(curve_arr) 
		curve = data{:, 2*i - 1:2*i}; % load x-y pairs
		curve = curve(sum(isnan(curve), 2)==0, :); % remove nan rows
		catchpole_fit_arr{i} = fit(curve(:, 1), curve(:, 2), 'smoothingspline');
        figure;
        plot(catchpole_fit_arr{i}, curve(:, 1), curve(:, 2));
	end

end
