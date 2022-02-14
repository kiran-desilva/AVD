%clear
%clc
close all

farrar_data_path = fullfile('..', 'digitised_plots', 'farrar_data', 'farrar_data');

farrar_data.f_vals = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5];
farrar_table = readtable(farrar_data_path);

R = @(theta) [cosd(theta) -sind(theta); sind(theta) cosd(theta)]';
middle_pt = [1.5347, 1.0642];

rotation_for_continuity = 150;
farrar_data.transform = @(d) (d - middle_pt)*R(rotation_for_continuity);
for i = 1:numel(farrar_data.f_vals)
	d = farrar_table{:, 2*i - 1:2*i};
	d = d(sum(isnan(d),2)==0, :) - middle_pt;
    d = d*R(rotation_for_continuity);
	[theta, rho] = cart2pol(d(:, 1), d(:, 2));
	farrar_data.data{i} = [theta, rho];
    
    d = farrar_data.data{i};
    f = fit(d(:, 1), d(:, 2), 'smoothingspline', 'SmoothingParam', 0.9999);
    %f = fit(d(:, 1), d(:, 2), 'linearinterp');
    %figure;
    %plot(f, d(:, 1), d(:, 2));
	farrar_data.fit{i} = @(theta) (pol_fit_to_cart(theta, f))*R(-rotation_for_continuity) + repmat(middle_pt, size(theta));
	farrar_data.x_range{i} = [min(d(:, 1)), max(d(:, 1))];
end

save('farrar_data')
index_into = @(A, idx) A(idx);  

if true
    figure;
    hold on;
    for i = 1:numel(farrar_data.f_vals)
		x_space = linspace(farrar_data.x_range{i}(1), farrar_data.x_range{i}(2), 100);
		% d = farrar_data.data{i};
		f = farrar_data.fit{i};
		XY = f(x_space');
		scatter(XY(:, 1), XY(:, 2),'r');
		% plot(x_space, farrar_data.fit{i}(x_space));
    end
    hold off;
	grid on;
end

function XY = pol_fit_to_cart(theta, fit)
	r = fit(theta);
	[x, y] = pol2cart(theta, r);
    XY = [x, y];
end
