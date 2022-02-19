%% Buckling stress function 
% fit: Fourier; 4 terms

% clear
% clc

% load ('a_b_table.mat') %load b greater than a data
% curve_arr=[1,1.5,2,3,10];
% 
% 	for i = 1:numel(curve_arr) 
% 		curve = a_b_final_tab{:, 2*i - 1:2*i}; % load x-y pairs
% 		curve = curve(sum(isnan(curve), 2)==0, :); % remove nan rows
% 		a_b_fit_arr{i} = fit(curve(:, 1), curve(:, 2), 'smoothingspline');
%         figure;
%         plot(a_b_fit_arr{i}, curve(:, 1), curve(:, 2));
% 	end
% 
% x_range=linspace(0,20,100000);
% a_b_range=[1,1.5,2,3,10];
% [x,a_b]=meshgrid(x_range,a_b_range);
% 
% for i=1:numel(a_b_range)
%     K_array_a_b(i,:) = feval(a_b_fit_arr{i},x_range);
% end
% figure
% subplot(1,2,1)
% mesh(x,a_b,K_array_a_b)
% xlabel('x')
% ylabel('b_a')
% zlabel('K_s')
% 
% catchpoleInt3D = griddedInterpolant(x',a_b',K_array_a_b','spline','spline');
% 
% save a_b_int
