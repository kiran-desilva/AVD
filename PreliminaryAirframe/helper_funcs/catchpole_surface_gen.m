clear
clc
load('catchpole_data.mat')

h_b_range = linspace(0.1,1,10000);
ts_t_range = [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2.0];
[h_b,ts_t] = meshgrid(h_b_range,ts_t_range);


for i = 1:numel(ts_t_range)
    K_array(i,:) = feval(fit_arr{i},h_b_range);
end
figure
subplot(1,2,1)
mesh(h_b,ts_t,K_array)
xlabel('h_b')
ylabel('ts_t')
zlabel('K')

% [xData, yData, zData] = prepareSurfaceData( h_b, ts_t, K_array );

% catchpole_3d_fit = fit([xData,yData],zData,'linearinterp')
catchpoleInt3D = griddedInterpolant(h_b',ts_t',K_array','linear','linear');

%plot extrapolation behaviour
subplot(1,2,2)

xrange = linspace(0,2);
yrange = linspace(0,3);
[x,y] = ndgrid(xrange,yrange);
zData = catchpoleInt3D(x,y);

mesh(x,y,zData)
xlabel('h_b')
ylabel('ts_t')
zlabel('K')
xlim([0.1 1])
ylim([0.5 2])
zlim([0 10])

save catchpole_int