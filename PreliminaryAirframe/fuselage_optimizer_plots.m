% clear
% clc

addpath(fullfile('.','fuselage_analysis_functions'))

% grid_data = load('test_weights_3.mat');
% analysis = fuselage_search_analysis(grid_data,0);

%plot valid function space

figure
hold on
get_boundary_N(80,analysis,'red');
get_boundary_N(70,analysis,'blue');
get_boundary_N(60,analysis,'green');


xlabel('Stringer Thickness [m]')
ylabel('Stringer Height [m]')
zlabel('Fuselage Thickness [m]')
legend('80 Stringers','70 Stringers','60 Stringers')
grid on
view([-35 10])

% set(gcf, 'rend', 'painters', 'Units', 'pixels', 'pos', ...
%         [100 100 800 600]);

% axis_handles=findobj(gcf,'type','axe');

% for i = 1:length(axis_handles)
%     ax = axis_handles(i);

%     % Change default font size (tick labels, legend, etc.)
%     set(ax, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 1);
    
%     set(ax, 'Box', 'on');

%     % Change font size for axis text labels
%     set(get(ax, 'XLabel'),'FontSize', 15, 'FontWeight', 'Bold');
%     set(get(ax, 'YLabel'),'FontSize', 15, 'FontWeight', 'Bold');

% end




mins = get_mins(analysis);
fs = mins(:,1);
ts = mins(:,2);
hs = mins(:,3);
n = mins(:,4);
w = mins(:,5);
w_scaled = log(w);


fs_n_w_int = scatteredInterpolant(fs,n,w_scaled,'natural','linear');
ts_n_w_int = scatteredInterpolant(ts,n,w_scaled,'natural','linear');
hs_n_w_int = scatteredInterpolant(hs,n,w_scaled,'natural','linear');

get_range = @(array) linspace(min(array),max(array),100);

%surface plots
figure

colormap(turbo)

tiledlayout(2,2);
nexttile

hold on

create_mesh(get_range(fs),get_range(n),fs_n_w_int);
create_contour(get_range(fs),get_range(n),fs_n_w_int);

scatter(fs,n,'red','x')
scatter3(fs,n,w_scaled,'red','x')

xlabel("Fuselage Thickness [m]")
ylabel("Number of Stringers")
zlabel("Log(Weight)")

legend('','Minimum Weight Surface','Minimum Weight Solutions')


view([-130 30])

grid on

nexttile

hold on
create_mesh(get_range(ts),get_range(n),ts_n_w_int);
create_contour(get_range(ts),get_range(n),ts_n_w_int);

scatter(ts,n,'red','x')
scatter3(ts,n,w_scaled,'red','x')

xlabel("Stringer Thickness [m]")
ylabel("Number of Stringers")
zlabel("Log(Weight)")

legend('','Minimum Weight Surface','Minimum Weight Solutions')


view([-130 30])
grid on

nexttile

hold on
create_mesh(get_range(hs),get_range(n),hs_n_w_int);
create_contour(get_range(hs),get_range(n),hs_n_w_int);

scatter(hs,n,'red','x')
scatter3(hs,n,w_scaled,'red','x')

xlabel("Stringer Height [m]")
ylabel("Number of Stringers")
zlabel("Log(Weight)")

legend('','Minimum Weight Surface','Minimum Weight Solutions')


view([-130 30])
grid on


set(gcf, 'rend', 'painters', 'Units', 'pixels', 'pos', ...
        [100 100 1600 1200]);

axis_handles=findobj(gcf,'type','axe');

for i = 1:length(axis_handles)
    ax = axis_handles(i);

    % Change default font size (tick labels, legend, etc.)
    set(ax, 'FontSize', 15, 'FontName', 'Arial', 'LineWidth', 1);
    
    set(ax, 'Box', 'on');

    % Change font size for axis text labels
    set(get(ax, 'XLabel'),'FontSize', 15, 'FontWeight', 'Bold');
    set(get(ax, 'YLabel'),'FontSize', 15, 'FontWeight', 'Bold');

end




function [tsurf,boundary_idx] = get_boundary_N(N,data,color)

    flatten = @(array) array(:);

    Ns = data.full_domain.N(1,1,1,:);
    Ns = Ns(:);
    Nidx = (Ns == N);

    x = data.valid_domain.TS(:,:,:,Nidx);
    y = data.valid_domain.HS(:,:,:,Nidx);
    z = data.valid_domain.FS(:,:,:,Nidx);
    xnan = isnan(x);
    ynan = isnan(y);
    znan = isnan(z);
    nanidx = xnan | ynan | znan;
    %get non nan
    x = flatten(x(~nanidx));
    y = flatten(y(~nanidx));
    z = flatten(z(~nanidx));
    % get boundary
    [boundary_idx] = boundary(x(:),y(:),z(:));



    tsurf = trimesh(boundary_idx,x(:),y(:),z(:),'EdgeColor',[0.8 0.8 0.8],'FaceColor',color,'FaceAlpha',0.2,'LineWidth',0.1);
    % tsurf = triangulation(boundary_idx,x(:),y(:),z(:))

end

function [mins] = get_mins(data)
    [~,idx] = min(data.valid_domain.weights,[],[1 2 3],'linear');
    fs = data.valid_domain.FS(idx);
    ts = data.valid_domain.TS(idx);
    hs = data.valid_domain.HS(idx);
    N = data.valid_domain.N(idx);
    weight = data.valid_domain.weights(idx);
    mins = [fs(:),ts(:),hs(:),N(:),weight(:)];
end

function [m] = create_mesh(x,y,fit)
    [X,Y] = meshgrid(x,y);
    Z = fit(X,Y);
    m = surf(X,Y,Z,'EdgeColor','interp');
    % m = contourf(X,Y,Z);
end  


function [m] = create_contour(x,y,fit)
    [X,Y] = meshgrid(x,y);
    Z = fit(X,Y);
    % m = mesh(X,Y,Z);
    m = contourf(X,Y,Z);
end  


