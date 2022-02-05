clear
clc

f = stlread('Fuselage_body.stl');
f1.vertices = f.Points;
f1.faces = f.ConnectivityList;

start = 4.19;
stop = -7.54;
disc = 1000;
range = linspace(start,stop,disc);
h = (stop-start)/(disc-1);

areas = zeros(size(range));
areas_2 = areas;

figure
hold on

for i=1:length(range)
    [y,z] = meshgrid(-2:2);             % meshgrid for plane
    h = surf(y*0+range(i),y,z,'visible','off');                % crossection plane
    f2 = surf2patch(h,'triangles');     % convert to patch
    [~,ff] = SurfaceIntersection(f1,f2);% find crossection line
    ff = rmfield(ff,'edges');           % remove edges
    f_csa = ff.faces;
    csa_yy = ff.vertices(:,2);
    csa_zz = ff.vertices(:,3);
    [ix, v] = boundary(csa_yy,csa_zz);
    areas(i) = v;

    plot3(range(i)*ones(size(csa_yy(ix))),csa_yy(ix),csa_zz(ix))

end

axis equal
grid on


figure
hold on
plot(areas)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( range, areas );

% Set up fittype and options.
ft = 'linearinterp';
excludedPoints = excludedata( xData, yData, 'Indices', [42 43 44 388 460 742 743 858 871 887 964] );
opts = fitoptions( 'Method', 'LinearInterpolant' );
opts.Normalize = 'on';
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'areas vs. range', 'Excluded areas vs. range', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'range', 'Interpreter', 'none' );
ylabel( 'areas', 'Interpreter', 'none' );
grid on


integrand = integrate(fitresult,xData,start)

voldist.areas = areas;
voldist.fit = fitresult;
voldist.totalvol = abs(integrand);

save('voldist.mat','voldist');