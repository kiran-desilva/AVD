

function [Kna,Ina,area,Teff] = stringer_panel_geometry(b,t,ts,td,d,h,type,doplot)
 
    shp = polyshape();
    if strcmp(type,"Z_dh03")
        assert(abs((d/h)-0.3)<1e-3,'d/h needs to be 0.3');
        coords = [-b/2, 0;...
                  b/2, 0;...
                  b/2, t;...
                  ts/2, t;...
                  ts/2, t+h;...
                  d, t+h;...
                  d, t+(h+td);...
                  -ts/2, t+(h+td);...
                  -ts/2, t+td;...
                  -d, t+td;...
                  -d,t;...
                  -b/2,t];

        % coords = [
        %           ts/2, t;...
        %           ts/2, t+h;...
        %           d, t+h;...
        %           d, t+(h+td);...
        %           -ts/2, t+(h+td);...
        %           -ts/2, t+td;...
        %           -d, t+td;...
        %           -d,t];
        shp = polyshape(coords(:,1),coords(:,2));


        Ina = (((0.633*b*t) + (0.37*h*ts))/((1.6*h*ts)+(b*t)))*((h^3)*ts); %from Farrar, only applies to d_h = 0.3
        % Ina = (2.54e-2)^4*Ina;
          
    end

    area = shp.area;
    [cx,cy] = shp.centroid;
    Teff = area/b;
    Kna = sqrt(Ina/area);

    [~,i,~] = polygeom(coords(:,1),coords(:,2));
    c = (i(1) + i(2))/2;
    r = sqrt((((i(1) - i(2))/2)^2) + i(3)^2);
    imin = c-r;
    imax = c+r;

 
    if (doplot)
        figure
        hold on
        plot(shp) 
        plot(cx,cy,'x','color','red')
        axis equal

        figure
        hold on
        theta_range = linspace(0,2*pi);
        mohr_circle_func = @(theta) [c + r*cos(theta);r*sin(theta)];
        circle_points = mohr_circle_func(theta_range);
        plot(circle_points(1,:),circle_points(2,:));
        plot(imin,0,'x');
        plot(imax,0,'x');
        axis equal
        grid on
        figure
    end

end