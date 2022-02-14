

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
        shp = polyshape(coords(:,1),coords(:,2));
        Ina = (((0.633*b*t) + (0.37*h*ts))/((1.6*h*ts)+(b*t)))*((h^3)*ts); %from Farrar, only applies to d_h = 0.3
          
    end

    area = shp.area;
    [cx,cy] = shp.centroid;
    Teff = area/b;
    Kna = sqrt(Ina/area);
 



    if (doplot)
        figure
        hold on
        plot(shp) 
        plot(cx,cy,'x','color','red')
    end

end