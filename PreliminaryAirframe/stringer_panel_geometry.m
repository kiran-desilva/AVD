

function [area,effThickness] = stringer_panel_geometry(b,t,ts,td,d,h,type,doplot)
    shp = polyshape();
    if strcmp(type,"Z")
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
        
    end

    area = shp.area;
    [cx,cy] = shp.centroid
    effThickness = area/b;

    shifted_coords(:,1) = coords(:,1) - cx;
    shifted_coords(:,2) = coords(:,2) - cy;
    %get centroidal moments of intertia
    [geo,intertial,cpm] = polygeom(shifted_coords(:,1),shifted_coords(:,2))


    



    if (doplot)
        figure
        hold on
        plot(shp) 
        plot(cx,cy,'x','color','red')
    end

end