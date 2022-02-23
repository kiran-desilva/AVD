function [required_skin_thickness,max_shear_flow,total_shear_flow,q_func,fig] = fuselage_shear_flow_analysis(loadcase,material,doplot)
    fig = 0;
    
    df = 1.68; %fuselage
    rf = df/2;

    %no idea what alpha does here ecentricity??
    q_func = @(T,P,Q,phi,alpha) ((T+(P*rf))/(2*pi*(rf^2))) + ((P*cos(phi-alpha))/(pi*rf)) + ((Q*sin(phi-alpha))/(pi*rf));
    
    phi_range = linspace(0,2*pi,1000);
    max_shear = max(loadcase.shear);
    total_shear_flow =@(phi) q_func(0,0,max_shear,phi,0);
    total_shear_flow_dist = total_shear_flow(phi_range);

    [max_shear_flow,idx] = max(abs(total_shear_flow_dist));
    max_shear_flow_phi = phi_range(idx);
    required_skin_thickness = max_shear_flow/material.shear_yield;


    if doplot
        fig = figure('Name','skin_shear_flow_dist')

        hold on

        unit_shear_flow_dist = abs(total_shear_flow_dist)./max(total_shear_flow_dist);
        fuselage = ones(size(phi_range));
        fuselage_dist = fuselage + unit_shear_flow_dist;

        [fx,fy] = pol2cart(phi_range+(pi/2),fuselage);
        [fsx,fsy] = pol2cart(phi_range+(pi/2),fuselage_dist);
        
        fuse = plot(fx,fy,'color','red')
        dist = plot(fsx,fsy,'color','green')
        

        

        for i=1:5:length(phi_range)
            p = plot([fx(i) fsx(i)],[fy(i) fsy((i))],'color','green');

        end

        legend([fuse dist],'Fuselage','Scaled Shear Flow Distribution')
        
        

        
        xlabel('X coordinate')
        ylabel('Y coordinate')

        grid on
        axis equal
        xlim([-2.25 2.25])



    end