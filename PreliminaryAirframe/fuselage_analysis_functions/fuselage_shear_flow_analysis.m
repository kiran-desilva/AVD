function [required_skin_thickness,max_shear_flow,total_shear_flow,q_func,fig] = fuselage_shear_flow_analysis(loadcase,material,doplot)
    
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
        fig = figure
        polarplot(phi_range,abs(total_shear_flow_dist));
        ax1 = gca;
        ax1.ThetaZeroLocation = 'bottom';
        hold on
        polarplot(phi_range,7500*ones(size(phi_range)));
        legend('Shear Flow Distribution','Fuselage')
        grid on



    end