function [required_skin_thickness] = fuselage_shear_flow_analysis(loadcase,material,doplot)
    
    df = 1.68; %fuselage
    rf = df/2;

    %no idea what alpha does here ecentricity??
    q_func = @(T,P,Q,phi,alpha) ((T+(P*rf))/(2*pi*(rf^2))) + ((P*cos(phi-alpha))/(pi*rf))
    
    if doplot
        figure
    end