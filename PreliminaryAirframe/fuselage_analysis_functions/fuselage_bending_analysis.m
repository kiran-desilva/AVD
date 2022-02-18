%%stringer
%%stringer.thickness
%%stringer.flange_to_web_ratio
%%stringer.web_height
%%
%%
%%%%

%%assuming z stringers here
function [weight_per_length,booms,max_principle_stress,fig] = fuselage_bending_analysis(loadcase,fuselage_skin_thickness,n_stringers,material,stringer,doplot)

    

    df = 1.68;
    cf = pi*df; %fuselage circumference

    booms.phi = linspace(0,2*pi,n_stringers+1);
    booms.phi(end) = []; %remove last element as this overlaps with 0th element
    booms.x = (df/2) * cos(booms.phi);
    booms.y = (df/2) * sin(booms.phi);
    %s coordinate -> s along circumference ccw
    booms.s = linspace(0,cf,n_stringers+1);
    booms.s(end) = []; 

    booms.skin_thickness = fuselage_skin_thickness;
    booms.stringer_pitch = abs(booms.s(2) - booms.s(1)); 

    booms.stringer_area = stringer.web_height*stringer.thickness*(1+(2*stringer.flange_to_web_ratio));

    booms.skin_equivalent_boom_area = 15*fuselage_skin_thickness^2; %15t^2

    single_boom_area = booms.skin_equivalent_boom_area + booms.stringer_area;

    booms.area = single_boom_area*ones(size(booms.phi));
    booms.ixx = booms.area.*(booms.y.^2);
    booms.iyy = booms.area.*(booms.x.^2);
    booms.total_ixx = sum(booms.ixx);
    booms.total_iyy = sum(booms.iyy);

    %applying My only so looking for max x
    % max_stress_edge = abs(max(booms.x));
    My = max(abs(loadcase.bending_moment));
    Mx = 0;
    booms.stress = (My/booms.total_iyy)*(booms.x) + (Mx/booms.total_ixx)*(booms.y);
    [max_principle_stress,max_principle_stress_idx] = max(booms.stress);

    max_yield_stress = material.tensile_yield;

    %check max princple stress is below yield stress
    assert(max_principle_stress<=max_yield_stress,'tensile yield stress exceeded!!');

    total_csa = (n_stringers*booms.stringer_area) + (fuselage_skin_thickness*cf);
    weight_per_length = (total_csa * material.rho);



    if doplot
        fig = figure
        hold on
        phirange = linspace(0,2*pi,1000);
        plot3(0.5*df*cos(phirange),0.5*df*sin(phirange),zeros(size(phirange)),'color','blue');
        plot3(booms.x,booms.y,zeros(size(booms.x)),'o','color','blue','MarkerFaceColor','blue','markersize',6);
        quiver3(booms.x,booms.y,zeros(size(booms.x)),zeros(size(booms.x)),zeros(size(booms.x)),booms.stress,...
                                                                                                            'color','black',...
                                                                                                            'LineWidth',1.3,...
                                                                                                            'AutoScaleFactor',1.5,...
                                                                                                            'MaxHeadSize',.1);
        xlabel('X coordinate [m]')
        ylabel('y coordinate [m]')
        zlabel('Stress due to bending')
        legend('','Boom Area','Scaled stress')
        grid on
        view([330 10])

    end