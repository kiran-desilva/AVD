function [frames,fig] = fuselage_light_frames(frame_material,sectionType,booms,stringer,doPlot)
    addpath(fullfile('.','helper_funcs'))
    fig = 0;
    frames.L_conventional = 0.5; % this is for an airline, might be diff for buisness hjet
    df = 1.68;
    lf = 11.37;


    stringer.flange_width = stringer.flange_to_web_ratio*stringer.web_height;
    [Kna,~,stringer_panel_area,~] = stringer_panel_geometry(booms.stringer_pitch,...
                                                            booms.skin_thickness,...
                                                            stringer.thickness,...
                                                            stringer.thickness,...
                                                            stringer.flange_width,...
                                                            stringer.web_height,...
                                                            "Z_dh03",...
                                                            0);



    %% Catchpole diagram calculations
    t_s_over_t = stringer.thickness/booms.skin_thickness;
    h_over_b = stringer.web_height/booms.stringer_pitch;
    K_catchpole = catchpole_calculator(h_over_b, t_s_over_t, booms.material.v,1);
    sigma_cr = K_catchpole * booms.material.E*((booms.skin_thickness/booms.stringer_pitch)^2);

    %set euler critical buckling to occur at same stress for best efficency
    frames.L = Kna*pi*sqrt(booms.material.E/sigma_cr);
   

    cf = 1/16000;
    frames.required_stiffness = cf*booms.max_bending_moment*(df^2)/frames.L;
    frames.required_I = frames.required_stiffness/frame_material.E;

    if strcmp(sectionType,'C')
        t_eq = @(b,h) frames.required_I./(((h.^3)/12) + ((b.*(h.^2))/2));
        area_eq = @(t,b,h) (2*b + h).*t;
        obj = @(b,h) area_eq(t_eq(b,h),b,h);
    else
        error("Invalid sectionType");
    end

    %manufacturing constraints -> need to come back to this later
    % hmin = 1e-3;
    % hmax = 0.1;
    hmin = 1.2 * stringer.web_height; % constrained by stringer height -> at least 1.2* stringer height
    hmax = 5e-2; % od - id for fuselage

    bmin = 1e-3;
    bmax = 0.1;
    tmin = 1e-3;
    tmax = 0.1;

    function [c,ceq] = ix_constraint(x)
        c = [];
        ceq = [t_eq(x(2),x(3)) - x(1)];
    end

    min_area_result = fmincon(@(x) area_eq(x(1),x(2),x(3)),[.01 .01 .01],[],[],[],[],[hmin,bmin,tmin],[hmax,bmax,tmax],@ix_constraint);
    
    % min_area_result = fmincon(@(x) obj(x(1),x(2)),[.01 .01],[],[],[],[],[hmin bmin],[hmax,bmax]);
    frames.t = min_area_result(1);
    frames.b = min_area_result(2);
    frames.h = min_area_result(3);
    frames.csa = area_eq(frames.t,frames.b,frames.h);
    frames.weight_per_frame = frames.csa * (pi*df) * frame_material.rho;
    frames.number = floor(lf/frames.L); % estimate number of light frames
    frames.weight = frames.weight_per_frame * frames.number;

    if doPlot
        %plot function space
        fig = figure("Name",'light_frame_ix_optimization');
        b_range = linspace(bmin,bmax,1000);
        h_range = linspace(hmin,hmax,1000);
        [B,H] = meshgrid(b_range,h_range);
        T = t_eq(B,H);
        Area = area_eq(T,B,H);
        hold on
        s = surf(B,H,T,Area);
        s.EdgeColor = 'interp';
        xlabel('Section Width [m]')
        ylabel('Section Height [m]')
        zlabel('Section Thickness [m]')
        plot3(frames.b,frames.h,frames.t,'x','color','red','MarkerSize',10,'LineWidth',5)
        c = colorbar;
        c.Label.String = 'Section Area [m^2]';
        colormap(turbo)
        grid on
        legend("Solution Surface","Optimum Solution")
        xlim([0 0.02])
        ylim([0 0.02])
        view([130 20])
    end
end
