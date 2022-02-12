function [figlist] = fuselage_load_plots(loadcase)
    figlist = [];

    %%Sectional Loading

    figlist(1) = figure;

    subplot(3,1,1)
    hold on

    plot(loadcase.x,loadcase.load)
    % plot(loadcase.Fi_cg,0,'x','color','red')

    xlabel('X coordinate from Tip (m)')
    ylabel('Sectional Loading (N)')
    grid on
    grid minor

    %%Shear force

    % figlist(2) = figure;
    subplot(3,1,2)
    hold on

    plot(loadcase.x,loadcase.shear)

    xlabel('X coordinate from Tip (m)')
    ylabel('Shear Force (N)')
    grid on
    grid minor

    %%Bending Moment
    
    % figlist(3) = figure;
    subplot(3,1,3)
    hold on

    plot(loadcase.x,loadcase.bending_moment)

    xlabel('X coordinate from Tip (m)')
    ylabel('Bending Moment (Nm)')
    grid on
    grid minor