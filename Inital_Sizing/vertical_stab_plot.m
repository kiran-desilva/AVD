function [fig] = vertical_stab_plot(struct)
    fig = figure
    hold on
    % z x
    wing = [0,0;
            0,struct.Croot;
            struct.b,struct.Croot - (struct.b * tand(struct.sweepLE));
            struct.b,(struct.Croot - (struct.b * tand(struct.sweepLE)) - struct.Ctip);
            0,0;
    ]


    plot(wing(:,2),wing(:,1));
    xlabel('x coordinate (m)')
    ylabel('z coordinate (m)')
    title("Vertical Tailplane")
    grid on
    axis equal
    hold off
end