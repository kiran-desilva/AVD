function [fig] = hdl(struct)

    fig = figure
    hold on
    % x y
    rhs_wing = [0,0;
              0,struct.Croot;
              struct.b/2,struct.Croot - (0.5*struct.b * tand(struct.sweepLE));
              struct.b/2,(struct.Croot - (0.5*struct.b * tand(struct.sweepLE)) - struct.Ctip);
              0,0;
    ]
    lhs_wing = [0,0;
              0,struct.Croot;
              -struct.b/2,struct.Croot - (0.5*struct.b * tand(struct.sweepLE));
              -struct.b/2,(struct.Croot - (0.5*struct.b * tand(struct.sweepLE)) - struct.Ctip);
              0,0;
    ]

    plot(rhs_wing(:,1),rhs_wing(:,2));
    plot(lhs_wing(:,1),lhs_wing(:,2));
    % plot(lhs_wing,'x');
    xlabel('x coordinate (m)')
    ylabel('y coordinate (m)')
    title("Horizontal Tailplane")
  
    hold on
    % x y
    rhs_wing_hdl = [0,0;
              0,struct.HDL_Croot;
              struct.b/2,(struct.HDL_Croot) - (0.5*struct.b * tand(struct.sweepHDL));
              struct.b/2,(struct.HDL_Croot - (0.5*struct.b * tand(struct.sweepHDL)) - struct.HDL_Ctip);
              0,0;
    ]
    lhs_wing_hdl = [0,0;
              0,struct.HDL_Croot;
              -struct.b/2,struct.HDL_Croot - (0.5*struct.b * tand(struct.sweepHDL));
              -struct.b/2,(struct.HDL_Croot - (0.5*struct.b * tand(struct.sweepHDL)) - struct.HDL_Ctip);
              0,0;
    ]

    plot(rhs_wing_hdl(:,1),rhs_wing_hdl(:,2));
    plot(lhs_wing_hdl(:,1),lhs_wing_hdl(:,2));
    % plot(lhs_wing,'x');
    xlabel('x coordinate (m)')
    ylabel('y coordinate (m)')
    title("Wing")
    grid on
    axis equal
    hold off
end