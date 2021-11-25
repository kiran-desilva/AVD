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
 
  
    hold on
    % x y
 
    aileron_start=0.7*struct.b/2
    C_HDL_root=((2*(struct.Ctip-struct.Croot))/struct.b)*(struct.HDL_start) +struct.Croot
    
    
    rhs_wing_hdl = [struct.HDL_start,-struct.HDL_start*tand(struct.sweepTE); struct.HDL_start,(-struct.HDL_start*tand(struct.sweepTE))+struct.C_HDL_root;
              aileron_start,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot;
              aileron_start,-aileron_start*tand(struct.sweepTE); 
              0,0;
    ]
    lhs_wing_hdl = [-struct.HDL_start,-struct.HDL_start*tand(struct.sweepTE); -struct.HDL_start,(-struct.HDL_start*tand(struct.sweepTE))+struct.C_HDL_root;
              -aileron_start,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot;
              -aileron_start,-aileron_start*tand(struct.sweepTE); 
              0,0;
    ]

    plot(rhs_wing_hdl(:,1),rhs_wing_hdl(:,2));
    plot(lhs_wing_hdl(:,1),lhs_wing_hdl(:,2));
    % plot(lhs_wing,'x');
    
    
    
     hold on
    % x y
    
    
    rhs_wing_aileron = [aileron_start,-aileron_start*tand(struct.sweepTE); aileron_start,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot; struct.b/2,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot-(0.5*(struct.b-0.7*struct.b)*tand(struct.sweepHDL));
                        struct.b/2,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot-(0.5*(struct.b-0.7*struct.b)*tand(struct.sweepHDL))-struct.HDL_Ctip;
                        aileron_start,-aileron_start*tand(struct.sweepTE)
              
    ]
    lhs_wing_aileron = [-aileron_start,-aileron_start*tand(struct.sweepTE); -aileron_start,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot;-struct.b/2,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot-(0.5*(struct.b-0.7*struct.b)*tand(struct.sweepHDL));
                        -struct.b/2,(-aileron_start*tand(struct.sweepTE))+struct.aileron_Croot-(0.5*(struct.b-0.7*struct.b)*tand(struct.sweepHDL))-struct.HDL_Ctip;
                        -aileron_start,-aileron_start*tand(struct.sweepTE)
    ]

    plot(rhs_wing_aileron(:,1),rhs_wing_aileron(:,2));
    plot(lhs_wing_aileron(:,1),lhs_wing_aileron(:,2));
    % plot(lhs_wing,'x');
    
    fues=1.52/2 
    plot([fues,fues],[-2,3])
    plot([-fues,-fues],[-2,3])
    xlabel('x coordinate (m)')
    ylabel('y coordinate (m)')
    title("Wing")
    grid on
    axis equal
    hold off
   
end