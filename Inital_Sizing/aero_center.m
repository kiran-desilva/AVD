function [xac] = aero_center(Ybar,Cmac,sweepLE)
    xac = (Ybar*tand(sweepLE) + (0.25*Cmac));
end