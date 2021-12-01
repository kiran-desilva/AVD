function D = Drag_model(M,h,W)
%h in meters
    load('aero_analysis.mat');
    load('wing.mat');
    v = M*atmos(h,2);
    CD = aero_analysis.drag.cd0(1) + ((W / (0.5 * atmos(h,4) * v^2 * wing.Sref))^2) / (pi * wing.AR * aero_analysis.e);
end