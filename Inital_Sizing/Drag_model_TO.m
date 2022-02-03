function D = Drag_model_TO(v,h,W)
%h in meters
    load('aero_analysis.mat');
    load('wing.mat');
   CD = aero_analysis.drag.cd0(3) + (((W ./ (0.5 * 1.225 .* v.^2 * wing.Sref)).^2) ./ (pi .* wing.Ar .* aero_analysis.summary.e_wing));
    D = CD .* 0.5 .* atmos(h,4) .* v.^2 * wing.Sref;
end