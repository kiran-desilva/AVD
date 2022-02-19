function [sigma_crit_c,sigma_crit_s] = fuselage_skin_buckling_analysis(booms)
    effective_distance = booms.stringer_pitch - (15*booms.skin_thickness);
    assert(effective_distance>0,'EFFECTIVEDISTANCE',"Error effective distance negative");
    %these depend on stringer dimentions? -> might need to interpolate esdu graphs??
    % want to find actual esdu grpahs for these coefficents 
    %assuming all bending becomes a compressive stress as bending moment is from center of fuselage
    % so effectivley there is no bending moment applied directly to the plate
    Kc = 3.62;
    Ks = 5;
    sigma_crit_func = @(K) K*booms.material.E*((booms.skin_thickness/effective_distance)^2);
    sigma_crit_c = sigma_crit_func(Kc);
    sigma_crit_s = sigma_crit_func(Ks);

end