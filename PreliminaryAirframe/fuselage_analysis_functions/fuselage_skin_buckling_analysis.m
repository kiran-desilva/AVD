function [sigma_crit_c,sigma_crit_s] = fuselage_skin_buckling_analysis(booms,material)
    effective_distance = booms.stringer_pitch - (15*skin_thickness);
    assert(effective_distance>1,"Error effective distance negative");
    %these depend on stringer dimentions? -> might need to interpolate esdu graphs??
    Kc = 3.62;
    Ks = 5;
    sigma_crit_func = @(K) K*material.E*((booms.skin_thickness/effective_distance)^2);
    sigma_crit_c = sigma_crit_func(Kc);
    sigma_crit_s = sigma_crit_func(Ks);

end