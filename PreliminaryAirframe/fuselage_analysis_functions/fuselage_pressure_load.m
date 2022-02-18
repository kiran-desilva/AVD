%Calculated the required thicknesses to resist the pressure load at 45000ft with a cabin pressurization of 8000ft
function [t1_t2,t_min_hoop,t_min_long,t_min_cap] = pressure_load(material)

    abs_ceiling_alt = 45000;
    cabin_alt = 8000;
    abs_ceiling_alt_m = distdim(abs_ceiling_alt,'ft','m');
    cabin_alt_m = distdim(cabin_alt,'ft','m');

    [~,~,abs_ceiling_p,~] = atmosisa(abs_ceiling_alt_m);
    [~,~,cabin_p,~] = atmosisa(cabin_alt_m);

    delta_p = (cabin_p - abs_ceiling_p); %Pa

    df = 1.68 ; %m fuselage diameter

    t1_t2 = (2-material.v)/(1-material.v); %end cap skin thickness ratio

    sigma_allow = material.tensile_yield;

    t_min_hoop = (df/(sigma_allow/delta_p)) * 0.5;

    t_min_long = (df/(sigma_allow/delta_p)) * 0.25; % this will always be smaller than hoop stress thickness

    t_min_cap = (df/(sigma_allow/delta_p)) * 0.25;


    t_dome_compatibility = t_min_hoop/t1_t2 ; 


end