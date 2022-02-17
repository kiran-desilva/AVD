%%stringer
%%stringer.thickness
%%stringer.flange_to_web_ratio
%%stringer.web_height
%%
%%
%%%%
function [] = fuselage_bending_analysis(loadcase,fuselage_skin_thickness,n_stringers,material,stringer,doplot)

    df = 1.68;
    cf = pi*df; %fuselage circumference
    stringer_pitch = cf/(n_stringers+1); %n+1 to account for wrap around at phi = 0
    stringer_angular_pitch = (2*pi)/(n_stringers+1);
    
    skin_equivalent_boom_area = 15*fuselage_skin_thickness^2; %15t^2

    max_yield_stress = material.tensile_yield;



    if doplot
        figure

    end