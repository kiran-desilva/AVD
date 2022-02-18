function [required_area] = fuselage_bending_min_stringer_area(loadcase,fuselage_skin_thickness,n_stringers,material)

    df = 1.68;
    cf = pi*df; %fuselage circumference

    booms.phi = linspace(0,2*pi,n_stringers+1);
    booms.phi(end) = []; %remove last element as this overlaps with 0th element
    booms.x = (df/2) * cos(booms.phi);
    booms.y = (df/2) * sin(booms.phi);


    %solve to find the required area so we are at the max yield stress i.e the minimum area allowed
    max_yield_stress = material.tensile_yield;
    %assuming only my applied
    My = max(abs(loadcase.bending_moment));
    max_x = max(abs(booms.x));
    required_iy = (My*max_x)/max_yield_stress;
    
    skin_equivalent_boom_area = 15*fuselage_skin_thickness^2; %15t^2

    unit_iyy = sum(booms.x.^2);
    required_boom_area = required_iy/unit_iyy;
    required_area = required_boom_area - skin_equivalent_boom_area;
   
end