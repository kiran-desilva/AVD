%wraps fuselage_generate and changes an error weight to a big number so the correct gradient can be calculated
%hides plotting option too
%implements limits aswell


function [total_weight] = fminopt_fuselage_generate(material,loadcase,n_stringer,stringer,fuselage_thickness,limits)

    if (n_stringer < limits.n_stringer_min) || (stringer.thickness < limits.stringer_thickness_min) || (stringer.web_height < limits.stringer_height_min) || (fuselage_thickness < limits.skin_min)
        total_weight = 2000000000;
       disp("out of range")
       return;
    end
    %enforce integer of n_stringer
    n_stringer = round(n_stringer);


    [~,total_weight] = fuselage_generate(material,loadcase,n_stringer,stringer,fuselage_thickness,0);
    
    if total_weight < 0
        total_weight = 2000000000;
    end

end