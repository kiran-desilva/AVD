%currently assuming same material for whole fuselage -> this can be easily changed later
%%stringer.thickness
%%stringer.flange_to_web_ratio
%%stringer.web_height

function [fuselage,total_weight,figlist] = fuselage_generate(material,loadcase,n_stringer,stringer,fuselage_thickness,doplot)
    % get max shear flow
    figlist = [];
    [~,shear_flow,~,~,fig] = fuselage_shear_flow_analysis(loadcase,material,doplot);
    figlist = [fig figlist];

    fuselage.max_shear_stress = shear_flow/fuselage_thickness;

    try
        [fuselage.stringerpanel.weight,fuselage.stringerpanel.booms,fuselage.max_principle_stress,~,fig] = fuselage_bending_analysis(loadcase,fuselage_thickness,n_stringer,stringer,material,doplot);
        figlist = [fig figlist];
    catch e
        if strcmp(e.message,'TENSILEYIELD')
            disp('Fuselage Tensile Yield')
            fuselage.total_weight = -1;
            total_weight=fuselage.total_weight;
            return;
        else
            rethrow(e)
        end
    end

    try
        [fuselage.sigma_crit_c,fuselage.sigma_crit_s] = fuselage_skin_buckling_analysis(fuselage.stringerpanel.booms);
    catch exception
        if strcmp(exception.message,'EFFECTIVEDISTANCE')
            disp("effective distance negative")
            fuselage.total_weight = -2;
            total_weight=fuselage.total_weight;
            return;
        else
            rethrow(exception)
        end
    end

    % check mixed mode buckling criteron

    Rcs = (fuselage.max_principle_stress/fuselage.sigma_crit_c) + (fuselage.max_shear_stress/fuselage.sigma_crit_s)^2;
    if (Rcs > 1)
        disp("Fuselage Buckled eek")
        fuselage.total_weight = -3;
        total_weight=fuselage.total_weight;
        return;
    end
   
    try
        [fuselage.frames,fig] = fuselage_light_frames(material,'C',fuselage.stringerpanel.booms,stringer,doplot);
    catch e
        if contains(e.message,'Tried to sample catchpole diagram out of distribution') % hack
            disp('catchpole out of range')
            fuselage.total_weight = -4;
            total_weight=fuselage.total_weight;
            return;
        else
            rethrow(e)
        end
    end

    figlist = [fig figlist];
    
    
    fuselage.total_weight = fuselage.frames.weight + fuselage.stringerpanel.weight;
    total_weight=fuselage.total_weight;




