function [figs,fout] = fuselage_gen_figs(fuselage_struct,save)
    total_weight = fuselage_struct.total_weight;
    n_stringer = length(fuselage_struct.stringerpanel.booms.phi);
    stringer = fuselage_struct.stringer;
    fuselage_thickness = fuselage_struct.stringerpanel.booms.skin_thickness;
    material = fuselage_struct.stringerpanel.booms.material;
    loadcase = fuselage_struct.loadcase;
    %regenerate fuselage with parameters to get plots
    [fout,~,figs] = fuselage_generate(material,loadcase,n_stringer,stringer,fuselage_thickness,1);


    if save
        for f = figs
            improvePlot(f)
            saveas(f,"Figures/" + f.Name,"epsc")
        end
    end
    