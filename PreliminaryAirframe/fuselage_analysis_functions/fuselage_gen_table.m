function [StringerPanel] = fuselage_gen_table(fuselage)
    StringerPanel.Labels = {'Weight';
                            'Fuselage Skin Thickness';
                            'Number of Stringers';
                            'Stringer Section Type';
                            'Stringer Thickness';
                            'Stringer Web Height';
                            'Stringer Web to Flange Ratio'};
    StringerPanel.Units = {
        'N';
        'mm';
        '-';
        '-';
        'mm';
        'mm';
        '-';
    };
    StringerPanel.Data = {round(fuselage.stringerpanel.weight,1);
                        round(fuselage.stringerpanel.booms.skin_thickness*(1e3),1);
                        length(fuselage.stringerpanel.booms.phi);
                        'Z';
                        round(fuselage.stringer.thickness*(1e3),1);
                        round(fuselage.stringer.web_height*(1e3),1);
                        fuselage.stringer.flange_to_web_ratio};
    
    StringerPanel.table = table(StringerPanel.Units,StringerPanel.Data)
    StringerPanel.table.Properties.RowNames = StringerPanel.Labels;
    StringerPanel.table.Properties.VariableNames = {'Units','Value'};
    table2latex(StringerPanel.table,'Tables/fuselage_stringerpanel_table_raw.tex');


    LightFrames.Labels = {'Weight';
                        'Spacing';
                        'Number of Frames';
                        'Frame Section Type'
                        'Frame Thickness';
                        'Frame Height';
                        'Frame Width'};

    LightFrames.Units = {
        'N';
        'm';
        '-';
        '-';
        'mm';
        'mm';
        'mm';
    };
    LightFrames.Data = {
        round(fuselage.frames.weight,1);
        fuselage.frames.L;
        fuselage.frames.number;
        'C';
        round(fuselage.frames.t*(1e3),1);
        round(fuselage.frames.h*(1e3),1);
        round(fuselage.frames.b*(1e3),1);

    };
    
    LightFrames.table = table(LightFrames.Units,LightFrames.Data)
    LightFrames.table.Properties.RowNames = LightFrames.Labels;
    LightFrames.table.Properties.VariableNames = {'Units','Value'};
    table2latex(LightFrames.table,'Tables/fuselage_lightframes_table_raw.tex');
