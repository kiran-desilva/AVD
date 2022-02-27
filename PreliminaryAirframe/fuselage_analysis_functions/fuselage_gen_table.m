function [] = fuselage_gen_table(fuselage)
    StringerPanel.Labels = {'Weight [N]';
                            'Fuselage Skin Thickness [m]';
                            'Number of Stringers';
                            'Stringer Thickness [m]';
                            'Stringer Web Height [m]';
                            'Stringer Web to Flange Ratio'};
    StringerPanel.Data = [fuselage.stringerpanel.weight;
                        fuselage.stringerpanel.booms.skin_thickness;
                        length(fuselage.stringerpanel.booms.phi);
                        fuselage.stringer.thickness;
                        fuselage.stringer.web_height;
                        fuselage.stringer.flange_to_web_ratio];
    
    StringerPanel.table = table(StringerPanel.Data)
    StringerPanel.table.Properties.RowNames = StringerPanel.Labels;
    table2latex(StringerPanel.table,'Tables/fuselage_stringerpanel_table.tex');
