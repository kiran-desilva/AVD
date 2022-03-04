clear
clc

fig_path = fullfile('./Tables');

load('HeavyFrames.mat');
heavyframes.Labels = {'Heavy Frame','Mass (Kg)','t_{Web} (m)','t_{Flange}','l_{Web}(m)','l_{Flange} (m)'};

framenames = {'Wing Front'; 'Wing Rear'; 'Tail Front'; 'Tail Rear'};
masses = [round(HeavyFrames.wing.massfront, 1); round(HeavyFrames.wing.massrear,1); round(HeavyFrames.engine.mass,1); round(HeavyFrames.tail.mass,1)];
tweb = [round(HeavyFrames.wing.twFRONT,3); round(HeavyFrames.wing.twREAR,3); round(HeavyFrames.engine.t,3); round(HeavyFrames.tail.tvREAR,3)];
lflange = [round(HeavyFrames.wing.lfwFRONT,3); round(HeavyFrames.wing.lfwREAR,3); round(HeavyFrames.engine.lf,3); round(HeavyFrames.tail.lfvREAR,3)];
hweb = [round(HeavyFrames.wing.Hfront,3); round(HeavyFrames.wing.Hrear, 3); round(HeavyFrames.engine.H,3); round(HeavyFrames.tail.H,3)];
%tweb = [HeavyFrames.wing.twFRONT; HeavyFrames.wing.twREAR; HeavyFrames.engine.t; HeavyFrames.tail.tvREAR];
tflange = [round(HeavyFrames.wing.tffront,3); round(HeavyFrames.wing.tfrear,3); round(HeavyFrames.engine.tf,3); round(HeavyFrames.tail.tf,3)];



  %  heavyframes.Data = [framenames,masses,tweb,tflange,hweb,lflange];
    
heavyframes.table = table(framenames,masses,tweb,tflange,hweb,lflange)
heavyframes.table.Properties.VariableNames = heavyframes.Labels;
table2latex(heavyframes.table,'Tables/heavyframes.tex');
 
 