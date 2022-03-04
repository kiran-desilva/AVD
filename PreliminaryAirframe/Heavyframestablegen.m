clear
clc

fig_path = fullfile('./Tables');

load('HeavyFrames.mat');
heavyframes.Labels = {'Heavy Frame','Mass (Kg)','Web Thickness (m)','Flange Thickness','Web Height (m)','Flange Length (m)'};

framenames = {'Wing Front'; 'Wing Rear'; 'Tail Front'; 'Tail Rear'};
masses = [round(HeavyFrames.wing.massfront, 1); round(HeavyFrames.wing.massrear,1); round(HeavyFrames.engine.mass,1); round(HeavyFrames.tail.mass,1)];
tweb = [HeavyFrames.wing.twFRONT; HeavyFrames.wing.twREAR; HeavyFrames.engine.t; HeavyFrames.tail.tvREAR];
lflange = [round(HeavyFrames.wing.lfwFRONT,3); round(HeavyFrames.wing.lfwREAR,3); round(HeavyFrames.engine.lf,3); round(HeavyFrames.tail.lfvREAR,3)];
hweb = [HeavyFrames.wing.Hfront; HeavyFrames.wing.Hrear; HeavyFrames.engine.H; HeavyFrames.tail.H];
%tweb = [HeavyFrames.wing.twFRONT; HeavyFrames.wing.twREAR; HeavyFrames.engine.t; HeavyFrames.tail.tvREAR];
tflange = [HeavyFrames.wing.tffront; HeavyFrames.wing.tfrear; HeavyFrames.engine.tf; HeavyFrames.tail.tf];



  %  heavyframes.Data = [framenames,masses,tweb,tflange,hweb,lflange];
    
heavyframes.table = table(framenames,masses,tweb,tflange,hweb,lflange)
heavyframes.table.Properties.VariableNames = heavyframes.Labels;
table2latex(heavyframes.table,'Tables/heavyframes.tex');
 
 