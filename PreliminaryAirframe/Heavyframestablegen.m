clear
clc

fig_path = fullfile('./Tables');

load('HeavyFrames.mat');
heavyframes.Labels = {'Heavy Frame','Mass (Kg)','Web Thickness (m)','Flange Thickness','Web Height (m)','Flange Length (m)'};

framenames = {'Wing Front Spar'; 'Wing Rear Spar'; 'Tail Front Spar'; 'Tail Rear Spar'};
masses = [HeavyFrames.wing.massfront; HeavyFrames.wing.massrear; HeavyFrames.engine.mass; HeavyFrames.tail.mass];
tweb = [HeavyFrames.wing.twFRONT; HeavyFrames.wing.twREAR; HeavyFrames.engine.t; HeavyFrames.tail.tvREAR];
lflange = [HeavyFrames.wing.lfwFRONT; HeavyFrames.wing.lfwREAR; HeavyFrames.engine.lf; HeavyFrames.tail.lfvREAR];
hweb = [HeavyFrames.wing.Hfront; HeavyFrames.wing.Hrear; HeavyFrames.engine.H; HeavyFrames.tail.H];
%tweb = [HeavyFrames.wing.twFRONT; HeavyFrames.wing.twREAR; HeavyFrames.engine.t; HeavyFrames.tail.tvREAR];
tflange = [HeavyFrames.wing.tffront; HeavyFrames.wing.tfrear; HeavyFrames.engine.tf; HeavyFrames.tail.tf];



  %  heavyframes.Data = [framenames,masses,tweb,tflange,hweb,lflange];
    
heavyframes.table = table(framenames,masses,tweb,tflange,hweb,lflange)
heavyframes.table.Properties.VariableNames = heavyframes.Labels;
table2latex(heavyframes.table,'Tables/heavyframes.tex');
 
 