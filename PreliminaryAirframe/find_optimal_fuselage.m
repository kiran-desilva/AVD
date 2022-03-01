clear
clc

addpath(fullfile('.','fuselage_analysis_functions'))
addpath(fullfile('.','helper_funcs'))

load fuselageLoading
load materialLib

%get limiting loadcase
loadcase = fuselageLoading.Vd_flight;
material = materialLib{4};

%determine min skin thickness 

[t1_t2,t_min_hoop,t_min_long,t_min_cap] = fuselage_pressure_load(material);
pressure_min = max([t_min_hoop,t_min_long]);
%determine min thicknes due to shear flow yield
[shear_min] = fuselage_shear_flow_analysis(loadcase,material,0);

limits.skin_min = max([shear_min,pressure_min 1e-3]);

limits.stringer_thickness_min = 1e-3;

limits.stringer_height_min = 1e-3;

limits.n_stringer_min = 1;

create_stringer = @(ts,hs,dh) struct('thickness',ts,'web_height',hs,'flange_to_web_ratio',dh);

% options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-12,'PlotFcns',{@optimplotfval,@param_plot});
% options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-12,'PlotFcns',{@optimplotfval,@param_plot},'UseParallel',true);
% options = optimset('Display','iter','InitialPopulation',{});
options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20,'UseParallel',true);
% options = optimoptions('fminunc','Display','iter','PlotFcns',{@optimplotfval,@optimplotx},'UseParallel',true);
% options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20,'UseParallel',true);
% x0 = [50,1.1e-3,0.02,0.0013]; %set x0 to the optimal found during the grid search

% x0 = [6e-3,1.5e-3,1.5e-2,1.1e-3]; %set x0 to the optimal found during the grid search
x0 = [7e-3,1.5e-3,1.5e-2,1.1e-3]; % note stringer panel number has been scaled so is a similar magnitude to other parameters

% options = optimoptions('ga','Display','iter','InitialPopulation',x0,'UseParallel',true,'PlotFcn',@gaplotbestf,'MutationFcn','mutationgaussian');
% x0 = [6.2e-3,1e-3,1.05e-2,1.2e-3];
% x0 = [6.5e-3,1.2e-3,1.444e-2,1.1e-3]; 
res = [];

% ms = MultiStart('UseParallel',true,'Display','iter');
gs = GlobalSearch('Display','iter','PlotFcn',@gsplotbestf);
problem = createOptimProblem('fmincon','objective',...
                            @(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),...
                            'x0',x0,...
                            'lb',[10e-4,limits.stringer_thickness_min,limits.stringer_height_min,limits.skin_min],...
                            'options',options);

% problem = createOptimProblem('fminsearch','objective',...
%                             @(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),...
%                             'x0',x0,...
%                             'options',options);

[res,f,exitflag,output,solutions] = run(gs,problem)

% while true
%     % [res,fval,exit,out] = fminsearch(@(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),x0,options)
%     % [res,fval,exit,out] = fminunc(@(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),x0,options)
%     % [res,fval,exit,out] = ga(@(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),4,[],[],[],[],[],[],[],options)
%     [res,fval,exit,out] = fmincon(@(x) fminopt_fuselage_generate(material,loadcase,x(1)*(1e4),create_stringer(x(2),x(3),0.3),x(4),limits),x0,[],[],[],[],[10e-4,limits.stringer_thickness_min,limits.stringer_height_min,limits.skin_min],[],[],options)
%     fuselage_generate(material,loadcase,round(res(1)*(1e4)),...
%                                                       create_stringer(res(2),res(3),0.3),...
%                                                       res(4),...
%                                                       0);
%     ui = input("Press Enter to Continue, Press E to finish",'S');
%     if strcmp(ui,"E")
%         break;
%     else
%         x0 = res;
%     end
% end

optimal.n_stringer = round(res(1)*(1e4));
optimal.ts = res(2);
optimal.hs = res(3);
optimal.fs = res(4);
optimal.fuselage = fuselage_generate(material,loadcase,optimal.n_stringer,...
                                                      create_stringer(optimal.ts,optimal.hs,0.3),...
                                                      optimal.fs,...
                                                      0);
optimal.fuselage.frames
                
save fuselageOptimizerOutput

function stop = param_plot(x,optimValues,state,varargin)
    persistent lines;

    stop = false;
   switch state
    case 'init'
        hold on
        grid on
        iteration = optimValues.iteration;
        lines{1} = plot(iteration,x(1));
        lines{2} = plot(iteration,x(2));
        lines{3} = plot(iteration,x(3));
        lines{4} = plot(iteration,x(4));
    case 'iter'
        iteration = optimValues.iteration;
        lines{1}.XData(end+1) = iteration;
        lines{2}.XData(end+1) = iteration;
        lines{3}.XData(end+1) = iteration;
        lines{4}.XData(end+1) = iteration;

        lines{1}.YData(end+1) = x(1);
        lines{2}.YData(end+1) = x(2);
        lines{3}.YData(end+1) = x(3);
        lines{4}.YData(end+1) = x(4);
    end
end

        
