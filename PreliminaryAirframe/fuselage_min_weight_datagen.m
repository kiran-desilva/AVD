clear
clc


addpath(fullfile('.','fuselage_analysis_functions'))
addpath(fullfile('.','helper_funcs'))

load fuselageLoading
load materialLib

%get limiting loadcase
loadcase = fuselageLoading.Va_flight;
material = materialLib{1};

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

stringer_range = [10:2:100];
min_weight_results = zeros(length(stringer_range),5);

parfor i = 1:length(stringer_range)

    gs = GlobalSearch('Display','iter','PlotFcn',@gsplotbestf);
    options = optimoptions('fmincon','Display','iter','ConstraintTolerance',1e-20,'OptimalityTolerance',1e-20,'UseParallel',true);
    x0 = [1.5e-3,1.5e-2,1.1e-3]; 
    gs = GlobalSearch('Display','iter','PlotFcn',@gsplotbestf);
    problem = createOptimProblem('fmincon','objective',...
                                @(x) fminopt_fuselage_generate(material,loadcase,stringer_range(i),create_stringer(x(1),x(2),0.3),x(3),limits),...
                                'x0',x0,...
                                'lb',[limits.stringer_thickness_min,limits.stringer_height_min,limits.skin_min],...
                                'options',options);
    [res,weight,exitflag,output,solutions] = run(gs,problem);
    min_weight_results(i,:) = [res,stringer_range(i),weight];
end
    
         

save fuselageMinWeightData min_weight_results



        
