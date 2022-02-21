clear
clc

addpath(fullfile('.','fuselage_analysis_functions'))
addpath(fullfile('.','helper_funcs'))

load fuselageLoading
load materialLib

%get limiting loadcase
loadcase = fuselageLoading.Vd_flight;
material = materialLib{1};

%determine min skin thickness 

[t1_t2,t_min_hoop,t_min_long,t_min_cap] = fuselage_pressure_load(material);
pressure_min = max([t_min_hoop,t_min_long]);
%determine min thicknes due to shear flow yield
[shear_min] = fuselage_shear_flow_analysis(loadcase,material,0);

skin_min = max([shear_min,pressure_min]);
skin_max = 1e-2;
% skin_thickness_range = generate_range(skin_min,skin_max);

stringer_thickness_min = 1e-3;
stringer_thickness_max = 1e-2;
% stringer_thickness_range = generate_range(stringer_thickness_min,stringer_thickness_max);

stringer_height_min = 1e-2;
stringer_height_max = 1e-1;
% stringer_height_range = generate_range(stringer_height_min,stringer_height_max);

n_stringer_min = 10;
n_stringer_max = 100;
n_stringer_range = [4:2:100];

create_stringer = @(ts,hs) struct()