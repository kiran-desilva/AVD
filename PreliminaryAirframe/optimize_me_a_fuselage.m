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

%setting boundaries for search
generate_range = @(min,max) linspace(min,max,10);

skin_min = max([shear_min,pressure_min]);
skin_max = 1e-2;
skin_thickness_range = generate_range(skin_min,skin_max);

stringer_thickness_min = 1e-3;
stringer_thickness_max = 1e-2;
stringer_thickness_range = generate_range(stringer_thickness_min,stringer_thickness_max);

stringer_height_min = 1e-2;
stringer_height_max = 1e-1;
stringer_height_range = generate_range(stringer_height_min,stringer_height_max);

n_stringer_min = 10;
n_stringer_max = 100;
n_stringer_range = [4:2:50];

[FS,TS,HS,N] = ndgrid(skin_thickness_range,stringer_thickness_range,stringer_height_range,n_stringer_range); %create indivdual vectors for easy plotting
[X] = ndgrid(skin_thickness_range,stringer_thickness_range,stringer_height_range,n_stringer_range);
fuselages = num2cell(repmat(struct(),size(X)));
weights = zeros(size(X));

skin_thickness_range_length = numel(skin_thickness_range);
stringer_thickness_range_length = numel(stringer_thickness_range);
stringer_height_range_length = numel(stringer_height_range);
n_stringer_range_length = numel(n_stringer_range);

parfor i_fs = 1:skin_thickness_range_length
    for i_ts = 1:stringer_thickness_range_length
        for i_hs = 1:stringer_height_range_length
            for i_n = 1:n_stringer_range_length
                skin_thickness = skin_thickness_range(i_fs);
                stringer = struct();
                stringer.flange_to_web_ratio = 0.3;
                stringer.thickness = stringer_thickness_range(i_ts);
                stringer.web_height = stringer_height_range(i_hs);
                N = n_stringer_range(i_n);
                [fuse,weight] = fuselage_generate(material,loadcase,N,stringer,skin_thickness,0);
                weights(i_fs,i_ts,i_hs,i_n) = weight;
                fuselages{i_fs,i_ts,i_hs,i_n} = fuse;
            end
        end
    end
end








save fuselageStructure