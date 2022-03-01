clear
clc

%SI UNITS

%for upper wing skin-stringer, ribs, fuselage stringers and frames
materialLib{1}.name = 'Al 7075';
materialLib{1}.E = 69e9; %N/m^2 (69-76)
materialLib{1}.v = 0.325; %(0.325-0.335)
materialLib{1}.G = 26e9; %N/m^2 (26-28)
materialLib{1}.K = 67e9; %N/m^2 (67-74)
materialLib{1}.tensile_strength = 530e6; %Pa (530-580)
materialLib{1}.yield_strength = 460e6; %Pa (460-530)
materialLib{1}.shear_yield = 0.5*460e6; %Pa
materialLib{1}.rho = 2.77e3; %kg/m^3 (2.77-2.82)

%for lower wing skin-stringer
materialLib{2}.name = 'Al 7475';
materialLib{2}.E = 69e9; %N/m^2 (69-72.5)
materialLib{2}.v = 0.33; % (0.33-0.343)
materialLib{2}.G = 26e9; %N/m^2 (26-27.3)
materialLib{2}.K = 67.6e9; %N/m^2 (67.6-71.1)
materialLib{2}.tensile_yield =  %Pa
materialLib{2}.shear_yield =  %Pa
materialLib{2}.rho = 2.78e3; %kg/m^3 (2.78-2.81)

%for spars
materialLib{3}.name = 'Al 7175';
materialLib{3}.E = 70e9; %N/m^2 (70-73.6)
materialLib{3}.v = 0.33; % (0.33-0.343)
materialLib{3}.G = 27e9; %N/m^2 (27-28.4)
materialLib{3}.K = 69e9; %N/m^2 (69-72.5)
materialLib{3}.tensile_yield =  %Pa
materialLib{3}.shear_yield =  %Pa
materialLib{3}.rho = 2.78e3; %kg/m^3 (2.78-2.81)

%for fuselage skin
materialLib{4}.name = 'Al 2524 T3';
materialLib{4}.E = 71e9; %N/m^2 (71-74.6)
materialLib{4}.v = 0.35; %(0.35-0.364)
materialLib{4}.G = 28e9; %N/m^2 (28-29.4)
materialLib{4}.K = 78.9e9; %N/m^2 (78.9-82.9)
materialLib{4}.tensile_yield =  %Pa
materialLib{4}.shear_yield =  %Pa
materialLib{4}.rho = 2.75e3; %kg/m^3 (2.75-2.78)


save('materialLib','materialLib')