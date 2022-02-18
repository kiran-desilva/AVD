clear
clc

%SI UNITS

materialLib{1}.name = 'Al 2024';
materialLib{1}.E = 73e9; %N/m^2
materialLib{1}.v = 1/3;
materialLib{1}.G = 27.444e9; %N/m^2
materialLib{1}.K = 71.569e9; %N/m^2
materialLib{1}.tensile_yield = 324e6 %Pa
materialLib{1}.shear_yield = 187.1e6; %Pa
materialLib{1}.rho = 2.73e3; %kg/m^3



save('materialLib','materialLib')