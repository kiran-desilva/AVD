clear
clc

%SI UNITS

materialLib{1}.name = 'Al 2024';
materialLib{1}.E = 73000; %N/mm^2
materialLib{1}.v = 0.33;
materialLib{1}.G = 27444; %N/mm^2
materialLib{1}.K = 71569; %N/mm^2
materialLib{1}.tensile_yield = 324 %Mpa
materialLib{1}.shear_yield = 187.1; %Mpa
materialLib{1}.rho = 2.73e3; %kg/m^3

save('materialLib','materialLib')