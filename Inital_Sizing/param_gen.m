%key parameters


parameters.cruise_range_km = 2500;%km
parameters.cruise_mach = 0.75;
parameters.cruise_alt_ft = 40000;%ft amsl
parameters.alternate_range_km = 370;%km
parameters.loiter_alt_ft = 5000;
parameters.loiter_duration = 45; %minutes

parameters.runway_length = 1200; 
parameters.max_mach = 0.78;
parameters.absolute_Ceiling = 45000; %amsl
parameters.max_cabin_press = 8000;

parameters.engine_number = 2;
parameters.rho_0 = 1.225;
parameters.target_landing_weight = 0.85;



save("parameters","parameters");