%key parameters


parameters.cruise_range_km = 2500;%km
parameters.cruise_mach = 0.75;
parameters.cruise_alt_ft = 40000;%ft amsl
parameters.alternate_range_km = 370;%km
parameters.loiter_alt_ft = 5000;
parameters.loiter_duration = 45; %minutes

parameters.min_runway = 1.2; 
parameters.max_mach = 0.78;
parameters.max_alt_ft = 45000; %amsl
parameters.max_cabin_press = 8000;

save("parameters","parameters");