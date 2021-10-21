clear 
clc

p = open('parameters.mat');
r = open('roskamdata.mat');

roskam.c = 0.7; %ish? 1/hr
roskam.LoverD_cruise = 11; %ish?
roskam.LoverD_loiter = 13; %ish?
roskam.Wx_W0 = 1; %initialize
roskam.A = 0.2678;
roskam.B = 0.9979;

raymer.c = 0.7; %ish? 1/hr
raymer.LoverD_cruise = 11; %ish?
raymer.LoverD_loiter = 13; %ish?
raymer.Wx_W0 = 1; %initialize
raymer.A = 1.51;
raymer.C = -0.1;

W_crew = 100 * 6 * 9.81; %???
W_pld = 20 * 6 * 9.81; %???
initialW0 = 80000;

%fuel mass fractions

roskam.fuelfrac(1) = 0.99 * 0.995 * 0.995; %engine start + taxi + takeoff
roskam.fuelfrac(2) = 0.98;    %climb
roskam.fuelfrac(3) = exp(-(p.parameters.cruise_range_km * 1000 * roskam.c / 3600) / (p.parameters.cruise_mach * 295.07 * roskam.LoverD_cruise));  %cruise 1 using breguet range
roskam.fuelfrac(4) = 0.99;  %descent 1
roskam.fuelfrac(5) = 0.98;  %climb and accelerate note assumed same as prev climb
roskam.fuelfrac(6) = exp(-(p.parameters.alternate_range_km * 1000 * roskam.c / 3600) / (p.parameters.cruise_mach * 295.07 * roskam.LoverD_cruise));  %cruise 2
roskam.fuelfrac(7) = exp(-(p.parameters.loiter_duration * 60 * roskam.c / 3600) / (roskam.LoverD_loiter)); % loiter using endurance eqn
roskam.fuelfrac(8) = 0.99; %descent 2
roskam.fuelfrac(9) = 0.992;    %landing + taxi

raymer.fuelfrac(1) = 0.97; %engine start + taxi + takeoff
raymer.fuelfrac(2) = 0.985;    %climb
raymer.fuelfrac(3) = exp(-(p.parameters.cruise_range_km * 1000 * raymer.c / 3600) / (p.parameters.cruise_mach * 295.07 * raymer.LoverD_cruise));  %cruise 1 using breguet range
raymer.fuelfrac(4) = 0.99;  %descent 1
raymer.fuelfrac(5) = 0.985;  %climb and accelerate note assumed same as prev climb
raymer.fuelfrac(6) = exp(-(p.parameters.alternate_range_km * 1000 * raymer.c / 3600) / (p.parameters.cruise_mach * 295.07 * raymer.LoverD_cruise));  %cruise 2
raymer.fuelfrac(7) = exp(-(p.parameters.loiter_duration * 60 * raymer.c / 3600) / (raymer.LoverD_loiter)); % loiter using endurance eqn
raymer.fuelfrac(8) = 0.99; %descent 2
raymer.fuelfrac(9) = 0.995;    %landing + taxi

for i = 1:length(roskam.fuelfrac)
    roskam.Wx_W0 = roskam.Wx_W0 * roskam.fuelfrac(i); %multiply fuel mass fractions 
end

roskam.Wf_W0 = 1.01 * (1 - roskam.Wx_W0);   %calculate total fuel mass fraction for mission


for i = 1:length(raymer.fuelfrac)
    raymer.Wx_W0 = raymer.Wx_W0 * raymer.fuelfrac(i); %multiply fuel mass fractions 
end

raymer.Wf_W0 = 1.01 * (1 - raymer.Wx_W0);   %calculate total fuel mass fraction for mission

%Roskam W0

roskam.W0 = initialW0; %initial W0 guess (N)
roskam.W0prev = 0; %initialise
count = 0; 
figure 
hold on
error = 10; 

while abs(error) > 1e-10
    roskam.W0prev = roskam.W0; 
    We_W0_roskam_regress = (10^((log10(roskam.W0) - roskam.A) / roskam.B)) / roskam.W0;   %Roskam regression method
    roskam.W0 = (W_crew + W_pld) / (1 - roskam.Wf_W0 - (We_W0_roskam_regress));    %W0 calculation 
    count = count + 1;  %count number of iterations
    error = roskam.W0 - roskam.W0prev;    %Calculate difference between consecutive W0 values
    plot(count, error, 'b*') %plot error 
    grid on
    pause(0.05)
end  
hold off

roskam.We_W0 = We_W0_roskam_regress; 
roskam.W0 = roskam.W0 * 0.2248; %N to ib force

%Raymer W0

raymer.W0 = initialW0; %initial W0 guess (N)
raymer.W0prev = 0; %initialise
raymer.Kvs = 1; %fixed sweep

count = 0; 
figure 
hold on
error = 10; 

while abs(error) > 1e-10
    raymer.W0prev = raymer.W0; 
    We_W0_raymer_regress = raymer.A * (raymer.W0 ^ raymer.C) * raymer.Kvs;   %Raymer regression method
    raymer.W0 = (W_crew + W_pld) / (1 - raymer.Wf_W0 - (We_W0_raymer_regress));    %W0 calculation 
    count = count + 1;  %count number of iterations
    error = raymer.W0 - raymer.W0prev;    %Calculate difference between consecutive W0 values
    plot(count, error, 'b*') %plot error 
    grid on
    pause(0.05)
end  
hold off

raymer.We_W0 = We_W0_raymer_regress; 
raymer.W0 = raymer.W0 * 0.2248; %N to ib force

figure
plot(roskam.W0, roskam.We_W0,'r*')
hold on
plot(raymer.W0, raymer.We_W0,'b*')
hold on
plot(r.roskamdata.W0, r.roskamdata.We_W0,'x')
hold off
grid on
xlabel("$W_{0}$ Ibs", 'interpreter', 'Latex','FontSize', 15)
ylabel("$\frac{W_{e}}{W_{0}}$ Ibs", 'interpreter', 'Latex','FontSize', 15)
title("Comparison to Roskam data",'interpreter', 'Latex','FontSize', 15)
