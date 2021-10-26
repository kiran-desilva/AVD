clear 
clc

p = open('parameters.mat');
r = open('roskamdata.mat');

sizing.AR = 7.8;
sizing.lambdaLE = 20; %ish? used t/c of abt 0.12-0.14
sizing.e = 4.61 * (1 - 0.045 * sizing.AR ^ 0.68) * ((cosd(sizing.lambdaLE))^0.15) - 3.1; %Raymer eqn from Gud
sizing.CD0 = 0.02; 
sizing.L_Dmax = 0.5 * sqrt((pi * sizing.AR * sizing.e) / sizing.CD0); %LoverDmax

roskam.c_1=0.6197; %cruise
roskam.c_2=0.6197; %alternate cruise
roskam.c_3=0.5191; %loiter
roskam.LoverD_cruise(1) = 11; %ish?
roskam.LoverD_cruise(2) = 0.866 * sizing.L_Dmax; % Raymer 
roskam.LoverD_loiter(1) = 13; %ish?
roskam.LoverD_loiter(2) = sizing.L_Dmax; %raymer?
roskam.Wx_W0(1) = 1; %initialize
roskam.Wx_W0(2) = 1; %initialize
roskam.A = 0.2678;
roskam.B = 0.9979;

raymer.c_1=0.6197; %cruise
raymer.c_2=0.6197; %alternate
raymer.c_3=0.5194; %loDesciter
raymer.LoverD_cruise(1) = 11; %ish?
raymer.LoverD_cruise(2) = 0.866 * sizing.L_Dmax; % Raymer 
raymer.LoverD_loiter(1) = 13; %ish?
raymer.LoverD_loiter(2) = sizing.L_Dmax; %raymer
raymer.Wx_W0(1) = 1; %initialize
raymer.Wx_W0(2) = 1; %initialize
raymer.A = 1.51;
raymer.C = -0.1;

W_crew = 94 * 6 * 2.20462; % lbs
W_pld = 20 * 6 * 2.20462; % lbs
initialW0 = 200000;

%fuel mass fractions


for j = 1:2
        roskam.fuelfrac(1) = 0.99 * 0.995 * 0.995; %engine start + taxi + takeoff
        roskam.fuelfrac(2) = 0.98;    %climb
        roskam.fuelfrac(3) = exp(-(p.parameters.cruise_range_km * 1000 * roskam.c_1 / 3600) / (p.parameters.cruise_mach * 295.07 * roskam.LoverD_cruise(j)));  %cruise 1 using breguet range
        roskam.fuelfrac(4) = 0.99;  %descent 1
        roskam.fuelfrac(5) = 0.98;  %climb and accelerate note assumed same as prev climb
        roskam.fuelfrac(6) = exp(-(p.parameters.alternate_range_km * 1000 * roskam.c_2 / 3600) / (p.parameters.cruise_mach * 295.07 * roskam.LoverD_cruise(j)));  %cruise 2
        roskam.fuelfrac(7) = exp(-(p.parameters.loiter_duration * 60 * roskam.c_3 / 3600) / (roskam.LoverD_loiter(j))); % loiter using endurance eqn
        roskam.fuelfrac(8) = 0.99; %descent 2
        roskam.fuelfrac(9) = 0.992;    %landing + taxi

        raymer.fuelfrac(1) = 0.97; %engine start + taxi + takeoff
        raymer.fuelfrac(2) = 0.985;    %climb
        raymer.fuelfrac(3) = exp(-(p.parameters.cruise_range_km * 1000 * raymer.c_1 / 3600) / (p.parameters.cruise_mach * 295.07 * raymer.LoverD_cruise(j)));  %cruise 1 using breguet range
        raymer.fuelfrac(4) = 0.99;  %descent 1
        raymer.fuelfrac(5) = 0.985;  %climb and accelerate note assumed same as prev climb
        raymer.fuelfrac(6) = exp(-(p.parameters.alternate_range_km * 1000 * raymer.c_2 / 3600) / (p.parameters.cruise_mach * 295.07 * raymer.LoverD_cruise(j)));  %cruise 2
        raymer.fuelfrac(7) = exp(-(p.parameters.loiter_duration * 60 * raymer.c_3 / 3600) / (raymer.LoverD_loiter(j))); % loiter using endurance eqn
        raymer.fuelfrac(8) = 0.99; %descent 2
        raymer.fuelfrac(9) = 0.995;    %landing + taxi

        for i = 1:length(roskam.fuelfrac)
            roskam.Wx_W0(j) = roskam.Wx_W0(j) * roskam.fuelfrac(i); %multiply fuel mass fractions 
        end

        roskam.Wf_W0(j) = 1.03 * (1 - roskam.Wx_W0(j));   %calculate total fuel mass fraction for mission
                                                          %2% fuel reserve
                                                          % + 1% unuseable

        for i = 1:length(raymer.fuelfrac)
            raymer.Wx_W0(j) = raymer.Wx_W0(j) * raymer.fuelfrac(i); %multiply fuel mass fractions 
        end

        raymer.Wf_W0(j) = 1.03 * (1 - raymer.Wx_W0(j));   %calculate total fuel mass fraction for mission
                                                          %2% fuel reserve
                                                          % + 1% unuseable
        %Roskam W0

        roskam.W0(j) = initialW0; %initial W0 guess (N)
        roskam.W0prev = 0; %initialise
        count = 0; 
        figure 
        hold on
        error = 10; 

        while abs(error) > 1e-8
            roskam.W0prev = roskam.W0(j); 
            We_W0_roskam_regress = (10^((log10(roskam.W0(j)) - roskam.A) / roskam.B)) / roskam.W0(j);   %Roskam regression method
            roskam.W0(j) = (W_crew + W_pld) / (1 - roskam.Wf_W0(j) - (We_W0_roskam_regress));    %W0 calculation 
            count = count + 1;  %count number of iterations
            error = roskam.W0(j) - roskam.W0prev;    %Calculate difference between consecutive W0 values
            plot(count, error, 'r*') %plot error 
            grid on
            pause(0.05)
        end  
        hold off

        roskam.We_W0(j) = We_W0_roskam_regress; 
        roskam.W0(j) = roskam.W0(j) * 0.453592 * 9.81; %ibs to kg to N;

        %Raymer W0

        raymer.W0(j) = initialW0; %initial W0 guess (N)
        raymer.W0prev = 0; %initialise
        raymer.Kvs = 1; %fixed sweep

        count = 0; 
        figure 
        hold on
        error = 10; 

        while abs(error) > 1e-8
            raymer.W0prev = raymer.W0(j); 
            We_W0_raymer_regress = raymer.A * (raymer.W0(j) ^ raymer.C) * raymer.Kvs;   %Raymer regression method
            raymer.W0(j) = (W_crew + W_pld) / (1 - raymer.Wf_W0(j) - (We_W0_raymer_regress));    %W0 calculation 
            count = count + 1;  %count number of iterations
            error = raymer.W0(j) - raymer.W0prev;    %Calculate difference between consecutive W0 values
            plot(count, error, 'b*') %plot error 
            grid on
            %pause(0.05)
        end  
        hold off

        raymer.We_W0(j) = We_W0_raymer_regress; 
        raymer.W0(j) = raymer.W0(j) * 0.453592 * 9.81; %ibs to kg to N
end

%Roskam_fit = fit(r.roskamdata.W0',r.roskamdata.We_W0','poly1');

figure
plot(roskam.W0(1), roskam.We_W0(1),'bo')
hold on
plot(roskam.W0(2), roskam.We_W0(2),'r*')
hold on
plot(raymer.W0(1), raymer.We_W0(1),'bo')
hold on
plot(raymer.W0(2), raymer.We_W0(2),'b*')
hold on
plot(r.roskamdata.W0, r.roskamdata.We_W0,'kx')
hold off
grid on
xlabel("$W_{0}$ [N]", 'interpreter', 'Latex','FontSize', 15)
ylabel("$\frac{W_{e}}{W_{0}}$", 'interpreter', 'Latex','FontSize', 15)
%title("Comparison to Roskam data",'interpreter', 'Latex','FontSize', 15)
legend("Roskam regression, Roskam L/D","Roskam regression, calculated L/D","Raymer regression, Roskam L/D","Raymer regression, calculated L/D","Roskam data for Business jets", 'interpreter', 'Latex')

%% Constraint diagram
%Fractions 

sizing.fraction.before_take_off=0.99*0.995;
sizing.fraction.before_cruise=sizing.fraction.before_take_off*roskam.fuelfrac(2);
sizing.fraction.before_alternate_cruise=sizing.fraction.before_cruise*roskam.fuelfrac(3)*roskam.fuelfrac(4)*roskam.fuelfrac(5);
sizing.fraction.before_loiter=sizing.fraction.before_alternate_cruise*roskam.fuelfrac(6);
sizing.fraction.end=roskam.Wx_W0(2);

%% save sizing to 
save('sizing','sizing');
