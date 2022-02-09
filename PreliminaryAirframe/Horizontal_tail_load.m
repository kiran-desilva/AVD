clear
clc

%loading on horizontal tail as a factor of unit lift force

s_h = 2.3887; %m
U = 100;
AR_h = 3.9;
rho = 1.225; 
Sref_h = 1.9283; 
Sref_w = 10.32;

X_cg = 14.94 * 0.3048; %ft to m ref point nose
X_acHtail = 10.3474; 
X_acWing = 4.3620; 
WH = 17.5 * 0.453592 * 9.81; %horizontal tail weight
ULF = 3.75; % ultimate load factor
W0 = 3128.2 * 9.81; %takeoff weight
Lwing = ULF * W0;
Vdive = 137.25;
Cm0 = -0.0266;
M0w = 0.5 * rho * Vdive^2 * Sref_w * Cm0;

L_Htail = (Lwing * (X_cg - X_acWing) + M0w) / (X_acHtail - X_cg); 

y = linspace(-s_h/2,s_h/2,100);
gamma0 = (8 * s_h * 1) / (pi * AR_h * rho * U * Sref_h); %L = 1
gamma = @(y) gamma0 * sqrt(1 - (y ./ (s_h/2)).^2); 
TailLoadinit = rho * U .* gamma(y); 
TailLoadintegral = trapz(y,TailLoadinit); 
TailLoad = @(y) L_Htail * rho * U * gamma(y) / TailLoadintegral; %unit overall tail lift 

figure
plot(y, TailLoad(y))
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Load (N)", 'interpreter', 'Latex')
grid on

R = WH - L_Htail; %vertical reaction from vertical stab

[shearforce] = ShearBending(y,s_h,R,TailLoad,WH);
bendingmoment = cumtrapz(y,shearforce);

figure
plot(y, shearforce)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Shear Force (N)", 'interpreter', 'Latex')
grid on

figure
plot(y, bendingmoment)
xlabel("y (m)", 'interpreter', 'Latex')
ylabel("Horizontal Tail Bending Moments (Nm)", 'interpreter', 'Latex')
grid on


function [ShearforceH] = ShearBending(y,s_h,R,TailLoad,WH)
    for i = 1:length(y)
        Y = linspace(-s_h/2, y(i), 1000);
        %Moments = TailMoments(Y, (y(i) - y), TailLoad);
        if y(i) > 0 
             ShearforceH(i) = -(R - WH +  trapz(Y, TailLoad(Y))); 
             %BendingmomentH(i) = -((-R + WH)*y(i) +  trapz(Y, Moments(Y)));
        else
            ShearforceH(i) = -(trapz(Y, TailLoad(Y))); 
            %BendingmomentH(i) = -(trapz(Y, Moments(Y)));
        end
    end
     
end

%function [TailMoments] = TailMoments(y, momentarms, TailLoad)
 %   TailMoments = TailLoad(y) .* momentarms;
%end