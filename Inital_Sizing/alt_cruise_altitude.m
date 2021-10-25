%min drag with g and rho
clear
clc
close all

alt=[0:5E+03:40E+03]; %ft
g=[32.174,32.159,32.143,32.128,32.112,32.097,32.082, 32.066,32.051]; %ft/s^2
rho=[0.0765,0.0659,0.0565,0.0481,0.0408,0.0343,0.0287,0.0237,0.0189]; %lbs/ft^3
hold on
figure
funct_1=1./sqrt(rho);
yyaxis left
plot(alt,funct_1, '-xk')
axis ([0 40E+03 0 8])
xlabel 'Altitude (ft)'
ylabel '1/sqrt(rho)'

climb_rate=2500 %ft/min

climb_time=alt./climb_rate %min
climb_range=climb_time*464.1*60 %ft (141.47 m/s)
range_total=1214E+03*ones(1,9) %ft (370 km)
cruise_range=range_total-climb_range %ft


cruise_range_km=cruise_range./(3.208*1000)
yyaxis right
plot(alt, cruise_range_km, '-om')
xlabel 'Altitude (ft)'
ylabel 'Cruising range (total-climb) (km)'
grid on
grid minor
axis ([0 40E+03 0 380])
hold off

