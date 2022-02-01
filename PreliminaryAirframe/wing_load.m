clear
clc
g = 9.81;

%% Plot fuselage
fus.r = 1.68/2;
th = 0:pi/50:2*pi;
xunit = fus.r * cos(th);
yunit = fus.r * sin(th);
color = 'k';

plot(xunit, yunit, color);
axis equal;
grid on;


%% Underbelly
ubelly.x_at_r = 0.51;
ubelly.y_at_r = -sqrt(fus.r^2 - ubelly.x_at_r^2);

ubelly.x_at_c = 0;
ubelly.y_at_c = -0.14-fus.r;

syms A B x
f(x) = A*exp(15*x) + B;

[solA, solB] = solve(f(ubelly.x_at_r) == ubelly.y_at_r, f(ubelly.x_at_c) == ubelly.y_at_c);
f(x) = subs(f, [A B], [solA, solB]);
x_range = linspace(ubelly.x_at_c, ubelly.x_at_r, 100);
hold on;
plot(x_range, f(x_range), color);
plot(-x_range, f(x_range), color);
hold off

c_root = 1.767;
c_tip = 0.533;
span = 8.94;
semispan = span/2;
t_c_ratio = 0.12;

c(x) = -(c_root - c_tip)/semispan*x + c_root;
h(x) = t_c_ratio*c(x);

x_heights = linspace(0, semispan, 100);
top_edge(x) = h(x)/2 - fus.r;
bottom_edge(x) = -h(x)/2 - fus.r;

top_edge_intercept = vpasolve(top_edge == f, x, [0 semispan]);
bottom_edge_intercept = vpasolve(bottom_edge == f, x, [0 semispan]);

top_linspace = linspace(top_edge_intercept, semispan, 100);
bottom_linspace = linspace(bottom_edge_intercept, semispan, 100);

hold on;
plot(top_linspace, top_edge(top_linspace), color);
plot(bottom_linspace, bottom_edge(bottom_linspace), color);
plot([semispan semispan], [bottom_edge(semispan) top_edge(semispan)], color)
hold off

%% Loading scale factor
loading_scale_factor = 1e+03;

%% Undercarriage loading
uc.spanwise_start = 1.18 - 0.049;
uc.spanwise_end = 1.18 + 0.049 + 0.2;
uc.wing_loading = 80*0.454*g/abs(uc.spanwise_start - uc.spanwise_end);
uc.loading(x) = piecewise((x >= uc.spanwise_start) & (x <= uc.spanwise_end), -uc.wing_loading, 0)/loading_scale_factor;

%% Fuel loading
fuel.volume_in_wings_m3 = 0.5341;
fuel.density = 775; % kg/m^3
fuel.fuel_weight = g*fuel.density*fuel.volume_in_wings_m3;
fuel.loading(x) = piecewise((x >= bottom_edge_intercept) & (x <= 0.95*semispan), -fuel.fuel_weight, 0)/loading_scale_factor;

%% Lift loading
wing.L0 = 8000; %% TODO: GET FROM CHORLEY
wing.loading(x) = wing.L0*sqrt(1 - (x/semispan)^2)/loading_scale_factor;
hold on;
plot(bottom_linspace, fuel.loading(bottom_linspace) + wing.loading(bottom_linspace) + uc.loading(bottom_linspace) - fus.r);
hold off;
