
airfoilcoords.top = [1.00000     0.00000;
0.95033     0.00986;
0.90066     0.01979;
0.85090     0.02974;
0.80103     0.03935;
0.75107     0.04847;
0.70101     0.05686;
0.65086     0.06440;
0.60064     0.07085;
0.55035     0.07602;
0.50000     0.07963;
0.44962     0.08139;
0.39923     0.08139;
0.34884     0.07971;
0.29846     0.07658;
0.24881     0.07193;
0.19781     0.06562;
0.14757     0.05741;
0.09746     0.04672;
0.07247     0.04010;
0.04757     0.03227;
0.02283     0.02234;
0.01059     0.01588;
0.00580     0.01236;
0.00347     0.01010;
0.00000     0.00000];

airfoilcoords.bottom = [ 0 0;
0.00653     -0.00810;
0.00920     -0.00956;
0.01441     -0.01160;
0.02717     -0.01490;
0.05243     -0.01963;
0.07753     -0.02314;
0.10254     -0.02604;
0.15243     -0.03049;
0.20219     -0.03378;
0.25189     -0.03613;
0.30154     -0.03770;
0.35116     -0.03851;
0.40077     -0.03855;
0.45038     -0.03759;
0.50000     -0.03551;
0.54965     -0.03222;
0.59936     -0.02801;
0.64914     -0.02320;
0.69899     -0.01798;
0.74893     -0.01267;
0.79897     -0.00751;
0.84910     -0.00282;
0.89934     0.00089;
0.94967     0.00278;
1.00000     0.00000];




x_frontspar = 0.1;
x_rearspar = 0.74;

k = x_rearspar - x_frontspar;


figure
hold on
plot(airfoilcoords.top(:,1),airfoilcoords.top(:,2),'-x')
plot(airfoilcoords.bottom(:,1),airfoilcoords.bottom(:,2),'-x')
axis equal

topfit = fit(airfoilcoords.top(:,1),airfoilcoords.top(:,2),'smoothingspline')
bottomfit = fit(airfoilcoords.bottom(:,1),airfoilcoords.bottom(:,2),['smoothingspline'])

plot(topfit)
plot(bottomfit)

x_range = linspace(x_frontspar,x_rearspar);
x_points = [x_range flip(x_range)];

y_points = [topfit(x_range)' bottomfit(flip(x_range))'];

fuel_tank = polyshape(x_points,y_points);
plot(fuel_tank)


t_c_frontspar = topfit(x_frontspar) - bottomfit(x_frontspar)
t_c_rearspar = topfit(x_rearspar) - bottomfit(x_rearspar)

plot([x_frontspar, x_frontspar],[topfit(x_frontspar),bottomfit(x_frontspar)],'linewidth',5);
plot([x_rearspar, x_rearspar],[topfit(x_rearspar),bottomfit(x_rearspar)],'linewidth',5);

span = 8.9729;

syms y
croot = 1.7674;
ctip = 0.5333;
c(y) = (((2*(ctip-croot))/(span))*y) + croot;


rib1 = span * 0/2;  
rib2 = span * 0.95/2;

area = fuel_tank.area;

raw_vol = double(2*int(c^2*area, rib1, rib2))

corrected_vol = raw_vol*0.975*0.975 %raymer -> fire and leak protection using foam


