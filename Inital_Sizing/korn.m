clear
clc

syms sweep Cl_design Mdd Cl_cruise t_c


sweep_historical = 22;
design_cl = 0.33759

% t_c = 0.12;
ka = 0.87;
Mcruise = 0.75;

Mdd = Mcruise*sqrt(cosd(sweep))

Cl_airfoil = 10*(ka - t_c - Mdd)

Cl_cruise(sweep,t_c) = Cl_airfoil*(0.9*0.95*cosd(sweep))

simplify(Cl_cruise)

figure
hold on
fplot(@(x) Cl_cruise(x,0.08),[0,100])
fplot(@(x) Cl_cruise(x,0.09),[0,100])
fplot(@(x) Cl_cruise(x,0.10),[0,100])
fplot(@(x) Cl_cruise(x,0.12),[0,100])
fplot(@(x) Cl_cruise(x,0.15),[0,100])

% fcn = subs(Cl_cruise,t_c,0.12);
options = optimset('Display','off')
required_sweep_08 = fsolve(@(x) double(Cl_cruise(x,0.08)) - design_cl, 20,options)
required_sweep_09 = fsolve(@(x) double(Cl_cruise(x,0.09)) - design_cl, 20,options)
required_sweep_10 = fsolve(@(x) double(Cl_cruise(x,0.10)) - design_cl, 20,options)
required_sweep_12 = fsolve(@(x) double(Cl_cruise(x,0.12)) - design_cl, 20,options)
required_sweep_15 = fsolve(@(x) double(Cl_cruise(x,0.15)) - design_cl, 20,options)

Mdd_08 = double(subs(Mdd,sweep,required_sweep_08))
Mdd_09 = double(subs(Mdd,sweep,required_sweep_09))
Mdd_10 = double(subs(Mdd,sweep,required_sweep_10))
Mdd_12 = double(subs(Mdd,sweep,required_sweep_12))
Mdd_15 = double(subs(Mdd,sweep,required_sweep_15))

Cl_section_design_08 = design_cl/(0.9*0.95*cosd(required_sweep_08))
Cl_section_design_09 = design_cl/(0.9*0.95*cosd(required_sweep_09))
Cl_section_design_10 = design_cl/(0.9*0.95*cosd(required_sweep_10))
Cl_section_design_12 = design_cl/(0.9*0.95*cosd(required_sweep_12))
Cl_section_design_15 = design_cl/(0.9*0.95*cosd(required_sweep_15))