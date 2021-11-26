syms t_c phi_25 MDD Cl MDD_eff phi_min M_crit M_cruise MDD_airfoil

Km = 1.05;

M_desired = 0.75;

ka = 0.87;

MDD_airfoil(Cl,t_c) = (Cl/10) + t_c - ka; 

MDD_eff = MDD * sqrt(cosd(phi_25));
t_c(MDD,Cl,phi_25) = 0.3*cosd(phi_25) * (( ((1 - ( (5+(MDD_eff^2)) / (5 + (Km-((0.25*Cl)^2)) ) )^(3.5)) * ((sqrt(1-(MDD_eff^2)))/(MDD_eff^2)) ))^(2/3));

simplify(t_c)


figure
hold on 
t_c_func(phi_25) = subs(t_c,[MDD,Cl],[M_desired,0])

fplot(t_c_func,[0 50])


fsolve(@(x) double(t_c(M_desired,0,x)) - 0.15,20)