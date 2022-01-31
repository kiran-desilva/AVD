clear
clc

cl_max = 2.27;
rho_0 = 1.225;
Sref = 10;
Mtow = 3.0688e4;

n_max = 2.5;
n_min = -1;

syms v cl

n =  0.5 * rho_0 * cl * (v^2) * Sref * (1/Mtow);

Va = solve(n(cl_max) == 1)

figure
hold on
fplot(subs(n,cl,cl_max),[0 300])

