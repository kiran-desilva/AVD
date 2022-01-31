clear
clc

cl_max_pos = 2.27;
cl_max_neg = -2;

rho_0 = 1.225; 
Sref = 10.32; %m^2 
Mtow = 3128*9.81; %n

n_max = 2.5;
n_min = -1;

syms v cl

n =  0.5 * rho_0 * cl * (v^2) * Sref * (1/Mtow);

pos_solution = @(array) array(array>=0)

Va = pos_solution(double(solve(subs(n,cl,cl_max_pos) == n_max,v)))

Vf = pos_solution(double(solve(subs(n,cl,cl_max_neg)  == n_min,v)))

% Vne = 

figure
hold on
fplot(subs(n,cl,cl_max_pos),[0 Va])
fplot(subs(n,cl,cl_max_neg),[0 Vf])
% yline

