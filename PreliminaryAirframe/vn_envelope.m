clear
clc

cl_max_pos = 1.362;

% cl_2d_max_neg = -0.8;
% cl_max_neg = 0.9*cl_2d_max_neg*cosd(18.17);
cl_max_neg = -0.684;

[T_0,a_0,p_0,rho_0] = atmosisa(0); %sea-level conditions
Sref = 10.32; %m^2 
Mtow = 3128*9.81; %n
cruisealt_ft = 40000;
cruisealt_m = distdim(cruisealt_ft,'ft','m')

n_max = 2.5;
n_min = -1;

syms v cl

n =  0.5 * rho_0 * cl * (v^2) * Sref * (1/Mtow); %stall speed constraint

pos_lift = subs(n,cl,cl_max_pos);
% pos_lift_flapped = subs()
neg_lift = subs(n,cl,cl_max_neg);


pos_solution = @(array) array(array>=0);

Vs1 = pos_solution(double(solve(subs(n,cl,cl_max_pos) == 1,v)));

Vc = m_to_eas(0.75,cruisealt_m);

%equivalent to Vs1*sqrt(n_max)
Va = pos_solution(double(solve(subs(n,cl,cl_max_pos) == n_max,v)));

Vf = pos_solution(double(solve(subs(n,cl,cl_max_neg)  == n_min,v)));

Vd = Vc/0.8;
Md = 0.75/0.8;
% Vd_2 = m_to_eas(Md,cruisealt_m); gives the same result




% Vne = 

figure
hold on
grid on
grid minor

fplot(subs(n,cl,cl_max_pos),[0 Va],'color','blue')
fplot(subs(n,cl,cl_max_neg),[0 Vf],'color','blue')

plot([Va Vd],[n_max n_max],'color','blue')
plot([Vf Vd],[n_min n_min],'color','blue')
plot([Vd Vd],[n_max n_min],'color','blue')


xline(Va,'--')
xline(Vc,'--')
xline(Vs1,'--')
xline(Vd,'--')

xlabel("EAS (m/s)")
ylabel("Load Factor")

xlim([0 150])
ylim([-3 4])

