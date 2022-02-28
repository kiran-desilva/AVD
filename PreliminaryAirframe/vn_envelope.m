clear
clc

kts_to_ms = 0.514444;
fts_to_ms = 0.3048;

cl_max_pos = 1.362;
% cl_max_pos = 1.6;

% cl_2d_max_neg = -0.8;
% cl_max_neg = 0.9*cl_2d_max_neg*cosd(18.17);
cl_max_neg = -0.684;

%clean clalpha for wing at takeoff
cl_alpha_w_0 = 5.6007;
% cl_alpha_w_0 = 7;
%clean clalpha at cruise altitude
cl_alpha_w_cruise =  5.42;

cmac = 1.2607;

[T_0,a_0,p_0,rho_0] = atmosisa(0); %sea-level conditions
[T_c,a_c,p_c,rho_c] = atmosisa(distdim(40000,'ft','m')); %sea-level conditions

Sref = 10.32; %m^2 

Mtow = 3128*9.81; %n
wf_end_cruise = 0.8423;
ws_tow = Mtow/Sref;
ws_end_cruise = wf_end_cruise*ws_tow;

cruisealt_ft = 40000;
cruisealt_m = distdim(cruisealt_ft,'ft','m');

n_max = 2.5;
n_min = -1;

syms v cl 

%%STALL CONSTRAINT%%
n =  0.5 * rho_0 * cl * (v^2) * Sref * (1/Mtow); %stall speed constraint

%%%%%%%%%%%%%%%%%%

pos_solution = @(array) array(array>=0);

Vs1 = pos_solution(double(solve(subs(n,cl,cl_max_pos) == 1,v)));

Vc = m_to_eas(0.75,cruisealt_m);

%equivalent to Vs1*sqrt(n_max)
Va = pos_solution(double(solve(subs(n,cl,cl_max_pos) == n_max,v)));

Vf = pos_solution(double(solve(subs(n,cl,cl_max_neg)  == n_min,v)));

Vd = Vc/0.8;
Md = 0.75/0.8;
% Vd_2 = m_to_eas(Md,cruisealt_m); gives the same result

%%GUST FACTORS%%
syms v cl cl_alpha ws rho Ude
mu = (2*ws)/(rho*9.81*cmac*cl_alpha);
%subsonic allieviation factor
k_sub = (0.88*mu)/(5.3+mu);
delta_n(v,cl_alpha,ws,rho,Ude) = simplify((rho_0*k_sub*Ude*v*cl_alpha)/(2*ws));

%taken from FAR-25 Transport Category Aircraft at 40kft
% Ude_Vd = 16.3*fts_to_ms; 
Ude_Vd = 25*fts_to_ms; 
% Ude_Vc = 33*fts_to_ms;
Ude_Vc = 50*fts_to_ms;
% Ude_Vb = 47*fts_to_ms; 
Ude_Vb = 66*fts_to_ms; 
% Ude_Vb = 20;

% using uref, density and cl_alpha at sea level , idk if this is right ngl 
Vb = Vs1*sqrt(1+double(delta_n(Vc,cl_alpha_w_0,ws_tow,rho_0,uref_func(0))));
Vb = Vc - (1.32*uref_func(cruisealt_m));


% Vb = Vs1*sqrt(1+double(delta_n(Vc,cl_alpha_w_0,ws_tow,rho_0,Ude_Vb)));
% Vb = Vc - (1.32*uref_func(0));




% Vb = Vs1*sqrt(1+double(delta_n(Vc,cl_alpha_w_0,ws_tow,rho_0,uref_func(0))));
% Vb = Vs1*sqrt(1+double(delta_n(Vc,cl_alpha_w_0,ws_tow,rho_0,Ude_Vb)));
% Vb_sym = solve(n == 1 + delta_n(v,cl_alpha_w_0,ws_tow,rho_0,Ude_Vb_0),v);
% Vb = pos_solution(double(subs(Vb_sym,cl,cl_max_pos)));

% Vb_delta_n = double(delta_n(Vb,cl_alpha_w_cruise,ws_end_cruise,rho_c,Ude_Vb));
% Vb_delta_n  = double(delta_n(Vb,cl_alpha_w_cruise,ws_tow,rho_c,Ude_Vb));
% Vc_delta_n = double(delta_n(Vc,cl_alpha_w_cruise,ws_tow,rho_c,Ude_Vc));
% Vd_delta_n = double(delta_n(Vd,cl_alpha_w_cruise,ws_tow,rho_c,Ude_Vd));


Vb_delta_n  = double(delta_n(Vb,cl_alpha_w_0,ws_tow,rho_0,Ude_Vb));
Vc_delta_n = double(delta_n(Vc,cl_alpha_w_0,ws_tow,rho_0,Ude_Vc));
Vd_delta_n = double(delta_n(Vd,cl_alpha_w_0,ws_tow,rho_0,Ude_Vd));
%%%%%%%%%%%%%%%%%% %%


figure
hold on
grid on
grid minor

%stall constraints
SC_plot = fplot(subs(n,cl,cl_max_pos),'-.','color','black');
fplot(subs(n,cl,cl_max_neg),'-.','color','black')

%%Important Speeds
xline(Va,'--','Va','LabelOrientation','horizontal','linewidth',1,'fontsize',12)
xline(Vc,'--','Vc','LabelOrientation','horizontal','linewidth',1,'fontsize',12)
xline(Vs1,'--','Vs1','LabelOrientation','horizontal','linewidth',1,'fontsize',12)
xline(Vd,'--','Vd','LabelOrientation','horizontal','linewidth',1,'fontsize',12)
xline(Vb,'--','Vb','LabelOrientation','horizontal','linewidth',1,'fontsize',12)

%%Manouver Envelope - UltimateLoad
n_ul_max = 1.5*2.5;
n_ul_min = 1.5*-1;

% Va_ul = Vs1*sqrt(n_ul_max);
Va_ul = Va;
% Vf_ul = Vf*sqrt(abs(n_ul_min)); 
Vf_ul = Vf;

ULM_color = 'red';
ULM_linewidth = 2;
ULM_plot = fplot(subs(n,cl,cl_max_pos),[0 Va_ul],'color',ULM_color,'linewidth',ULM_linewidth);
fplot(subs(n,cl,cl_max_neg),[0 Vf_ul],'color',ULM_color,'linewidth',ULM_linewidth)

plot([Va_ul Vd],[n_ul_max n_ul_max],'color',ULM_color,'linewidth',ULM_linewidth)
plot([Vf_ul Vc],[n_ul_min n_ul_min],'color',ULM_color,'linewidth',ULM_linewidth)
plot([Vc Vd],[n_ul_min 0],'color',ULM_color,'linewidth',ULM_linewidth)% CS 25.337 C2

plot([Vd Vd],[n_ul_max 0],'color',ULM_color,'linewidth',ULM_linewidth)


%%Manouver Envelope - LimitLoad
LLM_color = 'blue';
LLM_linewidth = 2;
LLM_plot = fplot(subs(n,cl,cl_max_pos),[0 Va],'color',LLM_color,'linewidth',LLM_linewidth);
fplot(subs(n,cl,cl_max_neg),[0 Vf],'color',LLM_color,'linewidth',LLM_linewidth)

plot([Va Vd],[n_max n_max],'color',LLM_color,'linewidth',LLM_linewidth)
plot([Vf Vc],[n_min n_min],'color',LLM_color,'linewidth',LLM_linewidth)
plot([Vc Vd],[n_min 0],'color',LLM_color,'linewidth',LLM_linewidth)% CS 25.337 C2

plot([Vd Vd],[n_max 0],'color',LLM_color,'linewidth',LLM_linewidth)


%%Gust Envelope
GE_color = 'green';
GE_linewidth = 2;

plot([0 Vb],[1 1+Vb_delta_n],'--','color','green')
plot([0 Vb],[1 1-Vb_delta_n],'--','color','green')

plot([0 Vc],[1 1+Vc_delta_n],'--','color','green')
plot([0 Vc],[1 1-Vc_delta_n],'--','color','green')

plot([0 Vd],[1 1+Vd_delta_n],'--','color','green')
plot([0 Vd],[1 1-Vd_delta_n],'--','color','green')

GE_plot = plot([0 Vb Vc Vd Vd Vc Vb 0], [1 1+Vb_delta_n 1+Vc_delta_n 1+Vd_delta_n 1-Vd_delta_n 1-Vc_delta_n 1-Vb_delta_n 1],'color',GE_color,'linewidth',GE_linewidth);



xlabel("EAS (m/s)")
ylabel("Load Factor")

legend([ULM_plot LLM_plot GE_plot SC_plot],'Ultimate Load Manouever Envelope','Limit Load Manouever Envelope','Gust Envelope','Stall Load Factor Constraint','location','southwest')

%xlim([0 150])
ylim([-3 4])

improvePlotNOLINE(gcf);

saveas(gcf,"Figures/Vndiagram",'epsc')

vn.Va = Va;
vn.Vd = Vd;
vn.cl_alpha_w = cl_alpha_w_cruise;
vn.cl_max = cl_max_pos;
vn.Mtow = Mtow;
vn.cmac = cmac;
vn.Sref = Sref;
vn.ws_tow = ws_tow;
vn.cruisealt_m = cruisealt_m;
vn.n_max = n_max;

save('vn.mat','vn')