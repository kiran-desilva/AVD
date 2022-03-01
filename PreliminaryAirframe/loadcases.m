clear
clc

load('vn.mat')

%%important wing constants at 40000 ft at cruise conditions

alpha_0 = -3*(pi/180); % radians
cl_0 = vn.cl_alpha_w*alpha_0;% this is unused 
cl_w = @(aoa) vn.cl_alpha_w.*(aoa+alpha_0); % cl of wing

cm_alpha = 0;
cm_0 = -0.0266;
cm = @(aoa) cm_0; %looking at thoery of wing sections, cm changes negligibly with aoa as naca 65412 is a cambered airfoil 


%%LOAD CASE 1 -> SYMMETRIC AT VA%%
idx = 1;

loadcase{idx}.n = vn.Va_nmax*1.5;
loadcase{idx}.lift = vn.Va_nmax * vn.Mtow; % the limit lift load
loadcase{idx}.v = vn.Va; %eas

loadcase{idx}.alt = vn.cruisealt_m;
[~,a,~,rho] = atmosisa(loadcase{idx}.alt);
loadcase{idx}.a = a;
loadcase{idx}.rho = rho;

loadcase{idx}.cl = (loadcase{idx}.lift)/(0.5*1.225*(loadcase{idx}.v^2 * vn.Sref));
loadcase{idx}.aoa = fsolve(@(x) cl_w(x) - loadcase{idx}.cl,0); %solve for required aoa

loadcase{idx}.cm = cm(loadcase{idx}.aoa);
loadcase{idx}.m = loadcase{idx}.cm * (0.5*1.225 * (loadcase{idx}.v^2) * vn.Sref * vn.cmac)

%%LOAD CASE 2 -> SYMMETRIC AT VD%%
idx = 2; 

loadcase{idx}.n = vn.Vd_nmax * 1.5;
loadcase{idx}.lift = vn.Vd_nmax * vn.Mtow; %same as above
loadcase{idx}.v = vn.Vd; %eas

loadcase{idx}.alt = vn.cruisealt_m;
[~,a,~,rho] = atmosisa(loadcase{idx}.alt);
loadcase{idx}.a = a;
loadcase{idx}.rho = rho;

loadcase{idx}.cl = (loadcase{idx}.lift)/(0.5*1.225*(loadcase{idx}.v^2) * vn.Sref);
loadcase{idx}.aoa = fsolve(@(x) cl_w(x) - loadcase{idx}.cl,0); %solve for required aoa

loadcase{idx}.cm = cm(loadcase{idx}.aoa);
loadcase{idx}.m = loadcase{idx}.cm * (0.5*1.225* (loadcase{idx}.v^2) * vn.Sref * vn.cmac)


save('loadcase.mat','loadcase');