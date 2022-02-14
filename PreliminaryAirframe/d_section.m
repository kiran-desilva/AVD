%% D section of wing
% not the most efficient script in the world but seems the best way to
% input the K_s value manually in a short period of time. Not sure how
% efficient it is to digitise the plots.
clear
clc

load wing_ribs.mat %in PreliminaryAirframe file
%% Wing box parameters
b_ref=3.83; %wingbox span [m] -NEED TO VERIFY
a=s; %d-cell span length
b=0.12.*chord; %curved panel length - proportion of x/c. [NEED TO ESTABLISH THE PROPORTION]
t1=2e-03; %d-cell thickness (initial guess - 2mm)
R1= 1; %D-cell radius

frac_A1=b./a; %b/a if b>a
frac_B1=a./sqrt(R1*t1); %if b>a
frac_A2=a./b; %a/b if a>b
frac_B2=b./sqrt(R1*t1); %if a>b



for i=1:length(a)
    mess_1=['Iteration ', num2str(i)];
    disp(mess_1)
    mess_3=['Value of a ', num2str(a(i))];
    disp(mess_3)
    mess_4=['Value of b ', num2str(b(i))];
    disp(mess_4)
    if a(i)>b(i)
        disp('a is greater than b')
        mess_2=['a/b = ', num2str(frac_A2(i))];
        mess_5=['b/sqrt(t*R)= ', num2str(frac_B2(i))];
    else
        disp('b is greater than a')
        mess_2=['b/a = ', num2str(frac_A1(i))];
        mess_5=['a/sqrt(t*R)= ', num2str(frac_B1(i))];
    end
    disp (mess_2)
    disp (mess_5)
    prompt='What is the K_s value?';
    k_s(i)=input(prompt);
    clc
end
E_panel=72e+09; %Young's modulus of Panel material
tau_cr=k_s.*E_panel.*(t1./b).^2; %buckling stress of curved panel in shear

save ('d_section.mat', 'tau_cr','k_s','b','t1','E_panel')


