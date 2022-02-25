clear 
clc

%fmincon
%see lightframes for area

% tf = %flange thickness          //////////
% tw = %web thickness                 //          ^
% lw = %web height                    //          | lw
% lf = %flange length                 //          |
D = 1.58; %fuselage diameter  m   //////////

load("materialLib");
yielddirect = materialLib{1}.tensile_yield;
yieldshear = materialLib{1}.shear_yield;
yieldbending = 2000*10^6;%materialLib{1}.bending_yield; 

syms tf tw lf 
lw = 0.1;%not quite
Ixx = (lw^3)*tw/12 + 2*((tf^3)*lf/12 + tf*lf*(lw+tf)^2 /4); %second moment f=of area of I beam
A = 2*tf*lf + lw*tw;    %frame sectionalarea
yc = 0.05; %section thickness / 2 (defined by fuselage paramaters)

P = 10000;
Q = 10000;
T = 100000;

[N, M, S] = SectionalLoadCalc(P, Q, T, D);     

directstressmax = N / (A); %Pa
shearstressmax = S / (A); %Pa
bendingstressmax = (M*yc) / (Ixx); %Pa

f1 = -yielddirect + directstressmax == 0;

f2 = -yieldshear + shearstressmax == 0;
f3 = -yieldbending + bendingstressmax == 0;


eqns = [f1 f2 f3];
vars = [tf; tw; lf];
Sol = solve(eqns, vars)

function [Nmax, Mmax, Smax] = SectionalLoadCalc(P, Q, T, D)

    R = D/2;
    theta = linspace(0,360,13);

    %tangential load case 

    NT = (P/(2*pi)) .* ((sind(theta)./2) - (pi - theta) .* cosd(theta));
    MT = (P * R / (2 * pi)) .* (1.5 .* sind(theta) + (pi - theta).*(cosd(theta) - 1));
    ST = (P/(2*pi)) * ((pi - theta) .*sind(theta) - 1 - (cosd(theta)./2));

    %radial load case

    NR = (Q/(2*pi)) .* (1.5.*cosd(theta) + (pi - theta).*sind(theta));
    MR = (Q*R/(2*pi)) .* (0.5.*cosd(theta) - (pi - theta).*sind(theta) + 1);
    SR = (Q/(2*pi)).* ((pi - theta).*cosd(theta) - 0.5.*sind(theta));

    %moment case

    NM = (T/(2*pi*R)).*(1.5.*cosd(theta) + (pi - theta).*sind(theta));
    MM = (T/(2*pi)).*(pi - 2.*sind(theta) - theta);
    SM = (T/(2*pi*R)).*(1+2.*cosd(theta));

    N = NT + NR + NM;
    M = MT + MR + MM;
    S = ST + SR + SM;
 
    Nmax = max(N);
    Mmax = max(M);
    Smax = max(S);
end