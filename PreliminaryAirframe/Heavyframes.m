%% HEAVY FRAMES CALCULATOR
%WINGS FRONT AND REAR SPAR, VERTICAL TAIL FRON AND REAR SPAR, ENGINES
clear 
clc

%% WINGS FRONT SPAR
D_wing = 1.58;
Loads_wingf = ;
angles_wingf = ;
[HeavyFrames.wing.twFRONT, HeavyFrames.wing.lfwFRONT] = framedimensioncalc(Loads_wingf, angles_wingf, D_wing);

%% WINGS REAR SPAR
D_wing = 1.58;
Loads_wingr = ;
angles_wingr = ;
[HeavyFrames.wing.twREAR, HeavyFrames.wing.lfwREAR] = framedimensioncalc(Loads_wingr, angles_wingr, D_wing);

%% VERTICAL TAIL FRONT SPAR
D_VTF = ;
Loads_VTF = ;
angles_VTF = ;
[HeavyFrames.tail.tvFRONT, HeavyFrames.tail.lfvFRONT] = framedimensioncalc(Loads_VTF, angles_VTF, D_VTF);

%% VERTICAL TAIL REAR SPAR
D_VTR = ;
Loads_VTR = ;
angles_VTR = ;
[HeavyFrames.tail.tvREAR, HeavyFrames.tail.lfvREAR] = framedimensioncalc(Loads_VTF, angles_VTF, D_VTF);

%% 
function [t, lf] = framedimensioncalc(Loads, angles, D)

    load("materialLib");
    yielddirect = materialLib{1}.tensile_yield;
    yieldshear = materialLib{1}.shear_yield;
    yieldbending = 2000*10^6;%materialLib{1}.bending_yield; 

    syms t lf Ad As %thicknesses equal
    H = 0.05; %wall thickness
    Ixx = ((H-2*t)^3)*t/12 + 2*((t^3)*lf/12 + t*lf*((H-2*t)+t)^2 /4); %second moment f=of area of I beam
    %A = 2*tf*lf + lw*tw;    %frame sectionalarea
    yc = 0.025; %section thickness / 2 (defined by fuselage paramaters)
    
    for i = 1:length(Loads)
        P = L(i)*sind(angles(i));
        Q = L(i)*cosd(angles(i));
        T = ;
        [N(i), M(i), S(i)] = SectionalLoadCalc(P, Q, T, D); 
    end
    
    Ntot = sum(N);  %sum together total stresses due to all forces on frame
    Mtot = sum(M);
    Stot = sum(S);     
    
    Nmax = max(Ntot);   %find max vals of stresses
    Mmax = max(Mtot);
    Smax = max(Stot);

    directstressmax = Nmax / Ad; %Pa
    shearstressmax = Smax / As; %Pa
    bendingstressmax = (Mmax*yc) / Ixx; %Pa

    f1 = -yielddirect + directstressmax == 0;
    f2 = -yieldshear + shearstressmax == 0;
    f3 = -yieldbending + bendingstressmax == 0;

    Ad = double(solve(f1,Ad));
    As = double(solve(f2,As));
    areas = [Ad,As];
    A = max(areas);
    f4 = A == 2*t*lf + (H-2*t)*t;    %frame sectionalarea

    eqns = [f3 f4];
    vars = [t; lf];
    Sol = solve(eqns, vars);
    t = double(Sol.t((imag(Sol.t)==0)));
    lf = double(Sol.lf((imag(Sol.t)==0)));
end

function [N, M, S] = SectionalLoadCalc(P, Q, T, D)

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
end