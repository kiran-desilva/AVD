%% HEAVY FRAMES CALCULATOR
%WINGS FRONT AND REAR SPAR, VERTICAL TAIL FRON AND REAR SPAR, ENGINES
clear 
clc
close all

load("fuselageLoading.mat");


%% WINGS FRONT SPAR
% D_wing = 1.58;
% Loads_wingf = [fuselageLoading.Vd_flight.Ff, fuselageLoading.Vd_flight.Ff];
% angles_wingf = [50, 50]; %CHANGE THIS  
% Torque_wingf = []
% [HeavyFrames.wing.twFRONT, HeavyFrames.wing.lfwFRONT] = framedimensioncalc(Loads_wingf, Torque_wingf, angles_wingf, D_wing, "Wing Front Spar Heavy Frame");
% 
% %% WINGS REAR SPAR
% D_wing = 1.58;
% Loads_wingr = [fuselageLoading.Vd_flight.Fr, fuselageLoading.Vd_flight.Fr];
% angles_wingr = [50, 50]; %CHANGE THIS
% Torque_wingr = []
% [HeavyFrames.wing.twREAR, HeavyFrames.wing.lfwREAR] = framedimensioncalc(Loads_wingr, Torque_wingr, angles_wingr, D_wing, "Wing Rear Spar Heavy Frame");
% 
% %% VERTICAL TAIL FRONT SPAR
% D_VTF = ;
% Loads_VTF = ;
% angles_VTF = ;
% Torque_VTF = [];
% [HeavyFrames.tail.tvFRONT, HeavyFrames.tail.lfvFRONT] = framedimensioncalc(Loads_VTF, Torque_VTF, angles_VTF, D_VTF, "Vertical Tail Front Spar Heavy Frame");
% 
% %% VERTICAL TAIL REAR SPAR
% D_VTR = ;
% Loads_VTR = ;
% angles_VTR = ;
% Torque_VTr = [];
% [HeavyFrames.tail.tvREAR, HeavyFrames.tail.lfvREAR] = framedimensioncalc(Loads_VTR, Torque_VTr, angles_VTR, D_VTR, "Vertical Tail Rear Spar Heavy Frame");

%% ENGINES
D_E = 1.1426;
Loads_E = [132*9.81, 132*9.81];
angles_E = [90, 90];
Torque_E = [132*9.81*(1.1-0.5*D_E), 132*9.81*(2.2-0.5*D_E)]; %assumed weight acts at center of engine?
[HeavyFrames.engine.t, HeavyFrames.engine.lf] = framedimensioncalc(Loads_E, Torque_E, angles_E, D_E, "Engines Heavy Frame");

%%
function [t, lf] = framedimensioncalc(L, Torque, angles, D, LoadcaseSTR)
%LOADS: VERTICAL LOADS 
%ANGLES: DOWNWARDS = 0DEG
%TORQUE

    load("materialLib");
    yielddirect = materialLib{1}.tensile_yield;
    yieldshear = materialLib{1}.shear_yield;
    %yieldbending = ;SAME AS DIRECT?

    syms t lf Ad As %thicknesses equal
    H = 0.05; %wall-to-wall thickness                     
    Ixx = ((H-2*t)^3)*t/12 + 2*((t^3)*lf/12 + t*lf*((H-2*t)+t)^2 /4); %second moment f=of area of I beam
    %A = 2*tf*lf + lw*tw;    %frame sectionalarea
    yc = 0.025; %section thickness / 2 (defined by fuselage paramaters)
    R = D/2;
    
    for i = 1:length(L)
        P(i) = L(i)*sind(angles(i));
        Q(i) = L(i)*cosd(angles(i));
        T(i) = Torque(i);
        [N(:,i), M(:,i), S(:,i)] = SectionalLoadCalc(P(i), Q(i), T(i), D); 
    end
    
    Ntot = sum(N,2);  %sum together total stresses due to all forces on frame
    Mtot = sum(M,2);
    Stot = sum(S,2); 
    
    theta = linspace(0,360,360);
    figure
    plot(theta, (Ntot/sum(P)))
    hold on
    plot(theta, (Mtot/(R*sum(P))))
    hold on
    plot(theta, (Stot/sum(P)))
    hold off
    grid on
    xlabel("Theta (deg)", 'interpreter', 'Latex')
    ylabel("$\frac{N}{P}, \frac{M}{PR}, \frac{S}{P}$", 'interpreter', 'Latex')
    legend("$\frac{N}{P}$","$\frac{M}{PR}$","$\frac{S}{P}$", 'interpreter', 'Latex')
    title(LoadcaseSTR)
    
    figure
    plot(theta, (Ntot/sum(Q)))
    hold on
    plot(theta, (Mtot/(R*sum(Q))))
    hold on
    plot(theta, (Stot/sum(Q)))
    hold off
    grid on
    xlabel("Theta (deg)", 'interpreter', 'Latex')
    ylabel("$\frac{N}{Q}, \frac{M}{QR}, \frac{S}{Q}$", 'interpreter', 'Latex')
    legend("$\frac{N}{Q}$","$\frac{M}{QR}$","$\frac{S}{Q}$", 'interpreter', 'Latex')
    title(LoadcaseSTR)

    figure
    plot(theta, (Ntot*R/sum(T)))
    hold on
    plot(theta, (Mtot/(sum(T))))
    hold on
    plot(theta, (Stot*R/sum(T)))
    hold off
    grid on
    xlabel("Theta (deg)", 'interpreter', 'Latex')
    ylabel("$\frac{NR}{T}, \frac{M}{T}, \frac{SR}{T}$", 'interpreter', 'Latex')
    legend("$\frac{NR}{T}$","$\frac{M}{T}$","$\frac{SR}{T}$", 'interpreter', 'Latex')
    title(LoadcaseSTR)

    
    Nmax = max(Ntot);   %find max vals of stresses
    Mmax = max(Mtot);
    Smax = max(Stot);

    directstressmax = Nmax / Ad; %Pa
    shearstressmax = Smax / As; %Pa
    bendingstressmax = (Mmax*yc) / Ixx; %Pa

    f1 = -yielddirect + directstressmax == 0;
    f2 = -yieldshear + shearstressmax == 0;
    f3 = -yielddirect + bendingstressmax == 0;

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
    theta = linspace(0,360,360);

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