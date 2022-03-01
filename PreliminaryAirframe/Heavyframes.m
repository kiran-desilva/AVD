%% HEAVY FRAMES CALCULATOR
%WINGS FRONT AND REAR SPAR, VERTICAL TAIL FRON AND REAR SPAR, ENGINES
clear 
clc
close all

load("fuselageLoading.mat");
load materialLib;


%% WINGS FRONT SPAR
D_wing = 1.58;
Loads_wingf = [fuselageLoading.Vd_flight.Ff, fuselageLoading.Vd_flight.Ff];
angles_wingf = [50, 50]; %CHANGE THIS  
Torque_wingf = [0,0];
[HeavyFrames.wing.twFRONT, HeavyFrames.wing.lfwFRONT] = framedimensioncalc(Loads_wingf, Torque_wingf, angles_wingf, D_wing, "Wing Front Spar Heavy Frame");

%% WINGS REAR SPAR
D_wing = 1.58;
Loads_wingr = [fuselageLoading.Vd_flight.Fr, fuselageLoading.Vd_flight.Fr];
angles_wingr = [50, 50]; %CHANGE THIS
Torque_wingr = [0,0];
[HeavyFrames.wing.twREAR, HeavyFrames.wing.lfwREAR] = framedimensioncalc(Loads_wingr, Torque_wingr, angles_wingr, D_wing, "Wing Rear Spar Heavy Frame");

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
    materialLib = struct();
    load materialLib;
    
    yielddirect = materialLib{1}.tensile_yield;
    yieldshear = materialLib{1}.shear_yield;

       
    
    ixx =@(t,H,lf) ((H-2*t)^3)*t/12 + 2*((t^3)*lf/12 + t*lf*((H-2*t)+t)^2 /4); %second moment f=of area of I beam
    A =@(t,H,lf) 2*t*lf + (H-2*t)*t;    %frame sectionalarea

    % yc = 0.025; %section thickness / 2 (defined by fuselage paramaters)
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
    plot(theta, (Stot)*R/sum(T))
    hold off
    grid on
    xlabel("Theta (deg)", 'interpreter', 'Latex')
    ylabel("$\frac{NR}{T}, \frac{M}{T}, \frac{SR}{T}$", 'interpreter', 'Latex')
    legend("$\frac{NR}{T}$","$\frac{M}{T}$","$\frac{SR}{T}$", 'interpreter', 'Latex')
    title(LoadcaseSTR)
    
    figure
    plot(theta, Ntot)
    hold on
    plot(theta, Mtot/(sum(T)))
    hold on
    plot(theta, Stot)
    hold off
    grid on
    xlabel("Theta (deg)", 'interpreter', 'Latex')
    ylabel("$N, M, S$", 'interpreter', 'Latex')
    legend("$N$","$M$","$S$", 'interpreter', 'Latex')
    title(LoadcaseSTR)

    
    Nmax = max(abs(Ntot));   %find max vals of stresses
    Mmax = max(abs(Mtot));
    Smax = max(abs(Stot));

    Ad_req = Nmax/yielddirect;
    As_req = Smax/yieldshear;

    A_req = max([Ad_req,As_req]);
    ixx_req = @(H) (Mmax*0.5*H)/yielddirect;

    constraints.ixx_req = ixx_req;
    constraints.A_req = A_req;
    constraints.ixx = ixx;
    constraints.A = A;

    function [c,ceq] = cons(x,constraints)
        c = [constraints.ixx_req(x(2))-constraints.ixx(x(1),x(2),x(3));
        constraints.A_req - constraints.A(x(1),x(2),x(3));
        (2*x(1))-x(2)]; % constrain height to be at elast 2*thickness
        ceq = [];
    end

    H_max = 0.05; %wall-to-wall thickness 
    x0 = [0.5e-2,0.9*H_max,1e-2];

    lb = [1e-3,1e-3,1e-3];
    ub = [1e-1,H_max,0.5]; %lol a meter

    opts = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'UseParallel',true)
    [res,f,flag] = fmincon(@(x) A(x(1),x(2),x(3)),x0,[],[],[],[],lb,ub,@(x) cons(x,constraints),opts)
    res 
    f

    t = res(1)
    H = res(2)
    lf = res(3)
    sol_ixx_req = ixx_req(H)
    sol_ixx = ixx(t,H,lf)
    sol_A_req = A_req
    sol_A = A(t,H,lf)


    % %x = [t H lf]


    % % directstressmax = Nmax / Ad; %Pa
    % % shearstressmax = Smax / As; %Pa
    % bendingstressmax = (Mmax*yc) / Ixx; %Pa

    % f1 = -yielddirect + directstressmax == 0;
    % f2 = -yieldshear + shearstressmax == 0;
    % f3 = -yielddirect + bendingstressmax == 0;

    % Ad = double(solve(f1,Ad));
    % As = double(solve(f2,As));
    % areas = [Ad,As];
    % A = max(areas);
    % f4 = A == 2*t*lf + (H-2*t)*t;    %frame sectionalarea
    
%     area_eq = @(t,lf) 2*t*lf + (H-2*t)*t;
%     
%     ixx_eq = matlabFunction(Ixx) % function handle @(lf,t)
%     
%     
%     
%     fmincon(
    
    
    
    % eqns = [f3 f4];
    % vars = [t; lf];
    % Sol = solve(eqns, vars);
    % t = double(Sol.t((imag(Sol.t)==0)));
    % lf = double(Sol.lf((imag(Sol.t)==0)));
    

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