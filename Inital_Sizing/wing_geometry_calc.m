function [Croot,Ctip,Cmac,b] = wing_geometry_calc(Sref,Ar,lambda,sweep_25,n)
    b = sqrt(Ar*Sref);

    Croot = 2 * n * Sref / (n*b*(1+lambda));
    Ctip = lambda * Croot;
    Cmac = (2/3) * Croot * (1+lambda+(lambda^2)) / (1+lambda); 
    Ybar = (n*b*6)*((1+(2*lambda))/(1+lambda));
end