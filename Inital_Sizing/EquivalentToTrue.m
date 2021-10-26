function [TrueAirspeed] = EquivalentToTrue(Equivalent,Altitude)

    [T_SL, A_SL, P_SL, RHO_SL] = atmosisa(0);

    [T_H, A_H, P_H, RHO_H] = atmosisa(distdim(Altitude,'ft','m'));
    gamma = 1.4;




    % Now we have VimD varying with wing loading as required for absolute
    % ceiling. We can use the specific excess power equation to get the
    % required point performance, however, that equation uses true airspeed in
    % the calculations. To convert to true airspeed, this equation uses the
    % compressible flow relations.

    A = ((((gamma-1).*(Equivalent.^2))./(2*A_SL*A_SL)) + 1).^(gamma/(gamma-1));

    TrueAirspeed = A_H .* ((2/(gamma-1).*((1 + ((P_SL/P_H).*(A-1)) ).^((gamma-1)/gamma) - 1))).^(0.5);
    

end