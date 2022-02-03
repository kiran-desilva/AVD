function [Uref] = uref_func(alt)
    syms x
    uref_fit(x) = piecewise(x < 4572,poly2sym(polyfit([0 4572],[17.07 13.41],1)),x >= 4572,poly2sym(polyfit([4572 18288],[13.41 6.26],1)));
    Uref = double(uref_fit(alt));