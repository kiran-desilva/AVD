function [phi_m] = sweep_angle(phi_n,m,n,ar,lambda)
    phi_m = arctand(tand(phi_n) - ( (4/ar)*((m-n)/100)*((1-lambda)/(1+lambda)) ));
end