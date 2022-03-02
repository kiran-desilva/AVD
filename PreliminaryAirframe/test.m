clear
clc

syms H t lf M s


ixx = ((H-2*t)^3)*t/12 + 2*((t^3)*lf/12 + t*lf*((H-2*t)+t)^2 /4);
A = 2*t*lf + (H-2*t)*t;

ixx_req_constraint = simplify(((M*H)/s) - ixx)
solve(ixx_req_constraint,H)