span = 9.8062;

syms y

k = ;
c(y) = -1.4044*y + 1.9556;
a = ;
b = ;

rib1 = span * 0.35;
rib2 = span * 0.95;

vol = int(k*c^2*(a+b)/2, rib1, rib2)