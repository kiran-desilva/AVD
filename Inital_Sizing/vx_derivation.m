syms Tw ws k Vx alpha theta_max

[~,~,~,rho] = atmosisa(distdim(45000,'ft','m'))




e = .8
a = 0.8296
b = 0.25
n = 1
AR = 8 
CD_min = .02

k = 1/(pi*AR*e)

% theta_max = asind((Tw) - sqrt(4*CD_min*k))
% 
% figure
% fplot(theta_max,[0 10])

% Vx_calculated(ws) = sqrt( (2/rho) * (ws) * sqrt(k/CD_min) * sqrt(1-(alpha^2)))

Vx(ws) = sqrt( (2/rho) * (ws) * sqrt(k/CD_min) * 1)

abs_ceiling_constraint(ws) = simplify((a/b) * ( ( (0.5*rho*(Vx^2)*CD_min)/(ws) ) + ( ((n^2)*(ws))/(0.5*rho*(Vx^2)*pi*AR*e) ) ))

double(abs_ceiling_constraint)
