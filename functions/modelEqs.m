function [dX] = modelEqs(tspan, X, pars)
% Right-hand side of the donation-game host-phage coevolution ODE, for
% use with ode45 (or similar). See manuscript for state/parameter
% definitions and derivation.
%
% X    = [H_U, H_C, H_D, H_M, V_C, V_D]
% pars = [mu1, mu2, r, d, kappa, lambda, xi, rho]

mu1 = pars(1);
mu2 = pars(2);
r = pars(3);
d = pars(4);
kappa = pars(5);
lambda = pars(6);
xi = pars(7);
rho = pars(8);

alpha = 1;
beta = 0;
g = 1+rho;
delta = rho;

U = X(1);
C = X(2);
D =X(3);
M =X(4);
PC = X(5);
PD =X(6);

dX=zeros(6,1);
dX(1) =  r.*U - xi.*(U.^2) - U.*mu1.*(PC+PD) - U.*mu2.*(PC.^2 + 2.*PC.*PD + PD.^2) ;
dX(2)= U.*PC.*(mu1 + mu2.*PC) - d.*C;
dX(3)= U.*PD.*(mu1+mu2.*PD) - d.*D;
dX(4)= 2.*mu2.*U.*PC.*PD - d.*M;
dX(5)= lambda.*d.*(2.*alpha.*C + beta.*M) - d.*kappa.*PC;
dX(6)= lambda.*d.*(g.*M + 2.*delta.*D) - d.*kappa.*PD;
end