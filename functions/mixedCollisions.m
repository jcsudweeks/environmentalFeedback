function [ecCollisions, edCollisions] = mixedCollisions(params)
% Returns rho values of transcritical bifurcations between E_M and
% cooperator-only (E_C) or defector-only (E_D) equilibria.
%
% Two families (requires mu2 = mu1^2; mu = mu1 throughout):
%
%   E_M <-> E_C  (ecCollisions): roots of the cubic
%       d*kappa*xi*rho^3 - 2*lambda*mu*r*rho^2 - 2*lambda*mu*rho + 2*lambda*mu
%
%   E_M <-> E_D  (edCollisions): direct quadratic formula
%       (lambda*mu +/- sqrt(lambda*mu*(lambda*mu - 2*d*kappa*xi + 4*lambda*mu*r)))
%       / (d*kappa*xi - 2*lambda*mu*r)
%
% With two output arguments returns the two families separately (sorted).
% With one output argument returns all roots combined and sorted.
%
% Usage:
%   [ecColl, edColl] = mixedCollisions([mu1,mu2,r,d,kappa,lambda,xi])
%   allColl          = mixedCollisions([mu1,mu2,r,d,kappa,lambda,xi])

mu     = params(1);   % mu1; mu2 = mu1^2 is assumed, not used directly
r      = params(3);
d      = params(4);
kappa  = params(5);
lambda = params(6);
xi     = params(7);

%% E_M <-> E_C collisions: cubic polynomial
cr = roots([d*kappa*xi,  -2*lambda*mu*r,  -2*lambda*mu,  2*lambda*mu]);
ecCollisions = sort(real(cr(abs(imag(cr)) < 1e-10 & real(cr) > 0 & real(cr) < 1)));

%% E_M <-> E_D collisions: quadratic formula
disc  = lambda*mu * (lambda*mu - 2*d*kappa*xi + 4*lambda*mu*r);
denom = d*kappa*xi - 2*lambda*mu*r;
if disc >= 0 && abs(denom) > 1e-15
    qr = [(lambda*mu + sqrt(disc)); (lambda*mu - sqrt(disc))] / denom;
    edCollisions = sort(real(qr(real(qr) > 0 & real(qr) < 1)));
else
    edCollisions = zeros(0, 1);
end

%% Combine if single output requested
if nargout < 2
    ecCollisions = sort([ecCollisions(:); edCollisions(:)]);
end
end
