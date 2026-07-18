function defEq = defectorOnlyEquilibria(rho, params)
% Finds defector-only equilibria numerically.
%
% Analytically eliminates H_U and H_D from the equilibrium conditions,
% reducing to a single equation in V_D which is solved with fzero.
%
% At equilibrium (H_C = H_M = V_C = 0):
%   dV_D = 0  =>  H_D = kappa*V_D / (2*rho*lambda)
%   dH_D = 0  =>  H_U = d*kappa / (2*rho*lambda*(mu1 + mu2*V_D))
%   dH_U = 0  =>  r - xi*H_U - V_D*(mu1 + mu2*V_D) = 0
%
% Substituting gives f(V_D) = 0, solved numerically.
%
% Returns rows of [H_U, 0, H_D, 0, 0, V_D], sorted by V_D ascending.
% The stable equilibrium has the higher V_D (per user verification).
%
% Usage:
%   defEqs = defectorOnlyEquilibria(rho, [mu1,mu2,r,d,kappa,lambda,xi])
%   initDef = defEqs(end,:);   % stable equilibrium

mu1    = params(1);
mu2    = params(2);
r      = params(3);
d      = params(4);
kappa  = params(5);
lambda = params(6);
xi     = params(7);

% Reduced equation in V_D
f = @(VD) r - xi .* (d.*kappa ./ (2.*rho.*lambda .* (mu1 + mu2.*VD))) ...
             - VD .* (mu1 + mu2.*VD);

% Scan over plausible V_D range to find sign changes
VDscan = linspace(1e-3, 2*r/mu1, 2000);
fvals  = arrayfun(f, VDscan);

defEq = zeros(0, 6);
for i = 1:length(VDscan)-1
    if fvals(i) * fvals(i+1) < 0
        VD_root = fzero(f, [VDscan(i), VDscan(i+1)]);
        HU = d*kappa / (2*rho*lambda*(mu1 + mu2*VD_root));
        HD = kappa*VD_root / (2*rho*lambda);
        defEq(end+1, :) = [HU, 0, HD, 0, 0, VD_root]; %#ok<AGROW>
    end
end

% Sort by V_D ascending — stable equilibrium is last (highest V_D)
if ~isempty(defEq)
    [~, idx] = sort(defEq(:, 6));
    defEq    = defEq(idx, :);
end
end
