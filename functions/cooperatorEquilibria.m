function coopEq = cooperatorEquilibria(params)
% Finds cooperator-only equilibria numerically.
%
% Analytically eliminates H_U and H_C from the equilibrium conditions,
% reducing to a single equation in V_C which is solved with fzero.
%
% At equilibrium (H_D = H_M = V_D = 0):
%   dV_C = 0  =>  H_C = kappa*V_C / (2*lambda)
%   dH_C = 0  =>  H_U = d*kappa / (2*lambda*(mu1 + mu2*V_C))
%   dH_U = 0  =>  r - xi*H_U - V_C*(mu1 + mu2*V_C) = 0
%
% Note: rho does not appear in the cooperator equations, so these
% equilibria are independent of rho.
%
% Returns rows of [H_U, H_C, 0, 0, V_C, 0], sorted by V_C ascending.
% In the host rich small burst regime, two equilibria exist:
%   row 1 (lower V_C): verify against bifurcation diagram which is E_C^U
%   row 2 (higher V_C): verify against bifurcation diagram which is E_C^S
%
% Usage:
%   coopEqs = cooperatorEquilibria([mu1, mu2, r, d, kappa, lambda, xi])

mu1    = params(1);
mu2    = params(2);
r      = params(3);
d      = params(4);
kappa  = params(5);
lambda = params(6);
xi     = params(7);

% Reduced equation in V_C (identical to defector case with rho = 1)
f = @(VC) r - xi .* (d.*kappa ./ (2.*lambda .* (mu1 + mu2.*VC))) ...
             - VC .* (mu1 + mu2.*VC);

% Scan over plausible V_C range to find sign changes
VCscan = linspace(1e-3, 2*r/mu1, 2000);
fvals  = arrayfun(f, VCscan);

coopEq = zeros(0, 6);
for i = 1:length(VCscan)-1
    if fvals(i) * fvals(i+1) < 0
        VC_root = fzero(f, [VCscan(i), VCscan(i+1)]);
        HU = d*kappa / (2*lambda*(mu1 + mu2*VC_root));
        HC = kappa*VC_root / (2*lambda);
        coopEq(end+1, :) = [HU, HC, 0, 0, VC_root, 0]; %#ok<AGROW>
    end
end

% Sort by V_C ascending
if ~isempty(coopEq)
    [~, idx] = sort(coopEq(:, 5));
    coopEq   = coopEq(idx, :);
end
end
