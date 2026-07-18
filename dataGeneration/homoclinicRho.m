% homoclinicRho.m
% Locates rho_HC for the rg1_2CO parameter set via limit cycle continuation
% and bisection.
%
% Method:
%   1. Coarse rho scan with continuation to bracket the mixed LC boundary.
%   2. Bisect the bracketed interval to tolerance 1e-4.
%
% Detection: max(V_C) in the tail of each simulation. Above threshold =
% mixed LC exists; near zero = LC gone (extinction attractor for rho_HC < rho < rho_S).
% V_C is used rather than V_D because near rho_HC the limit cycle approaches
% E_C^U where V_D ~ 0; V_C is large there, so max(V_C) is reliable throughout
% the mixed LC existence range.

addpath('../functions')

load('../data/rg1_2CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi')
load('../data/mixedLimitCycles.mat', 'MlimCyclePoint')

opts      = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
VC_thresh = 1;
t_scan    = 100000;
t_bisect  = 200000;

%% 1. Coarse scan to bracket the transition
rhoScan  = linspace(0.60, 0.70, 11);
pop0     = MlimCyclePoint(1,:);
lcExists = false(size(rhoScan));

fprintf('Coarse scan:\n')
for i = 1:length(rhoScan)
    rho   = rhoScan(i);
    param = [mu1, mu2, r, d, kappa, lambda, xi, rho];
    [~, pop] = ode45(@modelEqs, [0, t_scan], pop0, opts, param);
    tail  = pop(round(0.8*size(pop,1)):end, :);
    vc_max = max(tail(:, 5));
    lcExists(i) = vc_max > VC_thresh;
    if lcExists(i)
        pop0 = pop(end, :);   % continuation only when LC exists
    end
    if lcExists(i), lcStr = 'exists'; else, lcStr = 'gone'; end
    fprintf('  rho = %.4f   max V_C (tail) = %.4g   LC: %s\n', ...
        rho, vc_max, lcStr)
end

% Find first false
transIdx = find(~lcExists, 1);
if isempty(transIdx) || transIdx == 1
    error('homoclinicRho: transition not bracketed in [%.2f, %.2f] — extend rhoScan', ...
        rhoScan(1), rhoScan(end))
end
rhoLo = rhoScan(transIdx - 1);
rhoHi = rhoScan(transIdx);
fprintf('\nBracket: [%.4f, %.4f]\n', rhoLo, rhoHi)

% pop0 is already the endpoint from rhoLo (last successful continuation)
pop0_lo = pop0;

%% 2. Bisection
tol = 1e-7;
fprintf('\nBisecting to tol = %.0e ...\n', tol)
while (rhoHi - rhoLo) > tol
    rhoMid = (rhoLo + rhoHi) / 2;
    param  = [mu1, mu2, r, d, kappa, lambda, xi, rhoMid];
    [~, pop] = ode45(@modelEqs, [0, t_bisect], pop0_lo, opts, param);
    tail   = pop(round(0.8*size(pop,1)):end, :);
    vc_max = max(tail(:, 5));
    if vc_max > VC_thresh
        rhoLo   = rhoMid;
        pop0_lo = pop(end, :);   % continuation in lower half only
    else
        rhoHi = rhoMid;
    end
    fprintf('  [%.6f, %.6f]   mid = %.6f   max V_C = %.4g\n', ...
        rhoLo, rhoHi, rhoMid, vc_max)
end

rhoHC = (rhoLo + rhoHi) / 2;
fprintf('\nrho_HC = %.6f\n', rhoHC)
save('../data/homoclinicRhoData.mat', 'rhoHC')
