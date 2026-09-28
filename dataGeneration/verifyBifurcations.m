% verifyBifurcations.m
% Computes and prints all bifurcation values, regime by regime.
%
% Analytical expressions used where available; eigenvalue bisection for
% stability-change bifurcations (rho_M, rho_HM, rho_HD).
%
% rho_LC is loaded from limitCycleStabilityData.mat (pre-computed
% by cooperatorLimitCycleStability.m).
% rho_HC is the value found by homoclinicRho.m (pre-computed).
%
% Required: ../data/rl1_run.mat, rg1_1CO_run.mat, rg1_2CO_run.mat,
% limitCycleStabilityData.mat, homoclinicRhoData.mat

addpath('../functions')

EIGEN_TOL = 1e-6;
SCAN_PTS  = 50;

% Symbolic Jacobian computed once; jacFun(eq, params, rho) used throughout.
fprintf('Computing symbolic Jacobian (one-time cost)...\n')
J_fun  = buildSymbolicJacobian();
jacFun = @(eq, params, rho) J_fun(eq(1), eq(2), eq(3), eq(4), eq(5), eq(6), ...
    params(1), params(2), params(3), params(4), params(5), params(6), params(7), rho);

%% ================================================================
%  REGIME: rl1
%% ================================================================
load('../data/rl1_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi','stableCell')

params = [mu1, mu2, r, d, kappa, lambda, xi];

rho_TUD = d*kappa*xi / (2*lambda*mu1*r);
[ecColl, edColl] = mixedCollisions(params);

rows = {};
rows{end+1} = {'rho_TUD',  rho_TUD, 'analytical'};
for k = 1:length(ecColl)
    rows{end+1} = {sprintf('rho_TEC(%d)',k), ecColl(k), 'mixedCollisions (cubic)'};
end
for k = 1:length(edColl)
    rows{end+1} = {sprintf('rho_TED(%d)',k), edColl(k), 'mixedCollisions (quadratic)'};
end

% E_M stability changes (if any)
br_M = findAllBrackets(@mixedEquilibrium, jacFun, params, 0.05, 0.99, SCAN_PTS);
rl1_rhoM_vals = zeros(size(br_M, 1), 1);
for k = 1:size(br_M, 1)
    val = bisectEigen(@mixedEquilibrium, jacFun, params, br_M(k,1), br_M(k,2), EIGEN_TOL);
    rows{end+1} = {sprintf('rho_M(%d)',k), val, 'eigenvalue bisection at E_M'};
    rl1_rhoM_vals(k) = val;
end

printRegimeTable('rl1', rows)

rl1_bifVec_all = [rho_TUD; ecColl(:); edColl(:); rl1_rhoM_vals(:)];
rl1_bifVec = sort(rl1_bifVec_all(~isnan(rl1_bifVec_all) & rl1_bifVec_all > 0 & rl1_bifVec_all <= 1))';

%% External stability condition rho < 1/(1+mu1*V_C) at the C-only equilibrium
VcDense = stableCell{3}(1,5)
extStabCond = 1/(1+(mu1*VcDense))

%% ================================================================
%  REGIME: rg1_1CO
%% ================================================================
load('../data/rg1_1CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi','unStableCell')

params = [mu1, mu2, r, d, kappa, lambda, xi];

% Analytical
rho_SD_1CO  = (-2*mu1^3 - 9*mu1*mu2*r + 2*sqrt((mu1^2 + 3*mu2*r)^3)) * ...
              (d*kappa*xi) / (2*mu2*r^2*(mu1^2 + 4*mu2*r)*lambda);

rho_TUD_1CO = d*kappa*xi / (2*lambda*mu1*r);

[ecColl_1CO, edColl_1CO] = mixedCollisions(params);

% rho_M (~0.3894): targeted scan
br_M_1CO = findAllBrackets(@mixedEquilibrium, jacFun, params, 0.36, 0.43, SCAN_PTS);
if size(br_M_1CO, 1) >= 1
    rho_M_1CO = bisectEigen(@mixedEquilibrium, jacFun, params, ...
                            br_M_1CO(1,1), br_M_1CO(1,2), EIGEN_TOL);
else
    rho_M_1CO = NaN;
    warning('verifyBifurcations: rho_M not found in [0.36, 0.43] for rg1_1CO')
end

% rho_HD: last sign change in E_D^S eigenvalues
br_HD_1CO = findAllBrackets(@defectorEqStable, jacFun, params, ...
                            rho_SD_1CO + 0.05, 0.99, SCAN_PTS);
if size(br_HD_1CO, 1) >= 1
    rho_HD_1CO = bisectEigen(@defectorEqStable, jacFun, params, ...
                             br_HD_1CO(end,1), br_HD_1CO(end,2), EIGEN_TOL);
else
    rho_HD_1CO = NaN;
    warning('verifyBifurcations: rho_HD not found for rg1_1CO')
end

% rho_LC: pre-computed
lcd = load('../data/limitCycleStabilityData.mat', 'rhoTransition');
rho_LC = lcd.rhoTransition;

rows = {};
rows{end+1} = {'rho_SD',    rho_SD_1CO,   'analytical'};
for k = 1:length(ecColl_1CO)
    rows{end+1} = {sprintf('rho_TEC(%d)',k), ecColl_1CO(k), 'mixedCollisions (cubic)'};
end
rows{end+1} = {'rho_M',     rho_M_1CO,    'eigenvalue bisection at E_M'};
rows{end+1} = {'rho_TUD',   rho_TUD_1CO,  'analytical'};
rows{end+1} = {'rho_LC',    rho_LC,       'cooperatorLimitCycleStability (pre-computed)'};
for k = 1:length(edColl_1CO)
    rows{end+1} = {sprintf('rho_TED(%d)',k), edColl_1CO(k), 'mixedCollisions (quadratic)'};
end
rows{end+1} = {'rho_HD',    rho_HD_1CO,   'eigenvalue bisection at E_D^S'};

printRegimeTable('rg1_1CO', rows)

rg11_bifVec_all = [rho_SD_1CO; rho_TUD_1CO; ecColl_1CO(:); edColl_1CO(:); rho_M_1CO; rho_HD_1CO; rho_LC];
rg11_bifVec = sort(rg11_bifVec_all(~isnan(rg11_bifVec_all) & rg11_bifVec_all > 0 & rg11_bifVec_all <= 1))';

%% External stability condition rho < 1/(1+mu1*V_C) at the C-only equilibrium
VcDense = unStableCell{4}(1,5)
extStabCond = 1/(1+(mu1*VcDense))

%% ================================================================
%  REGIME: rg1_2CO
%% ================================================================
load('../data/rg1_2CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi','stableCell')

params = [mu1, mu2, r, d, kappa, lambda, xi];

% Analytical
rho_SD_2CO  = (-2*mu1^3 - 9*mu1*mu2*r + 2*sqrt((mu1^2 + 3*mu2*r)^3)) * ...
              (d*kappa*xi) / (2*mu2*r^2*(mu1^2 + 4*mu2*r)*lambda);

rho_TUD_2CO = d*kappa*xi / (2*lambda*mu1*r);
[ecColl_2CO, edColl_2CO] = mixedCollisions(params);

% rho_M (~0.3894): targeted scan
br_M_2CO = findAllBrackets(@mixedEquilibrium, jacFun, params, 0.36, 0.43, SCAN_PTS);
if size(br_M_2CO, 1) >= 1
    rho_M_2CO = bisectEigen(@mixedEquilibrium, jacFun, params, ...
                            br_M_2CO(1,1), br_M_2CO(1,2), EIGEN_TOL);
else
    rho_M_2CO = NaN;
    warning('verifyBifurcations: rho_M not found in [0.36, 0.43] for rg1_2CO')
end

% rho_HM (~0.5940): targeted scan
br_HM_2CO = findAllBrackets(@mixedEquilibrium, jacFun, params, 0.56, 0.63, SCAN_PTS);
if size(br_HM_2CO, 1) >= 1
    rho_HM_2CO = bisectEigen(@mixedEquilibrium, jacFun, params, ...
                             br_HM_2CO(1,1), br_HM_2CO(1,2), EIGEN_TOL);
else
    rho_HM_2CO = NaN;
    warning('verifyBifurcations: rho_HM not found in [0.56, 0.63] for rg1_2CO')
end

% rho_HD: last sign change in E_D^S eigenvalues
br_HD_2CO = findAllBrackets(@defectorEqStable, jacFun, params, ...
                            rho_SD_2CO + 0.05, 0.99, SCAN_PTS);
if size(br_HD_2CO, 1) >= 1
    rho_HD_2CO = bisectEigen(@defectorEqStable, jacFun, params, ...
                             br_HD_2CO(end,1), br_HD_2CO(end,2), EIGEN_TOL);
else
    rho_HD_2CO = NaN;
    warning('verifyBifurcations: rho_HD not found for rg1_2CO')
end

% rho_HC: saved by homoclinicRho.m
tmp = load('../data/homoclinicRhoData.mat', 'rhoHC'); rho_HC = tmp.rhoHC;


rows = {};
rows{end+1} = {'rho_SD',    rho_SD_2CO,  'analytical'};
for k = 1:length(ecColl_2CO)
    rows{end+1} = {sprintf('rho_TEC(%d)',k), ecColl_2CO(k), 'mixedCollisions (cubic)'};
end
rows{end+1} = {'rho_M',     rho_M_2CO,   'eigenvalue bisection at E_M'};
rows{end+1} = {'rho_HM',    rho_HM_2CO,  'eigenvalue bisection at E_M'};
rows{end+1} = {'rho_TUD',   rho_TUD_2CO, 'analytical'};
for k = 1:length(edColl_2CO)
    rows{end+1} = {sprintf('rho_TED(%d)',k), edColl_2CO(k), 'mixedCollisions (quadratic)'};
end
rows{end+1} = {'rho_HC',    rho_HC,      'homoclinicRho (pre-computed)'};
rows{end+1} = {'rho_HD',    rho_HD_2CO,  'eigenvalue bisection at E_D^S'};

printRegimeTable('rg1_2CO', rows)

rg12_bifVec_all = [rho_SD_2CO; rho_TUD_2CO; ecColl_2CO(:); edColl_2CO(:); rho_M_2CO; rho_HM_2CO; rho_HD_2CO; rho_HC];
rg12_bifVec = sort(rg12_bifVec_all(~isnan(rg12_bifVec_all) & rg12_bifVec_all > 0 & rg12_bifVec_all <= 1))';

save('../data/verifiedBifurcations.mat', 'rl1_bifVec', 'rg11_bifVec', 'rg12_bifVec')
fprintf('Saved verifiedBifurcations.mat\n')

%% External stability condition rho < 1/(1+mu1*V_C) at the C-only equilibrium
VcDense = stableCell{4}(1,5)
extStabCond = 1/(1+(mu1*VcDense))


%% ================================================================
%  LOCAL FUNCTIONS
%% ================================================================

function printRegimeTable(regimeName, rows)
fprintf('\n==============================\n')
fprintf('REGIME: %s\n', regimeName)
fprintf('==============================\n')
vals = cellfun(@(r) r{2}, rows);
[~, idx] = sort(vals);
rows = rows(idx);
fprintf('%-22s  %-12s  %s\n', 'Bifurcation', 'rho', 'Method')
fprintf('%s\n', repmat('-', 1, 70))
for i = 1:length(rows)
    if isnan(rows{i}{2})
        fprintf('%-22s  %-12s  %s\n', rows{i}{1}, 'NOT FOUND', rows{i}{3})
    else
        fprintf('%-22s  %-12.8f  %s\n', rows{i}{1}, rows{i}{2}, rows{i}{3})
    end
end
end


function J_fun = buildSymbolicJacobian()
% Builds a fast matlabFunction handle for the Jacobian of modelEqs.
% Called once at script startup; evaluation is then cheap.
syms U C D M VC VD mu1 mu2 r d kappa lambda xi rho real
f = [r*U - xi*U^2 - U*mu1*(VC+VD) - U*mu2*(VC^2 + 2*VC*VD + VD^2);
     U*VC*(mu1 + mu2*VC) - d*C;
     U*VD*(mu1 + mu2*VD) - d*D;
     2*mu2*U*VC*VD - d*M;
     lambda*d*(2*C) - d*kappa*VC;
     lambda*d*((1+rho)*M + 2*rho*D) - d*kappa*VD];
J_sym = jacobian(f, [U, C, D, M, VC, VD]);
J_fun = matlabFunction(J_sym, 'Vars', {U, C, D, M, VC, VD, mu1, mu2, r, d, kappa, lambda, xi, rho});
end


function eq = mixedEquilibrium(params, rho)
% Analytical mixed equilibrium [HU, HC, HD, HM, VC, VD].
% Returns NaN(1,6) if the equilibrium does not exist at this rho.
% Derivation: equilibrium conditions reduce to
%   S = VC+VD = (1-rho)/(rho*mu1),  HU = (r - (1-rho)/rho^2)/xi,
%   VC from HU = d*kappa/(2*lambda*(mu1+mu2*VC)).
mu1 = params(1); mu2 = params(2); r = params(3); d  = params(4);
kappa = params(5); lambda = params(6); xi = params(7);

R = r - (1 - rho)/rho^2;
if R <= 0, eq = NaN(1,6); return, end
HU = R / xi;
VC = (d*kappa*xi/(2*lambda*R) - mu1) / mu2;
S  = (1 - rho) / (rho*mu1);
VD = S - VC;
if VC <= 0 || VD <= 0, eq = NaN(1,6); return, end
HC = kappa*VC / (2*lambda);
HD = HU*VD*(mu1 + mu2*VD) / d;
HM = 2*mu2*HU*VC*VD / d;
eq = [HU, HC, HD, HM, VC, VD];
end


function eq = defectorEqStable(params, rho)
% Stable defector-only equilibrium (highest V_D). NaN(1,6) if none exist.
defs = defectorOnlyEquilibria(rho, params);
if isempty(defs), eq = NaN(1,6); else, eq = defs(end,:); end
end


function brackets = findAllBrackets(eqFun, jacFun, params, rhoMin, rhoMax, nPts)
% Coarse scan; returns [rhoLo, rhoHi] rows for each sign change in
% max(real(eigenvalues)). eqFun(params,rho) returns eq (NaN if nonexistent).
rhoV = linspace(rhoMin, rhoMax, nPts);
evV  = NaN(1, nPts);
for k = 1:nPts
    eq = eqFun(params, rhoV(k));
    if ~any(isnan(eq))
        evV(k) = max(real(eig(jacFun(eq, params, rhoV(k)))));
    end
end
brackets = zeros(0, 2);
for k = 1:nPts-1
    if ~isnan(evV(k)) && ~isnan(evV(k+1)) && sign(evV(k)) ~= sign(evV(k+1))
        brackets(end+1, :) = [rhoV(k), rhoV(k+1)]; %#ok<AGROW>
    end
end
end


function rhoOut = bisectEigen(eqFun, jacFun, params, rhoLo, rhoHi, tol)
% Bisects on max(real(eigenvalues)) = 0 between rhoLo and rhoHi.
    function ev = maxEig(rho)
        eq = eqFun(params, rho);
        ev = max(real(eig(jacFun(eq, params, rho))));
    end
evLo = maxEig(rhoLo);
evHi = maxEig(rhoHi);
if sign(evLo) == sign(evHi)
    error('bisectEigen: no sign change in [%.6f, %.6f]', rhoLo, rhoHi)
end
while rhoHi - rhoLo > tol
    rhoMid = (rhoLo + rhoHi) / 2;
    evMid  = maxEig(rhoMid);
    if sign(evMid) == sign(evLo)
        rhoLo = rhoMid; evLo = evMid;
    else
        rhoHi = rhoMid;
    end
end
rhoOut = (rhoLo + rhoHi) / 2;
end
