% mixedLimitCycle.m
% Traces the mixed limit cycle (H_M > 0) across rho in the host-rich,
% small-burst-size regime (rg1_2CO), between the cooperator-limit-cycle
% Hopf bifurcation and the global (homoclinic) bifurcation.
%
% For each rho, starts from the stable cooperator-only equilibrium,
% introduces a small defector population, and integrates long enough to
% converge onto the limit cycle; records its min/max amplitude.
%
% Output: mixedLimitCycles.mat -- loaded by figure1.m
% and mixedLCPhasePlaneData.m.

addpath('../functions')

load('../data/rg1_2CO_run.mat')
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Cooperator-only equations don't depend on rho, so this equilibrium is
% the same for every rho tested below.
coopEqs = cooperatorEquilibria([mu1,mu2,r,d,kappa,lambda,xi]);
coopPop = coopEqs(2,:);   % E_C^S -- stable cooperator-only equilibrium

% Sample rho between the cooperator-limit-cycle Hopf bifurcation (0.5940)
% and the global bifurcation (0.6287), staying just inside both bounds.
testingRho = linspace(0.5940+.0001, 0.6287-.0001, 5);
tspan = 0:.1:20000;

storeMinVec    = -1.*ones(length(testingRho),6);
storeMaxVec    = -1.*ones(length(testingRho),6);
MlimCyclePoint = -1.*ones(length(testingRho),6);
storePopMinVec = -1.*ones(length(testingRho),2);
storePopMaxVec = -1.*ones(length(testingRho),2);

for i = 1:length(testingRho)
    params = [mu1,mu2,r, d, kappa, lambda, xi, testingRho(i)];

    % Start from the cooperator equilibrium with a small defector introduced
    pop = [coopPop(1:5), 1];
    [T, pop] = ode45(@modelEqs, tspan, pop, opts, params);

    % Take the second half of the trajectory as converged onto the limit cycle
    limitCycletime = pop(end-50000:end,:);

    storeMaxVec(i,:)    = max(limitCycletime);
    storeMinVec(i,:)    = min(limitCycletime);
    storePopMaxVec(i,1) = max(sum(limitCycletime(:,1:4),2));
    storePopMaxVec(i,2) = max(sum(limitCycletime(:,5:6),2));
    storePopMinVec(i,1) = min(sum(limitCycletime(:,1:4),2));
    storePopMinVec(i,2) = min(sum(limitCycletime(:,5:6),2));
    MlimCyclePoint(i,:) = pop(end,:);
end

save('../data/mixedLimitCycles.mat','storeMaxVec','storeMinVec', ...
    'storePopMaxVec','storePopMinVec','MlimCyclePoint','testingRho')
