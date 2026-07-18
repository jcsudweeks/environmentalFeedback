% rg1_1CO_bifurcation.m
% Generates bifurcation data for the host-rich, large-burst-size regime
% (rg1_1CO).
%
% Output: data/rg1_1CO_run.mat -- loaded by figure1_combinedBifurcation.m,
% verifyBifurcations.m, cooperatorLimitCycleStability.m, and
% defectorLimitCycle.m.

addpath('../functions')

r = 5;
d = 60/108.9;               % derived from latent period time
kappa = 1;
mu1 = 9.9*(10^(-6.5));      % derived from adsorption rate
mu2 = mu1^2;
xi = r/1000;                 % sets carrying capacity to 100
lambda = 210;                 % from Turner and Chao paper

% Generate bifurcation branches over rho
inputRhoVec = [.01 1.2];
[stableCell, unStableCell, rhoVec] = bifurcationGenerator([mu1,mu2,r,d,kappa,lambda,xi], inputRhoVec);

% Collision points where mixed-equilibrium branches meet the
% cooperator-only / defector-only branches
collisions = mixedCollisions([mu1,mu2,r,d,kappa,lambda,xi]);

% Transcritical threshold between the uninfected and defector-only branches
rhoTUD = (d*kappa*xi)./(2*lambda*mu1*r);

bifVec = sort([collisions', rhoTUD]);

save('../data/rg1_1CO_run.mat')
