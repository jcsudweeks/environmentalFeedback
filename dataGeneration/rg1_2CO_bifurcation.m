% rg1_2CO_bifurcation.m
% Generates bifurcation data for the host-rich, small-burst-size regime
% (rg1_2CO).
%
% Output: data/rg1_2CO_run.mat -- loaded by figure1.m,
% verifyBifurcations.m, and mixedLimitCycle.m.

addpath('../functions')

r = 5;                        % must be greater than 1
d = 60/108.9;                 % derived from latent period time
kappa = 1;
mu1 = 9.9*(10^(-6.5));        % derived from adsorption rate
mu2 = mu1^2;
xi = r/1000;                   % sets carrying capacity to 1000

% lambda must lie between a lower and upper threshold for this regime;
% take the midpoint.
ltNum = (-2*mu1^3 - 9*mu1*mu2*r + 2*sqrt((mu1^2+3*r*mu2)^3))*d*kappa*xi;
ltDen = 2*mu2*r^2*(mu1^2+4*mu2*r);
lowerThresh = ltNum/ltDen;

utNum = d*kappa*xi;
utDen = 2*mu1*r;
upperThresh = utNum/utDen;

lambda = floor((lowerThresh+upperThresh)/2);   % right in the middle

% Generate bifurcation branches over rho
inputRhoVec = [0.01, 2];
[stableCell, unStableCell, rhoVec] = bifurcationGenerator([mu1,mu2,r,d,kappa,lambda,xi], inputRhoVec);

save('../data/rg1_2CO_run.mat')