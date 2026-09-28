% defectorLimitCycle.m
% Traces the defector-only limit cycle across rho in the host-rich,
% large-burst-size regime (rg1_1CO), from just above the defector Hopf
% bifurcation (rho = 0.8355) up to rho = 1.
%
% Each rho is initialized from the previous rho's converged limit-cycle
% point (continuation), so the trajectory only needs to relax onto the
% (slightly shifted) new cycle rather than find it from scratch.
%
% Output: defectorLimitCycleAmps.mat -- loaded by
% figure1.m.

addpath('../functions')

load('../data/rg1_1CO_run.mat')
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

n = 10;
rhoVec = linspace(0.8355+.01,1-.01,n);

DLCMaxAmp     = -1.*ones(n,6);
DLCMinAmp     = -1.*ones(n,6);
DLCPoint      = -1.*ones(n,6);   % point on the limit cycle, used to seed the next rho
DLCMaxHostAmp = -1.*ones(n,1);
DLCMinHostAmp = -1.*ones(n,1);

% First rho run separately: no prior limit-cycle point to start from, so
% it needs a longer run to converge from an arbitrary initial condition.
param = [mu1,mu2,r, d, kappa, lambda, xi,rhoVec(1)];
popVec = [10,0,10,0,0,10];
[t,popVec] = ode45(@modelEqs,0:.1:1000, popVec,opts, param);
DLCMaxAmp(1,:) = max(popVec(end-500:end,:));
DLCMinAmp(1,:) = min(popVec(end-500:end,:));
DLCPoint(1,:) = popVec(end,:);
DLCMaxHostAmp(1)= max(sum(popVec(end-500:end,1:4),2));
DLCMinHostAmp(1)= min(sum(popVec(end-500:end,1:4),2));

for i =2:n
 param = [mu1,mu2,r, d, kappa, lambda, xi,rhoVec(i)];
 popVec = DLCPoint(i-1,:);   % continuation from previous rho's cycle
 [t,popVec] = ode45(@modelEqs, 0:.1:1000, popVec,opts, param);

DLCMaxAmp(i,:) = max(popVec(end-500:end,:));
DLCMinAmp(i,:) = min(popVec(end-500:end,:));
DLCMaxHostAmp(i) = max(sum(popVec(end-500:end,1:4),2));
DLCMinHostAmp(i) = min(sum(popVec(end-500:end,1:4),2));

DLCPoint(i,:) = popVec(end,:);

end

DLCMax = [DLCMaxAmp(:,1),DLCMaxAmp(:,3),DLCMaxAmp(:,6)];
DLCMin = [DLCMinAmp(:,1),DLCMinAmp(:,3),DLCMinAmp(:,6)];
defectorLimitCycleRhoVec = rhoVec;

save('../data/defectorLimitCycleAmps.mat',"DLCMax","DLCMin","DLCPoint","defectorLimitCycleRhoVec",...
    "DLCMaxHostAmp","DLCMinHostAmp")

