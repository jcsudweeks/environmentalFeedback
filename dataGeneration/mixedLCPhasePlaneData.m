% mixedLCPhasePlaneData.m
% Generates trajectory data for the mixed limit cycle phase plane in
% (V_C, V_D) space, for the rg1_2CO parameter set. After a long
% transient, one full loop of the limit cycle (from one V_C peak to the
% next) is kept for each rho.
%
% Output: mixedLCPhasePlaneData.mat -- loaded by figure2.m.

addpath('../functions')

load('../data/rg1_2CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi')
load('../data/mixedLimitCycles.mat', 'MlimCyclePoint')

opts   = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
t_plot = 200000;

rhoVec = [0.60, 0.61, 0.62];          % representative values inside LC regime

VCdata  = cell(length(rhoVec), 1);
VDdata  = cell(length(rhoVec), 1);
fCdata  = cell(length(rhoVec), 1);

for k = 1:length(rhoVec)
    rho   = rhoVec(k);
    param = [mu1, mu2, r, d, kappa, lambda, xi, rho];
    pop0  = MlimCyclePoint(1, :);

    [~, pop] = ode45(@modelEqs, [0, t_plot], pop0, opts, param);

    [~, locs]   = findpeaks(pop(:, 5));
    seg         = locs(end-1):locs(end);    % last full loop, V_C peak to V_C peak
    VCdata{k}   = pop(seg, 5);
    VDdata{k}   = pop(seg, 6);
    fCdata{k}   = VCdata{k} ./ (VCdata{k} + VDdata{k});
end

threshold = mu1 / mu2;

save('../data/mixedLCPhasePlaneData.mat', 'VCdata', 'VDdata', 'fCdata', 'rhoVec', 'threshold')
