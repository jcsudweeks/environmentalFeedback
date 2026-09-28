% viralExtinctionData.m
% Generates trajectory and basin-of-attraction data for the viral
% extinction scenario in the rg1_2CO parameter set at rho = 0.65.
% Starting from E_C^S, defectors are introduced; cooperators are driven
% into the E_U basin and virus goes extinct.
%
% Output: viralExtinctionData.mat -- loaded by figure4.m.

addpath('../functions')

load('../data/rg1_2CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi')

rho    = 0.65;
t1     = 50;     % cooperator-only phase
t2     = 1000;   % post-introduction phase
t_run  = 1000;   % integration time for convergence and basin tests
VD_ts  = 100;    % V_D introduced for time series / basin trajectory

param = [mu1, mu2, r, d, kappa, lambda, xi, rho];
opts  = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

coopEqs = cooperatorEquilibria(param(1:7));
pop0    = coopEqs(end, :);   % E_C^S

% Phase 1: cooperators alone
[t_c, pop_c] = ode45(@modelEqs, [0, t1], pop0, opts, param);

% Phase 2: introduce defectors
pop1    = pop_c(end, :);
pop1(6) = VD_ts;
[t_d, pop_d] = ode45(@modelEqs, [0, t2], pop1, opts, param);
t_d = t_d + t1;

% Reference axis limits for (V_C, V_D) phase plane figures
ref_xl = [0, max(pop_d(:,5)) * 1.05];
ref_yl = [0, max(pop_d(:,6)) * 1.05];

% Convergence trajectories
VD_conv = [1, 10, 100, 1000, 10000];
conv_trajs = cell(length(VD_conv), 1);
for k = 1:length(VD_conv)
    ic    = pop0;
    ic(6) = VD_conv(k);
    [~, pop_k] = ode45(@modelEqs, [0, t_run], ic, opts, param);
    conv_trajs{k} = pop_k;
end

% Basin of attraction test along convergent trajectory (VD_ts = 100)
n_pts    = 50;
idx      = round(linspace(1, size(pop_d,1), n_pts));
test_pts = pop_d(idx, :);

CS    = coopEqs(end, :);
EU    = zeros(1, 6);
EU(1) = r / xi;

outcome = false(n_pts, 1);
parfor k = 1:n_pts
    ic_coop    = test_pts(k, :);
    ic_coop(3) = 0;   % H_D = 0
    ic_coop(4) = 0;   % H_M = 0
    ic_coop(6) = 0;   % V_D = 0
    [~, pop_coop] = ode45(@modelEqs, [0, t_run], ic_coop, opts, param);
    final      = pop_coop(end, :);
    dist_CS    = norm(final([2,5]) - CS([2,5]));
    dist_EU    = norm(final([2,5]) - EU([2,5]));
    outcome(k) = dist_CS < dist_EU;
end

VC_test = test_pts(:, 5);
trans   = find(diff(double(outcome)) < 0, 1);
if ~isempty(trans)
    VC_boundary = (VC_test(trans) + VC_test(trans+1)) / 2;
else
    VC_boundary = coopEqs(1, 5);
end

save('../data/viralExtinctionData.mat', ...
    'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi', 'rho', ...
    'param', 'coopEqs', 'pop0', ...
    't_c', 'pop_c', 't_d', 'pop_d', ...
    'VD_conv', 'conv_trajs', ...
    'outcome', 'VC_boundary', ...
    'ref_xl', 'ref_yl')
