% Checks stability of the cooperator limit cycle over a range of rho values.
% This is specifically in the host rich, large burst size regime. 
% Strategy:
%   1. Find the reference cooperator limit cycle accurately. Because rho does
%      not appear in the cooperator equations, this cycle is the same for all
%      rho. We find it by running a long trajectory to convergence and then
%      extracting exactly one period using peak detection on V_C.
%
%   2. For each rho, perturb off the limit cycle in all directions by a small
%      fixed amount and run a trajectory. Stability is assessed by computing
%      the minimum phase-space distance from the final point to the reference
%      limit cycle (each variable normalised by its range on the cycle).
%
% Output variables (matches previous version for downstream compatibility):
%   allMins, allMaxes       - min/max of each variable over the limit cycle (1x6)
%   limCycleHostMax/Min     - max/min of total host density over the limit cycle
%   stabVec                 - 2 = stable, 1 = unstable (length = rhoVec)
%   rhoVec                  - rho values tested

addpath('../functions')

%% 1. Find the reference cooperator limit cycle

load('../data/rg1_1CO_run.mat')
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% rho does not appear in cooperator equations — any value gives the same cycle
param0 = [mu1, mu2, r, d, kappa, lambda, xi, 0.5];
pop0   = [10, 10, 0, 0, 10, 0];    % start without defectors

tspan_long = 0:0.1:20000;
[~, popLong] = ode45(@modelEqs, tspan_long, pop0, opts, param0);

% Find peaks in V_C (column 5) to identify complete oscillation periods
VC = popLong(:, 5);
[~, peakIdx] = findpeaks(VC, 'MinPeakProminence', 0.1 * max(VC));

if length(peakIdx) < 2
    error('cooperatorLimitCycleStability: fewer than two peaks found in V_C — check that the parameter set admits a limit cycle.')
end

% Extract one clean period between the last two consecutive peaks
cycleStart = peakIdx(end-1);
cycleEnd   = peakIdx(end);
limitCycle = popLong(cycleStart:cycleEnd, :);

% Reference statistics used by downstream files
allMaxes        = max(limitCycle);
allMins         = min(limitCycle);
limCycleHostMax = max(sum(limitCycle(:, 1:4), 2));
limCycleHostMin = min(sum(limitCycle(:, 1:4), 2));

% Normalise limit cycle for distance computations.
% Variables with zero range on the cooperator cycle (H_D, H_M, V_D) get
% normScale = 1 so any nonzero value in those variables contributes naturally.
normScale              = allMaxes - allMins;
normScale(normScale < 1e-12) = 1;
lcNorm = (limitCycle - allMins) ./ normScale;

perturbFrac = 0.01;   % 1% relative perturbation for cooperator variables

% Four starting points spread evenly around the limit cycle (0%, 25%, 50%, 75%)
testIdx   = round(linspace(1, size(limitCycle, 1), 5));
testConds = limitCycle(testIdx(1:4), :);

% Stability tolerance: distance in normalised phase space below which we
% consider the trajectory to have returned to the limit cycle
stabTol = 0.05;

%% 2. Check stability over a range of rho values

% Upper bound kept below rho_T_MD = 0.4389 to avoid confounding that bifurcation
rhoVec  = linspace(0.01, 0.437, 50);
stabVec = zeros(size(rhoVec));

parfor i = 1:length(rhoVec)
    rho   = rhoVec(i);
    param = [mu1, mu2, r, d, kappa, lambda, xi, rho];

    % Test from each of the four starting points. Classify as unstable if
    % any starting point fails to return to the limit cycle.
    tspan_check = 0:0.1:10000;
    isStable = true;

    for k = 1:size(testConds, 1)
        pop_init            = testConds(k, :) .* (1 + perturbFrac);
        pop_init([3, 4, 6]) = 0.01;

        [~, popTraj] = ode45(@modelEqs, tspan_check, pop_init, opts, param);

        finalPt = (popTraj(end, :) - allMins) ./ normScale;
        dists   = sqrt(sum((lcNorm - finalPt).^2, 2));

        if min(dists) >= stabTol
            isStable = false;
            break;
        end
    end

    % Stable = 2, unstable = 1 (convention matches previous version)
    stabVec(i) = isStable + 1;
end


%% 3. Bisect to find the stability transition precisely

rhoLo = rhoVec(find(stabVec == 2, 1, 'last'));   % last stable rho
rhoHi = rhoVec(find(stabVec == 1, 1, 'first'));  % first unstable rho

bisectTol = 1e-6;

while (rhoHi - rhoLo) > bisectTol
    rhoMid = (rhoLo + rhoHi) / 2;
    param  = [mu1, mu2, r, d, kappa, lambda, xi, rhoMid];

    midStable = true;
    for k = 1:size(testConds, 1)
        pop_init            = testConds(k, :) .* (1 + perturbFrac);
        pop_init([3, 4, 6]) = 0.01;

        [~, popTraj] = ode45(@modelEqs, 0:0.1:10000, pop_init, opts, param);

        finalPt = (popTraj(end, :) - allMins) ./ normScale;
        dists   = sqrt(sum((lcNorm - finalPt).^2, 2));

        if min(dists) >= stabTol
            midStable = false;
            break;
        end
    end

    if midStable
        rhoLo = rhoMid;   % midpoint is stable, transition is higher
    else
        rhoHi = rhoMid;   % midpoint is unstable, transition is lower
    end
end

rhoTransition = (rhoLo + rhoHi) / 2;
fprintf('Limit cycle stability transition: rho = %.8f\n', rhoTransition)

%% Plot
figure(1)
plot(rhoVec, stabVec)
xline(rhoTransition, 'r--', sprintf('\\rho_LC = %.4f', rhoTransition))
xlabel('\rho')
ylabel('stability (2=stable, 1=unstable)')

%% Save output
CLCPoint = limitCycle(1,:);   % a point on the reference limit cycle, used to
                               % seed trajectories elsewhere (e.g. figure3.m)
save("../data/limitCycleStabilityData.mat", "allMins", "allMaxes", "stabVec", "rhoVec", ...
    "limCycleHostMax", "limCycleHostMin", "rhoTransition","limitCycle","CLCPoint")
