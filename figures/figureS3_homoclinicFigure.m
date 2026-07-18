% figureS3_homoclinicFigure.m
% Evidence for homoclinic bifurcation in host rich & small burst size regime.
%
% Panel A: mixed limit cycles in (H_C, V_C) phase plane at several rho
%          values approaching rho_HC. E_C^U is marked.
% Panel B: minimum normalised distance from limit cycle to E_C^U vs rho,
%          approaching zero at rho_HC.

addpath('../functions')

load('../data/rg1_2CO_run.mat', 'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi')
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

params = [mu1, mu2, r, d, kappa, lambda, xi];

%% 1. Find cooperator equilibria (independent of rho)
coopEqs = cooperatorEquilibria(params);
% Sorted by V_C ascending. In host rich small burst regime, two equilibria
% exist. 
% E_C^U is the saddle that the limit cycle collides with at rho_HC.
ECU = coopEqs(1, :);   
stableCoop = coopEqs(2,:);

% Normalisation scale for distance: use E_C^U values where nonzero,
% else fall back to 1
normScale              = abs(ECU);
normScale(normScale < 1e-12) = 1;

%% 2. Mixed limit cycles at rho values approaching rho_HC
rhoHC  = 0.628637;
% Dense near rho_HC, sparse near 0.60. Ordered far→close so continuation
% starts from well-behaved dynamics and progresses toward the bifurcation.
                                                                                                     
 rhoDists = fliplr(logspace(log10(1e-6), log10(rhoHC - 0.6), 3));                                                             
  rhoVec = rhoHC - rhoDists;



limitCycles = cell(size(rhoVec));
minDists    = NaN(size(rhoVec));

% Arbitrary initial condition for first rho value
pop0 = stableCoop + [0,0,100,100,0,100];

for i = 1:length(rhoVec)
    rho   = rhoVec(i);
    param = [mu1, mu2, r, d, kappa, lambda, xi, rho];

    [~, pop] = ode45(@modelEqs, 0:0.1:50000, pop0, opts, param);

    HM = pop(:, 4);
    [~, peakIdx] = findpeaks(HM, 'MinPeakProminence', 0.1 * max(HM));

    if length(peakIdx) < 2
        warning('homoclinicFigure: fewer than 2 peaks at rho = %.4f', rho)
        continue
    end

    lc = pop(peakIdx(end-1):peakIdx(end), :);
    limitCycles{i} = lc;

    % Use a point on this limit cycle as starting point for the next rho
    pop0 = lc(1, :);

    % Minimum normalised distance from limit cycle to E_C^U
    diffs = (lc - ECU); %./ normScale;
    dists = sqrt(sum(diffs.^2, 2));
    minDists(i) = min(dists);
end




%% 3. Plot
colors = {'#0072B2', '#E69F00', '#D55E00', '#009E73', '#56B4E9', '#CC79A7'};
LW     = 3;
MS     = 8;

figure(1)
clf
set(gcf, 'Position', [50 50 900 400], 'renderer', 'Painters')
subplot(1, 2, 1)
hold on
for i = 1:length(rhoVec)
    if isempty(limitCycles{i}), continue; end
    lc = limitCycles{i};
    ls = '-'; if i == length(rhoVec), ls = ':'; end
    plot(lc(:, 2), lc(:, 5), 'Color', colors{i}, 'LineWidth', LW, ...
        'LineStyle', ls, 'DisplayName', sprintf('\\rho = %.4f', rhoVec(i)))
end
plot(ECU(2), ECU(5), 'kx', 'MarkerSize', MS, 'LineWidth', 2, ...
    'HandleVisibility', 'off')
text(ECU(2) * 0.9, ECU(5), 'E_C^U','FontSize', 11, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
xlabel('H_C', 'FontSize', 20)
ylabel('V_C', 'FontSize', 20)
legend('Location', 'best', 'FontSize', 20)
title('Limit cycle projections', 'FontSize', 20)

subplot(1, 2, 2)
validIdx = ~isnan(minDists);
hold on
for i = find(validIdx)
    plot(log10(rhoHC - rhoVec(i)), log(minDists(i)), 'o', ...
        'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', MS)
end
xlabel('log of distance from \rho_{HC}', 'FontSize', 20)
%xlabel('log_{10}(\rho_{HC} - \rho)', 'FontSize', 13)
ylabel('log of minimum distance to E_C^U', 'FontSize', 20)
title('Minimum distance to E_C^U', 'FontSize', 20)

%%
figure(2)
subplot(1,2,1)
hold on
for i = 1:length(rhoVec)
    if isempty(limitCycles{i}), continue; end
    lc = limitCycles{i};
    plot(lc(:, 3), lc(:, 6), 'Color', colors{i}, 'LineWidth', LW, ...
        'DisplayName', sprintf('\\rho = %.5f', rhoVec(i)))
end
plot(ECU(3), ECU(6), 'kx', 'MarkerSize', MS, 'LineWidth', 2, ...
    'DisplayName', 'E_C^U')
xlabel('H_D', 'FontSize', 13)
ylabel('V_D', 'FontSize', 13)
legend('Location', 'best', 'FontSize', 10)
title('Phase plane', 'FontSize', 13)
subplot(1,2,2)
hold on
for i = 1:length(rhoVec)
    if isempty(limitCycles{i}), continue; end
    lc = limitCycles{i};
    plot(lc(:, 3), lc(:, 4), 'Color', colors{i}, 'LineWidth', LW, ...
        'DisplayName', sprintf('\\rho = %.5f', rhoVec(i)))
end
plot(ECU(4), ECU(6), 'kx', 'MarkerSize', MS, 'LineWidth', 2, ...
    'DisplayName', 'E_C^U')
xlabel('H_M', 'FontSize', 13)
ylabel('V_D', 'FontSize', 13)
legend('Location', 'best', 'FontSize', 10)
title('Phase plane', 'FontSize', 13)
