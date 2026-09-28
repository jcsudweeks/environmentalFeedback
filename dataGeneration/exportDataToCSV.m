% exportDataToCSV.m
% Exports the contents of every .mat file in ../data to plain-text CSV
% files in ../data/csv, so the data can be read without MATLAB. The .mat
% files are left unchanged and remain the inputs to all figure scripts.
%
% Nested MATLAB structures (cell arrays of branches or trajectories) are
% written in "long" format: one row per observation, with identifier
% columns (e.g. branch, rho, VD0) in place of the nesting. See README.md
% for a column-by-column description of each file.
%
% Run from the dataGeneration folder.

inDir  = '../data';
outDir = '../data/csv';
if ~exist(outDir, 'dir'), mkdir(outDir); end

stateNames = {'H_U', 'H_C', 'H_D', 'H_M', 'V_C', 'V_D'};
regimes    = {'rl1', 'rg1_1CO', 'rg1_2CO'};

%% 1. Parameters for each regime
parNames = {'mu1', 'mu2', 'r', 'd', 'kappa', 'lambda', 'xi'};
parDesc  = {'single-infection adsorption rate', ...
            'double-infection adsorption rate', ...
            'host growth rate', ...
            'lysis rate of infected hosts (1/latent period)', ...
            'virus decay rate relative to d (viruses decay at rate d*kappa)', ...
            'burst size', ...
            'host density dependence (r / carrying capacity)'};

regimeCol = {}; nameCol = {}; valCol = []; descCol = {};
for k = 1:numel(regimes)
    S = load(fullfile(inDir, [regimes{k} '_run.mat']));
    for j = 1:numel(parNames)
        regimeCol{end+1,1} = regimes{k};
        nameCol{end+1,1}   = parNames{j};
        valCol(end+1,1)    = S.(parNames{j});
        descCol{end+1,1}   = parDesc{j};
    end
end
% Burst-size bounds defining the rg1_2CO regime (see rg1_2CO_bifurcation.m)
S = load(fullfile(inDir, 'rg1_2CO_run.mat'));
extra = {'lowerThresh', S.lowerThresh, 'lambda_S: lower bound on burst size for the host rich & small burst size regime (= ltNum/ltDen)'; ...
         'upperThresh', S.upperThresh, 'lambda_T: upper bound on burst size for the host rich & small burst size regime (= utNum/utDen)'; ...
         'ltNum',       S.ltNum,       'numerator of lowerThresh'; ...
         'ltDen',       S.ltDen,       'denominator of lowerThresh'; ...
         'utNum',       S.utNum,       'numerator of upperThresh'; ...
         'utDen',       S.utDen,       'denominator of upperThresh'};
for j = 1:size(extra, 1)
    regimeCol{end+1,1} = 'rg1_2CO';
    nameCol{end+1,1}   = extra{j,1};
    valCol(end+1,1)    = extra{j,2};
    descCol{end+1,1}   = extra{j,3};
end
writetable(table(regimeCol, nameCol, valCol, descCol, ...
    'VariableNames', {'regime', 'parameter', 'value', 'description'}), ...
    fullfile(outDir, 'parameters.csv'));

%% 2. Bifurcation values
% verifiedBifurcations.mat stores sorted values without labels; labels
% below follow verifyBifurcations.m and were matched against the
% analytical expressions.
V = load(fullfile(inDir, 'verifiedBifurcations.mat'));
bif = { ...
 'rl1',     V.rl1_bifVec,  {'rho_TUD', 'rho_TMC', 'rho_TMD'}, ...
            {'analytical', 'mixedCollisions (cubic)', 'mixedCollisions (quadratic)'}; ...
 'rg1_1CO', V.rg11_bifVec, {'rho_SD', 'rho_TMC', 'rho_M', 'rho_TDU', 'rho_LC', 'rho_TMD', 'rho_HD'}, ...
            {'analytical', 'mixedCollisions (cubic)', 'eigenvalue bisection at E_M', 'analytical', ...
             'cooperatorLimitCycleStability.m', 'mixedCollisions (quadratic)', 'eigenvalue bisection at E_D^S'}; ...
 'rg1_2CO', V.rg12_bifVec, {'rho_TMC1', 'rho_HM', 'rho_HC', 'rho_TMC2', 'rho_SD'}, ...
            {'mixedCollisions (cubic)', 'eigenvalue bisection at E_M', 'homoclinicRho.m', ...
             'mixedCollisions (cubic)', 'analytical'}};

regimeCol = {}; labelCol = {}; rhoCol = []; methodCol = {};
for k = 1:size(bif, 1)
    vals = bif{k,2};
    assert(numel(vals) == numel(bif{k,3}), 'Label count mismatch for %s', bif{k,1});
    for j = 1:numel(vals)
        regimeCol{end+1,1} = bif{k,1};
        labelCol{end+1,1}  = bif{k,3}{j};
        rhoCol(end+1,1)    = vals(j);
        methodCol{end+1,1} = bif{k,4}{j};
    end
end
writetable(table(regimeCol, labelCol, rhoCol, methodCol, ...
    'VariableNames', {'regime', 'bifurcation', 'rho', 'method'}), ...
    fullfile(outDir, 'bifurcationValues.csv'));

%% 3. Equilibrium branches (bifurcation diagrams)
% stableCell/unStableCell: 1x9 cells, one per branch; rows are
% [H_U..V_D, rho, eigenvalues(1:6), stability (1 = unstable, 2 = stable)]
branchNames = {'E', 'U', 'C1', 'C2', 'C3', 'D1', 'D2', 'D3', 'M'};
eigNames = {};
for j = 1:6
    eigNames = [eigNames, {sprintf('eig%d_re', j), sprintf('eig%d_im', j)}];
end

for k = 1:numel(regimes)
    S = load(fullfile(inDir, [regimes{k} '_run.mat']));
    rows = []; branchCol = {}; stabCol = {};
    for b = 1:9
        M = [S.stableCell{b}; S.unStableCell{b}];
        if isempty(M), continue, end
        % Equilibria were filtered to real solutions before eigenvalues
        % were appended, so any imaginary part here is round-off
        assert(max(abs(imag(M(:, [1:7, 14]))), [], 'all') < 1e-8, ...
            'Non-negligible imaginary part in %s branch %s', regimes{k}, branchNames{b});
        M = sortrows(M, 7);
        E = M(:, 8:13);
        eigCols = zeros(size(E, 1), 12);
        eigCols(:, 1:2:end) = real(E);
        eigCols(:, 2:2:end) = imag(E);
        rows = [rows; real(M(:, 7)), real(M(:, 1:6)), eigCols];
        branchCol = [branchCol; repmat(branchNames(b), size(M, 1), 1)];
        stab = repmat({'unstable'}, size(M, 1), 1);
        stab(real(M(:, 14)) == 2) = {'stable'};
        stabCol = [stabCol; stab];
    end
    T = [table(branchCol, stabCol, 'VariableNames', {'branch', 'stability'}), ...
         array2table(rows, 'VariableNames', [{'rho'}, stateNames, eigNames])];
    writetable(T, fullfile(outDir, ['equilibria_' regimes{k} '.csv']));
end

%% 4. Cooperator limit cycle (rg1_1CO)
L = load(fullfile(inDir, 'limitCycleStabilityData.mat'));
stab = repmat({'unstable'}, numel(L.rhoVec), 1);
stab(L.stabVec == 2) = {'stable'};
writetable(table(L.rhoVec(:), stab, 'VariableNames', {'rho', 'stability'}), ...
    fullfile(outDir, 'cooperatorLC_stability.csv'));

% One period, sampled on the fixed output grid (dt = 0.1) used by
% cooperatorLimitCycleStability.m; t is measured from the start of the period
tCycle = (0:size(L.limitCycle, 1)-1)' * 0.1;
writetable([table(tCycle, 'VariableNames', {'t'}), ...
            array2table(L.limitCycle, 'VariableNames', stateNames)], ...
    fullfile(outDir, 'cooperatorLC_trajectory.csv'));

writetable(table([stateNames, {'H_total'}]', ...
                 [L.allMins, L.limCycleHostMin]', [L.allMaxes, L.limCycleHostMax]', ...
                 'VariableNames', {'variable', 'min', 'max'}), ...
    fullfile(outDir, 'cooperatorLC_summary.csv'));

%% 5. Defector limit cycle amplitudes (rg1_1CO)
D = load(fullfile(inDir, 'defectorLimitCycleAmps.mat'));
T = [table(D.defectorLimitCycleRhoVec(:), 'VariableNames', {'rho'}), ...
     array2table(D.DLCMax, 'VariableNames', {'H_U_max', 'H_D_max', 'V_D_max'}), ...
     array2table(D.DLCMin, 'VariableNames', {'H_U_min', 'H_D_min', 'V_D_min'}), ...
     table(D.DLCMaxHostAmp, D.DLCMinHostAmp, 'VariableNames', {'H_total_max', 'H_total_min'}), ...
     array2table(D.DLCPoint, 'VariableNames', strcat(stateNames, '_point'))];
writetable(T, fullfile(outDir, 'defectorLC_amplitudes.csv'));

%% 6. Mixed limit cycle amplitudes (rg1_2CO)
M = load(fullfile(inDir, 'mixedLimitCycles.mat'));
T = [table(M.testingRho(:), 'VariableNames', {'rho'}), ...
     array2table(M.storeMaxVec, 'VariableNames', strcat(stateNames, '_max')), ...
     array2table(M.storeMinVec, 'VariableNames', strcat(stateNames, '_min')), ...
     array2table(M.storePopMaxVec, 'VariableNames', {'H_total_max', 'V_total_max'}), ...
     array2table(M.storePopMinVec, 'VariableNames', {'H_total_min', 'V_total_min'}), ...
     array2table(M.MlimCyclePoint, 'VariableNames', strcat(stateNames, '_point'))];
writetable(T, fullfile(outDir, 'mixedLC_amplitudes.csv'));

%% 7. Mixed limit cycle phase plane (rg1_2CO)
% One full loop of the limit cycle per rho, in order around the loop.
P = load(fullfile(inDir, 'mixedLCPhasePlaneData.mat'));
rows = [];
for k = 1:numel(P.rhoVec)
    rows = [rows; repmat(P.rhoVec(k), numel(P.VCdata{k}), 1), ...
            P.VCdata{k}, P.VDdata{k}, P.fCdata{k}];
end
writetable(array2table(rows, 'VariableNames', {'rho', 'V_C', 'V_D', 'f_C'}), ...
    fullfile(outDir, 'mixedLC_phasePlane.csv'));

%% 8. Viral extinction (rg1_2CO, rho = 0.65)
X = load(fullfile(inDir, 'viralExtinctionData.mat'));

% Time series: cooperators alone, then defectors introduced at t = 50
phase = [repmat({'cooperators_only'}, numel(X.t_c), 1); ...
         repmat({'defectors_introduced'}, numel(X.t_d), 1)];
writetable([table(phase, [X.t_c; X.t_d], 'VariableNames', {'phase', 't'}), ...
            array2table([X.pop_c; X.pop_d], 'VariableNames', stateNames)], ...
    fullfile(outDir, 'viralExtinction_timeseries.csv'));

% Convergence trajectories from E_C^S with different initial V_D (no time
% vector saved; each runs over t = 0 to 1000)
rows = [];
for k = 1:numel(X.VD_conv)
    traj = X.conv_trajs{k};
    rows = [rows; repmat(X.VD_conv(k), size(traj, 1), 1), (1:size(traj, 1))', traj];
end
writetable(array2table(rows, 'VariableNames', [{'VD0', 'point_index'}, stateNames]), ...
    fullfile(outDir, 'viralExtinction_convergence.csv'));

% Basin test: the 50 points along pop_d tested in viralExtinctionData.m
idx     = round(linspace(1, size(X.pop_d, 1), numel(X.outcome)));
outcome = repmat({'E_U'}, numel(X.outcome), 1);
outcome(X.outcome) = {'E_C_S'};
writetable([table(X.t_d(idx), 'VariableNames', {'t'}), ...
            array2table(X.pop_d(idx, :), 'VariableNames', stateNames), ...
            table(outcome, 'VariableNames', {'converges_to'})], ...
    fullfile(outDir, 'viralExtinction_basinTest.csv'));

% Cooperator-only equilibria and the basin boundary
writetable([table({'E_C_U'; 'E_C_S'}, 'VariableNames', {'equilibrium'}), ...
            array2table(X.coopEqs, 'VariableNames', stateNames)], ...
    fullfile(outDir, 'viralExtinction_cooperatorEquilibria.csv'));
writetable(table({'rho'; 'VC_boundary'}, [X.rho; X.VC_boundary], ...
                 'VariableNames', {'quantity', 'value'}), ...
    fullfile(outDir, 'viralExtinction_scalars.csv'));

fprintf('CSV files written to %s\n', outDir)
