---
editor_options: 
  markdown: 
    wrap: 72
---

# Environmental feedback maintains cooperation in viruses

Jaye Sudweeks and Christoph Hauert

26/09/2026

Code and data for reproducing the model, bifurcation analysis and
figures.

### Requirements

All code was written and run in MATLAB R2024b (version 24.2.0.2773142,
Update 2), with Symbolic Math Toolbox 24.2 and Signal Processing Toolbox
24.2.

## Contents

1.  [Repository structure](#1-repository-structure)
2.  [Repository contents](#2-repository-contents)
3.  [Reproducing the figures](#3-reproducing-the-figures)
4.  [Data](#4-data)
5.  [File naming conventions](#5-file-naming-conventions)

## 1. Repository structure

| Folder | Contents |
|----|----|
| `functions/` | functions shared across scripts |
| `dataGeneration/` | scripts that generate data and save results to data/ |
| `data/` | precomputed data as .mat files |
| `data/csv/` | data in CSV format, written by dataGeneration/exportDataToCSV.m |
| `figures/` | scripts that load data/ and generate the manuscript figures |

## 2. Repository contents

```         
README.md                 this file
.gitattributes            Git line-ending settings
.gitignore                files excluded from the repository

functions/                functions shared across scripts
  modelEqs.m                model equations
  bifurcationGenerator.m    solves symbolically all equilibria; for varying rho:
                            evaluates equilibria and classifies their stability
  cooperatorEquilibria.m    solves for cooperator-only equilibria 
  defectorOnlyEquilibria.m  solves for defector-only equilibria at a given rho
  mixedCollisions.m         determined rho values of the transcritical
                            bifurcations between E_M and E_C / E_D

dataGeneration/           scripts that generate data and save results to data/
  rl1Bifurcation.m          generates bifurcation data for the host poor regime
                            -> rl1_run.mat
  rg1_1CO_bifurcation.m     generates bifurcation data for the host rich & large
                            burst size regime -> rg1_1CO_run.mat
  rg1_2CO_bifurcation.m     generates bifurcation data for the host rich & small
                            burst size regime -> rg1_2CO_run.mat
  verifyBifurcations.m      computes and prints every bifurcation value
                            -> verifiedBifurcations.mat
  homoclinicRho.m           determines rho_HC by continuation and bisection
                            -> homoclinicRhoData.mat
  cooperatorLimitCycleStability.m
                            cooperator limit cycle and its stability to
                            defector invasion -> limitCycleStabilityData.mat
  defectorLimitCycle.m      defector limit cycle amplitudes
                            -> defectorLimitCycleAmps.mat
  mixedLimitCycle.m         mixed limit cycle amplitudes
                            -> mixedLimitCycles.mat
  mixedLCPhasePlaneData.m   one loop of the mixed limit cycle
                            -> mixedLCPhasePlaneData.mat
  viralExtinctionData.m     explores viral extinction in the host rich & small
                            burst size regime -> viralExtinctionData.mat
  exportDataToCSV.m         writes every .mat file in data/ to data/csv/

data/                     precomputed data as .mat files
  rl1_run.mat               data for bifurcation diagram: host poor regime
  rg1_1CO_run.mat           data for bifurcation diagram: host rich & large
                            burst size regime
  rg1_2CO_run.mat           data for bifurcation diagram: host rich & small
                            burst size regime
  verifiedBifurcations.mat  bifurcation values of rho for each regime
  homoclinicRhoData.mat     rho_HC, the homoclinic bifurcation
  limitCycleStabilityData.mat
                            cooperator limit cycle and its stability to defector
                            invasion
  defectorLimitCycleAmps.mat
                            defector limit cycle amplitudes over rho
  mixedLimitCycles.mat      mixed limit cycle amplitudes over rho
  mixedLCPhasePlaneData.mat
                            one loop of the mixed limit cycle at three rho
                            values
  viralExtinctionData.mat   viral extinction trajectories and basin of
                            attraction test

data/csv/                 data in CSV format
  parameters.csv            from rl1_run.mat, rg1_1CO_run.mat, rg1_2CO_run.mat
  bifurcationValues.csv     from verifiedBifurcations.mat
  equilibria_rl1.csv        from rl1_run.mat
  equilibria_rg1_1CO.csv    from rg1_1CO_run.mat
  equilibria_rg1_2CO.csv    from rg1_2CO_run.mat
  cooperatorLC_stability.csv
                            from limitCycleStabilityData.mat
  cooperatorLC_trajectory.csv
                            from limitCycleStabilityData.mat
  cooperatorLC_summary.csv  from limitCycleStabilityData.mat
  defectorLC_amplitudes.csv
                            from defectorLimitCycleAmps.mat
  mixedLC_amplitudes.csv    from mixedLimitCycles.mat
  mixedLC_phasePlane.csv    from mixedLCPhasePlaneData.mat
  viralExtinction_timeseries.csv
                            from viralExtinctionData.mat
  viralExtinction_convergence.csv
                            from viralExtinctionData.mat
  viralExtinction_basinTest.csv
                            from viralExtinctionData.mat
  viralExtinction_cooperatorEquilibria.csv
                            from viralExtinctionData.mat
  viralExtinction_scalars.csv
                            from viralExtinctionData.mat

figures/                  scripts that load data/ and generate the manuscript
                          figures
  figure1.m                 Figure 1: bifurcation diagrams for all three regimes
  figure2.m                 Figure 2: mixed limit cycle in the (V_C, V_D) plane
  figure3.m                 Figure 3: host extinction dynamics (host rich &
                            large burst size regime)
  figure4.m                 Figure 4: viral extinction dynamics (host rich &
                            small burst size regime)
  figureS3.m                Figure S3: evidence for the homoclinic bifurcation
```

## 3. Reproducing the figures

The `data/` folder contains precomputed data, so the scripts in
`figures/` can be run directly:

| Figure | Script | Requires (in `data/`) |
|----|----|----|
| 1 | figure1.m | verifiedBifurcations.mat, rl1_run.mat, rg1_1CO_run.mat, rg1_2CO_run.mat, limitCycleStabilityData.mat, defectorLimitCycleAmps.mat, mixedLimitCycles.mat, homoclinicRhoData.mat |
| 2 | figure2.m | mixedLCPhasePlaneData.mat |
| 3 | figure3.m | limitCycleStabilityData.mat |
| 4 | figure4.m | viralExtinctionData.mat |
| S3 | figureS3.m | rg1_2CO_run.mat |

## 4. Data

### 4.1 Regenerating the data

To regenerate `data/` from scratch, run the `dataGeneration/` scripts in
this order:

1.  `rl1Bifurcation.m`, `rg1_1CO_bifurcation.m`, `rg1_2CO_bifurcation.m`
2.  `cooperatorLimitCycleStability.m`, `defectorLimitCycle.m` (need
    rg1_1CO_run.mat); `mixedLimitCycle.m`, `viralExtinctionData.m` (need
    rg1_2CO_run.mat)
3.  `homoclinicRho.m`, `mixedLCPhasePlaneData.m` (need
    mixedLimitCycles.mat)
4.  `verifyBifurcations.m` (needs everything above except
    mixedLCPhasePlaneData.mat and viralExtinctionData.mat)
5.  `exportDataToCSV.m` (writes `data/csv/` from all of the .mat files
    above)

### 4.2 .mat files

#### data/rl1_run.mat

Generated by: `rl1Bifurcation.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| mu1 | scalar | single-infection adsorption rate | per density per hour |
| mu2 | scalar | double-infection adsorption rate | per density² per hour |
| r | scalar | host growth rate | per hour |
| d | scalar | lysis rate of infected hosts (1/latent period) | per hour |
| kappa | scalar | virus decay rate relative to d (viruses decay at rate d·kappa) | dimensionless |
| lambda | scalar | burst size | virions per infected host |
| xi | scalar | host density dependence (r / carrying capacity) | per density per hour |
| inputRhoVec | [1 2] | range of rho swept, [0.01, 1.2] | dimensionless |
| rhoVec | [1 2500] | rho grid used, refined near bifurcations | dimensionless |
| stableCell | cell [1 9] | stable equilibria. One cell per equilibrium branch, in the order E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only), M (mixed). Each cell is a matrix with one row per rho value at which that branch is biologically valid (real, non-negative) and stable. Columns: 1–6 state [H_U H_C H_D H_M V_C V_D]; 7 rho; 8–13 the six Jacobian eigenvalues (complex); 14 stability code (2 = stable). | columns 1–6 density; 7 dimensionless; 8–13 per hour; 14 – |
| unStableCell | cell [1 9] | unstable equilibria. Same layout as stableCell, with one row per rho value at which that branch is biologically valid and unstable; stability code 1 = unstable. | as stableCell |
| collisions | [2 1] | rho values of the transcritical bifurcations between E_M and E_C, and between E_M and E_D (from mixedCollisions.m) | dimensionless |
| rhoTUD | scalar | rho value of the transcritical bifurcation between E_U and E_D | dimensionless |
| bifVec | [1 3] | collisions and rhoTUD combined, sorted | dimensionless |

#### data/rg1_1CO_run.mat

Generated by: `rg1_1CO_bifurcation.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| mu1 | scalar | single-infection adsorption rate | per density per hour |
| mu2 | scalar | double-infection adsorption rate | per density² per hour |
| r | scalar | host growth rate | per hour |
| d | scalar | lysis rate of infected hosts (1/latent period) | per hour |
| kappa | scalar | virus decay rate relative to d (viruses decay at rate d·kappa) | dimensionless |
| lambda | scalar | burst size | virions per infected host |
| xi | scalar | host density dependence (r / carrying capacity) | per density per hour |
| inputRhoVec | [1 2] | range of rho swept, [0.01, 1.2] | dimensionless |
| rhoVec | [1 2500] | rho grid used, refined near bifurcations | dimensionless |
| stableCell | cell [1 9] | stable equilibria. One cell per equilibrium branch, in the order E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only), M (mixed). Each cell is a matrix with one row per rho value at which that branch is biologically valid (real, non-negative) and stable. Columns: 1–6 state [H_U H_C H_D H_M V_C V_D]; 7 rho; 8–13 the six Jacobian eigenvalues (complex); 14 stability code (2 = stable). | columns 1–6 density; 7 dimensionless; 8–13 per hour; 14 – |
| unStableCell | cell [1 9] | unstable equilibria. Same layout as stableCell, with one row per rho value at which that branch is biologically valid and unstable; stability code 1 = unstable. | as stableCell |
| collisions | [2 1] | rho values of the transcritical bifurcations between E_M and E_C, and between E_M and E_D (from mixedCollisions.m) | dimensionless |
| rhoTUD | scalar | rho value of the transcritical bifurcation between E_U and E_D | dimensionless |
| bifVec | [1 3] | collisions and rhoTUD combined, sorted | dimensionless |

#### data/rg1_2CO_run.mat

Generated by: `rg1_2CO_bifurcation.m`. Used in: Figures 1 and S3.

| Variable | Size | Description | Units |
|----|----|----|----|
| mu1 | scalar | single-infection adsorption rate | per density per hour |
| mu2 | scalar | double-infection adsorption rate | per density² per hour |
| r | scalar | host growth rate | per hour |
| d | scalar | lysis rate of infected hosts (1/latent period) | per hour |
| kappa | scalar | virus decay rate relative to d (viruses decay at rate d·kappa) | dimensionless |
| lambda | scalar | burst size | virions per infected host |
| xi | scalar | host density dependence (r / carrying capacity) | per density per hour |
| inputRhoVec | [1 2] | range of rho swept, [0.01, 2] | dimensionless |
| rhoVec | [1 2500] | rho grid used, refined near bifurcations | dimensionless |
| stableCell | cell [1 9] | stable equilibria. One cell per equilibrium branch, in the order E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only), M (mixed). Each cell is a matrix with one row per rho value at which that branch is biologically valid (real, non-negative) and stable. Columns: 1–6 state [H_U H_C H_D H_M V_C V_D]; 7 rho; 8–13 the six Jacobian eigenvalues (complex); 14 stability code (2 = stable). | columns 1–6 density; 7 dimensionless; 8–13 per hour; 14 – |
| unStableCell | cell [1 9] | unstable equilibria. Same layout as stableCell, with one row per rho value at which that branch is biologically valid and unstable; stability code 1 = unstable. | as stableCell |
| lowerThresh | scalar | lambda_S, the lower bound on burst size for the host rich & small burst size regime (lambda_S \< lambda \< lambda_T) | virions per infected host |
| upperThresh | scalar | lambda_T, the upper bound on burst size for the host rich & small burst size regime; lambda = 77 is the middle of this acceptable range | virions per infected host |
| ltNum, ltDen | scalars | numerator and denominator of lowerThresh (lowerThresh = ltNum/ltDen) | – |
| utNum, utDen | scalars | numerator and denominator of upperThresh (upperThresh = utNum/utDen) | – |

#### data/verifiedBifurcations.mat

Generated by: `verifyBifurcations.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| rl1_bifVec | [1 3] | bifurcation values of rho for the host poor regime, sorted | dimensionless |
| rg11_bifVec | [1 7] | bifurcation values of rho for the host rich & large burst size regime, sorted | dimensionless |
| rg12_bifVec | [1 5] | bifurcation values of rho for the host rich & small burst size regime, sorted | dimensionless |

#### data/homoclinicRhoData.mat

Generated by: `homoclinicRho.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| rhoHC | scalar | rho_HC, the homoclinic bifurcation for the host rich & small burst size regime | dimensionless |

#### data/limitCycleStabilityData.mat

Generated by: `cooperatorLimitCycleStability.m`. Used in: Figures 1 and
3.

| Variable | Size | Description | Units |
|----|----|----|----|
| limitCycle | [77 6] | one period of the cooperator limit cycle, from one V_C peak to the next, sampled every 0.1 hours | density |
| allMins, allMaxes | [1 6] | minimum and maximum of each state variable over the cycle | density |
| limCycleHostMin, limCycleHostMax | scalars | minimum and maximum of H_total over the cycle | density |
| CLCPoint | [1 6] | a point on the cycle (first row of limitCycle); initial condition for Figure 3 | density |
| rhoVec | [1 50] | rho values tested | dimensionless |
| stabVec | [1 50] | stability of the cycle to invasion by defectors at each rho (1 = unstable, 2 = stable) | – |
| rhoTransition | scalar | rho_LC, where the cycle loses stability (bisection) | dimensionless |

#### data/defectorLimitCycleAmps.mat

Generated by: `defectorLimitCycle.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| defectorLimitCycleRhoVec | [1 10] | rho values | dimensionless |
| DLCMax, DLCMin | [10 3] | maximum / minimum of H_U, H_D and V_D over the defector limit cycle at each rho | density |
| DLCMaxHostAmp, DLCMinHostAmp | [10 1] | maximum / minimum of H_total | density |
| DLCPoint | [10 6] | a point on the cycle at each rho | density |

#### data/mixedLimitCycles.mat

Generated by: `mixedLimitCycle.m`. Used in: Figure 1.

| Variable | Size | Description | Units |
|----|----|----|----|
| testingRho | [1 5] | rho values | dimensionless |
| storeMaxVec, storeMinVec | [5 6] | maximum / minimum of each state variable over the mixed limit cycle at each rho | density |
| storePopMaxVec, storePopMinVec | [5 2] | maximum / minimum of H_total (column 1) and V_total (column 2) | density |
| MlimCyclePoint | [5 6] | a point on the cycle at each rho | density |

#### data/mixedLCPhasePlaneData.mat

Generated by: `mixedLCPhasePlaneData.m`. Used in: Figure 2.

| Variable | Size | Description | Units |
|----|----|----|----|
| rhoVec | [1 3] | rho values 0.60, 0.61, 0.62 | dimensionless |
| VCdata, VDdata | cell [3 1] | one cell per rho, each a column vector of V_C or V_D along the loop | density |
| fCdata | cell [3 1] | as above, cooperator frequency f_C | dimensionless |
| threshold | scalar | mu1/mu2: total virus density above which double infection outpaces single infection (dashed line in Figure 2) | density |

#### data/viralExtinctionData.mat

Generated by: `viralExtinctionData.m`. Used in: Figure 4.

| Variable | Size | Description | Units |
|----|----|----|----|
| mu1 | scalar | single-infection adsorption rate | per density per hour |
| mu2 | scalar | double-infection adsorption rate | per density² per hour |
| r | scalar | host growth rate | per hour |
| d | scalar | lysis rate of infected hosts (1/latent period) | per hour |
| kappa | scalar | virus decay rate relative to d (viruses decay at rate d·kappa) | dimensionless |
| lambda | scalar | burst size | virions per infected host |
| xi | scalar | host density dependence (r / carrying capacity) | per density per hour |
| rho | scalar | cost to benefit ratio (the bifurcation parameter); here 0.65 | dimensionless |
| param | [1 8] | the eight values above, in the order [mu1 mu2 r d kappa lambda xi rho] used by modelEqs.m | as for each parameter above |
| coopEqs | [2 6] | cooperator-only equilibria: row 1 E_C\^U, row 2 E_C\^S | density |
| pop0 | [1 6] | initial condition, E_C\^S | density |
| t_c | [161 1] | time, cooperators alone from E_C\^S, t = 0 to 50 | hours |
| pop_c | [161 6] | state at times t_c | density |
| t_d | [7213 1] | time after V_D = 100 is added at t = 50, up to t = 1050 | hours |
| pop_d | [7213 6] | state at times t_d | density |
| VD_conv | [1 5] | initial V_D values added to E_C\^S: 1, 10, 100, 1000, 10000 | density |
| conv_trajs | cell [5 1] | one trajectory per VD_conv value, each [n 6], integrated from t = 0 to 1000 hours (time points not saved) | density |
| outcome | [50 1] | basin test for 50 points along pop_d: true = converges to E_C\^S, false = to E_U (see `viralExtinction_basinTest.csv`) | – |
| VC_boundary | scalar | V_C value at which the basin test outcome changes (dotted line in Figure 4) | density |
| ref_xl, ref_yl | [1 2] | reference axis limits for V_C and V_D; not used by the current figure scripts | density |

### 4.3 CSV files

Data in csv format. Where a .mat file nests data (cell arrays of
branches or trajectories), the CSV is in "long" format: one row per
observation, with identifier columns (for example `branch`, `rho`,
`VD0`).

#### data/csv/parameters.csv

Source: `rl1_run.mat`, `rg1_1CO_run.mat`, `rg1_2CO_run.mat`. One row per
parameter per regime.

| Column | Description | Units |
|----|----|----|
| regime | `rl1` (host poor regime), `rg1_1CO` (host rich & large burst size regime) or `rg1_2CO` (host rich & small burst size regime) | – |
| parameter | mu1, mu2, r, d, kappa, lambda, xi; for the host rich & small burst size regime also lowerThresh, upperThresh, ltNum, ltDen, utNum, utDen (see `rg1_2CO_run.mat` above) | – |
| value | parameter value | depends on the parameter; see the table below |
| description | short description of the parameter | – |

Units of `value`, by parameter:

| parameter | Description | Units |
|----|----|----|
| mu1 | single-infection adsorption rate | per density per hour |
| mu2 | double-infection adsorption rate | per density² per hour |
| r | host growth rate | per hour |
| d | lysis rate of infected hosts (1/latent period) | per hour |
| kappa | virus decay rate relative to d (viruses decay at rate d·kappa) | dimensionless |
| lambda | burst size | virions per infected host |
| xi | host density dependence (r / carrying capacity) | per density per hour |
| lowerThresh, upperThresh | lambda_S and lambda_T, the lower and upper bounds on burst size for the host rich & small burst size regime (lambda_S \< lambda \< lambda_T) | virions per infected host |
| ltNum, ltDen, utNum, utDen | numerators and denominators of lowerThresh and upperThresh | – |

#### data/csv/bifurcationValues.csv

Source: `verifiedBifurcations.mat`. One row per bifurcation.

| Column | Description | Units |
|----|----|----|
| regime | `rl1` (host poor regime), `rg1_1CO` (host rich & large burst size regime) or `rg1_2CO` (host rich & small burst size regime) | – |
| bifurcation | name of the bifurcation (see the list below) | – |
| rho | value of rho at the bifurcation | dimensionless |
| method | how the value computed | – |

Values of `bifurcation`:

| bifurcation | Meaning |
|----|----|
| rho_SD | saddle-node bifurcation of the defector-only equilibria |
| rho_TUD, rho_TDU | transcritical bifurcation of E_U and E_D |
| rho_TMC (rho_TMC1, rho_TMC2) | transcritical bifurcation of E_M and E_C (the host rich & small burst size regime has two) |
| rho_TMD | transcritical bifurcation of E_M and E_D |
| rho_M | E_M gains stability |
| rho_LC | cooperator limit cycle loses stability to defector invasion |
| rho_HD | Hopf bifurcation of E_D |
| rho_HM | Hopf bifurcation of E_M |
| rho_HC | homoclinic bifurcation (end of the mixed limit cycle) |

#### data/csv/equilibria_rl1.csv

Source: stableCell and unStableCell in `rl1_run.mat`. One row per
equilibrium per rho value. Branches that are not biologically relevant
have no entries.

| Column | Description | Units |
|----|----|----|
| branch | E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only) or M (mixed) | – |
| stability | `stable` (all eigenvalues have negative real part) or `unstable` | – |
| rho | value of rho | dimensionless |
| H_U, H_C, H_D, H_M | equilibrium densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | equilibrium densities of cooperator and defector viruses, respectively | density |
| eig1_re, eig1_im, …, eig6_re, eig6_im | real and imaginary parts of the six Jacobian eigenvalues | per hour |

#### data/csv/equilibria_rg1_1CO.csv

Source: stableCell and unStableCell in `rg1_1CO_run.mat`. One row per
equilibrium per rho value. Branches that are not biologically relevant
have no entries.

| Column | Description | Units |
|----|----|----|
| branch | E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only) or M (mixed) | – |
| stability | `stable` (all eigenvalues have negative real part) or `unstable` | – |
| rho | value of rho | dimensionless |
| H_U, H_C, H_D, H_M | equilibrium densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | equilibrium densities of cooperator and defector viruses, respectively | density |
| eig1_re, eig1_im, …, eig6_re, eig6_im | real and imaginary parts of the six Jacobian eigenvalues | per hour |

#### data/csv/equilibria_rg1_2CO.csv

Source: stableCell and unStableCell in `rg1_2CO_run.mat`. One row per
equilibrium per rho value. Branches that are not biologically relevant
have no entries.

| Column | Description | Units |
|----|----|----|
| branch | E (extinct), U (uninfected), C1–C3 (cooperator-only), D1–D3 (defector-only) or M (mixed) | – |
| stability | `stable` (all eigenvalues have negative real part) or `unstable` | – |
| rho | value of rho | dimensionless |
| H_U, H_C, H_D, H_M | equilibrium densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | equilibrium densities of cooperator and defector viruses, respectively | density |
| eig1_re, eig1_im, …, eig6_re, eig6_im | real and imaginary parts of the six Jacobian eigenvalues | per hour |

#### data/csv/cooperatorLC_stability.csv

Source: `limitCycleStabilityData.mat`.

| Column | Description | Units |
|----|----|----|
| rho | value of rho | dimensionless |
| stability | `stable` or `unstable`: whether the cooperator limit cycle resists invasion by a small defector population | – |

#### data/csv/cooperatorLC_trajectory.csv

Source: `limitCycleStabilityData.mat`.

| Column | Description | Units |
|----|----|----|
| t | time from the start of the period, step 0.1 | hours |
| H_U, H_C, H_D, H_M | densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | densities of cooperator and defector viruses, respectively | density |

#### data/csv/cooperatorLC_summary.csv

Source: `limitCycleStabilityData.mat`.

| Column   | Description                                              | Units   |
|----------|----------------------------------------------------------|---------|
| variable | H_U, H_C, H_D, H_M, V_C, V_D or H_total                  | –       |
| min      | minimum of that variable over the cooperator limit cycle | density |
| max      | maximum of that variable over the cycle                  | density |

#### data/csv/defectorLC_amplitudes.csv

Source: `defectorLimitCycleAmps.mat`. One row per rho value.

| Column | Description | Units |
|----|----|----|
| rho | value of rho | dimensionless |
| H_U_max, H_D_max, H_total_max | maxima of H_U, H_D and H_total over the defector limit cycle | density |
| V_D_max | maximum of V_D over the cycle | density |
| H_U_min, H_D_min, H_total_min | minima of H_U, H_D and H_total over the cycle | density |
| V_D_min | minimum of V_D over the cycle | density |
| H_U_point, H_C_point, H_D_point, H_M_point | host densities at a point on the cycle | density |
| V_C_point, V_D_point | virus densities at the same point | density |

#### data/csv/mixedLC_amplitudes.csv

Source: `mixedLimitCycles.mat`. One row per rho value.

| Column | Description | Units |
|----|----|----|
| rho | value of rho | dimensionless |
| H_U_max, H_C_max, H_D_max, H_M_max, H_total_max | maxima of the host densities over the mixed limit cycle | density |
| V_C_max, V_D_max, V_total_max | maxima of the virus densities over the cycle | density |
| H_U_min, H_C_min, H_D_min, H_M_min, H_total_min | minima of the host densities over the cycle | density |
| V_C_min, V_D_min, V_total_min | minima of the virus densities over the cycle | density |
| H_U_point, H_C_point, H_D_point, H_M_point | host densities at a point on the cycle | density |
| V_C_point, V_D_point | virus densities at the same point | density |

#### data/csv/mixedLC_phasePlane.csv

Source: `mixedLCPhasePlaneData.mat`.

| Column | Description | Units |
|----|----|----|
| rho | value of rho (0.60, 0.61 or 0.62) | dimensionless |
| V_C, V_D | densities of cooperator and defector viruses along the loop | density |
| f_C | cooperator frequency, V_C / (V_C + V_D) | dimensionless |

#### data/csv/viralExtinction_timeseries.csv

Source: `viralExtinctionData.mat.`

| Column | Description | Units |
|----|----|----|
| phase | `cooperators_only` (cooperators alone from E_C\^S, t = 0–50) or `defectors_introduced` (after V_D = 100 is added at t = 50, up to t = 1050) | – |
| t | time | hours |
| H_U, H_C, H_D, H_M | densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | densities of cooperator and defector viruses, respectively | density |

#### data/csv/viralExtinction_convergence.csv

Source: `viralExtinctionData.mat`.

| Column | Description | Units |
|----|----|----|
| VD0 | initial V_D added: 1, 10, 100, 1000 or 10000 | density |
| point_index | position of the row along its trajectory (1, 2, …) | – |
| H_U, H_C, H_D, H_M | densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | densities of cooperator and defector viruses, respectively | density |

#### data/csv/viralExtinction_basinTest.csv

Source: `viralExtinctionData.mat`.

| Column | Description | Units |
|----|----|----|
| t | time of the point on the `defectors_introduced` trajectory | hours |
| H_U, H_C, H_D, H_M | densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts at that point, before defectors were removed | density |
| V_C, V_D | densities of cooperator and defector viruses at that point, before defectors were removed | density |
| converges_to | `E_C_S` or `E_U`: the equilibrium the cooperator-only system reached (the closer of the two in H_C and V_C) | – |

#### data/csv/viralExtinction_cooperatorEquilibria.csv

Source: `viralExtinctionData.mat`.

| Column | Description | Units |
|----|----|----|
| equilibrium | `E_C_U` (unstable) or `E_C_S` (stable) cooperator-only equilibrium | – |
| H_U, H_C, H_D, H_M | densities of uninfected, cooperator-infected, defector-infected and mixed-infected hosts, respectively | density |
| V_C, V_D | densities of cooperator and defector viruses, respectively | density |

#### data/csv/viralExtinction_scalars.csv

Source: `viralExtinctionData.mat`.

| Column | Description | Units |
|----|----|----|
| quantity | `rho` or `VC_boundary` | – |
| value | rho = 0.65; VC_boundary, the V_C value at which the basin test outcome changes (dotted line in Figure 4) | rho dimensionless; VC_boundary density |

## 5. File naming conventions

functions/

-   files are named for what they compute

dataGeneraion/

-   file names indicate the data generated (for instance,
    rl1Bifurcation.m generates bifurcation data)

-   Some file names reference regimes:

    -   rl1: host poor regime

    -   rl1_1C0: host poor & large burst size regime

    -   rl1_2C0: host poor & small burst size regime

data/

-   file names generally indicate the data contained (for
    instance,verifiedBifurcations.mat lists bifurcations from each
    regime)

-   Some file names reference regimes:

    -   rl1: host poor regime

    -   rl1_1C0: host rich & large burst size regime

    -   rl1_2C0: host rich & small burst size regime

data/csv

-   file names indicate the data contained within

figures/

-   files are named for the figure that they recreate
