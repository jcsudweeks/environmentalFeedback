README
17/07/2026

Environmental feedback maintains cooperation in viruses
by: Jaye Sudweeks and Christoph Hauert

Code for reproducing the model, bifurcation analysis, and figures in
the above manuscript.

REQUIREMENTS

MATLAB with the Symbolic Math Toolbox (used for equilibrium derivation
and Jacobians in functions/bifurcationGenerator.m and
dataGeneration/verifyBifurcations.m) and the Signal Processing Toolbox
(findpeaks, used to extract limit cycles). The Parallel Computing
Toolbox speeds up several scripts (parfor) but is not required --
without it, parfor loops run serially.

REPOSITORY STRUCTURE

functions/       core model and equilibrium-finding functions, shared
                  across the scripts below.
dataGeneration/   scripts that run the model and save results to data/.
figures/          scripts that load from data/ and produce the
                  manuscript figures.
data/             pre-computed .mat files, so figures can be
                  reproduced without rerunning the (sometimes slow)
                  data generation.

Scripts in dataGeneration/ and figures/ reference functions/ and data/
with relative paths (addpath('../functions'), load('../data/...')), so
run them with the script's own folder as MATLAB's current folder.

REPRODUCING THE FIGURES

The data/ folder already contains the output of every script in
dataGeneration/, so the scripts in figures/ can be run directly:

  Figure 1   figure1_combinedBifurcation.m
             requires: verifiedBifurcations.mat, rl1_run.mat,
             rg1_1CO_run.mat, rg1_2CO_run.mat, limitCycleStabilityData.mat,
             defectorLimitCycleAmps.mat, mixedLimitCycles.mat,
             homoclinicRhoData.mat

  Figure 2   figure2_mixedLCPhasePlane.m
             requires: mixedLCPhasePlaneData.mat

  Figure 3   figure3_HostExtinction.m
             requires: limitCycleStabilityData.mat

  Figure 4   figure4_viralExtinction.m
             requires: viralExtinctionData.mat

  Figure S3  figureS3_homoclinicFigure.m
             requires: rg1_2CO_run.mat

REGENERATING THE DATA

To regenerate data/ from scratch, run the dataGeneration/ scripts in
this order (each depends on the .mat output of the ones before it):

  1. rl1Bifurcation.m, rg1_1CO_bifurcation.m, rg1_2CO_bifurcation.m
  2. cooperatorLimitCycleStability.m, defectorLimitCycle.m
     (need rg1_1CO_run.mat)
     mixedLimitCycle.m, viralExtinctionData.m
     (need rg1_2CO_run.mat)
  3. homoclinicRho.m, mixedLCPhasePlaneData.m
     (need mixedLimitCycles.mat)
  4. verifyBifurcations.m
     (needs everything above except mixedLCPhasePlaneData.mat and
     viralExtinctionData.mat)
