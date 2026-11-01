FINAL CORE-ABLATION EVIDENCE PACKAGE
====================================

This folder aggregates the completed simulations; it does not rerun or retune controllers.

Successful strictly paired Monte Carlo trials used: 100

Primary comparison groups:
  1. Full Proposed: K_mis=1; B_a, K_a, G_a adaptive.
  2. No Mismatch Compensation: K_mis=0; B_a, K_a, G_a adaptive.
  3. Simplified Adaptation (B/K): K_mis=1; B_a and K_a adaptive, G_a fixed.

Primary metrics (lower is better): tracking RMSE, jerk RMS, raw interaction-residual RMS.
Pairwise difference convention: Comparator - Full. Positive values mean Full is lower.
Statistics: two-sided paired t test, t-based 95% CI, Cohen d_z, and Holm adjustment across 18 endpoint tests.

Mismatch K_mis={0,0.25,0.5,0.75,1} data are retained only as mechanism/sensitivity evidence.
They are not used to tune K_mis or to redefine the Full Proposed method.

Key files:
  TableA_core_ablation_mc_mean_std.csv
  TableB_core_ablation_pairwise_statistics.csv
  TableC_core_ablation_nominal.csv
  TableD_mismatch_sensitivity_mechanism.csv
  TableE_mismatch_endpoint_mechanism.csv
  TableF_mismatch_nominal_endpoint_mechanism.csv
  figures/FigA_core_ablation_mc.*
  figures/FigB_mismatch_tradeoff_and_leakage.*
  figures/FigC_pairwise_effect_forest.*
  figure_data/FigA_core_ablation_mc.csv/.mat
  figure_data/FigB_mismatch_tradeoff_and_leakage.csv/.mat
  figure_data/FigC_pairwise_effect_forest.csv/.mat

Every ablation PDF/PNG/JPEG has an identically named CSV and MAT figure-data snapshot.
