# IE-AAC: Adaptive Admittance Control for Lower-Limb Rehabilitation Robots

**MATLAB research code accompanying the manuscript**

> **Adaptive Admittance Control for Lower-Limb Rehabilitation Robots with Nominal Interaction Mismatch Compensation**  
> Jiaqi Zhou, Yuanliang Zhang, and Mingda Miao  
> Manuscript prepared for submission to *Machines* (MDPI).

Repository: https://github.com/15775124913/IE-AAC-Lower-Limb-Rehabilitation

## Overview

This repository provides a MATLAB simulation implementation of **interaction-error-driven adaptive admittance control (IE-AAC)** for a supine, three-degree-of-freedom lower-limb rehabilitation robot. The modeled sagittal-plane joints are the **hip, knee, and ankle**.

The controller uses **nominal interaction mismatch compensation** to distinguish state-dependent nominal interaction changes induced by the commanded motion from the interaction residual used for admittance adaptation. Virtual inertia is fixed; virtual damping, stiffness, and interaction mapping gain are adapted subject to parameter and motion limits.

The study is **simulation-based**. Therapist demonstrations, patient–robot interaction, disturbances, and uncertainty realizations are simulated; the repository does not represent physical-robot or clinical validation.

## Study workflow

The source implements the following pipeline:

1. Generate a simulated therapist demonstration at **10 Hz**, including time warping and measurement noise.
2. Screen key states and construct a piecewise-quintic reference trajectory at **100 Hz**.
3. Repeat the shared reference over **14 rehabilitation cycles**.
4. Compare four closed-loop controllers using a common reference and paired disturbance realizations:
   - **Fixed** nominal admittance;
   - **Zhang-2018-derived** variable admittance;
   - **Cao-2024-derived** momentum-based variable admittance;
   - **Proposed IE-AAC**.
5. Perform **100 paired Monte Carlo trials**, core-module ablation, and nominal-interaction-mismatch sensitivity analysis.
6. Export simulation metrics, audits, figures, and figure-level numerical data.

The Zhang-2018 and Cao-2024 comparison controllers are **literature-derived implementations adapted to this simulation framework**, not claims of exact replication of the original authors' full systems.

## Requirements

- **MATLAB R2024b** (the MATLAB release reported by the author for this project).
- **Signal Processing Toolbox** (`sgolayfilt` in `TeachingDataPipeline.m`).
- **Statistics and Machine Learning Toolbox** (`boxplot` used for Monte Carlo and paper figures).

**Validation note:** The code structure and dependencies have been inspected, but the complete MATLAB simulations have **not** been executed independently as part of preparing this README. Before the archival release, the authors should run the complete manuscript workflow on **MATLAB R2024b** and verify that the outputs match the reported results.

## Quick start

1. Download this repository using **Code → Download ZIP**, then extract it, or clone the repository with Git.
2. Open MATLAB and set **Current Folder** to the directory containing `main_simulation.m` and `parameters.m`.
3. In the MATLAB Command Window run:

```matlab
main_simulation
```

This runs the project's full one-click workflow. **Expect a substantial computation time**, because the default includes the primary 100-trial Monte Carlo comparison, a 100-trial core ablation, a 100-trial mismatch-sensitivity study, and final evidence aggregation.

The relevant settings near the top of `main_simulation.m` are:

```matlab
RUN_ABLATION_WORKFLOW = true;
ABLATION_NUM_RUNS = 100;
MAKE_ABLATION_FIGURES = true;
KMIS_GRID = [0 0.25 0.5 0.75 1];

% Later in the same script:
RUN_MONTE_CARLO = true;
MC_NUM_RUNS = 100;
CLEAN_OLD_MC = false;
```

**For manuscript-level reproduction, retain the formal settings.** Do not change the seed or comparison-model parameters without documenting the change. Keep `CLEAN_OLD_MC = false` if previously generated results should be retained.

### Optional short smoke test

To test the core-ablation code path with a smaller number of trials, without running the full primary workflow:

```matlab
ablation_study('core', 'NumRuns', 5, 'MakeFigures', false);
```

This is only a **functional smoke test** and cannot reproduce the manuscript's 100-trial statistical results.

## Key source files

| File | Role |
| --- | --- |
| `main_simulation.m` | Primary one-click workflow: paired Monte Carlo, nominal comparison, figures, ablations, and final evidence |
| `main_ablation_study.m` | Optional wrapper for selectively running ablation stages |
| `parameters.m` | Simulation, controller, uncertainty, and reference parameters |
| `TeachingDataPipeline.m` | Simulated demonstration, key-state screening, and quintic reference reconstruction |
| `event_model.m` | Joint-specific transient events and event windows |
| `run_admittance_controller.m` | Closed-loop implementations of all four comparison controllers |
| `compute_trajectories.m` | Planar forward kinematics for joint and foot-point visualization |
| `compute_metrics.m` | Tracking, smoothness, interaction-residual, and other metrics |
| `nominal_diagnostics.m` | Structural, supervisory, and motion-limit audits |
| `mc_trial_utils.m` | Paired uncertainty sampling and pre-generated disturbance sequences |
| `monte_carlo_runner.m` | Primary four-controller Monte Carlo evaluation |
| `save_results.m` | Nominal simulation outputs and tables |
| `generate_paper_figures.m` | Manuscript figures and associated numerical data |
| `paper_figure_utils.m` | Shared figure export and formatting utilities |
| `ablation_study.m` | Core ablation, mismatch sensitivity, and final-evidence aggregation |

All 15 MATLAB files are intended to remain together in the repository root so the current entry points resolve correctly.

## Generated files and result organization

The scripts create output directories as needed. Typical output locations include:

```text
IE-AAC-Lower-Limb-Rehabilitation/
├── main_simulation.m
├── ... (14 other MATLAB source files)
├── data/
│   └── config_*.mat                # Generated parameter snapshots
└── results/
    ├── <timestamp>/               # Nominal results, audits, and paper figures
    │   ├── paper_figures/
    │   └── figure_data/
    ├── MC_<timestamp>/            # Four-controller paired Monte Carlo outputs
    ├── ABLATION_<timestamp>/      # Core-ablation outputs
    ├── MISMATCH_*/                # Mismatch-sensitivity outputs
    └── ...                        # Final evidence and supporting files
```

Representative files written by the current scripts include:

- `MC_trial_level_four_controller_*.csv` — paired trial-level controller metrics.
- `MC_summary_primary_four_controller_*.csv` — primary Monte Carlo summary.
- `MC_paired_change_four_controller_*.csv` — paired changes between controllers.
- `ablation_mc_trial_metrics_long.csv` and `ablation_mc_paired_full_vs_ablations.csv` — core-ablation results.
- `mismatch_mc_sensitivity_compact.csv` — mismatch-sensitivity summary.
- `TableA_core_ablation_mc_mean_std.csv` through `TableF_mismatch_nominal_endpoint_mechanism.csv` — aggregated final-evidence tables.

Folder names and timestamps may differ across runs. **Generated results are not automatically included in a source-code download**; consult the repository contents or an explicitly linked data archive for any precomputed datasets made public by the authors. The output files in a new run should be distinguished from those used for the manuscript's reported statistics.

## Reproducibility and interpretation

Within each paired trial, the formal comparison methods share the reconstructed reference, sampled model parameters, and pre-generated exogenous noise/event sequences. State-dependent interaction terms are recomputed from each controller's own commanded state. Primary reported outcomes include joint-reference deviation, motion smoothness (jerk), and raw interaction-residual measures.

The simulated demonstration and phenomenological interaction model are not recordings from patients. The results should therefore be interpreted as **numerical evidence under the stated modeling assumptions**, rather than experimental or clinical effectiveness.

## Citation

If you use this code, please cite the associated manuscript once a bibliographic record becomes available. A versioned software citation and archive DOI can be added after an official repository release.

**Software repository:** https://github.com/15775124913/IE-AAC-Lower-Limb-Rehabilitation

## License

Refer to the repository's `LICENSE` file **if one has been added**. Public accessibility by itself does not grant a software reuse license. The authors should confirm rights and the intended open-source license before the final release.
