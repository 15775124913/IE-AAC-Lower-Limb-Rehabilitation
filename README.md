# IE-AAC: Interaction-Error-Driven Adaptive Admittance Control for Lower-Limb Rehabilitation

This repository contains the MATLAB simulation implementation and associated numerical evidence for **interaction-error-driven adaptive admittance control (IE-AAC)** in a simulated supine lower-limb rehabilitation configuration.

> **Research scope:** This is a numerical simulation study. The human–robot interaction (HRI) model is phenomenological, and the nonideal execution layer is a numerical stress-test model. The repository does **not** provide evidence of physical robot, patient, or clinical validation. Fixed-parameter stability checks and numerical audits should not be interpreted as a proof of stability or passivity of the complete time-varying adaptive closed loop.

## Main features

- Four-controller comparison under a shared reference and paired external realizations: **Fixed**, **Zhang-2018-derived**, **Cao-2024-derived**, and **Proposed IE-AAC**.
- Reference trajectory reconstruction and HRI-model-based nominal interaction compensation.
- Bounded virtual damping (**B**), stiffness (**K**), and interaction gain (**G**) adaptation, with fixed virtual inertia.
- Main 100-trial paired Monte Carlo evaluation, joint-level and three-joint statistics, including Holm-adjusted comparisons.
- Ablation, interaction-model mismatch, joint-specific parameter, nonideal numerical execution, dead-zone/noise sensitivity, and conditional stability analyses.
- Source data for selected manuscript figures and tables, together with numerical diagnostics.

## Requirements

- MATLAB **R2024b** was used for the archived implementation and numerical workflow. Other releases have not been validated in this release.
- Sufficient disk space and computing time for the complete batch of simulations.

The repository preserves the MATLAB source and previously generated result files. The presence of archived outputs does not mean that downloading the repository automatically reruns them.

## Quick start

1. Clone or download this repository.
2. Open MATLAB and set the **Current Folder** to the repository root (where `main_simulation.m` and `parameters.m` are located).
3. Ensure that old copies of the same function names from other IE-AAC project versions are **not** on the MATLAB path.
4. Run:

   ```matlab
   main_simulation
   ```

`main_simulation.m` is the single official workflow entry point. It initializes parameters, constructs the shared reference, evaluates all four controllers, runs the enabled uncertainty and supplementary studies, and writes dated results under `results/`.

**Important:** As distributed, the main script enables the full research workflow, including several 100-trial batches. Running it may take substantial time and generate additional files. Its `CLEAN_OLD_MC` flag is set to `false` by default; preserve historical evidence and avoid overwriting prior results. If you change any simulation setting for exploratory runs, record that change and store the new outputs separately.

## Repository organization

```text
IE-AAC-Lower-Limb-Rehabilitation/
├── main_simulation.m                 # Single workflow entry point
├── parameters.m                      # Parameters and study settings
├── TeachingDataPipeline.m           # Reference reconstruction
├── run_admittance_controller.m       # Four controller variants
├── monte_carlo_runner.m             # Primary paired simulations
├── export_main_mc_pairwise.m        # Paired inferential statistics
├── ablation_study.m                  # Controller ablations
├── run_execution_dynamics_validation.m
├── run_model_mismatch_validation.m
├── run_joint_specific_validation.m
├── run_deadzone_sensitivity_validation.m
├── run_stability_validation.m
├── ...                               # Other MATLAB helper modules (25 .m files total)
├── results/                          # Archived outputs and new dated runs
├── data/                             # Preexisting repository data/configuration files
├── README.md
└── LICENSE
```

All **25 `.m` files are located at the repository root**. The obsolete `main_ablation_study.m` entry point from the earlier repository version is not used in this revision.

## Locating archived results

The structure of `results/` reflects the original saved experiment outputs. **Preserve the dated subdirectories and original filenames**; the paths below describe file patterns, not guaranteed single fixed paths:

| Evidence | Representative files or folder patterns |
| --- | --- |
| Main paired Monte Carlo | `MC_trial_level_four_controller_*.csv`, `MC_summary_primary_four_controller_*.csv`, `MC_paired_change_four_controller_*.csv` |
| Main paired inference | `MainMC_27_joint_pairwise_Holm.csv`, `MainMC_9_threejoint_pairwise_Holm.csv`; in a complete new run these may be inside an `MC_*/paired_statistics/` folder |
| Nominal four-controller comparison | `nominal_primary_metrics_*.csv`, `nominal_supervisory_diagnostics_*.csv` |
| Ablation evidence | `ablation_mc_trial_metrics_long.csv`, `ablation_mc_summary_mean_std.csv`, `TableA_core_ablation_mc_mean_std.csv`, `TableB_core_ablation_pairwise_statistics.csv` |
| Model mismatch | `model_mismatch_trial_joint_metrics.csv`, `model_mismatch_group_summary.csv`, `mismatch_mc_trial_metrics_long.csv` |
| Joint-specific parameters | `joint_specific_trial_metrics.csv`, `joint_specific_group_summary.csv` |
| Nonideal execution | `execution_dynamics/execution_dynamics_results.csv` and related numerical audit files |
| Dead-zone/noise sensitivity | `figure_data/P11_deadzone_noise_*.csv` |
| Stability audit | `stability_validation/` and `figure_data/Section5_3_1_stability_*.csv` |
| Manuscript plot/table source data | `figure_data/` |
| Output inventory | `all_experiment_output_index.csv`, where present |

**Statistical interpretation:** The primary Monte Carlo design uses paired trials. The three-joint summary is computed within each trial before averaging across trials; joint/time samples must not be treated as independent paired trials. Separate Holm families are used for joint-level and three-joint-level tests.

## Reproducibility and scope notes

- The Zhang-2018-derived and Cao-2024-derived implementations are **study comparators**, not complete reproductions of the original papers' hardware systems.
- Changes to the randomization, controller settings, or time horizon can lead to results different from the archived numerical evidence.
- Nonideal execution tests do not constitute hardware-in-the-loop or physical-robot validation.
- Model-mismatch results are limited to the tested uncertainty settings.
- The stability analysis concerns specified fixed-parameter discrete subsystems and conditional checks, **not** a general proof for the complete adaptive HRI closed loop.
- Consult the associated manuscript for metric definitions, experimental protocols, and limitations. Manuscript figure numbering can differ from historical MATLAB-generated filenames after editorial revision.

## Version history

This revised codebase has **25 MATLAB source modules** and supersedes the earlier **15-module** layout preserved in this repository's Git history. The earlier repository version remains accessible through historical commits/tags. The archived results should be identified by their own filenames, timestamps, and provenance information rather than assumed to belong to a new run.

## License and citation

See [`LICENSE`](LICENSE) for the repository license. If using this code or data in a publication, please cite the associated IE-AAC study once its final bibliographic record is available. No new archival DOI is claimed here.
