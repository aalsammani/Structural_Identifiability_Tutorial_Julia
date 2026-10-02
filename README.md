# Structural Identifiability Tutorial in Julia: companion code (version 3)

Companion material for **A Tutorial on Symbolic Structural Identifiability Analysis of ODE Models in Julia** (A. Alsammani, *Bulletin of Mathematical Biology*, manuscript BMAB-D-26-00484, second revision). Section, figure, and table numbers below refer to that revision.

## Contents

| Path | Purpose |
|---|---|
| `scripts/01_exponential_decay.jl` … `scripts/14_siwr_exchange.jl` | One script per example; see the table below. `scripts/common.jl` holds shared helpers. |
| `run_all.jl` | Runs every script in a fresh Julia process and writes its output to `expected_output/`. |
| `expected_output/` | Output of each script in the tested environment (reference for comparison). |
| `figures/make_figures.jl` | Regenerates Figures 3 and 4 (`fig_verdict_matrix.pdf`, `fig_nonidentifiability.pdf`). Figure 3 is computed with StructuralIdentifiability.jl. |
| `models/enzyme_kinetics.xml` | SBML Level 3 file of the Michaelis–Menten mechanism used in Section 3.4. |
| `Structural_Identifiability_Tutorial_Julia_v3.ipynb` | Interactive notebook covering the case studies (outputs cleared; regenerate with Kernel → Restart & Run All). |
| `Project.toml`, `Manifest.toml` | The tested environment (see below). |
| `LICENSE` | MIT. |

## Scripts and manuscript sections

| Script | Content | Section |
|---|---|---|
| `01_exponential_decay.jl` | Exponential decay | 4.1 |
| `02_sir.jl` | SIR model, incidence observed | 4.2 |
| `03_two_compartment_pk.jl` | Two-compartment model: peripheral output, central output with gain, no input | 4.3 |
| `04_viral_dynamics.jl` | Viral dynamics: V observed, known T(0) (`known_ic`), V and T observed | 4.4 |
| `05_siwr.jl` | SIWR: incidence; incidence and W; prevalence; prevalence and W | 4.5 |
| `06_seirh.jl` | SEIR-H: admissions; admissions and incidence; fixed sigma; occupancy | 4.6 |
| `07_bilinear.jl` | Bilinear model | 5.1 |
| `08_modelingtoolkit_catalyst.jl` | SIR model built with ModelingToolkit and Catalyst | 3.4 |
| `09_sbml_import.jl` | SBML → SBMLImporter.jl → Catalyst → StructuralIdentifiability.jl | 3.4 |
| `10_multiple_conditions.jl` | Pooling two experimental conditions by model augmentation | 6.8 |
| `11_rationalization.jl` | Two liftings of the generalized growth model | 6.8 |
| `12_timing.jl` | Runtimes of the local and global tests | 6.1 |
| `13_io_coefficients.jl` | Example 2.14 of Hong et al. (2020): non-identifiable input–output coefficient | 2.6 |
| `14_siwr_exchange.jl` | Simulation check of the SIWR exchange transformation | 4.5, Figure 4(D) |
| `figures/make_figures.jl` | Figures 3 and 4 | 4.7, 5.1 |

## Reproducing the results

1. Install Julia 1.12.2 (for example `juliaup add 1.12.2`).
2. In this directory, install the recorded package versions:
   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate()"
   ```
3. Run all scripts (about 15 minutes on a laptop, most of it compilation):
   ```bash
   julia --project=. run_all.jl
   ```
   and compare `expected_output/` with the files distributed here, or run a single script, for example `julia --project=. scripts/04_viral_dynamics.jl`.
4. Regenerate the data figures:
   ```bash
   julia --project=. figures/make_figures.jl
   ```

The global algorithms of StructuralIdentifiability.jl are randomized (default probability of correctness 0.99). Identifiability verdicts are reproducible; the printed form and order of generating sets of identifiable functions may differ between runs, and runtimes vary between machines.

## Tested environment

All scripts, the figure script, and the code cells of the notebook were run on 2026-10-01 under Windows 11 with Julia 1.12.2 in the environment recorded in `Manifest.toml`. Main package versions:

| Package | Version |
|---|---|
| StructuralIdentifiability.jl | 0.5.26 |
| ModelingToolkit.jl | 11.38.1 |
| Catalyst.jl | 16.4.0 |
| SBMLImporter.jl | 4.1.2 |
| ReactionNetworkImporters.jl | 1.5.0 |
| OrdinaryDiffEq.jl | 7.3.0 |
| Plots.jl | 1.41.6 |
| TimerOutputs.jl | 0.5.29 (restricted, see below) |

**Why TimerOutputs.jl is restricted.** The environment of the first revision (archive version 2.0.0) did not contain Catalyst.jl. Adding it required resolving the dependencies again, which updated TimerOutputs.jl, an indirect dependency of StructuralIdentifiability.jl, to version 1.2.2. StructuralIdentifiability.jl 0.5.26 declares compatibility with TimerOutputs 0.5 and 1.0, but with version 1.2.2 its global analysis fails (`MethodError: no method matching power_series_solution(...)` in the Wronskian step). Version 0.5.29 works. `Project.toml` therefore lists TimerOutputs.jl with the compatibility bound `0.5.29`. Use `Pkg.instantiate()` with the shipped `Manifest.toml` rather than adding packages afresh.

## Changes with respect to version 2.0.0

- Scripts, expected output, figure scripts, and the SBML model file added.
- Catalyst.jl, SBMLImporter.jl, and ReactionNetworkImporters.jl added to the environment; TimerOutputs.jl restricted as explained above.
- Notebook corrected: the first code cell did not parse (`? @__DIR__ : pwd()` must be written `? (@__DIR__) : pwd()`); the Catalyst cell expanded `@reaction_network` before loading Catalyst; the `known_ic` cell passed a state of a different model; section pointers updated to the second revision; outdated statements ("seven case studies", "essentially free") corrected.
- The README of version 2.0.0 referred to `Structural_Identifiability_Tutorial_Julia_v2.ipynb`, which was not in that archive; this version contains the notebook it describes.

## Citation

Please cite the manuscript and the Zenodo record of this repository (concept DOI 10.5281/zenodo.18684343, which resolves to the latest version).
