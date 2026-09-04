# Structural Identifiability Tutorial in Julia — companion code (v2)

Companion material for **A Tutorial on Symbolic Structural Identifiability Analysis of ODE Models in Julia** (A. Alsammani, *Bulletin of Mathematical Biology*, manuscript BMAB-D-26-00484, revised version).

## Contents

| File | Purpose |
|---|---|
| `Structural_Identifiability_Tutorial_Julia_v2.ipynb` | Corrected companion notebook for the revised manuscript (six case studies, Sections 5–6; ModelingToolkit/Catalyst routes, Section 3.4; bilinear example, Section 6; Fisher-information demonstration, Section 7.6). Outputs are cleared; regenerate them by running the notebook top to bottom. |
| `Project.toml` | Direct dependencies (with `[compat]` bounds for Julia, StructuralIdentifiability.jl, ModelingToolkit.jl). |
| `Manifest.toml` | Exact pinned versions of every package in the environment (Julia 1.12.2; StructuralIdentifiability.jl 0.5.26; ModelingToolkit.jl 11.38.1; OrdinaryDiffEq.jl 7.3.0; Plots.jl 1.41.6). |
| `LICENSE` | MIT. |
| `Structural_Identifiability_Tutorial_Julia.ipynb` | The notebook of the original submission (Zenodo v1.0.1), retained unchanged in the repository for the record (it is not part of this v2 package build and must be carried over from the v1.0.1 archive). Its two-compartment, viral-dynamics, SIWR and SEIR-H results were superseded by v2; see the revised manuscript. |

## Reproducing the results

1. Install Julia 1.12 (https://julialang.org/downloads/ or `juliaup add 1.12`).
2. Clone or download this repository and start Julia **in its directory**:
   ```bash
   julia --project=.
   ```
3. Instantiate the pinned environment once:
   ```julia
   using Pkg; Pkg.instantiate()
   ```
4. Launch Jupyter with the Julia kernel (`using IJulia; notebook(dir=pwd())`, or `jupyter lab`) and run
   `Structural_Identifiability_Tutorial_Julia_v2.ipynb` with **Kernel → Restart & Run All**.
   The first code cell activates this directory's environment; the second prints the Julia and package versions.

The global algorithms of StructuralIdentifiability.jl are randomized (default correctness probability 0.99). Identifiability verdicts are reproducible; the printed form of a generating set of identifiable functions may differ between runs (equivalent generators), and runtimes vary.

### Note on Catalyst

`Catalyst.jl` is required only for the Catalyst route in the "Building models with ModelingToolkit and Catalyst" cell. It is declared in `Project.toml`; the shipped `Manifest.toml` was generated before that declaration, so on first use `Pkg.instantiate()` (or `Pkg.resolve()` followed by `Pkg.instantiate()`) resolves and adds Catalyst and rewrites `Manifest.toml`; the updated manifest should then be committed so that later users get a fully pinned environment. The cell is guarded: if Catalyst cannot be loaded, the notebook reports this and continues, and every other result is unaffected.

## Correspondence with the manuscript

| Notebook cell(s) | Manuscript |
|---|---|
| ModelingToolkit / Catalyst routes | Section 3.4 |
| Case Study 1 exponential decay (+ scaled observation) | Section 5.1 |
| Case Study 2 SIR (+ prevalence, (I,R)) | Section 5.2 |
| Case Study 3 two-compartment PK, y = x₂; central-compartment variant; Figure 5(B) check | Section 5.3, Table 3, Figure 5 |
| Case Study 4 viral dynamics (V; V and T; known T(0); p·s fibre) | Section 5.4 |
| Case Study 5 SIWR (incidence; incidence + W; prevalence) | Section 5.5 |
| Case Study 6 SEIR-H (admissions; admissions + incidence; extended bonus model) | Section 5.6 |
| Bilinear model, Figure 5(A) check, summary chart | Sections 5.7, 6.1–6.2 |
| Fisher-information demonstration | Section 7.6 |

## Citation

Please cite the manuscript and the Zenodo record of this repository (DOI in the manuscript's Data and Code Availability statement).
