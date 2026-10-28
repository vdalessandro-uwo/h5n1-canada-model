# Quantifying Human Risk from H5N1 in Canada: model code

This repository contains the code for the multi-species SEIR(D) model with an environmental reservoir described in:

> D'Alessandro V, Poupliopoulou DV, Billias N, MacDermid JC, Miller E, McClinchey H, Quinn K, Stranges S, Weese S, Boutsikari E, Pereira TV, Johnstone J, Cipriano LE, Hopkins JP, Cini E, Wong JJ, Bobos P. *Quantifying Human Risk from H5N1 in Canada: A Multi-species Mathematical Model of Repeated Spillover Events and Public Health and Social Measures.* (manuscript under review)

The model represents wild birds, backyard poultry, large- and small-scale commercial poultry, reassortment hosts, poultry workers and the general population, linked through an environmental reservoir. It is calibrated by Bayesian MCMC to synthetic Canada-scale weekly case targets, attributes modelled human infections to source pathways, compares five public health and social measures (PHSMs), and adds a one-year extension with restocking of depopulated commercial premises and seasonal wild-bird importation.

Repository: https://github.com/vdalessandro-uwo/h5n1-canada-model  
Archived version: https://doi.org/10.5281/zenodo.23733343

## Contents

| File | What it does | Manuscript outputs |
|---|---|---|
| `R/01_main_analysis.R` | Model definition, synthetic targets, MCMC calibration, posterior predictive checks, reproduction numbers, source attribution, Spearman sensitivity analysis, PHSM scenarios, deaths, one-year extension, figures | Tables 2, 4, 5, 6; Figures 2–6; Supplemental Tables 1–3; Supplemental Figure 1; Technical Appendix Figures 2–3 |
| `R/02_sensitivity_analyses.R` | Scale of the scenario, decomposition of the one-year burden, Eq. 18 sensitivity, importation-rate sensitivity, prior sensitivity of source attribution (three re-calibrations) | Table 7; Supplemental Tables 4, 5, 6; one-year Results paragraph |
| `R/03_recovery_experiments.R` | Simulation study: all scenario outputs at the data-generating values, recovery of the generating values (new Table 2), the beta_HH experiment, recovery experiments for four key coefficients with replicate noise seeds, and a weak-prior fit | Table 2, Figure 7, and the true-world versions of Tables 4–6 and Figures 4–6 |
| `figures/figure1_schematic.py` | Model structure schematic | Figure 1 |

## Requirements

- R 4.5.1 (earlier 4.x versions should work)
- R packages: `deSolve`, `dplyr`, `tidyr`, `tibble`, `ggplot2`, `purrr`, `stringr`, `forcats`, `coda`, `scales`, `patchwork`
- Python 3 with `matplotlib`, for Figure 1 only

```r
install.packages(c("deSolve", "dplyr", "tidyr", "tibble", "ggplot2", "purrr",
                   "stringr", "forcats", "coda", "scales", "patchwork"))
```

## How to run

Run the R scripts in order, in the same R session, from the repository folder:

```r
setwd("path/to/this/repository")
source("R/01_main_analysis.R")        # about 1-2 hours
source("R/02_sensitivity_analyses.R") # about 3 hours; sections can be switched off at the top
source("R/03_recovery_experiments.R") # about 6 hours; resumes from saved fits if interrupted
```

Script 1 writes figures and tables to the working directory. Scripts 2 and 3 write to `outputs/` and `outputs/recovery/` and summarises its results in `outputs/sensitivity_summary.txt`. On macOS and Linux the MCMC chains in script 2 run in parallel; on Windows they run one after another.

Figure 1 is drawn separately:

```bash
cd figures
python3 figure1_schematic.py
```

## Reproducibility

All random steps use fixed seeds, so repeated runs give the same results on the same platform. The main settings are:

- **MCMC:** 4 chains × 30,000 iterations, burn-in 10,000, thinning 10 (8,000 retained draws). Prior set (iii) in script 2 uses 60,000 iterations with burn-in 20,000.
- **Posterior predictive checks, Table 4 and the PHSM comparison (Table 5):** the same 400 posterior draws.
- **One-year extension:** posterior-median parameters, 200 replicates. Replicate *r* uses seed `20261020 + r` in every scenario, so scenarios are compared on common random numbers. The burden-linked spillover term (Eq. 18) is off in the main analysis (`kappa_link = 0`) and on (`kappa_link = 1e-6`) in the sensitivity analysis.

The calibration targets are synthetic: they are generated inside script 1 from the model with a reporting probability of 0.6 and a negative-binomial dispersion of 50. No observed human case data are used or needed.

## Licence

Code is released under the MIT Licence (see `LICENSE`).

## Citation

If you use this code, please cite the article above and this repository (see `CITATION.cff`).

## Contact

Victoria D'Alessandro, School of Physical Therapy, University of Western Ontario.
