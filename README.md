# Mixed effects models for healthcare longitudinal data with an informative visiting process: a Monte Carlo simulation study

This repository contains the code required to simulate data, fit each model, and produce summary tables and figures for the manuscript _Mixed effects models for healthcare longitudinal data with an informative visiting process: a Monte Carlo simulation study_ by Gasparini _et al_. A preprint of the manuscript is available on [arXiv](https://arxiv.org/abs/1808.00419).

The following files are included in this repository:

1. `sim1_make-sim-data.do`, a Stata `.do` file used to simulate the data under each data-generating mechanism;
2. `sim1_run-sim.do`, a Stata `.do` file used to fit the models and save the results;
3. six datasets in Stata format from a single replication of each data-generating mechanism (`df_1_1.dta` for scenario 1, `df_2_1.dta` for scenario 2, and so on);
4. `sim1_mf.R`, an R script used to summarise results from the simulation study and produce tables and figures included in the manuscript.
