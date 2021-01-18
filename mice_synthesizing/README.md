---
title: "Creating synthetic data with mice"
output: html_document
---

This folder contains all files that are used in the current study to generate synthetic data with the R-package `mice`. First, there are the .Rmd files [synthesizing_partially_synthetic_data_rules.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/synthesizing_partially_synthetic_data_rules.html) containing simulations and information on getting inferences from synthetic data using what Drechsler calls *partially synthetic data*; [paper.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/paper.html) containing the lay-out for a paper on generating synthetic data with `mice`, including a study into suitable pooling rules; [mice_updates.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/mice_updates.html) containing a list with adjustments that must be made to `mice` before it properly handles the creation of synthetic data; [mice_synthesizing_partitioned_data.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/mice_synthesizing_partitioned_data.html) containing simulations with a partitioned dataset (where the current project started with), as an approach to overcome the complicated nature of federated imputation; [Mice_synthesizing.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/Mice_Synthesizing.html) containing regular simulations without any partitioning involved, to determine the most suitable `mice` parameters; [formulas.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/formulas.html) containing the formulas for the variance of all variants for fully and partially synthetic data; and [data_checks.Rmd](https://amices.org/Federated_imputation/mice_synthesizing/data_checks.html) containing some plots and numbers to assess the suitability of synthetic data to overcome confidentiality constraints. Additionally, all corresponding `.html`-files are included and linked through the `.Rmd` file links. 