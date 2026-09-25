# ETIS Analysis Scripts

These scripts contain the code used to run the latest ETIS analysis.

The scripts are divided into 3 stages. Each stage can be run from its wrapper script:

- 1-Data_Processing.R:
  - Convert data csvs from the ETIS Online database into usable seizures data with trade route information
  - Build weight models and perform weight estimation for records with missing weights
  - Calculate covariates for modelling seizure rate and reporting rate bias-adjustment

- 2-Trend_Analysis.R:
  - Fit Bayesian hierarchical trend models and perform model averaging
  - Fit weight models and simulate weights for the Weight Indices
  - Plot Transaction Indices and Weight Indices

- 3-Cluster_Analysis.R:
  - Perform agglomerative hierarchical clustering on country-level variables
  - Perform sensitivity analysis with respect to posterior uncertainty in clustering variables
