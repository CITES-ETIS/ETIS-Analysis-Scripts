################################################################################
# Wrapper script for trend analysis steps                                      #
#                                                                              #
# Graham Laidler, May 2025                                                     #
#                                                                              #
################################################################################


# Load functions
source("Scripts/functions.R")

# Analysis parameters
analysis_name <- "CoP20"     # Name to use in filenames of outputs

yearfrom <- 2008             # Period of trend analysis...
yearto <- 2023               # ...



#### Model fitting - Transaction Index -----------------------------------------

# Determine maximum number of knots
K <- 1
while(knot_segment_min(k = K) >= 4){
  K <- K + 1
}
K <- K - 1
print(paste0("Run models up to ", K, " knots"))

#### Fit models in script 2-01 ####
# Note - may need to add/remove models from script 2-01 depending on the number of models to be run, K. Currently set for K = 6.
# Running the models can be time-consuming, so may want to run them individually from script 2-01
source("Scripts/2-01-fit_JAGS_trend_models.R")


#### Check MCMC diagnostics 

# Generate Rhat plots for all parameters as an html with R markdown
rmarkdown::render("Scripts/Rhat.Rmd", output_dir = "JAGS MCMC outputs", output_file = paste0(analysis_name, "_Rhat.html"))
# Can also check denplot, traplot, caterplot from mcmcplots, and gelman.plot from coda


#### Model averaging - adjust depending of number of models used
models <- paste0(c("1knot_", "2knot_", "3knot_", "4knot_", "5knot_", "6knot_"), analysis_name)

source("Scripts/2-02-model_averaging.R")


#### Model fitting - Weight Index ----------------------------------------------

# Fit weight models
source("Scripts/2-03-fit_JAGS_weight_models.R")

# Check MCMC diagnostics

# Simulate weights for Weight Index
source("Scripts/2-04-simulate_weights.R")



#### Graphical outputs ---------------------------------------------------------
# See functions.R script for details of these plotting functions

# Transaction Index
TIplot(filename = paste0("ave_", analysis_name))
ggsave(paste0("Figures/TI_", analysis_name, ".png"), width = 7, height = 8, units = "in", dpi = 600, bg = "white")

# Posterior predictive plot
PPplot(filename = paste0("ave_", analysis_name), df.analysis = df.analysis)      # df.analysis is the analysis data from script 2-01
ggsave(paste0("Figures/PP_", analysis_name, ".png"), width = 7, height = 8, units = "in", dpi = 600, bg = "white")

# Weight Index
WIplot(filename = analysis_name)
ggsave(paste0("Figures/WI_", analysis_name, ".png"), width = 10.5, height = 4.5, units = "in", dpi = 600, bg = "white")


