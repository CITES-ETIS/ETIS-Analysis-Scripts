################################################################################
# Wrapper script for cluster analysis steps                                    #
#                                                                              #
# Graham Laidler, May 2025                                                     #
#                                                                              #
################################################################################


# Load functions
source("Scripts/functions.R")

# Analysis parameters
analysis_name <- "CoP20"                      # Name to use in filenames of outputs

yearfrom <- 2008                              # Period of trend analysis...
yearto <- 2023                                # ...
CoP_Period <- 2021:2023                       # Period of cluster analysis

filename = paste0("ave_", analysis_name)      # Filename for trend analysis model to use


# Prepare cluster analysis variables
source("Scripts/3-01-prepare_variables.R")

# Run cluster analysis and prepare figures
source("Scripts/3-02-clustering_and_plotting.R")

# Dendrogram and heatmap figure
ggarrange(dendrogram_plot, heatmap_plot, ncol = 1,
  heights = c(1, 1.5),  # Adjust the height ratio
  align = "v"
)
ggsave(paste0("Figures/cluster_analysis_", filename, ".png"), width = 10, height = 5, units = "in", dpi = 600, bg = "white")


# Sensitivity analysis - uses variables from scripts 3-01 and 3-02, so make sure these have been run first
source("Scripts/3-03-sensitivity_analysis.R")

# Sensitivity analysis heatmap
sensitivity_heatmap
ggsave(paste0("Figures/sensitivity_analysis_", filename, ".png"), width = 8, height = 7, units = "in", dpi = 600, bg = "white")


