################################################################################
# Wrapper script for data processing steps                                     #
#                                                                              #
# Graham Laidler, May 2025                                                     #
#                                                                              #
################################################################################


# load functions
source("Scripts/functions.R")

# Analysis parameters
analysis_name <- "CoP20"     # Name to use in filenames of outputs

yearfrom <- 2008             # Period of trend analysis...
yearto <- 2023               # ...


# Date of data validation notification..
# ..This is used in script 1-01 to exclude non-MA records submitted after this date,
# ..or MA records submitted after this date which implicate another Party 
data_val_notification <- "2024-05-31"

sz_ids_remove <- c()         # Set up vector to store seizure ids to be removed from analysis



#### Remove certain seizures ---------------------------------------------------
# Fill in and run this section as appropriate to remove certain seizures from analysis
# E.g....

# ...To remove seizures with unresolved review requests
seizures <- read_csv("Original Data/seizures.csv")

unresolved_requests <- seizures %>%
  filter(grepl("unresolved", status_note))

sz_ids_remove <- c(sz_ids_remove, unresolved_requests$id)


# ...To remove South Sudan's seizures
seizures <- read_csv("Original Data/seizures.csv")

countries <- read_csv("Original Data/countries.csv") %>%
  replace_na(list(code = "NAM")) # NA for Namibia gets read as missing, so replace with NAM

SS_id <- filter(countries, code == "SS")$id

SS_seizures <- seizures %>%
  filter(discovered_country_id == SS_id)

sz_ids_remove <- c(sz_ids_remove, SS_seizures$id)


#### Run data processing scripts -----------------------------------------------
source("Scripts/1-01-create_seizure_table.R")
source("Scripts/1-02-weight_estimation_models.R")
source("Scripts/1-03-estimate_missing_weights.R")
source("Scripts/1-04-calculate_covariates.R")
source("Scripts/1-05-set_up_data_for_analysis.R")

