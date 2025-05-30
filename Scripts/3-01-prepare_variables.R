################################################################################
# - Prepare variables for cluster analysis                                     #
#                                                                              #
################################################################################


#### Transaction Index variables -----------------------------------------------

# Get lambda from model (filename) specified in cluster analysis wrapper script 3
lambda <- Get_lambda(filename)

# Configure dimensions
n.yrs <- length(yearfrom:yearto)
n.sim <- dim(lambda)[1]
n.ctry <- dim(lambda)[2]/n.yrs
dim(lambda) <- c(n.sim, n.ctry, n.yrs, 5)
lambda <- aperm(lambda, c(2,3,4,1))           # make lambda to have dimensions c(n.ctry, n.yr, 5 (n.classes), n.sim)

# Get trend analysis inclusion countries
inclusion.cts <- read_rds(paste0("Processed Data/inclusioncts_", analysis_name, ".rds"))

# mean lambda for each country and year and ivory class
mean_lambdas <- apply(lambda, c(1,2,3), mean, na.rm = TRUE)

# reduce to relevant CoP period
mean_lambdas <- mean_lambdas[, which(yearfrom:yearto %in% CoP_Period), ]

# Sum lambdas over the CoP period years
mean_lambdas <- apply(mean_lambdas, c(1,3), sum, na.rm = TRUE)

# Format into dataframe
mean_lambdas <- as.data.frame(mean_lambdas) %>%
  mutate(country = inclusion.cts) %>%
  relocate(country) %>%
  rename_with(~ paste0("lambda_", seq_along(.)), -country)




#### Phi and theta bias-adjustment for all countries and years -----------------
# For each posterior sample of the b and g coefficients from the trend analysis model, 
# calculate the seizure (phi) and reporting rate (theta) for each country and year
# by using the country and year specific covariates.
# Then get the mean of the seizure and reporting rates over the posterior samples.

# Get seizure and reporting rate covariates 
covariates <- read_csv(paste0("Processed Data/covars_", analysis_name, ".csv"))

# Get posterior samples of b and g coefficients
b1 <- readRDS(paste('JAGS MCMC outputs/Processed outputs/b1_', filename, '.rds',sep=""))
b2 <- readRDS(paste('JAGS MCMC outputs/Processed outputs/b2_', filename, '.rds',sep=""))
g1 <- readRDS(paste('JAGS MCMC outputs/Processed outputs/g1_', filename, '.rds',sep=""))
g2 <- readRDS(paste('JAGS MCMC outputs/Processed outputs/g2_', filename, '.rds',sep=""))

# Filter covariates to CoP period years
covariates <- covariates %>%
  filter(year %in% CoP_Period)

# Get unique countries and years
countries <- unique(covariates$country)
years <- unique(covariates$year)

# Define dimensions
n.ctry <- length(countries)
n.yr <- length(years)
n.sim <- length(b1)  # Should be 40000

# Expand the covariates dataframe to be repeated for each of the n.sim samples of b1, b2, g1 and g2
covariates_expanded <- covariates %>%
  slice(rep(1:n(), each = n.sim)) %>%
  mutate(sample = rep(1:n.sim, times = nrow(covariates)),
         b1 = rep(b1, times = nrow(covariates)),
         b2 = rep(b2, times = nrow(covariates)),
         g1 = rep(g1, times = nrow(covariates)),
         g2 = rep(g2, times = nrow(covariates)))

# Compute phi and theta for each sample
covariates_expanded <- covariates_expanded %>%
  mutate(phi = 1 / (1 + exp(-(b1 * LE1 + b2 * TCI))),
         theta = 1 / (1 + exp(-(g1 * ETIS.rep + g2 * CITES.rep.logit))))

# Reshape phis and thetas into 3D arrays
phi_array <- array(covariates_expanded$phi,
                   dim = c(n.sim, n.ctry, n.yr))

theta_array <- array(covariates_expanded$theta,
                     dim = c(n.sim, n.ctry, n.yr))

# Permute dimensions to [countries, years, samples]
phi_array <- aperm(phi_array, c(2, 3, 1))
theta_array <- aperm(theta_array, c(2, 3, 1))

# Calculate posterior mean phis and thetas for each country and year
mean_phis <- apply(phi_array, c(1,2), mean, na.rm = TRUE)
mean_thetas <- apply(theta_array, c(1,2), mean, na.rm = TRUE)

# Reformat to dataframes
mean_phis <- as.data.frame(mean_phis)
mean_phis$country <- countries
colnames(mean_phis)[1:n.yr] <- CoP_Period

mean_thetas <- as.data.frame(mean_thetas)
mean_thetas$country <- countries
colnames(mean_thetas)[1:n.yr] <- CoP_Period

# Pivot to long formats and combine
mean_phis <- mean_phis %>%
  pivot_longer(-country, names_to = "year", values_to = "phi")

mean_adjs <- mean_thetas %>%
  pivot_longer(-country, names_to = "year", values_to = "theta") %>%
  left_join(mean_phis)
# mean_adjs contains the mean theta and phi for each country and year

# Make the year column numeric
mean_adjs$year <- as.numeric(mean_adjs$year)




#### Bias-adjusted seizures ----------------------------------------------------

# Get seizures table
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_with_est_weights_", analysis_name, ".csv"))

# Fill in missing proportions when multiple countries of origin
seizures <- seizures %>%
  group_by(seizure_id, ivory_type) %>%
  mutate(n_origins = n()) %>% 
  mutate(proportion = ifelse(is.na(proportion), 100/n_origins, proportion)) %>%
  ungroup()

# Filter years to CoP period and merge seizures with their phi and theta bias-adjustment
seizures <- seizures %>%
  filter(seizure_year %in% CoP_Period) %>%
  left_join(mean_adjs, join_by("seizure_year" == "year", "discovered_country_code" == "country"))

# Add RIE weight column (raw and worked combined) and calculate bias-adjusted number and weight of seizures
RIEs <- seizures %>%
  distinct(seizure_id, ivory_type, .keep_all = TRUE) %>%
  group_by(seizure_id) %>%
  mutate(RIE_weight = sum(raw_weight, worked_weight, na.rm = TRUE)) %>%
  select(seizure_id, RIE_weight) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  ungroup()
seizures <- seizures %>%
  left_join(RIEs, by = "seizure_id") %>%
  mutate(sz_adj = 1/(theta*phi),
         RIE_adj = RIE_weight/(theta*phi))
# Note - RIEs contains RIE weight of each unique seizure...
# ... seizures still contains multiple rows where a seizure has both raw and worked or multiple countries of origin



#### Bias-adjusted seizure variables -------------------------------------------
# For each country and year (in CoP period), calculate the following variables:
# - sz_out_1: bias-adjusted seizures-out which are less than 500 kg
# - sz_out_2: bias-adjusted seizures-out which are at least 500 kg
# - wt_in_1: bias-adjusted weight-in from seizures-in which are less than 500 kg
# - wt_in_2: bias-adjusted weight-in from seizures-in which are at least 500 kg
# - wt_out_1: bias-adjusted weight-out from seizures-out which are less than 500 kg
# - wt_out_2: bias-adjusted weight-out from seizures-out which are at least 500 kg

# Initialise a dataframe to store the country-year variables
vars <- mean_adjs %>%
  dplyr::select(country, year) %>%
  mutate(sz_out_1 = 0,
         sz_out_2 = 0,
         wt_in_1 = 0,
         wt_in_2 = 0,
         wt_out_1 = 0,
         wt_out_2 = 0)

# Calculate the six variables for each country and year 
for (i in 1:nrow(vars)){
  
  cty <- vars$country[i]
  yr <- vars$year[i]
  
  # Seizures-out
  temp <- seizures %>%
    filter(seizure_year == yr & discovered_country_code != cty) %>%
    dplyr::select(seizure_id, RIE_weight, sz_adj, origin_country_code:transit_4_country_code) %>%      # don't include destination country in seizures-out
    filter_all(any_vars(grepl(cty, .))) %>%
    distinct(seizure_id, .keep_all = TRUE)
  
  if (nrow(temp) > 0){
    vars$sz_out_1[i] <- sum(filter(temp, RIE_weight < 500)$sz_adj)
    vars$sz_out_2[i] <- sum(filter(temp, RIE_weight >= 500)$sz_adj)
  }
  
  
  # Weights-in
  temp <- seizures %>%
    filter(seizure_year == yr & discovered_country_code == cty) %>%
    distinct(seizure_id, .keep_all = TRUE)
  
  if (nrow(temp) > 0){
    vars$wt_in_1[i] <- sum(filter(temp, RIE_weight < 500)$RIE_adj)
    vars$wt_in_2[i] <- sum(filter(temp, RIE_weight >= 500)$RIE_adj)
  }
  
  # Weights-out... 
  # ...need to take care w.r.t. multiple and single country of origin seizures, and mixed raw and worked seizures:
  # ...recalculate RIE weights after adjusting weights for origin countries with proportion < 100
  temp <- seizures %>%
    filter(seizure_year == yr & discovered_country_code != cty) %>%
    select(seizure_id, ivory_type, raw_weight, worked_weight, proportion, theta, phi, origin_country_code:transit_4_country_code) %>%
    filter(if_any(-ivory_type, ~ grepl(cty, .))) %>%              # Exclude ivory_type column in case country code appears in "raw" or "worked" words
    mutate(raw_weight = ifelse(!is.na(raw_weight) & !is.na(origin_country_code) & origin_country_code == cty,
                               raw_weight * proportion / 100, raw_weight),
           worked_weight = ifelse(!is.na(worked_weight) & !is.na(proportion) & !is.na(origin_country_code) & origin_country_code == cty,
                                  worked_weight * proportion / 100, worked_weight))
  
  if (nrow(temp) > 0){
    RIEs <- temp %>%
      distinct(seizure_id, ivory_type, .keep_all = TRUE) %>%
      group_by(seizure_id) %>%
      mutate(RIE_weight = sum(raw_weight, worked_weight, na.rm = TRUE)) %>%
      dplyr::select(seizure_id, RIE_weight) %>%
      distinct(seizure_id, .keep_all = TRUE) %>%
      ungroup()
    temp <- temp %>%
      left_join(RIEs, by = "seizure_id") %>%
      mutate(RIE_adj = RIE_weight/(theta*phi))
    
    vars$wt_out_1[i] <- sum(filter(temp, RIE_weight < 500)$RIE_adj)
    vars$wt_out_2[i] <- sum(filter(temp, RIE_weight >= 500)$RIE_adj)
  }
  
}

# Sum variables over CoP Period
vars <- vars %>%
  group_by(country) %>%
  summarise(across(everything(), sum)) %>%
  select(-year)



#### Combine all cluster analysis variables together and save ------------------
cluster_vars <- mean_lambdas %>%
  left_join(vars, by = "country")

write_csv(cluster_vars, file = paste0("Processed Data/cluster_vars_", analysis_name, ".csv"))




