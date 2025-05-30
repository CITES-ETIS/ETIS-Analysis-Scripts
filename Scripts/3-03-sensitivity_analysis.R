################################################################################
# - Sensitivity analysis for pairwise clustering similarities:                 #
# - Runs the cluster analysis for N posterior samples of the lambdas and       #
#   bias-adjustment, stores the dendrogram merge heights of each pair of       #
#   countries from each run, and prepares a heatmap figure to show the mean    #
#   merge heights (lower merge height = greater similarity) between each pair  #
#   of countries.                                                              #
# - This script uses objects and variables established in scripts 3-01 and     #
#   3-02, so ensure that these have been run first.                            #
#                                                                              #
################################################################################


# Load seizures table
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_with_est_weights_", analysis_name, ".csv")) %>%
  filter(seizure_year %in% CoP_Period)

# Fill in missing proportions when multiple countries of origin
seizures <- seizures %>%
  group_by(seizure_id, ivory_type) %>%
  mutate(n_origins = n()) %>% 
  mutate(proportion = ifelse(is.na(proportion), 100/n_origins, proportion)) %>%
  ungroup()

# Create RIEs table
RIEs <- seizures %>%
  distinct(seizure_id, ivory_type, .keep_all = TRUE) %>%
  group_by(seizure_id) %>%
  mutate(RIE_weight = sum(raw_weight, worked_weight, na.rm = TRUE)) %>%
  select(seizure_id, RIE_weight) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  ungroup() 

# Get list of countries (in order of phi and theta arrays, same as script 3-01)
countries <- unique(covariates$country)

# Trim lambda (from script 3-01) to CoP_Period
lambda <- lambda[, which(yearfrom:yearto %in% CoP_Period), , ]

# Sum lambdas over the CoP period years
lambda <- apply(lambda, c(1,3,4), sum, na.rm = TRUE)

# Set up empty array to store the merge heights between each pair of countries
merge_heights_mat <- matrix(0, dim(X)[1], dim(X)[1])    # X is the cluster analysis data matrix from script 3-02. dim(X)[1] is the number of countries in the cluster analysis

# Set number of cluster analysis iterations to run
N <- 1000 

# Run cluster analysis iterations
for (i in 1:N){  
  
  # Print progress statement every 10 iterations
  if (i %% 10 == 0){
    print(paste0("Completed ", i, " out of ", N, " iterations"))
  }
  
  #### Prepare variables
  
  # Format lambdas into dataframe
  lambdas.i <- as.data.frame(lambda[, , i]) %>%
    mutate(country = inclusion.cts) %>%
    relocate(country) %>%
    rename_with(~ paste0("lambda_", seq_along(.)), -country)
  
  
  # Phi and theta for iteration i
  phi.i <- as.data.frame(phi_array[, , i]) %>%
    mutate(country = countries)
  colnames(phi.i)[1:n.yr] <- CoP_Period
  
  theta.i <- as.data.frame(theta_array[, , i]) %>%
    mutate(country = countries)
  colnames(theta.i)[1:n.yr] <- CoP_Period
  
  # Pivot to long format
  adjs <- pivot_longer(theta.i, -country, names_to = "year", values_to = "theta") %>%
    left_join(pivot_longer(phi.i, -country, names_to = "year", values_to = "phi")) %>%
    mutate(year = as.numeric(year))
  
  # Add adjs to seizures table and mutate bias-adjusted seizure and RIE weight columns
  seizures.i <- seizures %>%
    left_join(adjs, join_by("seizure_year" == "year", "discovered_country_code" == "country")) %>%
    left_join(RIEs, by = "seizure_id") %>%
    mutate(sz_adj = 1/(theta*phi),
           RIE_adj = RIE_weight/(theta*phi))
  
  # Initialise df to store country-year variables
  vars.i <- mean_adjs %>%
    dplyr::select(country, year) %>%
    mutate(sz_out_1 = 0,
           sz_out_2 = 0,
           wt_in_1 = 0,
           wt_in_2 = 0,
           wt_out_1 = 0,
           wt_out_2 = 0)
  
  # Calculate the six bias-adjusted seizure variables for each country and year
  for (j in 1:nrow(vars.i)){
    cty <- vars.i$country[j]
    yr <- vars.i$year[j]
    
    # Seizures-out
    temp <- seizures.i %>%
      filter(seizure_year == yr & discovered_country_code != cty) %>%
      select(seizure_id, RIE_weight, sz_adj, origin_country_code:transit_4_country_code) %>%
      filter_all(any_vars(grepl(cty, .))) %>%
      distinct(seizure_id, .keep_all = TRUE)
    
    if (nrow(temp) > 0){
      vars.i$sz_out_1[j] <- sum(filter(temp, RIE_weight < 500)$sz_adj)
      vars.i$sz_out_2[j] <- sum(filter(temp, RIE_weight >= 500)$sz_adj)
    }
    
    
    # Weights-in
    temp <- seizures.i %>%
      filter(seizure_year == yr & discovered_country_code == cty) %>%
      distinct(seizure_id, .keep_all = TRUE)
    
    if (nrow(temp) > 0){
      vars.i$wt_in_1[j] <- sum(filter(temp, RIE_weight < 500)$RIE_adj)
      vars.i$wt_in_2[j] <- sum(filter(temp, RIE_weight >= 500)$RIE_adj)
    }
    
    # Weights-out... 
    # ...need to take care w.r.t. multiple and single country of origin seizures, and mixed raw and worked seizures:
    # ...recalculate RIE weights after adjusting weights for origin countries with proportion < 100
    temp <- seizures.i %>%
      filter(seizure_year == yr & discovered_country_code != cty) %>%
      select(seizure_id, ivory_type, raw_weight, worked_weight, proportion, theta, phi, origin_country_code:transit_4_country_code) %>%
      filter(if_any(-ivory_type, ~ grepl(cty, .))) %>%
      mutate(raw_weight = ifelse(!is.na(raw_weight) & !is.na(proportion) & origin_country_code == cty,
                                 raw_weight * proportion / 100, raw_weight),
             worked_weight = ifelse(!is.na(worked_weight) & !is.na(proportion) & origin_country_code == cty,
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
      
      vars.i$wt_out_1[j] <- sum(filter(temp, RIE_weight < 500)$RIE_adj)
      vars.i$wt_out_2[j] <- sum(filter(temp, RIE_weight >= 500)$RIE_adj)
    }
    
  }
  
  # Sum variables over CoP Period
  vars.i <- vars.i %>%
    group_by(country) %>%
    summarise(across(everything(), sum)) %>%
    dplyr::select(-year)
  
  # Combine all cluster variables together
  cluster_vars.i <- lambdas.i %>%
    left_join(vars.i, by = "country")
  
  
  #### Cluster analysis
  
  # Remove excl_ctys, and log transform all variables
  cluster_vars.i <- cluster_vars.i %>%
    filter(!(country %in% excl_ctys)) %>%
    mutate(across(where(is.numeric), ~ log(.x + 1)))
  
  # Get variables for clustering
  X.i <- select(cluster_vars.i, -country)
  
  # Run agglomerative hierarchical clustering
  cl.i <- agnes(X.i, stand = TRUE, method = "ward")
  
  # Add pairwise merge heights to the array
  merge_heights_mat <- merge_heights_mat + as.matrix(cophenetic(cl.i))
  
}

# Divide by number of iterations to get average merge height
merge_heights_ave <- merge_heights_mat/N

# Save as csv
write_csv(as.data.frame(merge_heights_ave), file = paste0("Processed Data/merge_heights_ave_", analysis_name, ".csv"))

# Reorder to the country order from the main cluster analysis (cty_order from script 3-02)
merge_heights_ave <- merge_heights_ave[cty_order, cty_order]


# Reformat data for heatmap
sensitivity_df <- as.data.frame(merge_heights_ave) %>%
  mutate(country.x = cluster_vars$country[cty_order]) %>%
  relocate(country.x)

colnames(sensitivity_df)[2:ncol(sensitivity_df)] <- as.character(cluster_vars$country[cty_order])

sensitivity_df <- sensitivity_df %>%
  pivot_longer(cols = -country.x, names_to = "country.y", values_to = "value")

sensitivity_df$country.y <- factor(sensitivity_df$country.y, levels = levels(factor(sensitivity_df$country.y))[cty_order])


# Heatmap plot
sensitivity_heatmap <- ggplot(sensitivity_df, aes(x = country.x, y = country.y, fill = value)) +
  geom_tile(color = "#CECACD", linewidth=0.1) +
  scale_fill_gradientn(colors = c("black","#7E777C","white"), name = "Average\ncophenetic\ndistance") +  
  ylab(" ") + xlab("") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8, hjust=0.5, vjust=0.5 ,angle=90),  
        plot.title = element_text(size=12, hjust = 0.5, face="bold"),
        axis.line = element_line(colour = "grey"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank())
