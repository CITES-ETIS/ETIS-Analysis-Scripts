################################################################################
# - Set up seizures data (weight classes) and covariates for analysis          #
# - Inclusion criteria to select countries in trend analysis                   #
#                                                                              #
################################################################################



# Load seizures data with estimated weights
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_with_est_weights_", analysis_name, ".csv"))

# Load covariates data
covariates <- read_csv(paste0("Processed Data/covars_", analysis_name, ".csv"))


#### Country and year weight class seizure totals ------------------------------
seizure_totals <- covariates %>%
  select(year, country) %>%                         # Get country and year columns
  filter(year >= yearfrom & year <= yearto) %>%     # Restrict to years of trend analysis
  rowwise() %>%                                     # For each country and year...
  # ...count seizures in raw small ivory class
  mutate(r1 = nrow(seizures %>%
                     filter(seizure_year == year & discovered_country_code == country) %>%
                     filter(ivory_type == "raw") %>%
                     filter(raw_weight < 10) %>%
                     distinct(seizure_id))) %>%
  # ...count seizures in raw medium ivory class
  mutate(r2 = nrow(seizures %>%
                     filter(seizure_year == year & discovered_country_code == country) %>%
                     filter(ivory_type == "raw") %>%
                     filter(raw_weight >= 10 & raw_weight < 100) %>%
                     distinct(seizure_id))) %>%
  # ...count seizures in raw large ivory class
  mutate(r3 = nrow(seizures %>%
                     filter(seizure_year == year & discovered_country_code == country) %>%
                     filter(ivory_type == "raw") %>%
                     filter(raw_weight >= 100) %>%
                     distinct(seizure_id))) %>%
  # ...count seizures in worked small ivory class
  mutate(w1 = nrow(seizures %>%
                     filter(seizure_year == year & discovered_country_code == country) %>%
                     filter(ivory_type == "worked") %>%
                     filter(worked_weight < 1) %>%
                     distinct(seizure_id))) %>%
  # ...count seizures in worked large ivory class
  mutate(w2 = nrow(seizures %>%
                     filter(seizure_year == year & discovered_country_code == country) %>%
                     filter(ivory_type == "worked") %>%
                     filter(worked_weight >= 1) %>%
                     distinct(seizure_id))) %>%
  ungroup()

#### Country inclusion criteria ------------------------------------------------

# Inclusion criteria considers for each country the counts of seizures -in and -out <10kg, 10-100kg, and >100kg (raw and worked combined, RIE)

# First calculate RIE weight for each seizure (some seizures have a raw and worked component)
RIEs <- seizures %>%
  distinct(seizure_id, ivory_type, .keep_all = TRUE) %>%
  group_by(seizure_id) %>%
  mutate(RIE_weight = sum(raw_weight, worked_weight, na.rm = TRUE)) %>%
  select(seizure_id, RIE_weight) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  ungroup()

# Add RIEs to seizure table
seizures <- seizures %>%
  left_join(RIEs, by = "seizure_id")

# Count small, medium and large seizures -in and -out for each country and year
inclusion_criteria_szs <- covariates %>%
  select(year, country) %>%              # Get country and year columns
  filter(year > yearto - 10) %>%         # Inclusion criteria considers the last 10 years
  rowwise() %>%                          # For each country and year...
  # ...count small seizures-in
  mutate(gp1in = nrow(seizures %>%
                        filter(seizure_year == year & discovered_country_code == country) %>%
                        filter(RIE_weight < 10) %>%
                        distinct(seizure_id))) %>%
  # ...count medium seizures-in
  mutate(gp2in = nrow(seizures %>%
                        filter(seizure_year == year & discovered_country_code == country) %>%
                        filter(RIE_weight >= 10 & RIE_weight < 100) %>%
                        distinct(seizure_id))) %>%
  # ...count large seizures-in
  mutate(gp3in = nrow(seizures %>%
                        filter(seizure_year == year & discovered_country_code == country) %>%
                        filter(RIE_weight >= 100) %>%
                        distinct(seizure_id))) %>%
  # ...count small seizures-out
  mutate(gp1out = nrow(seizures %>%
                        filter(seizure_year == year & discovered_country_code != country) %>%
                        filter(rowSums(select(., origin_country_code:destination_country_code) == country, na.rm = TRUE) >= 1) %>%
                        filter(RIE_weight < 10) %>%
                        distinct(seizure_id))) %>%
  # ...count medium seizures-out
  mutate(gp2out = nrow(seizures %>%
                         filter(seizure_year == year & discovered_country_code != country) %>%
                         filter(rowSums(select(., origin_country_code:destination_country_code) == country, na.rm = TRUE) >= 1) %>%
                         filter(RIE_weight >= 10 & RIE_weight < 100) %>%
                         distinct(seizure_id))) %>%
  # ...count large seizures-out
  mutate(gp3out = nrow(seizures %>%
                         filter(seizure_year == year & discovered_country_code != country) %>%
                         filter(rowSums(select(., origin_country_code:destination_country_code) == country, na.rm = TRUE) >= 1) %>%
                         filter(RIE_weight >= 100) %>%
                         distinct(seizure_id)))

# Calculate inclusion criteria score
inclusion_criteria <- inclusion_criteria_szs %>%
  # Score is 1 * small szs + 10 * medium szs + 100 * large szs
  mutate(score = 1*(gp1in + gp1out) + 10*(gp2in + gp2out) + 100*(gp3in + gp3out)) %>%
  select(-year) %>%
  group_by(country) %>%                                      # For each country,...
  summarise(across(everything(), sum), .groups = 'drop')     # ...sum the score over the 10 years


# Identify inclusion countries as those with a score of at least 100, and save
inclusion.cts <- filter(inclusion_criteria, score >= 100)$country

saveRDS(inclusion.cts, file=paste0('Processed Data/inclusioncts_', analysis_name, '.rds'))



#### Combine seizure totals and covariates for analysis ------------------------

# Filter seizure totals and covariates to countries and years in analysis
df.szs <- seizure_totals %>%
  filter(country %in% inclusion.cts)      # seizure_totals is already filtered to years in analysis

df.covars <- covariates %>%
  filter(year >= yearfrom & year <= yearto) %>%
  filter(country %in% inclusion.cts)

# Combine into one data frame for analysis
df.analysis <- df.szs %>%
  left_join(df.covars, by = c("year", "country")) %>%
  arrange(year, country)

# Save analysis data frame
write_csv(df.analysis, file=paste0("Processed Data/dfanalysis_", analysis_name, ".csv"))

