################################################################################
# - Set up covariates for seizure and reporting rate modelling                 #
#                                                                              #
################################################################################

# Load seizures table (don't need table with estimated weights for this script)
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_", analysis_name, ".csv"))

countries <- read_csv("Original Data/countries.csv") %>%
  replace_na(list(code = "NAM"))       # replace NA with NAM for Namibia

ids_and_codes <- countries %>%
  select(id, code)

# Create a dataframe with all combinations of years and countries - covariates will be calculated for each combination
all_combinations <- expand_grid(seizure_year = unique(seizures$seizure_year), discovered_country_code = unique(countries$code))


#### Reporting rate covariate columns ------------------------------------------

# By country and year, count total seizures and total MA-submitted seizures
ctry_year_covars <- seizures %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  group_by(seizure_year, discovered_country_code) %>%
  summarise(total_seizures = n(), total_ma_submitted = sum(ma_submitted), .groups = "drop") %>%
  right_join(all_combinations, by = c("seizure_year", "discovered_country_code")) %>%
  mutate(total_seizures = replace_na(total_seizures, 0),
         total_ma_submitted = replace_na(total_ma_submitted, 0)) %>%
  rename(year = seizure_year, country = discovered_country_code)
  

# Get CITES reporting scores - from processed table created from CITES annual reports file
CITESscores <- read_csv("Processed Data/CITESscores.csv")

ctry_year_covars <- ctry_year_covars %>%
  left_join(CITESscores, join_by(year, country))


# Calculate ETIS reporting score and CITES reporting score 
ctry_year_covars <- ctry_year_covars %>%
  mutate(ETIS.rep = total_ma_submitted / total_seizures) %>%
  mutate(CITES.rep.logit = emp.logit(CITES_reports, CITES_years)) %>%
  replace_na(list(ETIS.rep = 0)) %>%
  select(year, country, ETIS.rep, CITES.rep.logit) 


#### Seizure rate covariate columns --------------------------------------------

# Get seizure-in and seizure-out variables for LE ratio
LE_df <- seizures %>%
  # Seizure-in
  distinct(seizure_id, .keep_all = TRUE) %>%            # distinct seizure id, so as not to double count a seizure if it contains raw and worked, or has multiple countries of origin
  group_by(seizure_year, discovered_country_code) %>%   # group by year and discovery country
  summarise(sz_in = n(), .groups = "drop") %>%          # count seizures-in
  right_join(all_combinations, by = c("seizure_year", "discovered_country_code")) %>%
  replace_na(list(sz_in = 0)) %>%                       # replace NAs (for year-country combinations from all_combinations which didn't appear in seizures table) with 0
  rename(year = seizure_year, country = discovered_country_code) %>%
  arrange(year, country) %>%
  # Seizure-out
  rowwise() %>%                                         # for each country and year combination...
  mutate(sz_out = nrow(seizures %>%                     # ...filter seizures table to that year and not that discovery country, and count occurrences of that country in the other (non-destination) trade route columns
                         filter(seizure_year == year & discovered_country_code != country) %>%
                         select(seizure_id, origin_country_code:transit_4_country_code) %>%
                         filter_all(any_vars(grepl(country, .))) %>%
                         distinct(seizure_id)))


# Trade Chain Index (TCI):
# TCI = log((dest.score + 1) / (nondest.score + 1))
#
# - dest.score for country x in year y is the number of times country x is listed as a
#   country of destination divided by the number of seizures that list a country of destination.
#   dest.score is calculated separately for raw and worked ivory and then averaged.
#
# - nondest.score is the average of the origin, export and transit scores, which are 
#   calculated in the same way as the dest.score.
#

# Add dest and nondest scores to LE_df table
LE_df <- LE_df %>%
  rowwise() %>%
  # numerator of destination score for raw
  mutate(dest.num.r = nrow(seizures %>%
                             filter(ivory_type == "raw") %>%
                             filter(seizure_year == year & destination_country_code == country) %>%
                             distinct(seizure_id))) %>%
  # denominator of destination score for raw
  mutate(dest.denom.r = nrow(seizures %>%
                             filter(ivory_type == "raw") %>%
                             filter(seizure_year == year & !is.na(destination_country_code)) %>%
                             distinct(seizure_id))) %>%
  # destination score for raw
  mutate(dest.score.r = dest.num.r / dest.denom.r) %>%
  # numerator of destination score for worked
  mutate(dest.num.w = nrow(seizures %>%
                             filter(ivory_type == "worked") %>%
                             filter(seizure_year == year & destination_country_code == country) %>%
                             distinct(seizure_id))) %>%
  # denominator of destination score for worked
  mutate(dest.denom.w = nrow(seizures %>%
                               filter(ivory_type == "worked") %>%
                               filter(seizure_year == year & !is.na(destination_country_code)) %>%
                               distinct(seizure_id))) %>%
  # destination score for worked
  mutate(dest.score.w = dest.num.w / dest.denom.w) %>%
  # numerator of origin score for raw
  mutate(org.num.r = nrow(seizures %>%
                             filter(ivory_type == "raw") %>%
                             filter(seizure_year == year & origin_country_code == country) %>%
                             distinct(seizure_id))) %>%
  # denominator of origin score for raw
  mutate(org.denom.r = nrow(seizures %>%
                               filter(ivory_type == "raw") %>%
                               filter(seizure_year == year & !is.na(origin_country_code)) %>%
                               distinct(seizure_id))) %>%
  # origin score for raw
  mutate(org.score.r = org.num.r / org.denom.r) %>%
  # numerator of origin score for worked
  mutate(org.num.w = nrow(seizures %>%
                            filter(ivory_type == "worked") %>%
                            filter(seizure_year == year & origin_country_code == country) %>%
                            distinct(seizure_id))) %>%
  # denominator of origin score for worked
  mutate(org.denom.w = nrow(seizures %>%
                              filter(ivory_type == "worked") %>%
                              filter(seizure_year == year & !is.na(origin_country_code)) %>%
                              distinct(seizure_id))) %>%
  # origin score for worked
  mutate(org.score.w = org.num.w / org.denom.w) %>%
  # numerator of export score for raw
  mutate(exp.num.r = nrow(seizures %>%
                            filter(ivory_type == "raw") %>%
                            filter(seizure_year == year & (export_1_country_code == country | export_2_country_code == country)) %>%
                            distinct(seizure_id))) %>%
  # denominator of export score for raw
  mutate(exp.denom.r = nrow(seizures %>%
                              filter(ivory_type == "raw") %>%
                              filter(seizure_year == year & (!is.na(export_1_country_code) | !is.na(export_1_country_code))) %>%
                              distinct(seizure_id))) %>%
  # export score for raw
  mutate(exp.score.r = exp.num.r / exp.denom.r) %>%
  # numerator of export score for worked
  mutate(exp.num.w = nrow(seizures %>%
                            filter(ivory_type == "worked") %>%
                            filter(seizure_year == year & (export_1_country_code == country | export_2_country_code == country)) %>%
                            distinct(seizure_id))) %>%
  # denominator of export score for worked
  mutate(exp.denom.w = nrow(seizures %>%
                              filter(ivory_type == "worked") %>%
                              filter(seizure_year == year & (!is.na(export_1_country_code) | !is.na(export_1_country_code))) %>%
                              distinct(seizure_id))) %>%
  # export score for worked
  mutate(exp.score.w = exp.num.w / exp.denom.w) %>%
  # numerator of transit score for raw
  mutate(tra.num.r = nrow(seizures %>%
                            filter(ivory_type == "raw") %>%
                            filter(seizure_year == year & (transit_1_country_code == country | transit_2_country_code == country | transit_3_country_code == country | transit_4_country_code == country)) %>%
                            distinct(seizure_id))) %>%
  # denominator of transit score for raw
  mutate(tra.denom.r = nrow(seizures %>%
                              filter(ivory_type == "raw") %>%
                              filter(seizure_year == year & !is.na(transit_1_country_code)) %>%
                              distinct(seizure_id))) %>%
  # transit score for raw
  mutate(tra.score.r = tra.num.r / tra.denom.r) %>%
  # numerator of transit score for worked
  mutate(tra.num.w = nrow(seizures %>%
                            filter(ivory_type == "worked") %>%
                            filter(seizure_year == year & (transit_1_country_code == country | transit_2_country_code == country | transit_3_country_code == country | transit_4_country_code == country)) %>%
                            distinct(seizure_id))) %>%
  # denominator of transit score for worked
  mutate(tra.denom.w = nrow(seizures %>%
                              filter(ivory_type == "worked") %>%
                              filter(seizure_year == year & !is.na(transit_1_country_code)) %>%
                              distinct(seizure_id))) %>%
  # transit score for worked
  mutate(tra.score.w = tra.num.w / tra.denom.w) %>%
  # final destination, origin, export and transit scores
  mutate(dest.score = mean(c(dest.score.r, dest.denom.w))) %>%
  mutate(org.score = mean(c(org.score.r, org.denom.w))) %>%
  mutate(exp.score = mean(c(exp.score.r, exp.denom.w))) %>%
  mutate(tra.score = mean(c(tra.score.r, tra.denom.w))) %>%
  mutate(nondest.score = mean(c(org.score, exp.score, tra.score)))

# Generate LE1 (1-year lagged LE ratio) and TCI and remove other columns
LE_df <- LE_df %>%
  mutate(LE_ratio = sz_in / (sz_in + sz_out)) %>%
  replace_na(list(LE_ratio = 0)) %>%
  mutate(TCI = log((dest.score + 1) / (nondest.score + 1))) %>%
  group_by(country) %>%
  mutate(LE1 = lag(LE_ratio)) %>%       # Note LE_df is arranged by year
  select(year, country, LE1, TCI)




#### Combine reporting and seizure rate covariates and save --------------------
df_covars <- ctry_year_covars %>%
  left_join(LE_df, join_by(year, country))

write_csv(df_covars, file=paste0("Processed Data/covars_", analysis_name, ".csv"))

