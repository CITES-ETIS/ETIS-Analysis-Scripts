################################################################################
# - Generates seizures table with necessary columns                            #
# - Excludes non-ivory seizures                                                #
# - Amends and removes seizures as necessary                                   #
#                                                                              #
################################################################################


# Import data from csvs downloaded from ETIS Online
seizures <- read_csv("Original Data/seizures.csv")

seizure_origin_countries <- read_csv("Original Data/seizure_origin_countries.csv")
seizure_raw_origin_countries <- read_csv("Original Data/seizure_raw_origin_countries.csv")
seizure_worked_origin_countries <- read_csv("Original Data/seizure_worked_origin_countries.csv")

seizure_export_countries <- read_csv("Original Data/seizure_export_countries.csv")

seizure_transit_countries <- read_csv("Original Data/seizure_transit_countries.csv")

countries <- read_csv("Original Data/countries.csv") %>%
  replace_na(list(code = "NAM"))       # replace NA with NAM for Namibia

# Additional MA-submitted TRUEs - from ma_submitted_override and ma_approved status note
seizures <- seizures %>%
  mutate(ma_submitted = ifelse(!is.na(ma_submitted_override), TRUE, ma_submitted)) %>%
  mutate(ma_submitted = ifelse(!is.na(status_note) & status_note == "ma_approved", TRUE, ma_submitted))


#### Set up seizures lists -----------------------------------------------------
seizures_wide <- seizures %>%
  filter(status_id == 4) %>%           # Only include verified seizures
  select(seizure_id = id, seizure_year, 
         discovered_country_id, destination_country_id,
         raw_pieces, raw_weight, raw_present_amount_unknown, 
         worked_pieces, worked_weight, worked_present_amount_unknown,
         ivory_comment, ma_submitted)

# Split into raw and worked
seizures_raw <- seizures_wide %>%
  filter(!is.na(raw_weight) | !is.na(raw_pieces) | raw_present_amount_unknown == TRUE) %>%
  select(-worked_pieces, -worked_weight, -worked_present_amount_unknown) %>%
  mutate(ivory_type = "raw")

seizures_worked <- seizures_wide %>%
  filter(!is.na(worked_weight) | !is.na(worked_pieces) | worked_present_amount_unknown == TRUE) %>%
  select(-raw_pieces, -raw_weight, -raw_present_amount_unknown) %>%
  mutate(ivory_type = "worked")


#### Origin countries ----------------------------------------------------------
# Check in specific raw/worked origin spreadsheets first. If not found, check in combined origins spreadsheet

# Raw
seizures_raw <- seizures_raw %>%
  left_join(seizure_raw_origin_countries, by="seizure_id", multiple="all") %>%
  rename(origin_country_id = country_id) %>%
  select(-id, -created_at, -updated_at)

seizures_raw_missing_origin <- seizures_raw %>%
  filter(is.na(origin_country_id)) %>%
  select(-origin_country_id, -proportion) %>%
  left_join(seizure_origin_countries, by="seizure_id", multiple="all") %>%
  rename(origin_country_id = country_id) %>%
  select(-id, -created_at, -updated_at)

seizures_raw <- seizures_raw %>%
  filter(!is.na(origin_country_id)) %>%
  bind_rows(seizures_raw_missing_origin)

# Worked
seizures_worked <- seizures_worked %>%
  left_join(seizure_worked_origin_countries, by="seizure_id", multiple="all") %>%
  rename(origin_country_id = country_id) %>%
  select(-id, -created_at, -updated_at)

seizures_worked_missing_origin <- seizures_worked %>%
  filter(is.na(origin_country_id)) %>%
  select(-origin_country_id, -proportion) %>%
  left_join(seizure_origin_countries, by="seizure_id", multiple="all") %>%
  rename(origin_country_id = country_id) %>%
  select(-id, -created_at, -updated_at)

seizures_worked <- seizures_worked %>%
  filter(!is.na(origin_country_id)) %>%
  bind_rows(seizures_worked_missing_origin)

# Recombine - now excludes non-ivory as all seizures are either raw or worked
seizures_wide <- bind_rows(seizures_raw, seizures_worked)


#### Export countries ----------------------------------------------------------
# Widen into export1, export2 etc.
export_countries_wide <- seizure_export_countries %>%
  group_by(seizure_id) %>% 
  mutate(position = paste0("export_", row_number(), "_country_id")) %>%
  pivot_wider(id_cols=seizure_id, names_from=position, values_from=country_id)

# Add export countries to seizures table
seizures_wide <- seizures_wide %>%
  left_join(export_countries_wide, by="seizure_id")


#### Transit countries ---------------------------------------------------------
# Widen into transit1, transit2 etc.
transit_countries_wide <- seizure_transit_countries %>%
  group_by(seizure_id) %>%
  mutate(position = paste0("transit_", row_number(), "_country_id")) %>%
  pivot_wider(id_cols=seizure_id, names_from=position, values_from=country_id)

# Add transit countries to seizure table
seizures_wide <- seizures_wide %>%
  left_join(transit_countries_wide, by="seizure_id") 


#### TIDYING UP ----------------------------------------------------------------
# Reorder columns and rows
seizures_wide <- seizures_wide %>%
  relocate(destination_country_id, .after = last_col()) %>%
  relocate(ivory_type, .after = worked_present_amount_unknown) %>%
  relocate(discovered_country_id, .after = ivory_type) %>%
  relocate(origin_country_id, .after = discovered_country_id) %>%
  arrange(seizure_year, seizure_id, ivory_type) %>%
  distinct()

# Replace country ids with codes
ids_and_codes <- countries %>%
  select(id, code)

col_list <- seizures_wide %>%
  select(discovered_country_id:destination_country_id) %>%
  colnames()

for(col in col_list){
  seizures_wide <- seizures_wide %>%
    left_join(ids_and_codes, join_by(!!sym(col) == id)) %>%
    rename_with(~str_replace(col, "id", "code"), code) %>%
    select(-all_of(col))
}

# Remove 'region' countries (codes beginning with X)
seizures_wide <- seizures_wide %>%
  mutate(across(discovered_country_code:destination_country_code, ~ ifelse(str_detect(.x, "^X"), NA, .x)))


# Need to exclude records which are (not MA-submitted or which implicate another Party) and haven't been through a data validation cycle: 
# Records created after last data validation cycle not submitted by MA
# - need to go through validation before inclusion in analysis
sz_rm_validation <- seizures %>%
  filter(created_at > data_val_notification & ma_submitted == FALSE) %>%
  select(id, status_id)

# MA-submitted records created after last data validation cycle which implicate another party
# - can't be used in analysis as the implicated party needs a chance to validate
seizures <- seizures %>%
  left_join(select(seizures_wide, seizure_id, discovered_country_code:destination_country_code),
            join_by("id" == "seizure_id"))

sz_rm_validation_implicated <- seizures %>%
  filter(created_at > data_val_notification & ma_submitted == TRUE) %>%
  filter((!is.na(origin_country_code) & origin_country_code != discovered_country_code) |
           (!is.na(export_1_country_code) & export_1_country_code != discovered_country_code) |
           (!is.na(export_2_country_code) & export_2_country_code != discovered_country_code) |
           (!is.na(transit_1_country_code) & transit_1_country_code != discovered_country_code) |
           (!is.na(transit_2_country_code) & transit_2_country_code != discovered_country_code) |
           (!is.na(transit_3_country_code) & transit_3_country_code != discovered_country_code) |
           (!is.na(transit_4_country_code) & transit_4_country_code != discovered_country_code) |
           (!is.na(destination_country_code) & destination_country_code != discovered_country_code)) %>%
  select(id, status_id)

# Seizure ids to exclude. sz_ids_remove is created in script 1 and can be specific to the current analysis.
# ..the data validation-related ids will always be excluded
sz_rm <- c(sz_ids_remove, sz_rm_validation$id, sz_rm_validation_implicated$id)

# Remove seizure ids from seizures table
seizures_wide <- seizures_wide %>%
  filter(!(seizure_id %in% sz_rm))


#### Amend certain seizures as needed
# seizure id 115041 - year to change from 2020 to 2014
# seizure id 111801 worked component - year to change from 2016 to 2023
# seizure id 107677 - worked component, weight is only given for one piece although quoted as if for 71. Remove weight value so that this will instead be estimated from the weight model
# seizure id 118270 - confirmed not to include KM in the trade route
# seizure id 117854 - 200 raw pieces, but likely small teeth so weight estimation would result in overestimate. Remove from analysis. See correspondence with IN attached to record on ETIS Online
seizures_wide <- seizures_wide %>%
  mutate(seizure_year = ifelse(seizure_id == 115041, 2014, seizure_year)) %>%
  mutate(seizure_year = ifelse(seizure_id == 111801 & ivory_type == "worked", 2023, seizure_year)) %>%
  mutate(worked_weight = ifelse(seizure_id == 107677, NA, worked_weight)) %>%
  mutate(transit_1_country_code = ifelse(seizure_id == 118270, NA, transit_1_country_code)) %>%
  filter(seizure_id != 117854)


#### Save seizures table -------------------------------------------------------
write_csv(seizures_wide, file=paste0("Processed Data/seizures_trade_routes_", analysis_name, ".csv"))

