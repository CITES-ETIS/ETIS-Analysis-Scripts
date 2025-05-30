################################################################################
# - Estimate missing weights from saved weight estimation model                #
# - Convert worked weights to RIE                                              #
#                                                                              #
################################################################################


# Load seizures table
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_", analysis_name, ".csv"))

# Load weight models
load(paste0("Processed Data/wt est models_", analysis_name, ".Rdata"))


### Raw weight estimation ------------------------------------------------------

raw_estimation <- seizures %>%
  filter(ivory_type == "raw") %>%   
  distinct(seizure_id, .keep_all = TRUE) %>%
  filter(!is.na(raw_pieces) & is.na(raw_weight)) %>%
  select(seizure_id, pieces = raw_pieces, weight = raw_weight) %>%
  mutate(x = log(pieces + 1))

raw_preds <- predict(raw.ns, newdata = raw_estimation, se.fit = TRUE)
raw_estimation$y <- raw_preds$fit
raw_estimation$sigma <- raw_preds$se.fit
raw_estimation <- raw_estimation %>%
  mutate(weight = ((lambda.r*y + 1)^(1/lambda.r)) * (1 + (((sigma^2) * (1 - sigma))/(2 * ((lambda.r*y + 1)^2)))))

#### Worked weight estimation --------------------------------------------------

jewellery_strings <- c("jewellery", "jewelry", "bracelet", "necklace", " ring", "earring", "pendant", "bangle")

wkd_estimation <- seizures %>%
  filter(ivory_type == "worked") %>%  
  distinct(seizure_id, .keep_all = TRUE) %>%
  filter(!is.na(worked_pieces) & is.na(worked_weight)) %>%
  mutate(jewellery = ifelse(str_detect(ivory_comment, regex(paste(jewellery_strings, collapse = '|'), ignore_case = T)) == TRUE, 1, 0)) %>%
  replace_na(list(jewellery = 0)) %>%
  select(seizure_id, pieces = worked_pieces, weight = worked_weight, jewellery) %>%
  mutate(x = log(pieces + 1))

wkd_preds <- predict(wkd.ns, newdata = wkd_estimation, se.fit = TRUE)
wkd_estimation$y <- wkd_preds$fit
wkd_estimation$sigma <- wkd_preds$se.fit
wkd_estimation <- wkd_estimation %>%
  mutate(weight = ((lambda.w*y + 1)^(1/lambda.w)) * (1 + (((sigma^2) * (1 - sigma))/(2 * ((lambda.w*y + 1)^2)))))


#### Fill in seizures tables with estimated weights ----------------------------

seizures_raw <- seizures %>%
  filter(ivory_type == "raw") %>%
  left_join(select(raw_estimation, seizure_id, weight), by = "seizure_id") %>%
  mutate(raw_weight = coalesce(raw_weight, weight)) %>%
  select(-weight)

seizures_wkd <- seizures %>%
  filter(ivory_type == "worked") %>%
  left_join(select(wkd_estimation, seizure_id, weight), by = "seizure_id") %>%
  mutate(worked_weight = coalesce(worked_weight, weight)) %>%
  mutate(worked_weight = worked_weight/0.7) %>%          # RIE for worked weight
  select(-weight)

seizures <- seizures_raw %>%
  bind_rows(seizures_wkd) %>%
  arrange(seizure_year, seizure_id, ivory_type)

#### Save seizures table -------------------------------------------------------
write_csv(seizures, file=paste0("Processed Data/seizures_trade_routes_with_est_weights_", analysis_name, ".csv"))

