################################################################################
# - Build and save weight estimation models                                    #
#                                                                              #
################################################################################


# Load seizures table
seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_", analysis_name , ".csv"))

# Seizures to remove from weight model training sets
sz.rm <- c( 
  111934, # Number of pieces (1) seems too small for weight (116). Worked ivory. Remove from weight model training set
  31881, # Ivory comment indicated that weights and pieces are for raw and worked combined. Remove from weight model training set
  103825 # Weight and pieces inputted by ETIS for raw ivory, although documentation to support this is not found. Remove from weight model training set.
)

# Set knots for spline regression weight models
knots <- c(2,4,6)


#### Raw weight model ----------------------------------------------------------

raw_estimation <- seizures %>%
  filter(ivory_type == "raw" & seizure_year >= 2008) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  filter(!is.na(raw_pieces) & is.na(raw_weight))

raw_estimation_max <- range(raw_estimation$raw_pieces)[2]

raw_training <- seizures %>%
  filter(ivory_type == "raw") %>%
  filter(seizure_year >= 2008 & seizure_year <= yearto) %>%
  filter(!(seizure_id %in% sz.rm)) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  select(pieces = raw_pieces, weight = raw_weight) %>%
  filter(!is.na(pieces) & !is.na(weight)) %>%
  filter(pieces <= raw_estimation_max) %>%
  mutate(x = log(pieces + 1))

r0 <- lm(weight ~ ns(x, knots = knots), data = raw_training)
xx <- boxcox(r0, lambda = seq(0.01, 0.15, length = 50), plotit = FALSE)
lambda.r <- xx$x[xx$y == max(xx$y)]  
raw_training <- raw_training %>%
  mutate(y = ((weight ^ lambda.r) - 1) / lambda.r)

raw.ns <- lm(y ~ ns(x, knots = knots), data = raw_training)




#### Worked weight model -------------------------------------------------------

wkd_estimation <- seizures %>%
  filter(ivory_type == "worked" & seizure_year >= 2008) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  filter(!is.na(worked_pieces) & is.na(worked_weight))

wkd_estimation_max <- range(wkd_estimation$worked_pieces)[2]

jewellery_strings <- c("jewellery", "jewelry", "bracelet", "necklace", " ring", "earring", "pendant", "bangle")

wkd_training <- seizures %>%
  filter(ivory_type == "worked") %>%
  filter(seizure_year >= 2008 & seizure_year <= yearto) %>%
  filter(!(seizure_id %in% sz.rm)) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  filter(!is.na(worked_pieces) & !is.na(worked_weight)) %>%
  filter(worked_pieces <= wkd_estimation_max) %>%
  mutate(jewellery = ifelse(str_detect(ivory_comment, regex(paste(jewellery_strings, collapse = '|'), ignore_case = T)) == TRUE, 1, 0)) %>%
  replace_na(list(jewellery = 0)) %>%
  select(pieces = worked_pieces, weight = worked_weight, jewellery) %>%
  mutate(x = log(pieces + 1))
 
w0 <- lm(weight ~ ns(x, knots = knots) + as.factor(jewellery), data = wkd_training)
xx <- boxcox(w0, lambda = seq(0.01, 0.15, length = 50), plotit = FALSE)
lambda.w <- xx$x[xx$y == max(xx$y)]
wkd_training <- wkd_training %>%
  mutate(y = ((weight ^ lambda.w) - 1) / lambda.w)

wkd.ns <- lm(y ~ ns(x, knots = knots) + as.factor(jewellery), data = wkd_training)



#### Save weight estimation models ---------------------------------------------
save(raw.ns, wkd.ns, lambda.r, lambda.w, file = paste0("Processed Data/wt est models_", analysis_name, ".Rdata"))


