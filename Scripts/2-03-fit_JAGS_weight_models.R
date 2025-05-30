################################################################################
# - Fit models to obtain distribution of raw and worked weights                #
# - Only uses seizures for which raw weight is actually given                  # 
# - Includes a constraint that weight can't be more than 10 tonnes             #
#                                                                              #
################################################################################



# Upper weight limit
max.wt <- 10000


seizures <- read_csv(paste0("Processed Data/seizures_trade_routes_", analysis_name, ".csv"))

# These seizures are listed as worked ivory, but known to be a total weight of raw and worked, so remove from worked weight model
sz.rm <- c(106627, 109674)


# Define the model
mod.txt <- "model {
  for(i in 1:N) {
    l.wgt[i] ~ dt(mu[i], tau.y, nu)T(,max.wt)
    mu[i] <- a0 
  }
  a0 ~ dnorm(0, 1.0E-4)
  nu ~ dunif(1,30)
  tau.y ~ dgamma(0.001, 0.001)
  sigma2.y <- 1/tau.y
} 
"

# Write model to file
writeLines(mod.txt, paste('JAGS MCMC outputs/wt_model_', analysis_name, '.txt', sep = ''))


#### Raw weight model ----------------------------------------------------------

# Filter seizures
df.raw <- seizures %>%
  filter(seizure_year >= yearfrom & seizure_year <= yearto) %>%
  filter(ivory_type == "raw") %>%
  filter(!is.na(raw_weight)) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  dplyr::select(seizure_id, seizure_year, country = discovered_country_code, weight = raw_weight) %>%
  mutate(class = ifelse(weight < 10, 'r1', ifelse(weight < 100, 'r2', 'r3')))

# Model data
mod.raw.dat <- 
  list(N      = dim(df.raw)[1],
       l.wgt  = log(df.raw$weight),
       max.wt = log(max.wt)  
  )

# Inital parameters
mod.inits <- function(){
  list(a0    = rnorm(1, 0, 1),
       tau.y = runif(1, 0, 1),
       nu    = runif(1, 1, 5)
  )
}  

# Model parameters to output
pms <- c('a0', 'sigma2.y', 'nu') 

# Fit model
t1<- Sys.time()
mod.raw <- jags(data = mod.raw.dat,
                inits = mod.inits,
                parameters.to.save = pms,
                model.file = paste('JAGS MCMC outputs/wt_model_', analysis_name, '.txt', sep = ''),
                parallel = T,
                n.cores = 4,
                n.chains = 4,
                n.adapt = 1000,
                n.iter = 20000,
                n.burnin = 10000,
                n.thin = 2)
t2 <- Sys.time()
difftime(t2, t1, units = "hours")

# Save jags model output
saveRDS(mod.raw, file=paste0("JAGS MCMC outputs/raw_wt_model_", analysis_name,  ".rds"))



#### Worked weight model -------------------------------------------------------

# Filter seizures
df.wkd <- seizures %>%
  filter(seizure_year >= yearfrom & seizure_year <= yearto) %>%
  filter(ivory_type == "worked") %>%
  filter(!is.na(worked_weight)) %>%
  filter(!(seizure_id %in% sz.rm)) %>%
  distinct(seizure_id, .keep_all = TRUE) %>%
  mutate(worked_weight = worked_weight / 0.7) %>%       # convert weights to RIE
  dplyr::select(seizure_id, seizure_year, country = discovered_country_code, weight = worked_weight) %>%
  mutate(class = ifelse(weight < 1, 'w1', 'w2'))

# Model data
mod.wkd.dat <- 
  list(N      = dim(df.wkd)[1],
       l.wgt  = log(df.wkd$weight),
       max.wt = log(max.wt)  
  )

# Inital parameters
mod.inits <- function(){
  list(a0    = rnorm(1, 0, 1),
       tau.y = runif(1, 0, 1),
       nu    = runif(1, 1, 5)
  )
}  

# Model parameters to output
pms <- c('a0', 'sigma2.y', 'nu') 

# Fit model
t1<- Sys.time()
mod.wkd <- jags(data = mod.wkd.dat,
                inits = mod.inits,
                parameters.to.save = pms,
                model.file = paste('JAGS MCMC outputs/wt_model_', analysis_name, '.txt', sep = ''),
                parallel = T,
                n.cores = 4,
                n.chains = 4,
                n.adapt = 1000,
                n.iter = 20000,
                n.burnin = 10000,
                n.thin = 2)
t2 <- Sys.time()
difftime(t2, t1, units = "hours")

# Save jags model output
saveRDS(mod.wkd, file=paste0("JAGS MCMC outputs/wkd_wt_model_", analysis_name,  ".rds"))




##### MCMC Diagnostics ---------------------------------------------------------

# mcmcplot(mod.raw$samples)
# mcmcplot(mod.wkd$samples)


