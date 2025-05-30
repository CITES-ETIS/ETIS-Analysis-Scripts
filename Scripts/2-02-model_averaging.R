################################################################################
# - Model averaging of fitted Bayesian hierarchical models                     #
#                                                                              #
################################################################################


#### Get model weights ---------------------------------------------------------

# Get n values
load(paste0("JAGS MCMC outputs/n values_", analysis_name, ".Rdata"))

# Calculate LOO IC for each model (list of models is specified in wrapper script 2-Trend_Analysis.R)
LOOs <- list(NULL)
for (m in 1:length(models)){
  filename <- models[m]
  mod <- readRDS(file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))
  loglik <- mod$sims.list$LogLik
  dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
  LOOs[[m]] <- loo(loglik)
  rm(mod)
}

# Find model weights based on stacking method
weights <- loo_model_weights(LOOs, method = "stacking")

# Save model weights
saveRDS(weights, file=paste('JAGS MCMC outputs/weights_', analysis_name, '.rds',sep=""))



#### Draw samples from individual models according to model weights ------------

# Set parameters
N <- 40000                            # Choose number of samples for the averaged model
N_models <- round(N * weights)        # Number of samples to draw from each model

# Set up arrays to store samples from averaged model
lambda_ave <- array(NA, dim=c(0, n.ctry * n.yr, 5))
theta_ave <- array(NA, dim=c(0, n.ctry * n.yr, 5))
phi_ave <- array(NA, dim=c(0, n.ctry * n.yr, 5))
r_ave <- array(NA, dim=c(0, 5))
pi_ave <- array(NA, dim=c(0, n.ctry, 5))
b1_ave <- c()                         # b and g coefficients are needed to calculate... 
b2_ave <- c()                         # ...seizure and reporting rate bias-adjustment...
g1_ave <- c()                         # ...for all countries beyond just the trend analysis...
g2_ave <- c()                         # ...inclusion countries during the cluster analysis


# Populate arrays from model samples
set.seed(12345)
for (m in 1:length(models)){
  if (N_models[m] > 0){
    
    mod <- readRDS(paste('JAGS MCMC outputs/', models[m], '.rds',sep=""))
    lambda <- mod$sims.list$lambda
    theta <- mod$sims.list$theta
    phi <- mod$sims.list$phi
    r.cont <- mod$sims.list$r.cont
    pi <- mod$sims.list$pi
    b1 <- mod$sims.list$b1   
    b2 <- mod$sims.list$b2
    g1 <- mod$sims.list$g1
    g2 <- mod$sims.list$g2
    

    sample_m <- sample(1:n.sim, N_models[m], replace = TRUE)           # n.sim is no. samples available from each model
    
    lambda_ave <- abind(lambda_ave, lambda[sample_m, , ], along = 1)
    theta_ave <- abind(theta_ave, theta[sample_m, , ], along = 1)
    phi_ave <- abind(phi_ave, phi[sample_m, , ], along = 1)
    r_ave <- abind(r_ave, r.cont[sample_m, ], along = 1)
    pi_ave <- abind(pi_ave, pi[sample_m, , ], along = 1)
    b1_ave <- c(b1_ave, b1[sample_m])
    b2_ave <- c(b2_ave, b2[sample_m])
    g1_ave <- c(g1_ave, g1[sample_m])
    g2_ave <- c(g2_ave, g2[sample_m])
    
    
    rm(mod)
  }
}


#### Save posterior distributions for averaged model ---------------------------
saveRDS(lambda_ave, file=paste('JAGS MCMC outputs/Processed outputs/lambda_ave_', analysis_name, '.rds', sep=''))
saveRDS(theta_ave, file=paste('JAGS MCMC outputs/Processed outputs/theta_ave_', analysis_name, '.rds', sep=''))
saveRDS(phi_ave, file=paste('JAGS MCMC outputs/Processed outputs/phi_ave_', analysis_name, '.rds', sep=''))
saveRDS(r_ave, file=paste('JAGS MCMC outputs/Processed outputs/r.cont_ave_', analysis_name, '.rds', sep=''))
saveRDS(pi_ave, file=paste('JAGS MCMC outputs/Processed outputs/pi_ave_', analysis_name, '.rds', sep=''))
saveRDS(b1_ave, file=paste('JAGS MCMC outputs/Processed outputs/b1_ave_', analysis_name, '.rds', sep=''))
saveRDS(b2_ave, file=paste('JAGS MCMC outputs/Processed outputs/b2_ave_', analysis_name, '.rds', sep=''))
saveRDS(g1_ave, file=paste('JAGS MCMC outputs/Processed outputs/g1_ave_', analysis_name, '.rds', sep=''))
saveRDS(g2_ave, file=paste('JAGS MCMC outputs/Processed outputs/g2_ave_', analysis_name, '.rds', sep=''))




