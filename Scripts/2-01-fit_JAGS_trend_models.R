################################################################################
# - Fit Bayesian hierarchical models to seizures data                          #
#                                                                              #
################################################################################


# Get seizures and covariates data for analysis
df.analysis <- data.frame(read_csv(paste0("Processed Data/dfanalysis_", analysis_name, ".csv")))



#### Models --------------------------------------------------------------------

# 1 knot -----------------------------------------------------------------------
K <- 1
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod.txt <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
  }


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data

 
} #end model
"

# Write model to file
writeLines(mod.txt, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))


# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv', 
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
             inits = mod.inits,
             parameters.to.save = pms,
             model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
             parallel = TRUE,
             n.cores = 4,
             n.chains = 4,
             n.adapt = 50000,
             n.iter = 200000,
             n.burnin = 100000,
             n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))
saveRDS(mod, file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# 2 knot -----------------------------------------------------------------------
K <- 2
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod.txt <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5] + a6[ctry[i], k] * X[i,6]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
  }


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    a6[j, 1:5] ~ dmnorm(mu6[], Omega6.inv[,])
    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  mu6[1:5] ~ dmnorm(mn[], prec[,])
  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data
  Omega6.inv[1:5, 1:5] ~ dwish(R6[,], 5) # Rs in data
 
} #end model
"

# Write model to file
writeLines(mod.txt, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))


# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       R6     = diag.cov.a6,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a6      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       mu6     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       Omega6.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv', 'Omega6.inv', 
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
             inits = mod.inits,
             parameters.to.save = pms,
             model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
             parallel = TRUE,
             n.cores = 4,
             n.chains = 4,
             n.adapt = 50000,
             n.iter = 200000,
             n.burnin = 100000,
             n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))
saveRDS(mod,file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# 3 knot -----------------------------------------------------------------------
K <- 3
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod0 <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5] + a6[ctry[i], k] * X[i,6]
                         + a7[ctry[i], k] * X[i,7]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
}


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    a6[j, 1:5] ~ dmnorm(mu6[], Omega6.inv[,])
    a7[j, 1:5] ~ dmnorm(mu7[], Omega7.inv[,])
    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  mu6[1:5] ~ dmnorm(mn[], prec[,])
  mu7[1:5] ~ dmnorm(mn[], prec[,])
  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data
  Omega6.inv[1:5, 1:5] ~ dwish(R6[,], 5) # Rs in data
  Omega7.inv[1:5, 1:5] ~ dwish(R7[,], 5) # Rs in data
 
} #end model
"

# Write model to file
writeLines(mod0, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))

# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       R6     = diag.cov.a6,
       R7     = diag.cov.a7,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a6      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a7      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       mu6     = runif(5, -1, 1),
       mu7     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       Omega6.inv = diag(runif(5, 0, 1), 5),
       Omega7.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv', 'Omega6.inv', 'Omega7.inv',
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'mu7', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
             inits = mod.inits,
             parameters.to.save = pms,
             model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
             parallel = TRUE,
             n.cores = 4,
             n.chains = 4,
             n.adapt = 50000,
             n.iter = 200000,
             n.burnin = 100000,
             n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))

saveRDS(mod,file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# 4 knot -----------------------------------------------------------------------
K <- 4
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod0 <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5] + a6[ctry[i], k] * X[i,6]
                         + a7[ctry[i], k] * X[i,7] + a8[ctry[i], k] * X[i,8]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
  }


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    a6[j, 1:5] ~ dmnorm(mu6[], Omega6.inv[,])
    a7[j, 1:5] ~ dmnorm(mu7[], Omega7.inv[,])
    a8[j, 1:5] ~ dmnorm(mu8[], Omega8.inv[,])
    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  mu6[1:5] ~ dmnorm(mn[], prec[,])
  mu7[1:5] ~ dmnorm(mn[], prec[,])
  mu8[1:5] ~ dmnorm(mn[], prec[,])

  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data
  Omega6.inv[1:5, 1:5] ~ dwish(R6[,], 5) # Rs in data
  Omega7.inv[1:5, 1:5] ~ dwish(R7[,], 5) # Rs in data
  Omega8.inv[1:5, 1:5] ~ dwish(R8[,], 5) # Rs in data

 
} #end model
"

# Write model to file
writeLines(mod0, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))


# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       R6     = diag.cov.a6,
       R7     = diag.cov.a7,
       R8     = diag.cov.a8,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a6      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a7      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a8      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       mu6     = runif(5, -1, 1),
       mu7     = runif(5, -1, 1),
       mu8     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       Omega6.inv = diag(runif(5, 0, 1), 5),
       Omega7.inv = diag(runif(5, 0, 1), 5),
       Omega8.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv', 'Omega6.inv', 'Omega7.inv', 'Omega8.inv',
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'mu7', 'mu8', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
            inits = mod.inits,
            parameters.to.save = pms,
            model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
            parallel = TRUE,
            n.cores = 4,
            n.chains = 4,
            n.adapt = 50000,
            n.iter = 200000,
            n.burnin = 100000,
            n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))

saveRDS(mod,file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# 5 knot -----------------------------------------------------------------------
K <- 5
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod0 <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5] + a6[ctry[i], k] * X[i,6]
                         + a7[ctry[i], k] * X[i,7] + a8[ctry[i], k] * X[i,8]
                         + a9[ctry[i], k] * X[i,9]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
  }


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    a6[j, 1:5] ~ dmnorm(mu6[], Omega6.inv[,])
    a7[j, 1:5] ~ dmnorm(mu7[], Omega7.inv[,])
    a8[j, 1:5] ~ dmnorm(mu8[], Omega8.inv[,])
    a9[j, 1:5] ~ dmnorm(mu9[], Omega9.inv[,])
    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  mu6[1:5] ~ dmnorm(mn[], prec[,])
  mu7[1:5] ~ dmnorm(mn[], prec[,])
  mu8[1:5] ~ dmnorm(mn[], prec[,])
  mu9[1:5] ~ dmnorm(mn[], prec[,])

  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data
  Omega6.inv[1:5, 1:5] ~ dwish(R6[,], 5) # Rs in data
  Omega7.inv[1:5, 1:5] ~ dwish(R7[,], 5) # Rs in data
  Omega8.inv[1:5, 1:5] ~ dwish(R8[,], 5) # Rs in data
  Omega9.inv[1:5, 1:5] ~ dwish(R9[,], 5) # Rs in data

 
} #end model
"

# Write model to file
writeLines(mod0, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))


# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       R6     = diag.cov.a6,
       R7     = diag.cov.a7,
       R8     = diag.cov.a8,
       R9     = diag.cov.a9,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a6      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a7      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a8      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a9      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       mu6     = runif(5, -1, 1),
       mu7     = runif(5, -1, 1),
       mu8     = runif(5, -1, 1),
       mu9     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       Omega6.inv = diag(runif(5, 0, 1), 5),
       Omega7.inv = diag(runif(5, 0, 1), 5),
       Omega8.inv = diag(runif(5, 0, 1), 5),
       Omega9.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv',
         'Omega6.inv', 'Omega7.inv', 'Omega8.inv', 'Omega9.inv',
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'mu7', 'mu8', 'mu9', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
            inits = mod.inits,
            parameters.to.save = pms,
            model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
            parallel = TRUE,
            n.cores = 4,
            n.chains = 4,
            n.adapt = 50000,
            n.iter = 200000,
            n.burnin = 100000,
            n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))

saveRDS(mod,file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# 6 knot -----------------------------------------------------------------------
K <- 6
filename <- paste0(K, "knot_", analysis_name)


ctrys <- levels(as.factor(df.analysis$country)) 
n.ctry <- nlevels(as.factor(df.analysis$country)) 
n.yr <- length(unique(df.analysis$year))
classes <- c("r1", "r2", "r3", "w1", "w2")


# Define the model
mod0 <- "model {
  for (i in 1:N) {
    for (k in 1:5) {
      n.sz[i, k] ~ dnegbin(p[i,k], r[k])
      p[i, k] <- r[k]/(mu[i, k] + r[k])
      mu[i, k] <- phi[i,k]*theta[i,k] * lambda[i,k] * z[i, k] + 0.00001
      log(lambda[i, k])  <- a1[ctry[i], k] * X[i,1] + a2[ctry[i], k] * X[i,2]
                         + a3[ctry[i], k] * X[i,3] + a4[ctry[i], k] * X[i,4]  
                         + a5[ctry[i], k] * X[i,5] + a6[ctry[i], k] * X[i,6]
                         + a7[ctry[i], k] * X[i,7] + a8[ctry[i], k] * X[i,8]
                         + a9[ctry[i], k] * X[i,9] + a10[ctry[i], k] * X[i,10]
      logit(phi[i, k])   <- b1 * LE1[i]  + b2 * TCI[i] 
                          
      logit(theta[i, k]) <- g1 * ETIS.rep[i] + g2 * CITES.rep[i]
      
      z[i, k] ~ dbern(pi[ctry[i], k])
      
      LogLik[i, k] <- log((1-pi[ctry[i], k])*(n.sz[i,k] == 0) + pi[ctry[i], k]*dnegbin(n.sz[i,k],p[i,k],r[k]))
    } 
  }


  for (k in 1:5) {
    lgr[k] ~ dunif(0, 5)
    log(r.cont[k]) <- lgr[k]
    r[k] <- round(r.cont[k])
  }
  for (j in 1:N.c) {
    a1[j, 1:5] ~ dmnorm(mu1[], Omega1.inv[,])
    a2[j, 1:5] ~ dmnorm(mu2[], Omega2.inv[,])
    a3[j, 1:5] ~ dmnorm(mu3[], Omega3.inv[,])
    a4[j, 1:5] ~ dmnorm(mu4[], Omega4.inv[,])
    a5[j, 1:5] ~ dmnorm(mu5[], Omega5.inv[,])
    a6[j, 1:5] ~ dmnorm(mu6[], Omega6.inv[,])
    a7[j, 1:5] ~ dmnorm(mu7[], Omega7.inv[,])
    a8[j, 1:5] ~ dmnorm(mu8[], Omega8.inv[,])
    a9[j, 1:5] ~ dmnorm(mu9[], Omega9.inv[,])
    a10[j, 1:5] ~ dmnorm(mu10[], Omega10.inv[,])

    for (k in 1:5) {
      pi[j, k] ~ dunif(0,1)
    }
  }

  b1 ~ dnorm(0, 1.0E-04)
  b2 ~ dnorm(0, 1.0E-04)
  g1 ~ dnorm(0, 1.0E-04)
  g2 ~ dnorm(0, 1.0E-04)

  mu1[1:5] ~ dmnorm(mn[], prec[,])
  mu2[1:5] ~ dmnorm(mn[], prec[,])
  mu3[1:5] ~ dmnorm(mn[], prec[,])
  mu4[1:5] ~ dmnorm(mn[], prec[,])
  mu5[1:5] ~ dmnorm(mn[], prec[,])
  mu6[1:5] ~ dmnorm(mn[], prec[,])
  mu7[1:5] ~ dmnorm(mn[], prec[,])
  mu8[1:5] ~ dmnorm(mn[], prec[,])
  mu9[1:5] ~ dmnorm(mn[], prec[,])
  mu10[1:5] ~ dmnorm(mn[], prec[,])


  
  Omega1.inv[1:5, 1:5] ~ dwish(R1[,], 5) # Rs in data
  Omega2.inv[1:5, 1:5] ~ dwish(R2[,], 5) # Rs in data
  Omega3.inv[1:5, 1:5] ~ dwish(R3[,], 5) # Rs in data
  Omega4.inv[1:5, 1:5] ~ dwish(R4[,], 5) # Rs in data
  Omega5.inv[1:5, 1:5] ~ dwish(R5[,], 5) # Rs in data
  Omega6.inv[1:5, 1:5] ~ dwish(R6[,], 5) # Rs in data
  Omega7.inv[1:5, 1:5] ~ dwish(R7[,], 5) # Rs in data
  Omega8.inv[1:5, 1:5] ~ dwish(R8[,], 5) # Rs in data
  Omega9.inv[1:5, 1:5] ~ dwish(R9[,], 5) # Rs in data
  Omega10.inv[1:5, 1:5] ~ dwish(R10[,], 5) # Rs in data

 
} #end model
"

# Write model to file
writeLines(mod0, paste('JAGS MCMC outputs/', filename, '.txt', sep = ''))


# Spline basis functions
knots <- seq(yearfrom, yearto, length.out=K+2)[-c(1,K+2)]
X.bs <- bs(df.analysis$year, knots = knots, degree = 3, intercept = TRUE)     # column k of X.bs is the kth basis function (k = 1,...,3+K+1) at the year value of row i


# Estimate the covariance terms for the priors on Omega.inv 
# - Covariance matrices for the basis function coefficients for each country and seizure type on log scale
log.szs <- df.analysis %>%
  dplyr::select(year, country, r1, r2, r3, w1, w2) %>%
  mutate_at(vars(r1, r2, r3, w1, w2), function(x) (log(x+1)))
coefs <- array(0, dim = c(n.ctry, 5, 3+K+1))            # 3D array to store coefficients

for (i in 1:n.ctry){
  ctry.szs <- log.szs %>%
    filter(country == ctrys[i])
  X.yrs <- bs(ctry.szs$year, knots = knots, degree = 3, intercept = TRUE)
  for (k in 1:length(classes)){
    
    lmctry <- lm(ctry.szs[,classes[k]] ~ 0 + X.yrs)     # remove intercept term - B-spline basis covers this
    coefs[i, k, ] <- lmctry$coefficients
  }  
}

for (j in 1:(3+K+1)){
  assign(paste0("diag.cov.a",j), 5*diag(diag(cov(coefs[, , j]))))
}


# Model data
mod.dat <- 
  list(N      = dim(df.analysis)[1], 
       N.c    = n.ctry,
       n.sz   = as.matrix(df.analysis[,classes]),
       ctry   = as.numeric(as.factor(df.analysis$country)), 
       X      = X.bs,
       LE1    = zsc.fn(df.analysis$LE1),
       TCI =  zsc.fn(df.analysis$TCI), 
       ETIS.rep     = zsc.fn(df.analysis$ETIS.rep),
       CITES.rep  = zsc.fn(df.analysis$CITES.rep.logit),
       R1     = diag.cov.a1,
       R2     = diag.cov.a2,
       R3     = diag.cov.a3,
       R4     = diag.cov.a4,
       R5     = diag.cov.a5,
       R6     = diag.cov.a6,
       R7     = diag.cov.a7,
       R8     = diag.cov.a8,
       R9     = diag.cov.a9,
       R10     = diag.cov.a10,
       mn     = c(0,0,0,0,0),
       prec   = diag(0.0001,5,5)
  )

# Initial parameters
mod.inits <- function(){  
  list(a1      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a2      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a3      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a4      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)), 
       a5      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a6      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a7      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a8      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a9      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       a10      = matrix(ncol = 5, data = runif(n.ctry * 5, -1, 1)),
       b1      = runif(1, -1, 1),
       b2      = runif(1, -1, 1),
       g1      = runif(1, -1, 1),
       g2      = runif(1, -1, 1),
       lgr     = runif(5, 0, 1),
       mu1     = runif(5, -1, 1),
       mu2     = runif(5, -1, 1),
       mu3     = runif(5, -1, 1),
       mu4     = runif(5, -1, 1),
       mu5     = runif(5, -1, 1),
       mu6     = runif(5, -1, 1),
       mu7     = runif(5, -1, 1),
       mu8     = runif(5, -1, 1),
       mu9     = runif(5, -1, 1),
       mu10     = runif(5, -1, 1),
       Omega1.inv = diag(runif(5, 0, 1), 5),
       Omega2.inv = diag(runif(5, 0, 1), 5),
       Omega3.inv = diag(runif(5, 0, 1), 5),
       Omega4.inv = diag(runif(5, 0, 1), 5),
       Omega5.inv = diag(runif(5, 0, 1), 5),
       Omega6.inv = diag(runif(5, 0, 1), 5),
       Omega7.inv = diag(runif(5, 0, 1), 5),
       Omega8.inv = diag(runif(5, 0, 1), 5),
       Omega9.inv = diag(runif(5, 0, 1), 5),
       Omega10.inv = diag(runif(5, 0, 1), 5),
       pi     = matrix(ncol = 5, nrow = n.ctry, data = runif(n.ctry * 5, 0, 1))
  )  
}

# Model parameters to output
pms <- c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'b1', 'b2', 'g1', 'g2',
         'r.cont', 'Omega1.inv', 'Omega2.inv', 'Omega3.inv', 'Omega4.inv', 'Omega5.inv',
         'Omega6.inv', 'Omega7.inv', 'Omega8.inv', 'Omega9.inv', 'Omega10.inv',
         'mu1', 'mu2', 'mu3', 'mu4', 'mu5', 'mu6', 'mu7', 'mu8', 'mu9', 'mu10', 'pi',
         'LogLik', 'lambda', 'theta', 'phi')


# Fit model
t1 <- Sys.time()
set.seed(12345)
mod <- jags(data = mod.dat,
            inits = mod.inits,
            parameters.to.save = pms,
            model.file = paste('JAGS MCMC outputs/', filename, '.txt', sep = ''),
            parallel = TRUE,
            n.cores = 4,
            n.chains = 4,
            n.adapt = 50000,
            n.iter = 200000,
            n.burnin = 100000,
            n.thin = 10)
t2 <- Sys.time()          
difftime(t2, t1, units = "hours")

# Save files - csv of parameter summaries and jags model output
write_csv(rownames_to_column(as.data.frame(mod$summary), "param"), file=paste0("JAGS MCMC outputs/", filename, ".csv"))

saveRDS(mod,file=paste('JAGS MCMC outputs/', filename, '.rds',sep=""))


# Create Information Criterion and convergence summaries and save as csv
loglik <- mod$sims.list$LogLik
n.sim <- dim(loglik)[1]
dim(loglik) <- c(n.sim, n.ctry*n.yr*5)
WAIC <- waic(loglik)
LOO <- loo(loglik)
ICSummary <- data.frame(DIC = mod$DIC, pD = mod$pD,
                        WAIC = as.numeric(WAIC$estimates["waic","Estimate"]),
                        WAIC_SE = as.numeric(WAIC$estimates["waic","SE"]),
                        pWAIC = as.numeric(WAIC$estimates["p_waic","Estimate"]),
                        pWAIC_SE = as.numeric(WAIC$estimates["p_waic","SE"]),
                        LOOIC = as.numeric(LOO$estimates["looic","Estimate"]),
                        LOOIC_SE = as.numeric(LOO$estimates["looic","SE"]),
                        pLOO = as.numeric(LOO$estimates["p_loo","Estimate"]),
                        pLOO_SE = as.numeric(LOO$estimates["p_loo","SE"]),
                        meanDev = mean(mod$sims.list$deviance),
                        mean_Rhat = mean(unlist(mod$Rhat)),
                        max_Rhat = max(unlist(mod$Rhat)),
                        mean_neff = mean(unlist(mod$n.eff)),
                        min_neff = min(unlist(mod$n.eff)))
write_csv(ICSummary, file=paste0("JAGS MCMC outputs/", filename, "_ICs.csv"))




# Save analysis parameters -----------------------------------------------------
save(n.ctry, n.yr, n.sim, file = paste0("JAGS MCMC outputs/n values_", analysis_name, ".Rdata"))






