################################################################################
# - Simulate weights for Weight Index using lambdas from Transaction Index     #
#                                                                              #
# - Weights of ivory transactions for each year, country and ivory class are   #
# - simulated using the raw and worked weight distribution models (script 2-03)# 
# - and posterior samples of the lambdas (script 2-02) - the numbers of        #
# - transactions in each year, country and ivory class                         #
################################################################################


#### Get parameters ------------------------------------------------------------

# Raw weight model parameters
filename_wt <- paste0("raw_wt_model_", analysis_name)
mod <- readRDS(paste('JAGS MCMC outputs/', filename_wt, '.rds',sep=""))

a0_raw <- mod$sims.list$a0
sigma2.y_raw <- mod$sims.list$sigma2.y
nu_raw <- mod$sims.list$nu

# Worked weight model parameters
filename_wt <- paste0("wkd_wt_model_", analysis_name)
mod <- readRDS(paste('JAGS MCMC outputs/', filename_wt, '.rds',sep=""))

a0_wkd <- mod$sims.list$a0
sigma2.y_wkd <- mod$sims.list$sigma2.y
nu_wkd <- mod$sims.list$nu

# Lambda (Transaction Index samples)
lambda <- Get_lambda(paste0("ave_", analysis_name))   #This is for the averaged model. Can be changed to an individual model if desired

# Analysis parameters
n.yrs <- yearto - yearfrom + 1
n.ctry <- dim(lambda)[2]/n.yrs
n.sim <- dim(lambda)[1]
n.sim.wt <- length(a0_raw)



#### Simulate weights ----------------------------------------------------------

# Initialise array to store the weights
w.sum <- array(0, dim=dim(lambda))


# Loop to simulate weights
for (i in 1:n.sim){
  
  # loop indicator
  if(((i/1000) %% 1)==0){
    print(paste0("i=",i, " out of ", n.sim))
  }	
  
    
  nums.lambda <- list(round(lambda[i, , 1]), round(lambda[i, , 2]), round(lambda[i, , 3]),
                      round(lambda[i, , 4]), round(lambda[i, , 5]))
  total.lambda <- unlist(lapply(nums.lambda, sum, na.rm = TRUE))
  
  w.sim.gp <- list(NULL)
  
  #### Simulate raw weights 
  i.wt <- mod_n(i, n.sim.wt)           # if less simulations of weight model parameters than lambda, loop back round
  mu.i <- a0_raw[i.wt]   
  sigma.i <- sqrt(sigma2.y_raw[i.wt])
  df.i <- nu_raw[i.wt]
  
  cut.1 <- (log(10) - mu.i)/sigma.i
  cut.2 <- (log(100) - mu.i)/sigma.i
  cut.3 <- (log(10000) - mu.i)/sigma.i 
  
  # Estimate how many simulations we need from Student-t distribution to satisfy required totals in each class
  N.sim <- max(round(total.lambda[1] /pt(cut.1, df = df.i)),                    # pt gives P(X < cut.1), i.e., probability a simulated weight will fall into class r1
               round(total.lambda[2] /(pt(cut.2, df = df.i) - pt(cut.1, df = df.i))), 
               round(total.lambda[3] /(pt(cut.3, df = df.i) - pt(cut.2, df = df.i))))  

  w.sim <- rt(n = N.sim, df = df.i)
  w.r1 <- w.sim[which(w.sim < cut.1)]
  w.r2 <- w.sim[which(w.sim >= cut.1 & w.sim < cut.2)]
  w.r3 <- w.sim[which(w.sim >= cut.2 & w.sim < cut.3)]
  
  # Add more simulations if there are too few in any of the classes
  while(length(w.r1) < total.lambda[1]){
    w.extra <- rt(n = 1000, df = df.i)
    w.r1 <- c(w.r1, w.extra[which(w.extra < cut.1)])
  }
  while(length(w.r2) < total.lambda[2]){
    w.extra <- rt(n = 1000, df = df.i)
    w.r2 <- c(w.r2, w.extra[which(w.extra >= cut.1 & w.extra < cut.2)])
  }
  while(length(w.r3) < total.lambda[3]){
    w.extra <- rt(n = 1000, df = df.i)
    w.r3 <- c(w.r3, w.extra[which(w.extra >= cut.2 & w.extra < cut.3)])
  }
  
  # Add simulated weight values in each category to list
  w.sim.gp[[1]] <- exp(w.r1[1:total.lambda[1]] * sigma.i + mu.i)
  w.sim.gp[[2]] <- exp(w.r2[1:total.lambda[2]] * sigma.i + mu.i)
  w.sim.gp[[3]] <- exp(w.r3[1:total.lambda[3]] * sigma.i + mu.i)
  
  
  
  #### Simulate worked weights
  mu.i <- a0_wkd[i.wt]
  sigma.i <- sqrt(sigma2.y_wkd[i.wt])
  df.i <- nu_wkd[i.wt]
  
  cut.1 <- (log(1) - mu.i)/sigma.i
  cut.2 <- (log(10000) - mu.i)/sigma.i
  
  N.sim <- max(round(total.lambda[4] /pt(cut.1, df = df.i)), 
               round(total.lambda[5] /(pt(cut.2, df = df.i) - pt(cut.1, df = df.i))))  
  
  y.sim <- rt(n = N.sim, df = df.i)
  w.w1 <- y.sim[which(y.sim < cut.1)]
  w.w2 <- y.sim[which(y.sim >= cut.1 & y.sim < cut.2)]

  # Add more simulations if there are too few in any of the classes
  while(length(w.w1) < total.lambda[4]){
    w.extra <- rt(n = 1000, df = df.i)
    w.w1 <- c(w.w1, w.extra[which(w.extra < cut.1)])
  }
  while(length(w.w2) < total.lambda[5]){
    w.extra <- rt(n = 1000, df = df.i)
    w.w2 <- c(w.w2, w.extra[which(w.extra >= cut.1 & w.extra < cut.2)])
  }
  
  # Add simulated weight values in each category to list
  w.sim.gp[[4]] <- exp(w.w1[1:total.lambda[4]] * sigma.i + mu.i)
  w.sim.gp[[5]] <- exp(w.w2[1:total.lambda[5]] * sigma.i + mu.i)
  
  
  #### Allocate the weights to each country and year so that we can do the correct sums
  id.ctry.yr <- list(NULL)
  for (k in 1:5){
    id.ctry.yr[[k]] <- rep(1:(n.ctry*n.yrs), nums.lambda[[k]])
  }
  # id.ctry.yr shows which ctry*year number each simulated weight will be assigned to
  
  for (k in 1:5){
    if (length(id.ctry.yr[[k]]) > 0){
      w.sum[i, nums.lambda[[k]] > 0, k] <- tapply(w.sim.gp[[k]], id.ctry.yr[[k]], sum)
    }
  }
}


#### Save results --------------------------------------------------------------
saveRDS(w.sum, file=paste('JAGS MCMC outputs/Processed outputs/weights_', analysis_name, '.rds',sep=""))
