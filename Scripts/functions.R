################################################################################
# This script contains purpose-built functions used in the data processing,    #
# analysis, and graphical outputs throughout this repository                   #
#                                                                              #
# Graham Laidler, May 2025                                                     #
################################################################################

# Load packages
library(MASS)
library(tidyverse)
library(splines)
library(jagsUI)
library(mcmcplots)
library(loo)
library(abind)
library(ggpubr)
library(bayestestR)
library(cluster)
library(ggdendro)
library(dendsort)




# To avoid clashes with the select function between MASS and tidyverse packages
select <- dplyr::select

# Empirical logit function
emp.logit <- function(y,n){
  return(log((y+0.5)/(n-y+0.5)))
}

# Standardization function
zsc.fn <- function(x){
  zz <- (x-mean(x))/sd(x)
  return(zz)
}

# For a given analysis year range and number of knots, identify the minimum number of points in a segment (year data points and segment endpoints)
knot_segment_min <- function(year_from = yearfrom, year_to = yearto, k){
  knots <- seq(year_from, year_to, length.out = k+2)
  years <- year_from:year_to
  
  points <- sort(unique(c(knots, years)))
  
  segments <- 1:(length(knots)-1)
  n <- rep(NA, length(segments))
  
  for (i in segments){
    n[i] <- length(which(points >= knots[i] & points <= knots[i+1]))
  }
  return(min(n))
}

# Modular arithmetic function
mod_n <- function(x, n) { 
  res <- x %% n
  if (res == 0) {
    return(n)
  } 
  else {
    return(res)
  }
}

# Get lambda - for when the file path depends on whether using an individual or an averaged model
Get_lambda <- function(filename){ 
  if (file.exists(paste('JAGS MCMC outputs/Processed outputs/lambda_', filename, '.rds', sep=''))){
    lambda <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/lambda_', filename, '.rds', sep=''))
  }
  else{
    mod <- readRDS(paste('JAGS MCMC outputs/', filename, '.rds',sep=""))
    lambda <- mod$sims.list$lambda
  }
  return(lambda)
}

# Function for plotting R.hat from MCMC parameters - used in Rmd file
Rhatplot <- function(rhat.dat, rhat.guide = 1.1, title = NULL){
  
  if (is.null(dim(rhat.dat)) | length(dim(rhat.dat)) == 1){
    if (is.null(title)){
      stop("\n ************************************************
        \n Function will produce a single plot. Include a title 
        \n ************************************************")
    }
    else{
      df <- data.frame(index = 1:length(rhat.dat), R.hat=rhat.dat, title = title)
      ggplot(df, aes(x=index,y=R.hat)) + 
        geom_point() + 
        geom_hline(yintercept=rhat.guide,colour="mediumpurple",linetype=2) + 
        xlab(ifelse(length(unique(df$index)) == 5, "class", "")) +
        theme_bw() +
        facet_grid(. ~ title)
    }
  }
  else if (length(dim(rhat.dat)) == 2){
    df <- data.frame(rhat.dat)
    names(df) <- c('r1', 'r2', 'r3', 'w1', 'w2')
    df$index <- 1:nrow(df)
    plot_df <- gather(df, key = "class", value = "R.hat", -index)
    ggplot(plot_df, aes(x=index,y=R.hat)) + 
      geom_point() + 
      xlab(ifelse(length(unique(plot_df$index)) > 100, "country*year", ifelse(length(unique(plot_df$index)) > 5, "country", "class"))) +
      facet_wrap(~class,nrow=1) +
      geom_hline(yintercept=rhat.guide,colour="mediumpurple",linetype=2) + 
      ggtitle(title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face="bold", size=12))
  }
}





# Plotting functions -----------------------------------------------------------

#### Transaction Index
################################################################################
# Inputs:
# - years: years to plot trend - default is yearfrom:yearto
# - filename: vector of filenames (single for single trend or multiple (up to 3) for comparison plot)
# - CI: quantiles or hpd
# - mid: when using CI=quantiles, choose median (default) or mean
# - interval: size of credible interval (under either CI method)
# - country: plot TI over all countries (default) or specify a particular country (iso code)
# - inclusion.cts: if plotting TI for a single country, you must provide inclusion.cts (as saved in script 1-05)
# - worked_threshold: Weight threshold for worked ivory classes (shown in plot title)
  
TIplot <- function(years=yearfrom:yearto, 
                   filename, 
                   CI="quantiles", 
                   mid="median", 
                   interval=0.9, 
                   country="all", 
                   inclusion.cts=NULL,
                   worked_threshold=1){
  
  for (i in 1:length(filename)) {
    
    if (file.exists(paste('JAGS MCMC outputs/', filename[i], '.rds',sep=""))){
      mod <- readRDS(paste('JAGS MCMC outputs/', filename[i], '.rds',sep=""))
      
      lambda <- mod$sims.list$lambda
    }
    else{
      lambda <- readRDS(paste('JAGS MCMC outputs/Processed outputs/lambda_', filename[i], '.rds',sep=""))
    }
    
    # Configure dimensions: lambda will be dim c(n.ctry, n.yr, 5 (n.classes), n.sim)
    n.yr <- length(years)
    n.sim <- dim(lambda)[1]
    n.ctry <- dim(lambda)[2]/n.yr
    dim(lambda) <- c(n.sim, n.ctry, n.yr, 5)
    lambda <- aperm(lambda, c(2,3,4,1))
    
    # Year summaries
    if (country == "all"){
      x.yr.sim <- apply(lambda, c(2,3,4), sum, na.rm = TRUE)   # sum over countries: x.yr.sim is n.yr x 5 x n.sim
    }
    else if (country %in% inclusion.cts){
      x.yr.sim <- lambda[which(inclusion.cts == country), , , ]   # select relevant country: x.yr.sim is n.yr x 5 x n.sim
    }
    else{
      stop("\n ************************************************
      \n country argument must be 'all' or an element of inclusion.cts
      \n ************************************************")
    }
    
    # Add sum over classes (2nd dimension) for composite plot
    x.yr.sim <- abind(x.yr.sim, apply(x.yr.sim, c(1,3), sum), along = 2)
    
    if (CI == "quantiles"){
      x.yr.lo <- apply(x.yr.sim, c(1, 2), quantile, (1-interval)/2, na.rm=T)
      x.yr.hi <- apply(x.yr.sim, c(1, 2), quantile, 1-((1-interval)/2), na.rm=T)
      if (mid == "median"){
        x.yr.mid <- apply(x.yr.sim, c(1, 2), median, na.rm=T)
      }
      else if (mid == "mean"){
        x.yr.mid <- apply(x.yr.sim, c(1, 2) , mean, na.rm=T)
      }
      else {
        stop("\n ************************************************
      \n mid argument must be 'mean' or 'median'
      \n ************************************************")
      }
    }
    else if (CI == "hpd"){
      x.yr.lo <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(hdi(x, ci=interval)[2]))
      x.yr.mid <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(map_estimate(x)))
      x.yr.hi <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(hdi(x, ci=interval)[3]))
    }
    else{
      stop("\n ************************************************
      \n CI argument must be 'quantiles' or 'hpd'
      \n ************************************************")
    }
    
    assign(paste0("x.yr.lo.",i), x.yr.lo)
    assign(paste0("x.yr.mid.",i), x.yr.mid)
    assign(paste0("x.yr.hi.",i), x.yr.hi)
    
  }
  
  # Objects x.yr.lo, x.yr.mid, x.yr.hi: each is an array of dims n.yr x 6 (5 classes + composite)
  
  #### ggplot
  gp <- c('r1', 'r2', 'r3', 'w1', 'w2', 'comp')
  titles <- c('Raw  < 10 kg','Raw  10 - 100 kg', 'Raw  \u2265 100 kg', 
              paste0('Worked  < ',worked_threshold,' kg'), paste0('Worked  \u2265 ',worked_threshold,' kg'),
              'Composite')

  for (k in 1:6){
    tempdf <- data.frame(year = years, lo = x.yr.lo.1[, k], hi = x.yr.hi.1[, k], mid = x.yr.mid.1[, k]) %>%
      mutate_at(vars(lo:mid), ~ .*100/first(mid))   # relativize lambdas to first year
    
    tempplot <- ggplot(tempdf, aes(x = year, y = mid)) +
      geom_point() +
      geom_linerange(aes(ymin=lo,ymax=hi), linewidth=0.25) +
      scale_x_continuous(breaks = seq(years[1],years[length(years)], by=4)) +
      xlab("") + ylab("") +
      ggtitle(titles[k]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
    
    if (length(filename) == 2){
      tempdf <- data.frame(year = years, lo = x.yr.lo.2[, k], hi = x.yr.hi.2[, k], mid = x.yr.mid.2[, k]) %>%
        mutate_at(vars(lo:mid), ~ .*100/first(mid))
      
      tempplot <- tempplot +
        geom_point(tempdf, mapping=aes(x = year+0.25, y = mid), shape=23, color="grey70", fill="grey70") +
        geom_linerange(tempdf, mapping=aes(x = year+0.25, ymin = lo, ymax = hi), linewidth=0.5, color="grey70")
      
      tempdf <- data.frame(year = years, lo = x.yr.lo.1[, k], hi = x.yr.hi.1[, k], mid = x.yr.mid.1[, k]) %>%
        mutate_at(vars(lo:mid), ~ .*100/first(mid))
      
      tempplot <- tempplot +
        geom_point(tempdf, mapping=aes(x = year, y = mid)) +
        geom_linerange(tempdf, mapping=aes(x = year, ymin = lo, ymax = hi), linewidth=0.5)
      
    }
    
    else if(length(filename) == 3){
      tempdf <- data.frame(year = years, lo = x.yr.lo.3[, k], hi = x.yr.hi.3[, k], mid = x.yr.mid.3[, k]) %>%
        mutate_at(vars(lo:mid), ~ .*100/first(mid))
      
      tempplot <- tempplot +
        geom_point(tempdf, mapping=aes(x = year+0.5, y = mid), shape=23, color="grey85", fill="grey85") +
        geom_linerange(tempdf, mapping=aes(x = year+0.5, ymin = lo, ymax = hi), linewidth=0.5, color="grey70")
      
      tempdf <- data.frame(year = years, lo = x.yr.lo.2[, k], hi = x.yr.hi.2[, k], mid = x.yr.mid.2[, k]) %>%
        mutate_at(vars(lo:mid), ~ .*100/first(mid))
      
      tempplot <- tempplot +
        geom_point(tempdf, mapping=aes(x = year+0.25, y = mid), shape=22, color="grey50", fill="grey50") +
        geom_linerange(tempdf, mapping=aes(x = year+0.25, ymin = lo, ymax = hi), linewidth=0.5, color="grey50")
      
      tempdf <- data.frame(year = years, lo = x.yr.lo.1[, k], hi = x.yr.hi.1[, k], mid = x.yr.mid.1[, k]) %>%
        mutate_at(vars(lo:mid), ~ .*100/first(mid))
      
      tempplot <- tempplot +
        geom_point(tempdf, mapping=aes(x = year, y = mid)) +
        geom_linerange(tempdf, mapping=aes(x = year, ymin = lo, ymax = hi), linewidth=0.5)
    }
    
    assign(paste0("TI.",gp[k]), tempplot)
  }
  
  if (country %in% inclusion.cts){
    TIfigure <- ggarrange(TI.r1, TI.w1, TI.r2, TI.w2, TI.r3, TI.comp, ncol = 2, nrow=3)
    annotate_figure(TIfigure, left = text_grob("Relative number of transactions", rot = 90), 
                    top = text_grob(country))
  }
  else if (country == "all"){
    TIfigure <- ggarrange(TI.r1, TI.w1, TI.r2, TI.w2, TI.r3, TI.comp, ncol = 2, nrow=3)
    annotate_figure(TIfigure, left = text_grob("Relative number of transactions", rot = 90))
  }
  
}


#### Posterior predictive distributions
################################################################################
# Looks at posterior predictive distributions of the ys and compares to the data
# - This can be done for each ivory class, year or country
# - Default uses 50% intervals as a way of clearly seeing how the model is working
#
# Inputs:
# - years: years to plot - default is yearfrom:yearto
# - filename: single filename of model to plot
# - df.analysis: provide data used in model fitting (to plot observed data points)
# - interval: size of quantile interval to show
# - ctry: plot over all countries (default) or specify a particular country (iso code)
# - inclusion.cts: if plotting for a single country, you must provide inclusion.cts (as saved in script 1-05)
# - type: options are "by year" (default - year on x axis, y summed over countries on y axis) or "by party" (countries on x axis, y summed over years on y axis)
# - linecol: graphical parameter for line colour of median and quantile lines

PPplot <- function(years=yearfrom:yearto, 
                   filename, 
                   df.analysis, 
                   interval=0.5, 
                   ctry="all", 
                   inclusion.cts=NULL, 
                   type = "by year",
                   linecol = "black"){
  
  # Read in files
  if (file.exists(paste('JAGS MCMC outputs/Processed outputs/lambda_', filename, '.rds', sep=''))){
    lambda <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/lambda_', filename, '.rds', sep=''))
    theta <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/theta_', filename, '.rds', sep=''))
    phi <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/phi_', filename, '.rds', sep=''))
    r.cont <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/r.cont_', filename, '.rds', sep=''))
    if (file.exists(paste('JAGS MCMC outputs/Processed outputs/pi_', filename, '.rds', sep=''))){
      pi <- readRDS(file=paste('JAGS MCMC outputs/Processed outputs/pi_', filename, '.rds', sep=''))
    }
  }
  else{
    mod <- readRDS(paste('JAGS MCMC outputs/', filename, '.rds',sep=""))
    
    lambda <- mod$sims.list$lambda
    theta <- mod$sims.list$theta
    phi <- mod$sims.list$phi
    r.cont <- mod$sims.list$r.cont
    pi <- mod$sims.list$pi
  }
  
  # Generate mu
  mu <- lambda * theta * phi
  
  # Configure dimensions: mu is dim c(n.ctry, n.yr, 5, n.sim)
  n.sim <- dim(lambda)[1]
  n.yr <- length(years)
  n.ctry <- dim(lambda)[2]/n.yr
  dim(mu) <- c(n.sim, n.ctry, n.yr, 5)
  mu <- aperm(mu, c(2,3,4,1))
  
  # Simulate y (zero-inflated) by simulating yneg from negbinom, z from bernoulli, and multiplying
  yneg <- array(NA, dim = dim(mu))
  set.seed(12345)
  for (s in 1:n.sim){
    # for (i in 1:n.ctry){
    #   for (t in 1:n.yr){
    #     yneg[i, t, , s] <- rnbinom(n = 5, size = r.cont[s, ], mu = mu[i, t, , s])
    #   }
    # }
    yneg[, , , s] <- rnbinom(n = n.ctry * n.yr * 5, size = r.cont[s, ], mu = mu[, , , s])
    # Note this does the same as using the 3 for loops (commented out) but is faster
  }
  
  if (!is.null(pi)){
    z <- array(NA, dim = dim(mu))
    set.seed(12345)
    for (s in 1:n.sim){
      for (k in 1:5){
        z[, , k, s] <- rbinom(n = n.ctry * n.yr, size = 1, prob = pi[s, , k])
      }
    }
    
    y <- yneg * z
  }
  else{
    y <- yneg
  }
  
  # Set labels for plots
  gp <- c('r1', 'r2', 'r3', 'w1', 'w2', 'comp')
  titles <- c('Raw  < 10 kg','Raw  10 - 100 kg', 'Raw  \u2265 100 kg', 'Worked  < 1 kg','Worked  \u2265 1 kg','Composite')
  
  if (type == "by year"){
    
    if (ctry == "all"){
      y.yrs.sum <- apply(y,c(2,3,4),sum,na.rm=T)   # sum over countries: y.yrs.sum is n.yr x 5 x n.sim
      
      df.obs <- df.analysis
    }
    else if (ctry %in% inclusion.cts){
      y.yrs.sum <- y[which(inclusion.cts == country), , , ]   # select relevant country: y.yrs.sum is n.yr x 5 x n.sim
      
      # Filter df.obs here for correct obs.yrs.sum later
      df.obs <- df.analysis %>%
        filter(country == ctry)
    }
    else{
      stop("\n ************************************************
        \n ctry argument must be 'all' or an element of inclusion.cts
        \n ************************************************")
    }
    
    # Add sum over classes (2nd dimension) for composite plot
    y.yrs.sum <- abind(y.yrs.sum, apply(y.yrs.sum, c(1,3), sum), along = 2)   # now y.yrs.sum is n.yr x 6 x n.sim
    
    # Calculate quantiles and median 
    y.yrs.lo <- apply(y.yrs.sum, c(1, 2), quantile, (1-interval)/2, na.rm=TRUE)
    y.yrs.hi <- apply(y.yrs.sum, c(1, 2), quantile, 1-((1-interval)/2), na.rm=TRUE)
    y.yrs.md <- apply(y.yrs.sum, c(1, 2), median, na.rm=TRUE)
    
    # Calculate relevant summary values for observed data
    obs.yrs.sum <- matrix(nrow=n.yr, ncol=6, data=0)
    for (k in 1:5){
      obs.yrs.sum[, k] <- tapply(df.obs[, gp[k]], df.obs$year, sum)
    }
    obs.yrs.sum[, 6] <- rowSums(obs.yrs.sum[,1:5])
    
    # ggplot
    for (k in 1:6){
      tempdf <- data.frame(year = years, lo = y.yrs.lo[,k], hi = y.yrs.hi[,k], mid = y.yrs.md[,k], obs = obs.yrs.sum[,k])
      
      tempplot <- ggplot(tempdf, aes(x = year, y = mid)) +
        geom_line(colour = linecol) +
        geom_line(aes(y = lo), linetype="dashed", colour = linecol) +
        geom_line(aes(y = hi), linetype="dashed", colour = linecol) +
        geom_point(aes(y = obs), color = "#d55e00") +
        xlab("") + ylab("") +
        ggtitle(titles[k]) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
      
      assign(paste0("PP.",gp[k]), tempplot)
      
    }
    
    if (ctry %in% inclusion.cts){
      PPfigure <- ggarrange(PP.r1, PP.w1, PP.r2, PP.w2, PP.r3, PP.comp, ncol = 2, nrow=3)
      annotate_figure(PPfigure, left = text_grob("Reported seizures", rot = 90), 
                      top = text_grob(country))
    }
    else if (ctry == "all"){
      PPfigure <- ggarrange(PP.r1, PP.w1, PP.r2, PP.w2, PP.r3, PP.comp, ncol = 2, nrow=3)
      annotate_figure(PPfigure, left = text_grob("Reported seizures", rot = 90))
    }
    else{
      stop("\n ************************************************
      \n ctry argument must be 'all' or an element of inclusion.cts
      \n ************************************************")
    }
    
  }
  else if (type == "by Party"){
    # For country PP plot (countries on x axes, y summed over years on y axes)
    
    y.ctry.sum <- apply(y, c(1,3,4), sum, na.rm = T)
    y.ctry.sum <- abind(y.ctry.sum, apply(y.ctry.sum, c(1,3), sum), along = 2)
    
    y.ctry.lo <- apply(y.ctry.sum, c(1, 2), quantile, (1-interval)/2, na.rm=TRUE)
    y.ctry.hi <- apply(y.ctry.sum, c(1, 2), quantile, 1-((1-interval)/2), na.rm=TRUE)
    y.ctry.md <- apply(y.ctry.sum, c(1, 2), median, na.rm=TRUE)
    
    obs.ctry.sum <- matrix(nrow=n.ctry, ncol=6, data=0)
    for (k in 1:5){
      obs.ctry.sum[, k] <- tapply(df.analysis[, gp[k]], df.analysis$country, sum)
    }
    obs.ctry.sum[, 6] <- rowSums(obs.ctry.sum[,1:5])
    
    #### ggplot
    for (k in 1:6){
      
      tempdf <- data.frame(country = inclusion.cts, lo = y.ctry.lo[,k], hi = y.ctry.hi[,k], mid = y.ctry.md[,k], obs = obs.ctry.sum[,k])
      
      tempplot <- ggplot(tempdf, aes(x = country, y = mid)) +
        geom_linerange(aes(ymin=lo,ymax=hi), color = linecol) +
        geom_point(aes(y = obs), color = "#d55e00") +
        xlab("") + ylab("") +
        ggtitle(titles[k]) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=10),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))
      
      assign(paste0("PP.",gp[k]), tempplot)
    }
    
    PPfigure <- ggarrange(PP.r1, PP.w1, PP.r2, PP.w2, PP.r3, PP.comp, ncol = 2, nrow=3)
    annotate_figure(PPfigure, left = text_grob("Reported seizures", rot = 90))
    
  }
  else{
    stop("\n ************************************************
        \n type argument must be 'by year' (default) or 'by Party'
        \n ************************************************")
  }
  
}


#### Weight Index
################################################################################

# WIplot allows plotting of WI as points with CIs next to stacked bars (if plotting a single model),
# or as points with CIs for the composite WI if comparing two models. 
# Inputs:
# - years: years to plot - default is yearfrom:yearto
# - filename: vector of filenames (single for composite WI plot with credible intervals next to a stacked (median WI by ivory class) bar chart or multiple (up to 2) for composite WI with credible intervals comparison plot)
# - CI: quantiles (default) or hpd
# - mid: when using CI=quantiles, choose median (default) or mean
# - interval: size of credible interval (under either CI method)
# - country: plot TI over all countries (default) or specify a particular country (iso code)
# - inclusion.cts: if plotting TI for a single country, you must provide inclusion.cts (as saved in script 1-05)

WIplot <- function(years=yearfrom:yearto, 
                   filename, 
                   CI="quantiles", 
                   mid="median", 
                   interval=0.9, 
                   country="all", 
                   inclusion.cts=NULL){
  
  for (i in 1:length(filename)) {
    
    w.sum <- readRDS(paste('JAGS MCMC outputs/Processed outputs/weights_', filename[i], '.rds',sep=""))
    
    # Configure dimensions: w.sum will be dim c(n.ctry, n.yr, 5 (n.classes), n.sim)
    n.yr <- length(years)
    n.sim <- dim(w.sum)[1]
    n.ctry <- dim(w.sum)[2]/n.yr
    dim(w.sum) <- c(n.sim, n.ctry, n.yr, 5)
    w.sum <- aperm(w.sum, c(2,3,4,1))
    
    # Year summaries
    if (country == "all"){
      x.yr.sim <- apply(w.sum,c(2,3,4),sum,na.rm=T)   # sum over countries: x.yr.sim is n.yr x 5 x n.sim
    }
    else if (country %in% inclusion.cts){
      x.yr.sim <- w.sum[which(inclusion.cts == country), , , ]   # select relevant country: x.yr.sim is n.yr x 5 x n.sim
    }
    else{
      stop("\n ************************************************
      \n country argument must be 'all' or an element of inclusion.cts
      \n ************************************************")
    }
    
    # Add sum over classes (2nd dimension) for composite plot
    x.yr.sim <- abind(x.yr.sim, apply(x.yr.sim, c(1,3), sum), along = 2)
    
    if (CI == "quantiles"){
      x.yr.lo <- apply(x.yr.sim, c(1, 2), quantile, (1-interval)/2, na.rm=T)
      x.yr.hi <- apply(x.yr.sim, c(1, 2), quantile, 1-((1-interval)/2), na.rm=T)
      if (mid == "median"){
        x.yr.mid <- apply(x.yr.sim, c(1, 2), median, na.rm=T)
      }
      else if (mid == "mean"){
        x.yr.mid <- apply(x.yr.sim, c(1, 2) , mean, na.rm=T)
      }
      else {
        stop("\n ************************************************
      \n mid argument must be 'mean' or 'median'
      \n ************************************************")
      }
    }
    else if (CI == "hpd"){
      x.yr.lo <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(hdi(x, ci=interval)[2]))
      x.yr.mid <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(map_estimate(x)))
      x.yr.hi <- apply(x.yr.sim, c(1, 2), function(x) as.numeric(hdi(x, ci=interval)[3]))
    }
    else{
      stop("\n ************************************************
      \n CI argument must be 'quantiles' or 'hpd'
      \n ************************************************")
    }
    
    assign(paste0("x.yr.lo.",i), x.yr.lo)
    assign(paste0("x.yr.mid.",i), x.yr.mid)
    assign(paste0("x.yr.hi.",i), x.yr.hi)
    
  }
  
  # Objects x.yr.lo, x.yr.mid, x.yr.hi: each is an array of dims n.yr x 6 (5 classes + composite)
  
  #### ggplot - Points for composite WI with CI bars
  tempdf <- data.frame(year = years, lo = x.yr.lo.1[, 6], hi = x.yr.hi.1[, 6], mid = x.yr.mid.1[, 6]) %>%
    mutate_at(vars(lo:mid), ~ .*100/first(mid))   # relativize lambdas to first year
  
  WIcomp <- ggplot(tempdf, aes(x = year, y = mid)) +
    geom_point() +
    geom_linerange(aes(ymin=lo,ymax=hi), linewidth=0.25) +
    xlab("") + ylab("") +
    scale_x_continuous(breaks = seq(years[1],years[length(years)], by=2)) +
    ggtitle("Composite Weight Index") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
  
  if (length(filename) == 2){
    tempdf <- data.frame(year = years, lo = x.yr.lo.2[, 6], hi = x.yr.hi.2[, 6], mid = x.yr.mid.2[, 6]) %>%
      mutate_at(vars(lo:mid), ~ .*100/first(mid))
    
    WIcomp <- WIcomp +
      geom_point(tempdf, mapping=aes(x = year+0.25, y = mid), color="#d55e00", shape = 18) +
      geom_linerange(tempdf, mapping=aes(x = year+0.25, ymin = lo, ymax = hi), color="#d55e00")
  }
    
  
  if (country %in% inclusion.cts){
    WIcomp <- ggarrange(WI.comp, ncol = 1, nrow=1)
    WIcomp <- annotate_figure(WIcomp, left = text_grob("Relative weight", rot = 90), 
                              top = text_grob(country))
  }
  else if (country == "all"){
    WIcomp <- ggarrange(WI.comp, ncol = 1, nrow=1)
    WIcomp <- annotate_figure(WIcomp, left = text_grob("Relative weight", rot = 90))
  }
  
  if (length(filename) == 2){
    WIcomp
  }
  else{
    yupper <- ceiling(max(tempdf$hi)/100)*100
    WI.comp1 <- ggplot(tempdf, aes(x = year, y = mid)) +
      geom_point() +
      geom_linerange(aes(ymin=lo,ymax=hi), linewidth=0.25) +
      xlab("") + ylab("") +
      scale_x_continuous(breaks = seq(years[1],years[length(years)], by=2)) + 
      ylim(0,yupper) +
      ggtitle(titles[k]) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))
    
    #### ggplot - Stacked bar chart with median WI per class
    tempdf <- data.frame(year = years, r1 = x.yr.mid.1[, 1], r2 = x.yr.mid.1[, 2], r3 = x.yr.mid.1[, 3],
                         w1 = x.yr.mid.1[, 4], w2 = x.yr.mid.1[, 5]) %>%
      mutate(sum = r1+r2+r3+w1+w2) %>%
      mutate_at(vars(r1:w2), ~ .*100/first(sum)) %>%
      select(-sum) %>%
      gather("class", "weight", -year)
    
    labels <- c('Worked  \u2265 1 kg', 'Worked  < 1 kg', 'Raw  \u2265 100 kg', 'Raw  10 - 100 kg', 'Raw  < 10 kg')
    
    WIBar <- ggplot(tempdf, aes(x=year, y=weight)) +
      geom_bar(aes(fill=forcats::fct_rev(class)), colour="black", position = "stack", stat="identity",width = 0.8,linewidth=0.25) + 
      scale_fill_manual(name = "",
                        label=labels,
                        values=  c("#F5ECE0", "#EBD9C1", "#9E9FA9", "#BEBFC5", "#DFDFE2")) + 
      xlab("") + ylab("Relative weight") +
      scale_x_continuous(breaks = seq(years[1],years[length(years)], by=2)) +
      ylim(0,yupper) +
      ggtitle("Weight Index by ivory class") +
      theme_bw() +
      theme(axis.title = element_text(size=12),
            plot.title = element_text(hjust = 0.5, face="bold", size=10),
            legend.title = element_blank(),
            legend.position = "inside",
            legend.justification.inside = c(0.99,0.99))
    
    ggarrange(WIBar, WI.comp1, ncol = 2, nrow=1, align = "h")
    
  }

}


