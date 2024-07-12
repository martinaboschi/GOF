#-------------------------------------------------------------------------------
# Script: Simulation_randomeffects_p50.R
# Author: Martina Boschi
# Date: May 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2. Setting number of simulations, events and actors 
n.iter <- 500
n <- 7000
p <- 50
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Simulating the data with one element driving the DGP
# (intrinsic non-linear reciprocity driving the DGP)
dat.gams <- sapply(1:n.iter, function(x)
  data.generation_random_effects(n = n, p = p, seed = x, act.sd = 0.5),
  simplify = FALSE)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 4. Fitting correctly specified and misspecified models
gam.fitsCS <- gam.fitsMS <- list()

for (iter in 1:n.iter){
  
  print(iter)
  
  # 4.1. Extrapolating the data and setting a fixed response equal to 1
  data <- dat.gams[[iter]]
  y <- rep(1, nrow(data))
  
  s1 <- data$s1
  s2 <- data$s2
  ss <- factor(c(s1,s2))
  dim(ss) <- c(length(s1),2)
  Ls <- Ls_star <- matrix(0,length(s1),2); 
  
  Ls_star[,1] <- data$stp
  Ls_star[,2] <- -data$stp
  
  Ls[,1] <- 1
  Ls[,2] <- -1
  
  # 4.1. Correctly specified model: random intercept for sender
  gam.fitsCS[[iter]] <- gam(y~ -1 + s(ss, by=Ls, bs="re"),
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  # random slope for time depending on sender
  gam.fitsMS[[iter]] <- list()
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 5. Uniforming notations
objects <- convert.to.objects(dat.gams = dat.gams, 
                              gam.fitsCS = gam.fitsCS, 
                              gam.fitsMS = gam.fitsMS)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 6. Reproducing the asymptotic behaviour of the MR process
BB.stat <- BB.single(dim.k=p)[[2]]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 7. Testing the correctly specified model for coverage
# Testing via KS test with multivariate statistics
pvalue.CS <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsCS[[x]], 
                                          index = 1:p, 
                                          BB.stat = BB.stat)[[1]])
plot(ecdf(pvalue.CS))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power
# Testing via Kolmogorov test
pvalue.MS <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsMS[[x]][[1]], 
                                          index = 1:p, 
                                          BB.stat = BB.stat)[[1]])
plot(ecdf(pvalue.MS))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, 
     BB.stat, 
     pvalue.CS,
     pvalue.MS,
     file="01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p50.RData")
#-------------------------------------------------------------------------------
