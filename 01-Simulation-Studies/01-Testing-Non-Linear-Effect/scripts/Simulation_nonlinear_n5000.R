#-------------------------------------------------------------------------------
# Script: Simulation_nonlinear_n5000.R
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
n <- 5000
p <- 30
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Simulating the data with one element driving the DGP
# (intrinsic non-linear reciprocity driving the DGP)
dat.gams <- sapply(1:n.iter, function(x)
  data.generation_nonlinear_reciprocity(n = n, p = p, seed = x),
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
  
  W <- T <- cbind(data$recExpTime1,
                  data$recExpTime2)
  W[,1] <- data$recId1; 
  W[,2] <- -data$recId2;
  
  # 4.2. Correctly specified model: non-linear function of exp. reciprocity time
  gam.fitsCS[[iter]] <- gam(y~ -1 + s(T, by=W),
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  # exp. reciprocity time with FLE
  gam.fitsMS[[iter]] <- list()
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + recExpTime,
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
BB.stat_10 <- BB.single(dim.k=10)[[2]]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 7. Testing the correctly specified model for coverage
# Testing via KS test with multivariate statistics
pvalue.CS <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsCS[[x]], 
                                          index = 1:10, 
                                          BB.stat = BB.stat_10)[[1]])
plot(ecdf(pvalue.CS))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power
# Testing via Kolmogorov test
pvalue.MS <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsMS[[x]][[1]], 
                                          index = 1)[[1]])
plot(ecdf(pvalue.MS))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, 
     BB.stat_10, 
     pvalue.CS,
     pvalue.MS,
     file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n5000.RData")
#-------------------------------------------------------------------------------
