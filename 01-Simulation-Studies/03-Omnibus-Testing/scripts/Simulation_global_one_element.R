#-------------------------------------------------------------------------------
# Script: Simulation_global_one_element.R
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
n <- 2000
p <- 50
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Simulating the data with one element driving the DGP
# (time-related definition of reciprocity with FLE)
dat.gams <- sapply(1:n.iter, function(x)
  data.simulation(n = n, p = p, seed = x, L=1),
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
  
  # 4.2. Correctly specified model: time-related reciprocity with FLE
  gam.fitsCS[[iter]] <- gam(y~ -1 + last.rec,
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  # misspecified functional form for the computation of the covariate
  gam.fitsMS[[iter]] <- list()
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + last.rec_mis,
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
BB.stat <- BB.simulator(set_covariates = list(1))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 7. Testing the correctly specified model for coverage
# Testing via global test
pvalue.CS <- sapply(1:n.iter, function(x) GOF_pvalue_global(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsCS[[x]], 
                                          set_covariates = list(1), 
                                          BB.stat = BB.stat))
plot(ecdf(pvalue.CS))
abline(0,1)

# Testing via Kolmogorov test
pvalue.CS_alt <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                              dat.gams[[x]],
                                              gam.fit = gam.fitsCS[[x]], 
                                              index = 1)[[1]])
plot(ecdf(pvalue.CS_alt))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power
# Testing via global test
pvalue.MS <- sapply(1:n.iter, function(x) GOF_pvalue_global(data = 
                                          dat.gams[[x]],
                                          gam.fit = gam.fitsMS[[x]][[1]], 
                                          set_covariates = list(1), 
                                          BB.stat = BB.stat))

plot(ecdf(pvalue.MS))
abline(0,1)

# Testing via Kolmogorov test
pvalue.MS_alt <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                              dat.gams[[x]],
                                              gam.fit = gam.fitsMS[[x]][[1]], 
                                              index = 1)[[1]])
plot(ecdf(pvalue.MS_alt))
abline(0,1)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, 
     BB.stat, 
     pvalue.CS, pvalue.CS_alt, 
     pvalue.MS, pvalue.MS_alt, 
     file="01-Simulation-Studies/03-Omnibus-Testing/output/global_one_element.RData")
#-------------------------------------------------------------------------------
