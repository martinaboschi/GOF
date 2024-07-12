#-------------------------------------------------------------------------------
# Script: Simulation_global_two_elements.R
# Author: Martina Boschi
# Date: July 2024
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
# 3. Simulating the data with two elements driving the DGP
# (time-related definition of reciprocity with FLE, 
# covariate sampled from an exponential distribution with FLE)
dat.gams <- sapply(1:n.iter, function(x)
  data.simulation(n = n, p = p, seed = x+1000, L=2),
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
  
  # 4.2. Correctly specified model: time-related reciprocity with FLE + x1
  gam.fitsCS[[iter]] <- gam(y~ -1 + last.rec + x1,
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  # misspecified functional form for the computation of the covariate
  gam.fitsMS[[iter]] <- list()
  
  # I mispecification :
  # miscomputed reciprocity with FLE
  # x1 correctly specified
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + last.rec_mis + x1,
                                 family = binomial,
                                 data=data)
  
  # II mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  gam.fitsMS[[iter]][[2]] <- gam(y~ -1 + last.rec + x4,
                                 family = binomial,
                                 data=data)
  
  # III mispecification :
  # miscomputed reciprocity with FLE
  # x1 misspecified (multiplied by time)
  gam.fitsMS[[iter]][[3]] <- gam(y~ -1 + last.rec_mis + x4,
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
BB.stat_1 <- BB.simulator(set_covariates = list(1))
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 7. Testing the correctly specified model for coverage

# Inspecting single components
pvalue.CS_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                dat.gams[[x]],
                                                gam.fit = gam.fitsCS[[x]], 
                                                index = 1)[[1]])
plot(ecdf(pvalue.CS_comp1))
abline(0,1)

pvalue.CS_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                dat.gams[[x]],
                                                gam.fit = gam.fitsCS[[x]], 
                                                index = 2)[[1]])
plot(ecdf(pvalue.CS_comp2))
abline(0,1)

# Testing via global test

T_g_CS <- (tan(pi*(0.5-pvalue.CS_comp1)) + 
           tan(pi*(0.5-pvalue.CS_comp2)))/2
pvalue_Cauchy_CS <- 1/2 - atan(T_g_CS)/pi
plot(ecdf(pvalue_Cauchy_CS), cex=0.5,
     main="2 components - correctly specified model")
abline(0,1)

pvalue.CS <- pvalue_Cauchy_CS
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power

# I mispecification :
# Indicator for reciprocity with FLE
# x1 correctly specified

# Inspecting single components
pvalue.MS1_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[1]], 
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS1_comp1))
abline(0,1)

pvalue.MS1_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[1]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS1_comp2))
abline(0,1)

# Testing via global test

T_g_MS1 <- (tan(pi*(0.5-pvalue.MS1_comp1)) + tan(pi*(0.5-pvalue.MS1_comp2)))/2
pvalue_Cauchy_MS1 <- 1/2 - atan(T_g_MS1)/pi
plot(ecdf(pvalue_Cauchy_MS1))
abline(0,1)

pvalue.MS1 <- pvalue_Cauchy_MS1

# II mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)

# Inspecting single components
pvalue.MS2_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[2]], 
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS2_comp1))
abline(0,1)

pvalue.MS2_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[2]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS2_comp2))
abline(0,1)

# Testing via global test

T_g_MS2 <- (tan(pi*(0.5-pvalue.MS2_comp1)) + tan(pi*(0.5-pvalue.MS2_comp2)))/2
pvalue_Cauchy_MS2 <- 1/2 - atan(T_g_MS2)/pi
plot(ecdf(pvalue_Cauchy_MS2))
abline(0,1)

pvalue.MS2 <- pvalue_Cauchy_MS2

# III mispecification :
# Indicator for reciprocity with FLE
# x1 misspecified (multiplied by time)

# Inspecting single components
pvalue.MS3_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[3]], 
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS3_comp1))
abline(0,1)

pvalue.MS3_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[3]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS3_comp2))
abline(0,1)

# Testing via global test

T_g_MS3 <- (tan(pi*(0.5-pvalue.MS3_comp1)) + tan(pi*(0.5-pvalue.MS3_comp2)))/2
pvalue_Cauchy_MS3 <- 1/2 - atan(T_g_MS3)/pi
plot(ecdf(pvalue_Cauchy_MS3))
abline(0,1)

pvalue.MS3 <- pvalue_Cauchy_MS3

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, 
     BB.stat_1, 
     pvalue.CS, pvalue.CS_comp1, pvalue.CS_comp2, 
     pvalue.MS1, pvalue.MS1_comp1, pvalue.MS1_comp2, 
     pvalue.MS2, pvalue.MS2_comp1, pvalue.MS2_comp2, 
     pvalue.MS3, pvalue.MS3_comp1, pvalue.MS3_comp2, 
     file="01-Simulation-Studies/03-Omnibus-Testing/output/global_two_elements.RData")
#-------------------------------------------------------------------------------