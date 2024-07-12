#-------------------------------------------------------------------------------
# Script: Simulation_global_three_elements.R
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
# 3. Simulating the data with three elements driving the DGP
# (time-related definition of reciprocity with FLE, 
# covariate sampled from an exponential distribution with FLE,
# time-transformation of a covariate sampled form a normal with NLE)
dat.gams <- sapply(1:n.iter, function(x)
  data.simulation(n = n, p = p, seed = x+2000, L=3),
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
  
  S <- U <- cbind(data$x21,
                  data$x22)
  S[,1] <- 1
  S[,2] <- -1
  
  # 4.2. Correctly specified model: time-related reciprocity with FLE + x1
  #                                 + x2 with NLE
  gam.fitsCS[[iter]] <- gam(y~ -1 + last.rec + x1 + s(U, by=S),
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  
  gam.fitsMS[[iter]] <- list()
  
  # I mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 correctly specified
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + last.rec_mis + x1 + s(U, by=S),
                                 family = binomial,
                                 data=data)
  
  # II mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  gam.fitsMS[[iter]][[2]] <- gam(y~ -1 + last.rec + x4 + s(U, by=S),
                                 family = binomial,
                                 data=data)
  
  # III mispecification :
  # Time-related reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  gam.fitsMS[[iter]][[3]] <- gam(y~ -1 + last.rec + x1 + x2,
                                 family = binomial,
                                 data=data)
  
  # IV mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  gam.fitsMS[[iter]][[4]] <- gam(y~ -1 + last.rec_mis + x4 + s(U, by=S),
                                 family = binomial,
                                 data=data)
  
  # V mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  gam.fitsMS[[iter]][[5]] <- gam(y~ -1 + last.rec_mis + x1 + x2,
                                 family = binomial,
                                 data=data)
  # VI mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  gam.fitsMS[[iter]][[6]] <- gam(y~ -1 + last.rec + x4 + x2,
                                 family = binomial,
                                 data=data)
  
  # VII mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  gam.fitsMS[[iter]][[7]] <- gam(y~ -1 + last.rec_mis + x4 + x2,
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
BB.stat_9 <- BB.single(dim.k=9)[[2]]
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

pvalue.CS_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                                dat.gams[[x]],
                                                gam.fit = gam.fitsCS[[x]], 
                                                index = 3:11,
                                                BB.stat = BB.stat_9)[[1]])
plot(ecdf(pvalue.CS_comp3))
abline(0,1)

# Testing via global test

T_g_CS <- (tan(pi*(0.5-pvalue.CS_comp1)) + 
           tan(pi*(0.5-pvalue.CS_comp2)) +
           tan(pi*(0.5-pvalue.CS_comp3)))/3
pvalue_Cauchy_CS <- 1/2 - atan(T_g_CS)/pi
plot(ecdf(pvalue_Cauchy_CS), cex=0.5,
     main="3 components - correctly specified model")
abline(0,1)

pvalue.CS <- pvalue_Cauchy_CS

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power

# I mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 correctly specified

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

pvalue.MS1_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[1]], 
                                                 index = 3:11, 
                                                 BB.stat = BB.stat_9)[[1]])
plot(ecdf(pvalue.MS1_comp3))
abline(0,1)

# Testing via global test

T_g_MS1 <- (tan(pi*(0.5-pvalue.MS1_comp1)) + 
            tan(pi*(0.5-pvalue.MS1_comp2)) +
            tan(pi*(0.5-pvalue.MS1_comp3)))/3
pvalue_Cauchy_MS1 <- 1/2 - atan(T_g_MS1)/pi
plot(ecdf(pvalue_Cauchy_MS1), cex=0.5,
     main="3 components - I misspecified model")
abline(0,1)

pvalue.MS1 <- pvalue_Cauchy_MS1

# II mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified

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

pvalue.MS2_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[2]], 
                                                 index = 3:11, 
                                                 BB.stat = BB.stat_9)[[1]])
plot(ecdf(pvalue.MS2_comp3))
abline(0,1)

# Testing via global test

T_g_MS2 <- (tan(pi*(0.5-pvalue.MS2_comp1)) + 
            tan(pi*(0.5-pvalue.MS2_comp2)) +
            tan(pi*(0.5-pvalue.MS2_comp3)))/3
pvalue_Cauchy_MS2 <- 1/2 - atan(T_g_MS2)/pi
plot(ecdf(pvalue_Cauchy_MS2), cex=0.5,
     main="3 components - II misspecified model")
abline(0,1)

pvalue.MS2 <- pvalue_Cauchy_MS2

# III mispecification :
# Time-related reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)

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

pvalue.MS3_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[3]], 
                                                 index = 3)[[1]])
plot(ecdf(pvalue.MS3_comp3))
abline(0,1)

# Testing via global test

T_g_MS3 <- (tan(pi*(0.5-pvalue.MS3_comp1)) + 
            tan(pi*(0.5-pvalue.MS3_comp2)) +
            tan(pi*(0.5-pvalue.MS3_comp3)))/3
pvalue_Cauchy_MS3 <- 1/2 - atan(T_g_MS3)/pi
plot(ecdf(pvalue_Cauchy_MS3), cex=0.5,
     main="3 components - III misspecified model")
abline(0,1)

pvalue.MS3 <- pvalue_Cauchy_MS3

# IV mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified

# Inspecting single components

pvalue.MS4_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[4]], 
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS4_comp1))
abline(0,1)

pvalue.MS4_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[4]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS4_comp2))
abline(0,1)

pvalue.MS4_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[4]], 
                                                 index = 3:11, 
                                                 BB.stat = BB.stat_9)[[1]])
plot(ecdf(pvalue.MS4_comp3))
abline(0,1)

# Testing via global test

T_g_MS4 <- (tan(pi*(0.5-pvalue.MS4_comp1)) + 
            tan(pi*(0.5-pvalue.MS4_comp2)) +
            tan(pi*(0.5-pvalue.MS4_comp3)))/3
pvalue_Cauchy_MS4 <- 1/2 - atan(T_g_MS4)/pi
plot(ecdf(pvalue_Cauchy_MS4), cex=0.5,
     main="3 components - IV misspecified model")
abline(0,1)

pvalue.MS4 <- pvalue_Cauchy_MS4

# V mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)

# Inspecting single components
pvalue.MS5_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[5]],
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS5_comp1))
abline(0,1)

pvalue.MS5_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[5]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS5_comp2))
abline(0,1)

pvalue.MS5_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[5]], 
                                                 index = 3)[[1]])
plot(ecdf(pvalue.MS5_comp3))
abline(0,1)

# Testing via global test

T_g_MS5 <- (tan(pi*(0.5-pvalue.MS5_comp1)) + 
              tan(pi*(0.5-pvalue.MS5_comp2)) +
              tan(pi*(0.5-pvalue.MS5_comp3)))/3
pvalue_Cauchy_MS5 <- 1/2 - atan(T_g_MS5)/pi
plot(ecdf(pvalue_Cauchy_MS5), cex=0.5,
     main="3 components - V misspecified model")
abline(0,1)

pvalue.MS5 <- pvalue_Cauchy_MS5

# VI mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)

# Inspecting single components
pvalue.MS6_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[6]],
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS6_comp1))
abline(0,1)

pvalue.MS6_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[6]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS6_comp2))
abline(0,1)

pvalue.MS6_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[6]], 
                                                 index = 3)[[1]])
plot(ecdf(pvalue.MS6_comp3))
abline(0,1)

# Testing via global test

T_g_MS6 <- (tan(pi*(0.5-pvalue.MS6_comp1)) + 
            tan(pi*(0.5-pvalue.MS6_comp2)) +
            tan(pi*(0.5-pvalue.MS6_comp3)))/3
pvalue_Cauchy_MS6 <- 1/2 - atan(T_g_MS6)/pi
plot(ecdf(pvalue_Cauchy_MS6), cex=0.5,
     main="3 components - VI misspecified model")
abline(0,1)

pvalue.MS6 <- pvalue_Cauchy_MS6

# VII mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)

# Inspecting single components
pvalue.MS7_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[7]],
                                                 index = 1)[[1]])
plot(ecdf(pvalue.MS7_comp1))
abline(0,1)

pvalue.MS7_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[7]], 
                                                 index = 2)[[1]])
plot(ecdf(pvalue.MS7_comp2))
abline(0,1)

pvalue.MS7_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = 
                                                 dat.gams[[x]],
                                                 gam.fit = gam.fitsMS[[x]][[7]], 
                                                 index = 3)[[1]])
plot(ecdf(pvalue.MS7_comp3))
abline(0,1)


# Testing via global test

T_g_MS7 <- (tan(pi*(0.5-pvalue.MS7_comp1)) + 
            tan(pi*(0.5-pvalue.MS7_comp2)) +
            tan(pi*(0.5-pvalue.MS7_comp3)))/3
pvalue_Cauchy_MS7 <- 1/2 - atan(T_g_MS7)/pi
plot(ecdf(pvalue_Cauchy_MS7), cex=0.5,
     main="3 components - VII misspecified model")
abline(0,1)

pvalue.MS7 <- pvalue_Cauchy_MS7

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, 
     BB.stat_1, BB.stat_9, 
     pvalue.CS, pvalue.CS_comp1, pvalue.CS_comp2, pvalue.CS_comp3,
     pvalue.MS1, pvalue.MS1_comp1, pvalue.MS1_comp2, pvalue.MS1_comp3, 
     pvalue.MS2, pvalue.MS2_comp1, pvalue.MS2_comp2, pvalue.MS2_comp3,
     pvalue.MS3, pvalue.MS3_comp1, pvalue.MS3_comp2, pvalue.MS3_comp3,
     pvalue.MS4, pvalue.MS4_comp1, pvalue.MS4_comp2, pvalue.MS4_comp3,
     pvalue.MS5, pvalue.MS5_comp1, pvalue.MS5_comp2, pvalue.MS5_comp3,
     pvalue.MS6, pvalue.MS6_comp1, pvalue.MS6_comp2, pvalue.MS6_comp3,
     pvalue.MS7, pvalue.MS7_comp1, pvalue.MS7_comp2, pvalue.MS7_comp3,
     file="01-Simulation-Studies/03-Omnibus-Testing/output/global_three_elements.RData")
#-------------------------------------------------------------------------------