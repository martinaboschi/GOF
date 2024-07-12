#-------------------------------------------------------------------------------
# Script: Simulation_global_four_elements.R
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
# 3. Simulating the data with four elements driving the DGP
# (time-related definition of reciprocity with FLE, 
# covariate sampled from an exponential distribution with FLE,
# time-transformation of a covariate sampled form a normal with NLE, 
# random effect for sender activity)
dat.gams <- sapply(1:n.iter, function(x)
  data.simulation(n = n, p = p, seed = x+3000, L=4),
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
  
  s1 <- data$s1
  s2 <- data$s2
  ss <- factor(c(s1,s2))
  dim(ss) <- c(length(s1),2)
  Ls <- Ls_star <- matrix(0,length(s1),2); 
  
  Ls_star[,1] <- data$stp
  Ls_star[,2] <- -data$stp
  
  Ls[,1] <- 1
  Ls[,2] <- -1
  
  # 4.2. Correctly specified model: time-related reciprocity with FLE + x1
  #                                 + x2 with NLE
  gam.fitsCS[[iter]] <- gam(y~ -1 + last.rec + x1 + 
                              s(U, by=S) + s(ss, by=Ls, bs="re"),
                            family = binomial,
                            data=data)
  
  # 4.2. Misspecified model: 
  
  gam.fitsMS[[iter]] <- list()
  
  # I mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 correctly specified
  # re properly specified
  gam.fitsMS[[iter]][[1]] <- gam(y~ -1 + last.rec_mis + x1 + 
                                   s(U, by=S) + s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # II mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  # re properly specified
  gam.fitsMS[[iter]][[2]] <- gam(y~ -1 + last.rec + x4 + 
                                   s(U, by=S) + s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # III mispecification :
  # Time-related reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  # re properly specified
  gam.fitsMS[[iter]][[3]] <- gam(y~ -1 + last.rec + x1 + x2 + 
                                   s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # IV mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  # re properly specified
  gam.fitsMS[[iter]][[4]] <- gam(y~ -1 + last.rec_mis + x4 + 
                                   s(U, by=S) + s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # V mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  # re properly specified
  gam.fitsMS[[iter]][[5]] <- gam(y~ -1 + last.rec_mis + x1 + x2
                                 + s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  # VI mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  # re properly specified
  gam.fitsMS[[iter]][[6]] <- gam(y~ -1 + last.rec + x4 + x2 + 
                                   s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # VII mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  # re properly specified
  gam.fitsMS[[iter]][[7]] <- gam(y~ -1 + last.rec_mis + x4 + x2 + 
                                   s(ss, by=Ls, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # VIII mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 correctly specified
  # re misspecified
  gam.fitsMS[[iter]][[8]] <- gam(y~ -1 + last.rec_mis + x1 + 
                                   s(U, by=S) + s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # II mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  # re misspecified
  gam.fitsMS[[iter]][[9]] <- gam(y~ -1 + last.rec + x4 + 
                                   s(U, by=S) + s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # III mispecification :
  # Time-related reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  # re misspecified
  gam.fitsMS[[iter]][[10]] <- gam(y~ -1 + last.rec + x1 + x2 + 
                                   s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # IV mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 correctly specified
  # re misspecified
  gam.fitsMS[[iter]][[11]] <- gam(y~ -1 + last.rec_mis + x4 + 
                                   s(U, by=S) + s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # V mispecification :
  # misspecified reciprocity with FLE
  # x1 correctly specified
  # x2 misspecified (FLE instead of NLE)
  # re misspecified
  gam.fitsMS[[iter]][[12]] <- gam(y~ -1 + last.rec_mis + x1 + x2
                                 + s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  # VI mispecification :
  # Time-related reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  # re misspecified
  gam.fitsMS[[iter]][[13]] <- gam(y~ -1 + last.rec + x4 + x2 + 
                                   s(ss, by=Ls_star, bs="re"),
                                 family = binomial,
                                 data=data)
  
  # VII mispecification :
  # misspecified reciprocity with FLE
  # x1 misspecified (multiplied by time)
  # x2 misspecified (FLE instead of NLE)
  # re misspecified
  gam.fitsMS[[iter]][[14]] <- gam(y~ -1 + last.rec_mis + x4 + x2 + 
                                   s(ss, by=Ls_star, bs="re"),
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
BB.stat_50 <- BB.single(dim.k=50)[[2]]
BB.stat_9 <- BB.single(dim.k=9)[[2]]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 7. Testing the correctly specified model for coverage

# Inspecting single components

coefficients(objects[[1]][[2]])[1]
pvalue.CS_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                               gam.fit = objects[[x]][[2]], 
                                                               index = 1)[[1]])

coefficients(objects[[1]][[2]])[2]
pvalue.CS_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                               gam.fit = objects[[x]][[2]], 
                                                               index = 2)[[1]])

coefficients(objects[[1]][[2]])[3:11]
pvalue.CS_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[2]], 
                                                                 index = 3:11,
                                                                 BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[2]])[12:61]
pvalue.CS_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[2]],
                                                                 index = 12:61,
                                                                 BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_CS <- (tan(pi*(0.5-pvalue.CS_comp1)) + 
           tan(pi*(0.5-pvalue.CS_comp2)) + 
           tan(pi*(0.5-pvalue.CS_comp3)) +
           tan(pi*(0.5-pvalue.CS_comp4)))/4

pvalue_Cauchy_CS <- 1/2 - atan(T_g_CS)/pi
plot(ecdf(pvalue_Cauchy_CS), cex=0.5,
     main="4 components - correctly specified model")
abline(0,1)
pvalue.CS <- pvalue_Cauchy_CS

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. Testing the misspecified model for power

# I mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 correctly specified
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[1]])[1]
pvalue.MS1_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[1]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[1]])[2]
pvalue.MS1_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[1]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[1]])[3:11]
pvalue.MS1_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[1]],
                                                                  index = 3:11,
                                                                  BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[1]])[12:61]
pvalue.MS1_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[1]],
                                                                  index = 12:61,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS1 <- (tan(pi*(0.5-pvalue.MS1_comp1)) + 
            tan(pi*(0.5-pvalue.MS1_comp2)) + 
            tan(pi*(0.5-pvalue.MS1_comp3)) +
            tan(pi*(0.5-pvalue.MS1_comp4)))/4

pvalue_Cauchy_MS1 <- 1/2 - atan(T_g_MS1)/pi
plot(ecdf(pvalue_Cauchy_MS1), cex=0.5,
     main="4 components - misspecified reciprocity with FLE")

pvalue.MS1 <- pvalue_Cauchy_MS1

# II mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[2]])[1]
pvalue.MS2_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[2]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[2]])[2]
pvalue.MS2_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[2]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[2]])[3:11]
pvalue.MS2_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[2]],
                                                                  index = 3:11,
                                                                  BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[2]])[12:61]
pvalue.MS2_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[2]],
                                                                  index = 12:61,
                                                                  BB.stat = BB.stat_50)[[1]])
# Testing via global test

T_g_MS2 <- (tan(pi*(0.5-pvalue.MS2_comp1)) + 
            tan(pi*(0.5-pvalue.MS2_comp2)) + 
            tan(pi*(0.5-pvalue.MS2_comp3)) +
            tan(pi*(0.5-pvalue.MS2_comp4)))/4

pvalue_Cauchy_MS2 <- 1/2 - atan(T_g_MS2)/pi
plot(ecdf(pvalue_Cauchy_MS2), cex=0.5,
     main="4 components - x1 misspecified (multiplied by time)")
abline(0,1)

pvalue.MS2 <- pvalue_Cauchy_MS2

# III mispecification :
# Time-related reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[3]])[1]
pvalue.MS3_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[10]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[3]])[2]
pvalue.MS3_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[10]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[3]])[3]
pvalue.MS3_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[10]],
                                                                index = 3)[[1]])

coefficients(objects[[1]][[3]][[3]])[4:53]
pvalue.MS3_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[10]],
                                                                  index = 4:53,
                                                                  BB.stat = BB.stat_50)[[1]])
# Testing via global test

T_g_MS3 <- (tan(pi*(0.5-pvalue.MS3_comp1)) + 
            tan(pi*(0.5-pvalue.MS3_comp2)) + 
            tan(pi*(0.5-pvalue.MS3_comp3)) +
            tan(pi*(0.5-pvalue.MS3_comp4)))/4

pvalue_Cauchy_MS3 <- 1/2 - atan(T_g_MS3)/pi
plot(ecdf(pvalue_Cauchy_MS3), cex=0.5,
     main="4 components - x2 misspecified (FLE instead of NLE)")
abline(0,1)

pvalue.MS3 <- pvalue_Cauchy_MS3

# IV mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[4]])[1]
pvalue.MS4_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[4]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[4]])[2]
pvalue.MS4_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[4]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[4]])[3:11]
pvalue.MS4_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[4]],
                                                                  index = 3:11,
                                                                  BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[4]])[12:61]
pvalue.MS4_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[4]],
                                                                  index = 12:61,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS4 <- (tan(pi*(0.5-pvalue.MS4_comp1)) + 
            tan(pi*(0.5-pvalue.MS4_comp2)) + 
            tan(pi*(0.5-pvalue.MS4_comp3)) +
            tan(pi*(0.5-pvalue.MS4_comp4)))/4

pvalue_Cauchy_MS4 <- 1/2 - atan(T_g_MS4)/pi
plot(ecdf(pvalue_Cauchy_MS4), cex=0.5,
     main="4 components - misspecified reciprocity with FLE and x1")
abline(0,1)

pvalue.MS4 <- pvalue_Cauchy_MS4

# V mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[5]])[1]
pvalue.MS5_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[5]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[5]])[2]
pvalue.MS5_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[5]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[5]])[3]
pvalue.MS5_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[5]],
                                                                index = 3)[[1]])

coefficients(objects[[1]][[3]][[5]])[4:53]
pvalue.MS5_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[5]],
                                                                  index = 4:53,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS5 <- (tan(pi*(0.5-pvalue.MS5_comp1)) + 
            tan(pi*(0.5-pvalue.MS5_comp2)) + 
            tan(pi*(0.5-pvalue.MS5_comp3)) +
            tan(pi*(0.5-pvalue.MS5_comp4)))/4

pvalue_Cauchy_MS5 <- 1/2 - atan(T_g_MS5)/pi
plot(ecdf(pvalue_Cauchy_MS5), cex=0.5,
     main="4 components - misspecified reciprocity with FLE and x2")
abline(0,1)

pvalue.MS5 <- pvalue_Cauchy_MS5

# VI mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[6]])[1]
pvalue.MS6_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[6]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[6]])[2]
pvalue.MS6_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[6]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[6]])[3]
pvalue.MS6_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[6]],
                                                                index = 3)[[1]])

coefficients(objects[[1]][[3]][[6]])[4:53]
pvalue.MS6_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[6]],
                                                                  index = 4:53,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS6 <- (tan(pi*(0.5-pvalue.MS6_comp1)) + 
            tan(pi*(0.5-pvalue.MS6_comp2)) + 
            tan(pi*(0.5-pvalue.MS6_comp3)) +
            tan(pi*(0.5-pvalue.MS6_comp4)))/4

pvalue_Cauchy_MS6 <- 1/2 - atan(T_g_MS6)/pi
plot(ecdf(pvalue_Cauchy_MS6), cex=0.5,
     main="4 components - misspecified x1 and x2")
abline(0,1)

pvalue.MS6 <- pvalue_Cauchy_MS6

# VII mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)
# re properly specified

# Inspecting single components

coefficients(objects[[1]][[3]][[7]])[1]
pvalue.MS7_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[7]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[7]])[2]
pvalue.MS7_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[7]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[7]])[3]
pvalue.MS7_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[7]],
                                                                index = 3)[[1]])

coefficients(objects[[1]][[3]][[7]])[4:53]
pvalue.MS7_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[7]],
                                                                  index = 4:53,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS7 <- (tan(pi*(0.5-pvalue.MS7_comp1)) + 
            tan(pi*(0.5-pvalue.MS7_comp2)) + 
            tan(pi*(0.5-pvalue.MS7_comp3)) +
            tan(pi*(0.5-pvalue.MS7_comp4)))/4

pvalue_Cauchy_MS7 <- 1/2 - atan(T_g_MS7)/pi
plot(ecdf(pvalue_Cauchy_MS7), cex=0.5,
     main="4 components - misspecified reciprocity with FLE, x1, and x2")
abline(0,1)

pvalue.MS7 <- pvalue_Cauchy_MS7

# VIII mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 correctly specified
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[8]])[1]
pvalue.MS8_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[8]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[8]])[2]
pvalue.MS8_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[8]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[8]])[3:11]
pvalue.MS8_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[8]],
                                                                  index = 3:11,
                                                                  BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[8]])[12:61]
pvalue.MS8_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[8]],
                                                                  index = 12:61,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS8 <- (tan(pi*(0.5-pvalue.MS8_comp1)) + 
            tan(pi*(0.5-pvalue.MS8_comp2)) + 
            tan(pi*(0.5-pvalue.MS8_comp3)) +
            tan(pi*(0.5-pvalue.MS8_comp4)))/4

pvalue_Cauchy_MS8 <- 1/2 - atan(T_g_MS8)/pi
plot(ecdf(pvalue_Cauchy_MS8), cex=0.5,
     main="4 components - misspecified reciprocity with FLE")
abline(0,1)

pvalue.MS8 <- pvalue_Cauchy_MS8

# IX mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[9]])[1]
pvalue.MS9_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[9]], 
                                                                index = 1)[[1]])

coefficients(objects[[1]][[3]][[9]])[2]
pvalue.MS9_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                gam.fit = objects[[x]][[3]][[9]],
                                                                index = 2)[[1]])

coefficients(objects[[1]][[3]][[9]])[3:11]
pvalue.MS9_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[9]],
                                                                  index = 3:11,
                                                                  BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[9]])[12:61]
pvalue.MS9_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                  gam.fit = objects[[x]][[3]][[9]],
                                                                  index = 12:61,
                                                                  BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS9 <- (tan(pi*(0.5-pvalue.MS9_comp1)) + 
            tan(pi*(0.5-pvalue.MS9_comp2)) + 
            tan(pi*(0.5-pvalue.MS9_comp3)) +
            tan(pi*(0.5-pvalue.MS9_comp4)))/4

pvalue_Cauchy_MS9 <- 1/2 - atan(T_g_MS9)/pi
plot(ecdf(pvalue_Cauchy_MS9), cex=0.5,
     main="4 components - x1 misspecified (multiplied by time)")
abline(0,1)

pvalue.MS9 <- pvalue_Cauchy_MS9

# X mispecification :
# Time-related reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)
# misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[10]])[1]
pvalue.MS10_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[3]], 
                                                                 index = 1)[[1]])

coefficients(objects[[1]][[3]][[10]])[2]
pvalue.MS10_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[3]],
                                                                 index = 2)[[1]])

coefficients(objects[[1]][[3]][[10]])[3]
pvalue.MS10_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[3]],
                                                                 index = 3)[[1]])

coefficients(objects[[1]][[3]][[10]])[4:53]
pvalue.MS10_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[3]],
                                                                   index = 4:53,
                                                                   BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS10 <- (tan(pi*(0.5-pvalue.MS10_comp1)) + 
             tan(pi*(0.5-pvalue.MS10_comp2)) + 
             tan(pi*(0.5-pvalue.MS10_comp3)) +
             tan(pi*(0.5-pvalue.MS10_comp4)))/4

pvalue_Cauchy_MS10 <- 1/2 - atan(T_g_MS10)/pi
plot(ecdf(pvalue_Cauchy_MS10), cex=0.5,
     main="4 components - x2 misspecified (FLE instead of NLE)")
abline(0,1)

pvalue.MS10 <- pvalue_Cauchy_MS10

# XI mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 correctly specified
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[11]])[1]
pvalue.MS11_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[11]], 
                                                                 index = 1)[[1]])

coefficients(objects[[1]][[3]][[11]])[2]
pvalue.MS11_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[11]],
                                                                 index = 2)[[1]])

coefficients(objects[[1]][[3]][[11]])[3:11]
pvalue.MS11_comp3 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[11]],
                                                                   index = 3:11,
                                                                   BB.stat = BB.stat_9)[[1]])

coefficients(objects[[1]][[3]][[11]])[12:61]
pvalue.MS11_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[11]],
                                                                   index = 12:61,
                                                                   BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS11 <- (tan(pi*(0.5-pvalue.MS11_comp1)) + 
             tan(pi*(0.5-pvalue.MS11_comp2)) + 
             tan(pi*(0.5-pvalue.MS11_comp3)) +
             tan(pi*(0.5-pvalue.MS11_comp4)))/4

pvalue_Cauchy_MS11 <- 1/2 - atan(T_g_MS11)/pi
plot(ecdf(pvalue_Cauchy_MS11), cex=0.5,
     main="4 components - misspecified reciprocity with FLE and x1")
abline(0,1)

pvalue.MS11 <- pvalue_Cauchy_MS11

# XII mispecification :
# misspecified reciprocity with FLE
# x1 correctly specified
# x2 misspecified (FLE instead of NLE)
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[12]])[1]
pvalue.MS12_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[12]], 
                                                                 index = 1)[[1]])

coefficients(objects[[1]][[3]][[12]])[2]
pvalue.MS12_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[12]],
                                                                 index = 2)[[1]])

coefficients(objects[[1]][[3]][[12]])[3]
pvalue.MS12_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[12]],
                                                                 index = 3)[[1]])

coefficients(objects[[1]][[3]][[12]])[4:53]
pvalue.MS12_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[12]],
                                                                   index = 4:53,
                                                                   BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS12 <- (tan(pi*(0.5-pvalue.MS12_comp1)) + 
             tan(pi*(0.5-pvalue.MS12_comp2)) + 
             tan(pi*(0.5-pvalue.MS12_comp3)) +
             tan(pi*(0.5-pvalue.MS12_comp4)))/4

pvalue_Cauchy_MS12 <- 1/2 - atan(T_g_MS12)/pi
plot(ecdf(pvalue_Cauchy_MS12), cex=0.5,
     main="4 components - misspecified reciprocity with FLE and x2")
abline(0,1)

pvalue.MS12 <- pvalue_Cauchy_MS12

# XIII mispecification :
# Time-related reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[13]])[1]
pvalue.MS13_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[13]], 
                                                                 index = 1)[[1]])

coefficients(objects[[1]][[3]][[13]])[2]
pvalue.MS13_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[13]],
                                                                 index = 2)[[1]])

coefficients(objects[[1]][[3]][[13]])[3]
pvalue.MS13_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[13]],
                                                                 index = 3)[[1]])

coefficients(objects[[1]][[3]][[13]])[4:53]
pvalue.MS13_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[13]],
                                                                   index = 4:53,
                                                                   BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS13 <- (tan(pi*(0.5-pvalue.MS13_comp1)) + 
             tan(pi*(0.5-pvalue.MS13_comp2)) + 
             tan(pi*(0.5-pvalue.MS13_comp3)) +
             tan(pi*(0.5-pvalue.MS13_comp4)))/4

pvalue_Cauchy_MS13 <- 1/2 - atan(T_g_MS13)/pi
plot(ecdf(pvalue_Cauchy_MS13), cex=0.5,
     main="4 components - misspecified x1 and x2")
abline(0,1)

pvalue.MS13 <- pvalue_Cauchy_MS13

# XIV mispecification :
# misspecified reciprocity with FLE
# x1 misspecified (multiplied by time)
# x2 misspecified (FLE instead of NLE)
# re misspecified

# Inspecting single components

coefficients(objects[[1]][[3]][[14]])[1]
pvalue.MS14_comp1 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[14]], 
                                                                 index = 1)[[1]])

coefficients(objects[[1]][[3]][[14]])[2]
pvalue.MS14_comp2 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[14]],
                                                                 index = 2)[[1]])

coefficients(objects[[1]][[3]][[14]])[3]
pvalue.MS14_comp3 <- sapply(1:n.iter, function(x) GOF_univariate(data = objects[[x]][[1]],
                                                                 gam.fit = objects[[x]][[3]][[14]],
                                                                 index = 3)[[1]])

coefficients(objects[[1]][[3]][[14]])[4:53]
pvalue.MS14_comp4 <- sapply(1:n.iter, function(x) GOF_multivariate(data = objects[[x]][[1]],
                                                                   gam.fit = objects[[x]][[3]][[14]],
                                                                   index = 4:53,
                                                                   BB.stat = BB.stat_50)[[1]])

# Testing via global test

T_g_MS14 <- (tan(pi*(0.5-pvalue.MS14_comp1)) + 
             tan(pi*(0.5-pvalue.MS14_comp2)) + 
             tan(pi*(0.5-pvalue.MS14_comp3)) +
             tan(pi*(0.5-pvalue.MS14_comp4)))/4

pvalue_Cauchy_MS14 <- 1/2 - atan(T_g_MS14)/pi
plot(ecdf(pvalue_Cauchy_MS14), cex=0.5,
     main="4 components - misspecified reciprocity with FLE, x1, and x2")
abline(0,1)

pvalue.MS14 <- pvalue_Cauchy_MS14

#-------------------------------------------------------------------------------
# 9. Storing the results
save(objects, BB.stat_1, BB.stat_9, BB.stat_50,
     pvalue.CS,
     pvalue.MS1,
     pvalue.MS2,
     pvalue.MS3,
     pvalue.MS4,
     pvalue.MS5,
     pvalue.MS6,
     pvalue.MS7,
     pvalue.MS8,
     pvalue.MS9,
     pvalue.MS10,
     pvalue.MS11,
     pvalue.MS12,
     pvalue.MS13,
     pvalue.MS14,
     pvalue.CS_comp1,
     pvalue.CS_comp2,
     pvalue.CS_comp3,
     pvalue.CS_comp4,
     pvalue.MS1_comp1,
     pvalue.MS1_comp2,
     pvalue.MS1_comp3,
     pvalue.MS1_comp4,
     pvalue.MS2_comp1,
     pvalue.MS2_comp2,
     pvalue.MS2_comp3,
     pvalue.MS2_comp4,
     pvalue.MS3_comp1,
     pvalue.MS3_comp2,
     pvalue.MS3_comp3,
     pvalue.MS3_comp4,
     pvalue.MS4_comp1,
     pvalue.MS4_comp2,
     pvalue.MS4_comp3,
     pvalue.MS4_comp4,
     pvalue.MS5_comp1,
     pvalue.MS5_comp2,
     pvalue.MS5_comp3,
     pvalue.MS5_comp4,
     pvalue.MS6_comp1,
     pvalue.MS6_comp2,
     pvalue.MS6_comp3,
     pvalue.MS6_comp4,
     pvalue.MS7_comp1,
     pvalue.MS7_comp2,
     pvalue.MS7_comp3,
     pvalue.MS7_comp4,
     pvalue.MS8_comp1,
     pvalue.MS8_comp2,
     pvalue.MS8_comp3,
     pvalue.MS8_comp4,
     pvalue.MS9_comp1,
     pvalue.MS9_comp2,
     pvalue.MS9_comp3,
     pvalue.MS9_comp4,
     pvalue.MS10_comp1,
     pvalue.MS10_comp2,
     pvalue.MS10_comp3,
     pvalue.MS10_comp4,
     pvalue.MS11_comp1,
     pvalue.MS11_comp2,
     pvalue.MS11_comp3,
     pvalue.MS11_comp4,
     pvalue.MS12_comp1,
     pvalue.MS12_comp2,
     pvalue.MS12_comp3,
     pvalue.MS12_comp4,
     pvalue.MS13_comp1,
     pvalue.MS13_comp2,
     pvalue.MS13_comp3,
     pvalue.MS13_comp4,
     pvalue.MS14_comp1,
     pvalue.MS14_comp2,
     pvalue.MS14_comp3,
     pvalue.MS14_comp4,
     file="01-Simulation-Studies/03-Omnibus-Testing/output/global_four_elements.RData")
#-------------------------------------------------------------------------------