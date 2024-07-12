#-------------------------------------------------------------------------------
# Script: Functions.R
# Author: Martina Boschi
# Date: May 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 0. Loading required packages ####
library(mgcv)
library(mgcViz)
library(RColorBrewer)
library(gridExtra)
library(sde)
library(Matrix)
library(expm)
library(gratia)
library(e1071)
library(extraDistr)
library(lubridate)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Creating palettes for visualization ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2. Data-simulation functions for global testing ####
# Data simulations with up to 4 groups of elements involved in the model matrix:
# - Time-based covariate for reciprocity
#   (+ computation of a wrong variation of reciprocity);
# - x1 (exponential)
#   (+ computation of a variation of x1, multiplied by time, named x4);
# - Time-dependent x2 (obtained as a normal + current value of time) 
#   with a non-linear effect (log of absolute value);
# - Random effect for sender activity;

# 1. Function data.simulation_one_element
# Includes the computation of time-based covariates for reciprocity
# and the computation of a wrong variation of reciprocity;
# but only correct reciprocity is driving the actual data generating process
data.simulation_one_element <- function(n, p, seed){
  
  # Set number of senders equal to the number of receivers
  s = r = p
  
  # Set seed for reproducibility
  set.seed(seed)  
  
  # Initialize time
  t <- 0  
  
  # Initialize last interaction times for all pairs
  last_time <- last_time_mis <- rep(-Inf, p^2)  
  
  # Function to compute time-based covariate for reciprocity
  last.rec <- function(t){exp(-(t - as.vector(t(matrix(last_time, 
                                                       ncol=p, nrow=p)))))}
  
  # Function to compute time-based covariate for reciprocity (mispecified)
  last.rec_mis <- function(t){log(t/(t-as.vector(t(matrix(last_time_mis, 
                                                        ncol=p, nrow=p))))+1)}
  
  # Initial rate update based on reciprocity
  rate.upd <- 3 * last.rec(t)  
  
  # Initialize data storage matrix
  dat.gam <- matrix(NA, nrow=n, ncol=12)  
  
  # Initialize base rate
  rate <- -2
  
  # Calculate total rate
  tot.rate <- exp(rate + rate.upd)  
  tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
  # No self-interactions
  diag(tot.rate) <- 0  
  tot.rate <- as.vector(tot.rate)
  
  for (i in 1:n){
    
    # Draw interarrival time from exponential distribution
    dt <- rexp(1, sum(tot.rate))  
    
    # Update current time
    t <- t + dt  
    
    # Calculate probability of each dyad
    prob <- tot.rate / sum(tot.rate)  
    
    # Sample interacting dyad based on calculated probabilities
    sr.ev <- sample(s * r, 1, prob = prob)  
    s.ev <- (sr.ev - 1) %% s + 1  
    r.ev <- (sr.ev - 1) %/% s + 1 
    
    # Sample a non-interacting pair
    sr.nv <- sample(setdiff((1:(r * s)), sr.ev), 1)  
    s.nv <- (sr.nv - 1) %% s + 1 
    r.nv <- (sr.nv - 1) %/% s + 1 
    
    # Record event and non-event information in the matrix
    dat.gam[i,] <- c(1, t, s.ev, r.ev, s.nv, r.nv,
                     last.rec(t)[sr.ev],
                     last.rec(t)[sr.nv],
                     last.rec(t)[sr.ev] - last.rec(t)[sr.nv], 
                     last.rec_mis(t)[sr.ev],
                     last.rec_mis(t)[sr.nv],
                     last.rec_mis(t)[sr.ev] - last.rec_mis(t)[sr.nv])
    
    # Update last interaction time
    # Remark: if the event already occurred the clock is restarted 
    # Namely, if (s,r) was already occurred in the past
    # (r,s) is not anymore reciprocal until 
    # (s,r) is not anymore repeated until 
    # another (s,r) occurs
    last_time[sr.ev] <- ifelse(last_time[sr.ev] != -Inf, -Inf, t)
    last_time_mis[sr.ev] <- t
    
    # Update rate based on new last interaction time
    rate.upd <- 3 * last.rec(t)  
    tot.rate <- exp(rate + rate.upd)
    tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
    # Ensure no self-interactions
    diag(tot.rate) <- 0  
    tot.rate <- as.vector(tot.rate)
  }
  
  # Convert matrix to data frame with appropriate column names
  dat.gam <- data.frame(y = dat.gam[,1],
                        stp = dat.gam[,2],
                        s1 = dat.gam[,3],
                        r1 = dat.gam[,4],
                        s2 = dat.gam[,5],
                        r2 = dat.gam[,6],
                        last.rec1 = dat.gam[,7],
                        last.rec2 = dat.gam[,8],
                        last.rec = dat.gam[,9],
                        last.rec_mis1 = dat.gam[,10],
                        last.rec_mis2 = dat.gam[,11],
                        last.rec_mis = dat.gam[,12])
  
  # Return the simulated data frame
  return(dat.gam)  
}

# 2. Function data.simulation_two_elements
# Includes the computation of time-based cov. for (correct&wrong) reciprocity
# and two covariates: x1 sampled from an exponential
# and variation of x1, multiplied by time, named x4
# but only reciprocity and x1 are driving the actual data generating process
data.simulation_two_elements <- function(n, p, seed){
  
  # Set number of senders equal to the number of receivers
  s = r = p
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize time
  t <- 0
  
  # Initialize last interaction times for all pairs
  last_time <- last_time_mis <- rep(-Inf, p^2)  
  
  # Function to compute time-based covariate for reciprocity
  last.rec <- function(t){exp(-(t - as.vector(t(matrix(last_time, 
                                                       ncol=p, nrow=p)))))}
  
  # Function to compute time-based covariate for reciprocity (mispecified)
  last.rec_mis <- function(t){log(t/(t-as.vector(t(matrix(last_time_mis, 
                                                          ncol=p, nrow=p))))+1)}
  
  # Initial rate update based on reciprocity
  rate.upd <- 3 * last.rec(t)
  
  # Initialize covariates x1 and x4
  # Draw initial values for x1 from exponential distribution
  x1 <- rexp(p^2, 4)  
  # Copy x1 to x4 for initialization
  x4 <- x1  
  
  # Initialize rate with x1
  rate <- x1 - 1
  
  # Initialize data storage matrix
  dat.gam <- matrix(NA, nrow=n, ncol=18)
  
  # Calculate total rate
  tot.rate <- exp(rate + rate.upd)
  tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
  # No self-interactions
  diag(tot.rate) <- 0
  tot.rate <- as.vector(tot.rate)
  
  for (i in 1:n){
    
    # Draw interarrival time from exponential distribution
    dt <- rexp(1, sum(tot.rate))
    
    # Update current time
    t <- t + dt
    
    # Calculate probability of each dyad
    prob <- tot.rate / sum(tot.rate)
    
    # Sample interacting dyad based on calculated probabilities
    sr.ev <- sample(s * r, 1, prob = prob)
    s.ev <- (sr.ev - 1) %% s + 1
    r.ev <- (sr.ev - 1) %/% s + 1
    
    # Sample a non-interacting pair
    sr.nv <- sample(setdiff((1:(r * s)), sr.ev), 1)
    s.nv <- (sr.nv - 1) %% s + 1
    r.nv <- (sr.nv - 1) %/% s + 1
    
    # Record event and non-event information in the matrix
    dat.gam[i,] <- c(1, t, s.ev, r.ev, s.nv, r.nv,
                     last.rec(t)[sr.ev],
                     last.rec(t)[sr.nv],
                     last.rec(t)[sr.ev] - last.rec(t)[sr.nv],
                     last.rec_mis(t)[sr.ev],
                     last.rec_mis(t)[sr.nv],
                     last.rec_mis(t)[sr.ev] - last.rec_mis(t)[sr.nv],
                     x1[sr.ev], x1[sr.nv], 
                     x1[sr.ev] - x1[sr.nv], 
                     x4[sr.ev] * t, x4[sr.nv] * t, 
                     x4[sr.ev] * t - x4[sr.nv] * t)
    
    # Update last interaction time
    # Remark: if the event already occurred the clock is restarted 
    last_time[sr.ev] <- ifelse(last_time[sr.ev] != -Inf, -Inf, t)
    last_time_mis[sr.ev] <- t
    
    # Update rate based on new last interaction time
    rate.upd <- 3 * last.rec(t)
    tot.rate <- exp(rate + rate.upd)
    tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
    # Ensure no self-interactions
    diag(tot.rate) <- 0
    tot.rate <- as.vector(tot.rate)
  }
  
  # Convert matrix to data frame with appropriate column names
  dat.gam <- data.frame(y = dat.gam[,1],
                        stp = dat.gam[,2],
                        s1 = dat.gam[,3],
                        r1 = dat.gam[,4],
                        s2 = dat.gam[,5],
                        r2 = dat.gam[,6],
                        last.rec1 = dat.gam[,7],
                        last.rec2 = dat.gam[,8],
                        last.rec = dat.gam[,9],
                        last.rec_mis1 = dat.gam[,10],
                        last.rec_mis2 = dat.gam[,11],
                        last.rec_mis = dat.gam[,12],
                        x11 = dat.gam[,13],
                        x12 = dat.gam[,14],
                        x1 = dat.gam[,15], 
                        x41 = dat.gam[,16],
                        x42 = dat.gam[,17],
                        x4 = dat.gam[,18])
  
  # Return the simulated data frame
  return(dat.gam)  
}

# 3. Function data.simulation_three_elements
# Includes the computation of time-based covariates for reciprocity
# and three covariates: x1 sampled from an exponential,
# a variation of x1, multiplied by time, named x4,
# time-transformation of a covariate x2 sampled form a normal with NLE
# reciprocity, x1 and non-linear function of x2
# are driving the actual data generating process
data.simulation_three_elements <- function(n, p, seed){
  
  # Set number of senders equal to the number of receivers
  s = r = p
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize time
  t <- 0
  
  # Initialize last interaction times for all pairs
  last_time <- last_time_mis <- rep(-Inf, p^2)  
  
  # Function to compute time-based covariate for reciprocity
  last.rec <- function(t){exp(-(t - as.vector(t(matrix(last_time, 
                                                       ncol=p, nrow=p)))))}
  
  # Function to compute time-based covariate for reciprocity (mispecified)
  last.rec_mis <- function(t){log(t/(t-as.vector(t(matrix(last_time_mis, 
                                                          ncol=p, nrow=p))))+1)}
  
  # Initialize x2 from normal distribution
  x2 <- rnorm(p^2, 0.5, 1)
  
  # Function to compute x2 value over time
  x2.funct <- function(t, x2){
    t + t^2 +  as.vector(matrix(x2,ncol=p, nrow=p))
  }
  
  # Initial x2 values
  x2.value <- x2.funct(t, x2)
  
  # Initial rate update based on reciprocity
  rate.upd <- 3 * last.rec(t) + (cos(x2.value*2) + exp(abs(x2.value))/20)
  
  # Initialize covariates x1 and x4
  # Draw initial values for x1 from exponential distribution
  x1 <- rexp(p^2, 4)  
  
  # Copy x1 to x4 for initialization
  x4 <- x1
  
  # Initialize rate with x1
  rate <- 1.5*x1 - 2
  
  # Initialize data storage matrix
  dat.gam <- matrix(NA, nrow=n, ncol=21)
  
  # Calculate total rate
  tot.rate <- exp(rate + rate.upd)
  tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
  # No self-interactions
  diag(tot.rate) <- 0
  tot.rate <- as.vector(tot.rate)
  
  for (i in 1:n){

    # Draw interarrival time from exponential distribution
    dt <- rexp(1, sum(tot.rate))

    # Update current time
    t <- t + dt
    
    # Calculate probability of each dyad
    prob <- tot.rate / sum(tot.rate)
    
    # Sample interacting dyad based on calculated probabilities
    sr.ev <- sample(s * r, 1, prob = prob)
    s.ev <- (sr.ev - 1) %% s + 1
    r.ev <- (sr.ev - 1) %/% s + 1
    
    # Sample a non-interacting pair
    sr.nv <- sample(setdiff((1:(r * s)), sr.ev), 1)
    s.nv <- (sr.nv - 1) %% s + 1
    r.nv <- (sr.nv - 1) %/% s + 1
    
    # Record event and non-event information in the matrix
    dat.gam[i,] <- c(1, t, s.ev, r.ev, s.nv, r.nv,
                     last.rec(t)[sr.ev],
                     last.rec(t)[sr.nv],
                     last.rec(t)[sr.ev] - last.rec(t)[sr.nv],
                     last.rec_mis(t)[sr.ev],
                     last.rec_mis(t)[sr.nv],
                     last.rec_mis(t)[sr.ev] - last.rec_mis(t)[sr.nv],
                     x1[sr.ev], x1[sr.nv], 
                     x1[sr.ev] - x1[sr.nv], 
                     x4[sr.ev] * t, x4[sr.nv] * t, 
                     x4[sr.ev] * t - x4[sr.nv] * t,
                     x2.funct(t, x2)[sr.ev], x2.funct(t, x2)[sr.nv], 
                     x2.funct(t, x2)[sr.ev] - x2.funct(t, x2)[sr.nv])
    
    # Update last interaction time
    # Remark: if the event already occurred the clock is restarted 
    last_time[sr.ev] <- ifelse(last_time[sr.ev] != -Inf, -Inf, t)
    # For the misspecification version, reciprocity is not closed
    last_time_mis[sr.ev] <- t
    
    # Update x2 value
    x2.value <- x2.funct(t, x2)
    
    # Update rate based on new last interaction time and x2
    rate.upd <- 3 * last.rec(t) + (cos(x2.value*2) + exp(abs(x2.value))/20)
    
    tot.rate <- exp(rate + rate.upd)
    tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
    # Ensure no self-interactions
    diag(tot.rate) <- 0
    tot.rate <- as.vector(tot.rate)
  }
  
  # Convert matrix to data frame with appropriate column names
  dat.gam <- data.frame(y = dat.gam[,1],
                        stp = dat.gam[,2],
                        s1 = dat.gam[,3],
                        r1 = dat.gam[,4],
                        s2 = dat.gam[,5],
                        r2 = dat.gam[,6],
                        last.rec1 = dat.gam[,7],
                        last.rec2 = dat.gam[,8],
                        last.rec = dat.gam[,9],
                        last.rec_mis1 = dat.gam[,10],
                        last.rec_mis2 = dat.gam[,11],
                        last.rec_mis = dat.gam[,12],
                        x11 = dat.gam[,13],
                        x12 = dat.gam[,14],
                        x1 = dat.gam[,15], 
                        x41 = dat.gam[,16],
                        x42 = dat.gam[,17],
                        x4 = dat.gam[,18],
                        x21 = dat.gam[,19],
                        x22 = dat.gam[,20],
                        x2 = dat.gam[,21])
  
  # Return the simulated data frame
  return(dat.gam)  
}

# 4. Function data.simulation_four_elements
# Includes the computation of time-based covariates for reciprocity
# and three covariates: x1 sampled from an exponential,
# a variation of x1, multiplied by time, named x4,
# time-transformation of a covariate x2 sampled form a normal with NLE
# reciprocity, x1 and non-linear function of x2
# are driving the actual data generating process
# plus random effects for sender activity
# are driving the actual data generating process
data.simulation_four_elements <- function(n, p, seed, act.sd = 0.5){
  
  # Set number of senders equal to the number of receivers
  s = r = p
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize time
  t <- 0
  
  # Initialize last interaction times for all pairs
  last_time <- last_time_mis <- rep(-Inf, p^2)  
  
  # Function to compute time-based covariate for reciprocity
  last.rec <- function(t){exp(-(t - as.vector(t(matrix(last_time, 
                                                       ncol=p, nrow=p)))))}
  
  # Function to compute time-based covariate for reciprocity (mispecified)
  last.rec_mis <- function(t){log(t/(t-as.vector(t(matrix(last_time_mis, 
                                                          ncol=p, nrow=p))))+1)}
  
  # Initialize x2 from normal distribution
  x2 <- rnorm(p^2, 0.5, 1)
  
  # Function to compute x2 value over time
  x2.funct <- function(t, x2){
    t + t^2 +  as.vector(matrix(x2,ncol=p, nrow=p))
  }
  
  # Initial x2 values
  x2.value <- x2.funct(t, x2)
  
  # Initial rate update based on reciprocity
  rate.upd <- 3 * last.rec(t) + (cos(x2.value*2) + exp(abs(x2.value))/20)
  
  # Initialize covariates x1 and x4
  # Draw initial values for x1 from exponential distribution
  x1 <- rexp(p^2, 4)  
  
  # Copy x1 to x4 for initialization
  x4 <- x1
  
  # Initialize random effects for sender activity
  act.eff <- rep(rnorm(p, 0, act.sd), p)
  
  # Update rate with x1, x2, and random effects for sender activity
  rate <- 1.5*x1 + act.eff - 2
  
  # Initialize data storage matrix
  dat.gam <- matrix(NA, nrow=n, ncol=21)
  
  # Calculate total rate
  tot.rate <- exp(rate + rate.upd)
  tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
  # No self-interactions
  diag(tot.rate) <- 0
  tot.rate <- as.vector(tot.rate)
  
  for (i in 1:n){
  
    # Draw interarrival time from exponential distribution
    dt <- rexp(1, sum(tot.rate))
    
    # Update current time
    t <- t + dt
    
    # Calculate probability of each dyad
    prob <- tot.rate / sum(tot.rate)
    
    # Sample interacting dyad based on calculated probabilities
    sr.ev <- sample(s * r, 1, prob = prob)
    s.ev <- (sr.ev - 1) %% s + 1
    r.ev <- (sr.ev - 1) %/% s + 1
    
    # Sample a non-interacting pair
    sr.nv <- sample(setdiff((1:(r * s)), sr.ev), 1)
    s.nv <- (sr.nv - 1) %% s + 1
    r.nv <- (sr.nv - 1) %/% s + 1
    
    # Record event and non-event information in the matrix
    dat.gam[i,] <- c(1, t, s.ev, r.ev, s.nv, r.nv,
                     last.rec(t)[sr.ev],
                     last.rec(t)[sr.nv],
                     last.rec(t)[sr.ev] - last.rec(t)[sr.nv],
                     last.rec_mis(t)[sr.ev],
                     last.rec_mis(t)[sr.nv],
                     last.rec_mis(t)[sr.ev] - last.rec_mis(t)[sr.nv],
                     x1[sr.ev], x1[sr.nv], 
                     x1[sr.ev] - x1[sr.nv], 
                     x4[sr.ev] * t, x4[sr.nv] * t, 
                     x4[sr.ev] * t - x4[sr.nv] * t,
                     x2.funct(t, x2)[sr.ev], x2.funct(t, x2)[sr.nv], 
                     x2.funct(t, x2)[sr.ev] - x2.funct(t, x2)[sr.nv])
    
    # Update last interaction time
    # Remark: if the event already occurred the clock is restarted 
    last_time[sr.ev] <- ifelse(last_time[sr.ev] != -Inf, -Inf, t)
    last_time_mis[sr.ev] <- t
    
    # Update x2 value
    x2.value <- x2.funct(t, x2)
    
    # Update rate based on new last interaction time and x2
    rate.upd <- 3 * last.rec(t) + (cos(x2.value*2) + exp(abs(x2.value))/20)
    
    tot.rate <- exp(rate + rate.upd)
    tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
    # Ensure no self-interactions
    diag(tot.rate) <- 0
    tot.rate <- as.vector(tot.rate)
  }
  
  # Convert matrix to data frame with appropriate column names
  dat.gam <- data.frame(y = dat.gam[,1],
                        stp = dat.gam[,2],
                        s1 = dat.gam[,3],
                        r1 = dat.gam[,4],
                        s2 = dat.gam[,5],
                        r2 = dat.gam[,6],
                        last.rec1 = dat.gam[,7],
                        last.rec2 = dat.gam[,8],
                        last.rec = dat.gam[,9],
                        last.rec_mis1 = dat.gam[,10],
                        last.rec_mis2 = dat.gam[,11],
                        last.rec_mis = dat.gam[,12],
                        x11 = dat.gam[,13],
                        x12 = dat.gam[,14],
                        x1 = dat.gam[,15], 
                        x41 = dat.gam[,16],
                        x42 = dat.gam[,17],
                        x4 = dat.gam[,18],
                        x21 = dat.gam[,19],
                        x22 = dat.gam[,20],
                        x2 = dat.gam[,21])
  
  # Return the simulated data frame
  return(dat.gam)  
}

# Function data.simulation
# This function generates simulated data based on the number of elements (L) 
# involved in the model matrix.
# It calls different data simulation functions 
# based on the number of elements specified.
# Parameters:
#   - n: Number of observations to generate.
#   - p: Number of actors
#   - seed: Seed for reproducibility
#   - L: Number indicating the complexity of the model matrix:
#        1 for one element, 2 for two elements, 
#        3 for three elements, and 4 for four elements.
#   - act.sd: Standard deviation for the random effects of sender activity 
#     (default is 0.5).
# Returns:
#   - Simulated data frame generated based on the specified number of elements
data.simulation <- function(n, p, seed, L, act.sd=0.5){
  # Check the complexity of the model matrix
  if (L == 1){
    # Call the data simulation function for one element
    return(data.simulation_one_element(n, p, seed))
  } else {
    if (L == 2){
      # Call the data simulation function for two elements
      return(data.simulation_two_elements(n, p, seed))
    } else {
      if (L == 3){
        # Call the data simulation function for three elements
        return(data.simulation_three_elements(n, p, seed))
      } else {
        # Call the data simulation function for four elements
        return(data.simulation_four_elements(n, p, seed, act.sd))
      }
    }
  }
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Data-simulation function for one-covariate test with NLE ####
# Reciprocity is designed with an intrinsic non-linear behaviour
# and represents the only driver of the data generating process
data.generation_nonlinear_reciprocity <- function(n, p, 
                              # Weibull(scale=1, shape=2) ends in hazard h(t)=2t
                              l.wb = 1, k = 2, 
                              l.exp = 0.1,
                              seed = 1234, 
                              b = 1,
                              avoid.slf = TRUE){
  
  # Initialize variables for reciprocity
  # Array to store last reciprocity time for each pair
  rec.tm <- rep(-Inf, p^2)  
  # Array to store reciprocity status for each pair
  rec <- rep(0, p^2)        
  
  # Initialize starting time
  start <- 0
  
  # Generate initial spontaneous events
  set.seed(seed)
  tms <- start + rexp(p^2, rate = l.exp)
  
  # Avoid self-loops by setting times to infinity
  if (avoid.slf) {
    slf <- (0:(p - 1)) * p + (1:p)  # Indices representing self-loops
    tms[slf] <- Inf
  } else {
    slf <- NULL
  }
  
  # Initialize data storage
  dat <- NULL 
  
  # Number of non-events
  m <- 1
  
  for (i in 1:n) {
    
    # Identify the event with the minimum time
    event <- which.min(tms)
    # Time of the event
    tm <- min(tms)  
    # Reciprocity status of the event
    ev.rec <- rec[event]  
    # Time since last reciprocity for the event
    ev.rectm <- tm - rec.tm[event]
    # Sender of the event
    s.e <- (event - 1) %% p + 1  
    # Receiver of the event
    r.e <- (event - 1) %/% p + 1  
    
    # Randomly select non-events
    nonevents <- sample((1:p^2)[-c(event, slf)], m, replace = FALSE)  
    # Time since last reciprocity for non-events
    nonev.rectm <- tm - rec.tm[nonevents]  
    # Reciprocity status of non-events
    nonev.rec <- rec[nonevents]  
    # Senders of non-events
    s.n <- (nonevents - 1) %% p + 1  
    # Receivers of non-events
    r.n <- (nonevents - 1) %/% p + 1  
    
    # Update time for the current dyad
    tms[event] <- tm + rexp(1, l.exp)
    
    # Reset reciprocity status if the event was reciprocal
    if (ev.rec == 1) {
      rec[event] <- 0
      rec.tm[event] <- -Inf
    }
    
    # Set reciprocity for the event (r.e, s.e)
    new.rec <- (s.e - 1) * p + r.e
    # Remove self-loops
    new.rec <- setdiff(new.rec, slf)  
    # Set reciprocity status
    rec[new.rec] <- 1 
    # Update last reciprocity time
    rec.tm[new.rec] <- tm  
    # Generate new reciprocity times
    tms[new.rec] <- tm + rweibull(length(new.rec), k, l.wb)  
    
    # Add data entry to the dataset
    dat <- rbind(dat, cbind(rep(tm, m), rep(s.e, m), rep(r.e, m),
                            rep(ev.rectm, m), rep(ev.rec, m),
                            s.n, r.n, nonev.rectm, nonev.rec,
                            rep(event, m)))
  }
  
  # Add exponential decay for reciprocity times
  dat <- cbind(dat, exp(-b * dat[, 4]))
  dat <- cbind(dat, exp(-b * dat[, 8]))
  
  # Create the final dataset
  data <- data.frame(dat)
  # Remove the column representing event indices
  data <- data[, -10]  
  # Rename columns
  colnames(data) <- c("stp", "s1", "r1", "recTime1", "recId1", 
                      "s2", "r2", "recTime2", "recId2", 
                      "recExpTime1", "recExpTime2")  
  # Compute difference in reciprocity times
  data$recExpTime <- data$recExpTime1 - data$recExpTime2  
  
  # Return the final dataset
  return(data)  
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 4. Data-simulation function for random effects ####
# Sender activity represents the only driver of the data generating process
data.generation_random_effects <- function(n, p, 
                                           act.sd,
                                           seed=1234,
                                           avoid.slf=TRUE
){
  
  # Set number of senders equal to the number of receivers
  s = r = p
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Initialize time
  t <- 0
  
  # Initialize data storage matrix
  dat.gam <- matrix(NA, nrow=n, ncol=6)
  
  # Initialize random effects for sender activity
  act.eff <- rep(rnorm(p, 0, act.sd), p)
  
  # Update rate with sender activity
  rate <- act.eff
  
  # Calculate total rate
  tot.rate <- exp(rate)
  tot.rate <- matrix(tot.rate, ncol=p, nrow=p)
  # No self-interactions
  diag(tot.rate) <- 0
  tot.rate <- as.vector(tot.rate)
  
  # Simulation loop
  for (i in 1:n){
    
    # Draw interarrival time from exponential distribution
    dt <- rexp(1, sum(tot.rate))
    t <- t + dt
    
    # Calculate probability of each dyad
    prob <- tot.rate / sum(tot.rate)
    
    # Sample interacting dyad based on calculated probabilities
    sr.ev <- sample(s * r, 1, prob = prob)
    s.ev <- (sr.ev - 1) %% s + 1
    r.ev <- (sr.ev - 1) %/% s + 1
    
    # Sample a non-interacting pair
    sr.nv <- sample(setdiff((1:(r * s)), sr.ev), 1)
    s.nv <- (sr.nv - 1) %% s + 1
    r.nv <- (sr.nv - 1) %/% s + 1
    
    # Record event and non-event information in the matrix
    dat.gam[i,] <- c(1, t, s.ev, r.ev, s.nv, r.nv)
    
  }
  
  # Convert matrix to data frame with appropriate column names
  dat.gam <- data.frame(y = dat.gam[,1],
                        stp = dat.gam[,2],
                        s1 = dat.gam[,3],
                        r1 = dat.gam[,4],
                        s2 = dat.gam[,5],
                        r2 = dat.gam[,6])
  
  # Return the simulated data frame
  return(dat.gam)
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 5. Function aimed at uniforming different notations in a standard one ####
convert.to.objects <- function(dat.gams, gam.fitsCS, gam.fitsMS){
  objects <- list()
  for (iter in 1:n.iter){
    objects[[iter]] <- list()
    objects[[iter]][[1]] <- dat.gams[[iter]]
    objects[[iter]][[2]] <- gam.fitsCS[[iter]]
    objects[[iter]][[3]] <- gam.fitsMS[[iter]]
  }
  return(objects)
}

#-------------------------------------------------------------------------------
# 6. Function aimed at reproducing the asymptotic behaviour of MR process ####
BB.simulator <- function(set_covariates, n = 2000, nsim = 2000) {
  # Initialize an empty vector to store maximum values
  Bridge.data <- replicate(nsim, {
    max_s <- NULL
    # Loop through each set of covariates
    for (index in set_covariates) {
      
      # Ensure index is unlisted
      index = unlist(index)
      
      # Determine the dimension of the covariate set
      dim.k = length(index)
      
      # Generate n Brownian Bridge paths for each covariate
      BB.data <- replicate(dim.k, BBridge(0, 0, N = n - 1))
      
      # Calculate the process (squared norm) for each path and 
      # store the maximum value divided by the dimension
      process <- apply(BB.data, 1, function(x) crossprod(x, x))
      max_s <- c(max_s, max(process) / dim.k)
      
    }
    # Return the maximum value across all covariate sets
    return(max(max_s))
  })
  
  # Return the vector of maximum values
  return(Bridge.data)
}

BB.single <- function(dim.k, n.sim=5000, n=2000){
  Bridge.data <- replicate(n.sim, 
            {
              BB.data <- replicate(dim.k, BBridge(0,0,N=n-1))
              process <- apply(BB.data, 1, function(x) crossprod(x,x))
              return(process)
            })
  BB.stat <- apply(Bridge.data, 2, function(x) max(abs(x)))
  return(list(Bridge.data, BB.stat))
}

#-------------------------------------------------------------------------------
# 7. Approximating Kolmogorov Distribution ####
f <- function(x,i) {exp(-(2*i-1)^2*pi^2/(8*x^2)) }
kolm <- function(x) {sqrt(2*pi)/x*(f(x,1)+f(x,2)+f(x,3)+f(x,4)+
                                     f(x,5)+f(x,6)+f(x,7)+f(x,8)+
                                     f(x,9)+f(x,10))}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 8. GOF ####

## 8.1. GOF - Global Testing ####
GOF_pvalue_global <- function(data, gam.fit, set_covariates, BB.stat){
  
  # Get the number of observations
  n.e <- nrow(data)
  
  # Extract the design matrix from the GAM fit
  X <- model.matrix(gam.fit)
  # Calculate the mean response values
  mu <- gam.fit$fitted.values
  
  # Compute the working residuals
  w <- gam.fit$family$mu.eta(gam.fit$linear.predictors)*
    (gam.fit$y - mu)/(gam.fit$sig2*gam.fit$family$variance(mu))
  
  # Compute the weights for the process
  wt <- sqrt(mu*(1-mu))
  
  # Initialize an empty vector to store maximum statistics
  max_s <- NULL
  
  for (index in set_covariates){
    
    # Ensure index is unlisted
    index = unlist(index)
    
    # Determine the dimension of the covariate set
    dim.k = length(index)
    
    # Extract the columns corresponding to the covariate set from the design
    dim.spl = index
    
    if (dim.k == 1){
      # Psi for univariate covariates is a vector
      Psi <- X[,dim.spl] * w * wt
      Psi.last <- sum(Psi)
      Psi.tilde <- Psi - Psi.last/n.e
      W <- cumsum(Psi.tilde)/sqrt(n.e)
      # Inverse Variance of Psi
      Vs <- ginv(crossprod(Psi.tilde)/(n.e))
      # Squared Inverse Variance of Psi
      E <- eigen(Vs)
      P <- E$vectors
      D.val <- pmax(E$value, 0)
      D.sq <- (sqrt(D.val))
      B.inv.sq <- (P%*%D.sq%*%t(P))
      
    } else {
      # Psi for multivariate covariates is a matrix
      Psi <- X[,dim.spl] * w * wt
      Psi.last <- matrix(apply(Psi, 2, sum), ncol=dim.k)
      Psi.tilde <- Psi - apply(Psi.last, 2, 
                               function(x) x * rep(1/n.e, n.e))
      W <- apply(Psi.tilde, 2, cumsum)/sqrt(n.e)
      # Inverse Variance of Psi
      Vs <- ginv(crossprod(Psi.tilde)/(n.e))
      # Squared Inverse Variance of Psi
      E <- eigen(Vs)
      P <- E$vectors
      D.val <- pmax(E$value, 0)
      D.sq <- diag(sqrt(D.val))
      B.inv.sq <- (P%*%D.sq%*%t(P))
    }
    
    # Standardized Martingale-Residual Process
    efp <- t(B.inv.sq %*% t(W))
    # Compute the maximum statistics
    U.n <- apply(efp, 1, function(x) crossprod(x,x))
    max_s <- c(max_s, max(abs(U.n))/dim.k)
  }
  
  # Compute the test statistic
  stat = max(max_s)
  # Compute the p-value
  pvalue = mean(BB.stat>=stat)
  
  return(pvalue)
}

GOF_univariate <- function(data, gam.fit, index){
  
  # Get the number of observations
  n.e <- nrow(data)
  
  # Extract the design matrix from the GAM fit
  X <- model.matrix(gam.fit)
  # Calculate the mean response values
  mu <- gam.fit$fitted.values
  
  # Compute the working residuals
  w <- gam.fit$family$mu.eta(gam.fit$linear.predictors)*
    (gam.fit$y - mu)/(gam.fit$sig2*gam.fit$family$variance(mu))
  
  # Psi process 
  Psi <- X[,index] * w
  W <- cumsum(Psi)/sqrt(n.e)
  # Inverse Variance of Psi
  Vs <- ginv(crossprod(Psi)/(n.e))
  # Squared Inverse Variance of Psi
  E <- eigen(Vs)
  P <- E$vectors
  D.val <- pmax(E$value, 0)
  D.sq <- (sqrt(D.val))
  B.inv.sq <- (P%*%D.sq%*%t(P))
  
  # Standardized Martingale-Residual Process
  efp <- t(B.inv.sq %*% t(W))
  
  # Compute the p-value
  pvalue = 1 - kolm(max(abs(efp)))
  
  return(list(pvalue, efp))
}

GOF_multivariate <- function(data, gam.fit, index, BB.stat){
  
  # Get the number of observations
  n.e <- nrow(data)
  
  # Extract the design matrix from the GAM fit
  X <- model.matrix(gam.fit)
  # Calculate the mean response values
  mu <- gam.fit$fitted.values
  
  # Compute the working residuals
  w <- gam.fit$family$mu.eta(gam.fit$linear.predictors)*
    (gam.fit$y - mu)/(gam.fit$sig2*gam.fit$family$variance(mu))
  
  # Compute the weights for the process
  wt <- sqrt(mu*(1-mu))
  
  # Determine the dimension of the covariate set
  dim.k = length(index)
  
  # Extract the columns corresponding to the covariate set from the design
  dim.spl = index
  
  # Psi for multivariate covariates is a matrix
  Psi <- X[,dim.spl] * w * wt
  Psi.last <- matrix(apply(Psi, 2, sum), ncol=dim.k)
  Psi.tilde <- Psi - apply(Psi.last, 2, 
                           function(x) x * rep(1/n.e, n.e))
  W <- apply(Psi.tilde, 2, cumsum)/sqrt(n.e)
  # Inverse Variance of Psi
  Vs <- ginv(crossprod(Psi.tilde)/(n.e))
  # Squared Inverse Variance of Psi
  E <- eigen(Vs)
  P <- E$vectors
  D.val <- pmax(E$value, 0)
  D.sq <- diag(sqrt(D.val))
  B.inv.sq <- (P%*%D.sq%*%t(P))
  
  # Standardized Martingale-Residual Process
  efp <- t(B.inv.sq %*% t(W))
  
  # Compute the maximum statistics
  U.n <- apply(efp, 1, function(x) crossprod(x,x))
  stat <- max(abs(U.n))
  
  # Compute the p-value
  pvalue = mean(BB.stat>=stat)
  
  return(list(pvalue, efp))
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 9. Application functions ####

combinations <- function(covariates){
  response_variable <- "y"
  all_combinations <- unlist(lapply(1:length(covariates),
                                    function(r) combn(covariates, r,
                                                      simplify = FALSE)),
                             recursive = FALSE)
  regression_formulas <- lapply(all_combinations, function(combo) {
    if (length(which(grepl("time", combo)))==0) {
      formula <- paste(response_variable, "~", paste(combo, collapse = " + "), "-1")
    } else {
      if (sum(sapply(which(grepl("time", combo)), function(x)
        any(grepl(tolower(substr(combo[x], 5, nchar(combo[x]))), combo)))) == 0){
        formula <- paste(response_variable, "~", paste(combo, collapse = " + "), "-1")
        return(formula)
      } else {
        return(NA)
      }
    }
  })
  regression_formulas <- regression_formulas[!sapply(regression_formulas, is.na)]
  regression_formulas.with.time <- regression_formulas[sapply(regression_formulas,
                                   function(combo) any(grepl("time", combo)))]
  regression_formulas.without.time <- regression_formulas[!sapply(regression_formulas,
                                      function(combo) any(grepl("time", combo)))]
  return(list(regression_formulas.without.time, regression_formulas.with.time))
}

gam.fitting <- function(regression_formulas, dat.gam){
  gam.fits <- list()
  aics <- rep(0,length(regression_formulas))
  for (i in 1:length(regression_formulas)){
    print(paste0(i,"/",length(regression_formulas)))
    system.time(gam.fits[[i]] <-
                  gam(as.formula(regression_formulas[[i]]),
                      family="binomial"(link = 'logit'), 
                      method="REML",
                      data=dat.gam))
    aics[i] <- AIC(gam.fits[[i]])
  }
  return(list(gam.fits, aics))
}

#-------------------------------------------------------------------------------

# 10. Function data.simulation_example (for example in theoretical section)
# time-based reciprocity is driving the actual data generating process
data.simulation_example <- function(n, p, beta0, beta1, m.ncc, seed, 
                                    close.rec=TRUE){
  
  set.seed(seed)
  dat.gam <- NULL
  
  last_time <- first_time <- matrix(-Inf, ncol = p,nrow = p)
  
  last.rec <- function(t){exp(-(t-t(last_time)))}
  first.rec <- function(t){exp(-(t-t(first_time)))}
  rec <- function(t){t>t(last_time) & t(last_time)!=-Inf}
  
  t <- 0
  rate <- matrix(beta0, ncol = p,nrow = p)
  rate.upd <- beta1 * last.rec(t)
  tot.rate <- exp(rate + rate.upd)
  diag(tot.rate) <- 0
  
  for (i in 1:n){
    
    dt <- rexp(1,sum(tot.rate))
    t <- t + dt
    
    prob <- tot.rate/sum(tot.rate)
    sr.ev <- sample(p^2,1,prob = prob)
    
    s.ev<-(sr.ev-1)%%p+1
    r.ev<-(sr.ev-1)%/%p+1
    
    for (j in 1:(m.ncc-1)){
      sr.nv<-sample(setdiff((1:(p^2)),sr.ev),1)
      s.nv<-(sr.nv-1)%%p+1
      r.nv<-(sr.nv-1)%/%p+1
      dat.gam <- rbind(dat.gam, c(t, s.ev, r.ev, s.nv, r.nv,
                                  last.rec(t)[s.ev, r.ev],
                                  last.rec(t)[s.nv, r.nv],
                                  last.rec(t)[s.ev, r.ev] - 
                                    last.rec(t)[s.nv, r.nv],
                                  first.rec(t)[s.ev, r.ev],
                                  first.rec(t)[s.nv, r.nv],
                                  first.rec(t)[s.ev, r.ev] - 
                                    first.rec(t)[s.nv, r.nv], 
                                  rec(t)[s.ev, r.ev],
                                  rec(t)[s.nv, r.nv],
                                  rec(t)[s.ev, r.ev] - rec(t)[s.nv, r.nv], 
                                  j))
    }
    
    first_time[s.ev,r.ev] <- ifelse(first_time[s.ev,r.ev]==-Inf, t, 
                                    first_time[s.ev,r.ev])
    last_time[s.ev,r.ev] <- t
    
    rate.upd <- beta1 * last.rec(t)
    tot.rate <- exp(rate + rate.upd)
    diag(tot.rate) <- 0
    
    if(close.rec==TRUE){
      first_time[r.ev,s.ev] <- ifelse(rec(t)[s.ev,r.ev], -Inf, 
                                      first_time[r.ev,s.ev])
      last_time[r.ev,s.ev] <- ifelse(rec(t)[s.ev,r.ev], -Inf, 
                                     last_time[r.ev,s.ev])
    }
  }
  
  dat.gam <- data.frame(y = rep(1, nrow(dat.gam)), 
                        stp = dat.gam[,1],
                        s1 = dat.gam[,2],
                        r1 = dat.gam[,3],
                        s2 = dat.gam[,4],
                        r2 = dat.gam[,5],
                        last.rec1 = dat.gam[,6],
                        last.rec2 = dat.gam[,7],
                        last.rec = dat.gam[,8],
                        first.rec1 = dat.gam[,9],
                        first.rec2 = dat.gam[,10],
                        first.rec = dat.gam[,11],
                        rec1 = dat.gam[,12],
                        rec2 = dat.gam[,13],
                        rec = dat.gam[,14],
                        m.ncc = dat.gam[,15]) 
  return(dat.gam)
}
