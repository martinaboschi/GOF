# GOF
Goodness of fit of relational event models

This repository contains the codes for the paper **Goodness of fit of relational event models**. 

## Abstract

A type of dynamic networks involves temporally ordered interactions between actors, where past network configurations may influence future ones. The relational event model can be used to identify the underlying dynamics that drive interactions among system components. Despite the rapid development of this model over the past 15 years, an ongoing area of research revolves around evaluating the goodness of fit of this model, especially when it incorporates time-varying and random effects. Current methodologies often rely on comparing observed and simulated events using specific statistics, but this can be computationally intensive, and requires various assumptions.
		
We propose an additive mixed-effect relational event model estimated via case-control sampling, and introduce a versatile framework for testing the goodness of fit of such models using weighted martingale residuals. Our focus is on a Kolmogorov-Smirnov type test designed to assess if covariates are accurately modeled. Our approach can be easily extended to evaluate whether other features of network dynamics have been appropriately incorporated into the model. We assess the goodness of fit of various relational event models using synthetic data to evaluate the test's power and coverage. Furthermore, we apply the method to a social study involving 57,791 emails sent by 159 employees of a Polish manufacturing company in 2010.
		
The method is implemented in the R package `mgcv`.

## Structure of the Repository

First, we report the list of the necessary `R` packages to be installed: 
`mgcv`, `mgcViz`, `RColorBrewer`, `gridExtra`, `sde`, `Matrix`, `expm`, `gratia`, `e1071`, `lubridate`.

The reader will find four folders: 

A.	`00-Functions`: contains all the functions used in the simulation and application study. 
    **IMPORTANT**: `R` packages are loaded in the file `functions.R`. The user can directly call them using `source(00-Functions/functions.R)`.

B.	`01-Simulation-Studies`: 3 main simulations studies have been conducted in order to produce the relative section in the paper. 

    -	`01-Testing-Non-Linear-Effect`: this simulation study is aimed at verifying coverage and power of the test involving **covariates that may be represented by univariate or multivariate components of the model matrix**.     
    -	`02-Testing-Random-Effect`: this simulation study is aimed at verifying coverage and power of the test involving **components of the model effects referring to random effects fitted as splines**. 
    -	`03-Omnibus-Testing`: this simulation study is aimed at verifying coverage and power of the test whenever inspecting **global adequacy of the model**. 
    
C.	`02-Application`: testing goodness of fit of several models fitted to data involving 57,791 emails sent by 159 employees of a Polish manufacturing company in 2010.

D.	`03-Supplementary-Materials`: this folder contains the code which is necessary to generate the plots in the Supplementary Materials. 