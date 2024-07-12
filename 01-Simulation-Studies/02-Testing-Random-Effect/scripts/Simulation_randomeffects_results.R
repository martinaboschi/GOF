#-------------------------------------------------------------------------------
# Script: Simulation_randomeffects_results.R
# Author: Martina Boschi
# Date: July 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
#-------------------------------------------------------------------------------

# 2.    Coverage of the Statistical Test

## 2.1.   Sender Activity Random Intercept

# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p10.RData")
# pvalueCS.10 <- pvalue.CS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p50.RData")
# pvalueCS.50 <- pvalue.CS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p100.RData")
# pvalueCS.100 <- pvalue.CS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p150.RData")
# pvalueCS.150 <- pvalue.CS
# save(pvalueCS.10, pvalueCS.50, pvalueCS.100, pvalueCS.150,
#      file="01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_resultsCS.RData")

load(file="01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_resultsCS.RData")

pdf("01-Simulation-Studies/02-Testing-Random-Effect/pictures/randomeffects_coverage.pdf", width=15, height = 15)
plot(ecdf(pvalueCS.10), cex=0.2, col=pal.blue[4], lty=1,
     main="Testing Sender Activity - Increasing p", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.sub=2,
     cex.main=2.5, 
     cex.lab=2.2)  
lines(ecdf(pvalueCS.50), cex=0.2, col=pal.blue[8], lty=1)
lines(ecdf(pvalueCS.100), cex=0.2, col=pal.yellow[4], lty=1)
lines(ecdf(pvalueCS.150), cex=0.2, col=pal.yellow[8], lty=1)
abline(0,1,lwd=2, col=1)

legend("bottomright", 
       legend=c("Uniform Distribution", 
                paste("p = 10, prop_rej:", as.character(mean(pvalueCS.10<0.05))), 
                paste("p = 50, prop_rej:", as.character(mean(pvalueCS.50<0.05))),
                paste("p = 100, prop_rej:", as.character(mean(pvalueCS.100<0.05))), 
                paste("p = 150, prop_rej:", as.character(mean(pvalueCS.150<0.05)))),
       col=c(1,
             pal.blue[4], 
             pal.blue[8],
             pal.yellow[4], 
             pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
       cex=2)
dev.off()

rm(pvalueCS.10, pvalueCS.50, pvalueCS.100, pvalueCS.150)

# 2.    Power of the Statistical Test

## 2.1.   Sender Activity Random Slope for Time

# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p10.RData")
# pvalueMS.10 <- pvalue.MS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p50.RData")
# pvalueMS.50 <- pvalue.MS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p100.RData")
# pvalueMS.100 <- pvalue.MS
# load("01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_p150.RData")
# pvalueMS.150 <- pvalue.MS
# save(pvalueMS.10, pvalueMS.50, pvalueMS.100, pvalueMS.150,
#      file="01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_resultsMS.RData")

load(file="01-Simulation-Studies/02-Testing-Random-Effect/output/randomeffects_resultsMS.RData")

pdf("01-Simulation-Studies/02-Testing-Random-Effect/pictures/randomeffects_power.pdf", width=15, height = 15)
plot(ecdf(pvalueMS.10), cex=0.2, col=pal.blue[4], lty=1,
     main="Testing Sender Random Slope for Time - Increasing p", 
     xlab="pvalue", 
     ylab="ecdf", 
     xlim=c(0,1),
     cex.sub=2,
     cex.main=2.3, 
     cex.lab=2.2)  
lines(ecdf(pvalueMS.50), cex=0.2, col=pal.blue[8], lty=1)
lines(ecdf(pvalueMS.100), cex=0.2, col=pal.yellow[4], lty=1)
lines(ecdf(pvalueMS.150), cex=0.2, col=pal.yellow[8], lty=1)
abline(0,1,lwd=2, col=1)

legend("bottomright", 
       legend=c("Uniform Distribution", 
                paste("p = 10, prop_rej:", as.character(mean(pvalueMS.10<0.05))), 
                paste("p = 50, prop_rej:", as.character(mean(pvalueMS.50<0.05))),
                paste("p = 100, prop_rej:", as.character(mean(pvalueMS.100<0.05))), 
                paste("p = 150, prop_rej:", as.character(mean(pvalueMS.150<0.05)))),
       col=c(1,
             pal.blue[4], 
             pal.blue[8],
             pal.yellow[4], 
             pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
       cex=2)
dev.off()

rm(pvalueMS.10, pvalueMS.50, pvalueMS.100, pvalueMS.150)