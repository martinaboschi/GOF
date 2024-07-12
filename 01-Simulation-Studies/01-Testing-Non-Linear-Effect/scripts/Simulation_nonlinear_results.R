#-------------------------------------------------------------------------------
# Script: Simulation_nonlinear_results.R
# Author: Martina Boschi
# Date: July 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
#-------------------------------------------------------------------------------

# 2.    Coverage of the Statistical Test

## 2.1.   Non-Linear Effect of Reciprocity

# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n1000.RData")
# pvalueCS.1000 <- pvalue.CS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n5000.RData")
# pvalueCS.5000 <- pvalue.CS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n10000.RData")
# pvalueCS.10000 <- pvalue.CS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n50000.RData")
# pvalueCS.50000 <- pvalue.CS
# save(pvalueCS.1000, pvalueCS.5000, pvalueCS.10000, pvalueCS.50000,
#      file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_resultsCS.RData")

load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_resultsCS.RData")

pdf("01-Simulation-Studies/01-Testing-Non-Linear-Effect/pictures/nonlinear_coverage.pdf", width=15, height = 15)
plot(ecdf(pvalueCS.1000), cex=0.2, col=pal.blue[4], lty=1,
     main="Testing a 9-dimensional TPRS - Increasing n",
     xlab="pvalue", 
     ylab="ecdf",
     cex.sub=2,
     cex.main=2.5, 
     cex.lab=2.2)  
lines(ecdf(pvalueCS.5000), cex=0.2, col=pal.blue[8], lty=1)
lines(ecdf(pvalueCS.10000), cex=0.2, col=pal.yellow[4], lty=1)
lines(ecdf(pvalueCS.50000), cex=0.2, col=pal.yellow[8], lty=1)
abline(0,1,lwd=2, col=1)

legend("bottomright", 
       legend=c("Uniform Distribution", 
                paste("n = 1000, prop_rej:", as.character(mean(pvalueCS.1000<0.05))), 
                paste("n = 5000, prop_rej:", as.character(mean(pvalueCS.5000<0.05))),
                paste("n = 10000, prop_rej:", as.character(mean(pvalueCS.10000<0.05))), 
                paste("n = 50000, prop_rej:", as.character(mean(pvalueCS.50000<0.05)))),
       col=c(1,
             pal.blue[4], 
             pal.blue[8],
             pal.yellow[4], 
             pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
       cex=2)
dev.off()

rm(pvalueCS.1000, pvalueCS.5000, pvalueCS.10000, pvalueCS.50000)

# 2.    Power of the Statistical Test

## 2.1.   True Non-Linear Effect of Reciprocity fitted with a fixed effect

# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n1000.RData")
# pvalueMS.1000 <- pvalue.MS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n5000.RData")
# pvalueMS.5000 <- pvalue.MS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n10000.RData")
# pvalueMS.10000 <- pvalue.MS
# load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_n50000.RData")
# pvalueMS.50000 <- pvalue.MS
# save(pvalueMS.1000, pvalueMS.5000, pvalueMS.10000, pvalueMS.50000,
#      file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_resultsMS.RData")

load(file="01-Simulation-Studies/01-Testing-Non-Linear-Effect/output/nonlinear_resultsMS.RData")

pdf("01-Simulation-Studies/01-Testing-Non-Linear-Effect/pictures/nonlinear_power.pdf", width=15, height = 15)
plot(ecdf(pvalueMS.1000), cex=0.2, col=pal.blue[4], lty=1,
     main="Testing FLE for Reciprocity - Increasing n", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.sub=2,
     cex.main=2.5, 
     cex.lab=2)  
lines(ecdf(pvalueMS.5000), cex=0.2, col=pal.blue[8], lty=1)
lines(ecdf(pvalueMS.10000), cex=0.2, col=pal.yellow[4], lty=1, lwd=4)
lines(ecdf(pvalueMS.50000), cex=0.2, col=pal.yellow[8], lty=1)
abline(0,1,lwd=2, col=1)

legend("bottomright", 
       legend=c("Uniform Distribution", 
                paste("n = 1000, prop_rej:", as.character(mean(pvalueMS.1000<0.05))), 
                paste("n = 5000, prop_rej:", as.character(mean(pvalueMS.5000<0.05))),
                paste("n = 10000, prop_rej:", as.character(mean(pvalueMS.10000<0.05))), 
                paste("n = 50000, prop_rej:", as.character(mean(pvalueMS.50000<0.05)))),
       col=c(1,
             pal.blue[4], 
             pal.blue[8],
             pal.yellow[4], 
             pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
       cex=2)
dev.off()

rm(pvalueMS.1000, pvalueMS.5000, pvalueMS.10000, pvalueMS.50000)

