#-------------------------------------------------------------------------------
# Script: Simulation_global_results.R
# Author: Martina Boschi
# Date: July 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
#-------------------------------------------------------------------------------

# 2.    Coverage of the Statistical Test

# load("01-Simulation-Studies/03-Omnibus-Testing/output/global_one_element.RData")
# pvalueCS.L1 <- pvalue.CS
# load("01-Simulation-Studies/03-Omnibus-Testing/output/global_two_elements.RData")
# pvalueCS.L2 <- pvalue.CS
# load("01-Simulation-Studies/03-Omnibus-Testing/output/global_three_elements.RData")
# pvalueCS.L3 <- pvalue.CS
# load("01-Simulation-Studies/03-Omnibus-Testing/output/global_four_elements.RData")
# pvalueCS.L4 <- pvalue.CS
# save(pvalueCS.L1, pvalueCS.L2, pvalueCS.L3, pvalueCS.L4,
#      file="01-Simulation-Studies/03-Omnibus-Testing/output/global_resultsCS.RData")

load(file="01-Simulation-Studies/03-Omnibus-Testing/output/global_resultsCS.RData")

pdf("01-Simulation-Studies/03-Omnibus-Testing/pictures/global_coverage.pdf", width=15, height = 15)
plot(ecdf(pvalueCS.L1), cex=0.2, col=pal.blue[4], lty=1,
     main="Global Test of Correctly Specified Models - Increasing N. Components L",
     xlab="pvalue", 
     ylab="ecdf",
     cex.sub=2,
     cex.main=2, 
     cex.lab=2.2)  
lines(ecdf(pvalueCS.L2), cex=0.2, col=pal.blue[8], lty=1)
lines(ecdf(pvalueCS.L3), cex=0.2, col=pal.yellow[4], lty=1)
lines(ecdf(pvalueCS.L4), cex=0.2, col=pal.yellow[8], lty=1)
abline(0,1,lwd=2, col=1)

legend("bottomright", 
       legend=c("Uniform Distribution", 
                paste("L = 1, prop_rej:", as.character(mean(pvalueCS.L1<0.05))), 
                paste("L = 2, prop_rej:", as.character(mean(pvalueCS.L2<0.05))),
                paste("L = 3, prop_rej:", as.character(mean(pvalueCS.L3<0.05))), 
                paste("L = 4, prop_rej:", as.character(mean(pvalueCS.L4<0.05)))),
       col=c(1,
             pal.blue[4], 
             pal.blue[8],
             pal.yellow[4], 
             pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
       cex=2)
dev.off()

rm(pvalueCS.L1, pvalueCS.L2, pvalueCS.L3, pvalueCS.L4)

# 2. Power of the Statistical Test

load("01-Simulation-Studies/03-Omnibus-Testing/output/global_two_elements.RData")
pdf("01-Simulation-Studies/03-Omnibus-Testing/pictures/global_power.pdf", width=15, height = 15)
plot(ecdf(pvalue.CS), cex=0.2, col=pal.blue[8], 
     main="Global Test of Model including 2 covariates",
     xlab="pvalue", 
     ylab="ecdf",
     cex.sub=2,
     cex.main=2.3, 
     cex.lab=2.2) 
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS1), cex=0.2, col=pal.yellow[4])
lines(ecdf(pvalue.MS2), cex=0.2, col=pal.blue[4])
lines(ecdf(pvalue.MS3), cex=0.2, col=pal.yellow[8])
legend("bottomright", cex = 1.5,
       legend=c("Uniform Distribution", 
                paste("CS, prop_rej:", as.character(mean(pvalue.CS<0.05))),
                paste("MS - Mispecified Reciprocity, prop_rej:", 
                      as.character(mean(pvalue.MS1<0.05))),
                paste("MS - Mispecified Exponential Covariate, prop_rej:", 
                      as.character(mean(pvalue.MS2<0.05))),
                paste("MS - All elements misspecified, prop_rej:", 
                      as.character(mean(pvalue.MS3<0.05)))
                ),
       col=c(1, pal.blue[8], pal.yellow[4], pal.blue[4], pal.yellow[8]), 
       lwd=c(2,1,1,1,1), 
)
dev.off()

