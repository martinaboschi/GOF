#-------------------------------------------------------------------------------
# Script: Simulation_Supplementary.R
# Author: Martina Boschi
# Date: July 2024
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
# 2. Result visualization ####

#-------------------------------------------------------------------------------
## 2.1. One-element-included ####
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_one_element.RData")
pdf("03-Supplementary-Materials/pictures/global_power_one_element.pdf", width=15, height = 15)
plot(ecdf(pvalue.CS), cex=0.2, col=pal.blue[8], 
     main="Global Test of Model including 1 covariate", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=2,
     cex.lab=2) 
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS), cex=0.2, col=pal.yellow[4])
legend("bottomright", cex = 1.5,
       legend=c("Uniform Distribution", 
                paste("CS, prop_rej:", as.character(mean(pvalue.CS<0.05))),
                paste("MS - Mispecified Reciprocity, prop_rej:", 
                      as.character(mean(pvalue.MS<0.05)))
       ),
       col=c(1, pal.blue[8], pal.yellow[4]), 
       lwd=c(2,1,1), 
)
dev.off()
rm(list = ls())
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
## 2.2. Two-elements-included ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_two_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_two_elements.pdf", width=15, height = 15)
plot(ecdf(pvalue.CS), cex=0.2, col=pal.blue[8], 
     main="Global Test of Model including 2 covariates",
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=2,
     cex.lab=2) 
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS1), cex=0.2, col=pal.yellow[4])
lines(ecdf(pvalue.MS2), cex=0.2, col=pal.blue[4])
lines(ecdf(pvalue.MS3), cex=0.2, col=pal.yellow[8])
legend("bottomright", cex = 1,
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
rm(list = ls())
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
## 2.3. Three-elements-included ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_three_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_three_elements.pdf", width=15, height = 15)
plot(ecdf(pvalue.CS), cex=0.2, col=pal.blue[8], 
     main="Global Test of Model including 3 covariates", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=2,
     cex.lab=2)
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS1), cex=0.2, col=pal.yellow[4])
lines(ecdf(pvalue.MS2), cex=0.2, col=pal.blue[4])
lines(ecdf(pvalue.MS3), cex=1, col=pal.yellow[8])
lines(ecdf(pvalue.MS4), cex=0.2, col=pal.yellow[6])
lines(ecdf(pvalue.MS5), cex=0.7, col=pal.blue[6])
lines(ecdf(pvalue.MS6), cex=0.5, col=pal.greys[6])
lines(ecdf(pvalue.MS7), cex=0.3, col=pal.greys[8])
legend("bottomright", cex = 1,
       legend=c("Uniform Distribution", 
                paste("CS, prop_rej:", as.character(mean(pvalue.CS<0.05))),
                paste("MS - Mispecified Reciprocity, prop_rej:", 
                      as.character(mean(pvalue.MS1<0.05))),
                paste("MS - Mispecified Exponential Covariate, prop_rej:", 
                      as.character(mean(pvalue.MS2<0.05))),
                paste("MS - Mispecified Non-Linear Effect, prop_rej:", 
                      as.character(mean(pvalue.MS3<0.05))),
                paste("MS - Mispecified Reciprocity & Exp. Cov, prop_rej:", 
                      as.character(mean(pvalue.MS4<0.05))),
                paste("MS - Mispecified Reciprocity & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS5<0.05))),
                paste("MS - Mispecified Exp. Cov & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS6<0.05))),
                paste("MS - All elements misspecified, prop_rej:", 
                      as.character(mean(pvalue.MS7<0.05)))
       ),
       col=c(1, pal.blue[8], pal.yellow[4], pal.blue[4], pal.yellow[8], 
             pal.yellow[6], pal.blue[6], pal.greys[6], pal.greys[8]),
       lwd=c(1,1,1,1,2,1,1.8,1.5,1.3), 
)
dev.off()
rm(list = ls())
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_four_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_four_elements.pdf", width=30, height = 15)
par(mfrow=c(1,2))
plot(ecdf(pvalue.CS), cex=0.2, col=pal.blue[8], 
     main="Global Test of Model including 3 covariates and correctly specified random effects", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=2,
     cex.lab=2)
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS1), cex=0.2, col=pal.yellow[4])
lines(ecdf(pvalue.MS2), cex=0.2, col=pal.blue[4])
lines(ecdf(pvalue.MS3), cex=0.2, col=pal.yellow[8])
lines(ecdf(pvalue.MS4), cex=0.2, col=pal.yellow[6])
lines(ecdf(pvalue.MS5), cex=0.5, col=pal.blue[6])
lines(ecdf(pvalue.MS6), cex=0.2, col=pal.greys[6])
lines(ecdf(pvalue.MS7), cex=0.2, col=pal.greys[8])
legend("bottomright", cex = 1.6,
       legend=c("Uniform Distribution", 
                paste("CS, prop_rej:", as.character(mean(pvalue.CS<0.05))),
                paste("MS - Mispecified Reciprocity, prop_rej:", 
                      as.character(mean(pvalue.MS1<0.05))),
                paste("MS - Mispecified Exponential Covariate, prop_rej:", 
                      as.character(mean(pvalue.MS2<0.05))),
                paste("MS - Mispecified Non-Linear Effect, prop_rej:", 
                      as.character(mean(pvalue.MS3<0.05))),
                paste("MS - Mispecified Reciprocity & Exp. Cov, prop_rej:", 
                      as.character(mean(pvalue.MS4<0.05))),
                paste("MS - Mispecified Reciprocity & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS5<0.05))),
                paste("MS - Mispecified Exp. Cov & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS6<0.05))),
                paste("MS - All elements misspecified, prop_rej:", 
                      as.character(mean(pvalue.MS7<0.05)))
       ),
       col=c(1, pal.blue[8], pal.yellow[4], pal.blue[4], pal.yellow[8], 
             pal.yellow[6], pal.blue[6], pal.greys[6], pal.greys[8]),
       lwd=c(1,1,1,1,1,1,1.5,1,1), 
)
plot(ecdf(pvalue.MS8), cex=0.2, col=pal.yellow[4], 
     main="Global Test of Model including 3 covariates and misspecified random effects", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=2,
     cex.lab=2)
abline(0,1, lwd=2)
lines(ecdf(pvalue.MS9), cex=0.2, col=pal.blue[4])
lines(ecdf(pvalue.MS10), cex=0.2, col=pal.yellow[8])
lines(ecdf(pvalue.MS11), cex=0.2, col=pal.yellow[6])
lines(ecdf(pvalue.MS12), cex=0.5, col=pal.blue[6])
lines(ecdf(pvalue.MS13), cex=0.2, col=pal.greys[6])
lines(ecdf(pvalue.MS14), cex=0.2, col=pal.greys[8])
legend("bottomright", cex = 1.6,
       legend=c("Uniform Distribution",
                paste("MS - Mispecified Reciprocity, prop_rej:", 
                      as.character(mean(pvalue.MS8<0.05))),
                paste("MS - Mispecified Exponential Covariate, prop_rej:", 
                      as.character(mean(pvalue.MS9<0.05))),
                paste("MS - Mispecified Non-Linear Effect, prop_rej:", 
                      as.character(mean(pvalue.MS10<0.05))),
                paste("MS - Mispecified Reciprocity & Exp. Cov, prop_rej:", 
                      as.character(mean(pvalue.MS11<0.05))),
                paste("MS - Mispecified Reciprocity & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS12<0.05))),
                paste("MS - Mispecified Exp. Cov & NLE, prop_rej:", 
                      as.character(mean(pvalue.MS13<0.05))),
                paste("MS - All elements misspecified, prop_rej:", 
                      as.character(mean(pvalue.MS14<0.05)))
       ),
       col=c(1, pal.yellow[4], pal.blue[4], pal.yellow[8], 
             pal.yellow[6], pal.blue[6], pal.greys[6], pal.greys[8]),
       lwd=c(1,1,1,1,1,1.5,1,1), 
)
dev.off()
rm(list = ls())
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Single-plot evaluation

#-------------------------------------------------------------------------------
## 3.1. Two-elements-included ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_two_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_two_elements_single.pdf", width=20, height = 30)

par(mfrow=c(4,2))
plot(ecdf(pvalue.CS_comp1), cex=0.5, col=pal.yellow[4], 
     main="CS: Testing (correctly specified) reciprocity", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.CS_comp2), cex=0.5, col=pal.blue[8], 
     main="CS: Testing (correctly specified) x1", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS1_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS1: Testing (misspecified) reciprocity", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS1_comp2), cex=0.5, col=pal.blue[8], 
     main="MS1: Testing (correctly specified) x1", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS2_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS2: Testing (correctly specified) reciprocity", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS2_comp2), cex=0.5, col=pal.blue[8], 
     main="MS2: Testing (misspecified) x1", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS3_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS3: Testing (misspecified) reciprocity", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
plot(ecdf(pvalue.MS3_comp2), cex=0.5, col=pal.blue[8], 
     main="MS3: Testing (misspecified) x1", 
     sub="Empirical Distribution of the p-value", 
     xlab="pvalue", 
     ylab="ecdf")
abline(0,1)
dev.off()

#-------------------------------------------------------------------------------
## 3.3. Three-elements-included ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_three_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_three_elements_single.pdf", width=70, height = 50)

par(mfrow=c(4,6))
plot(ecdf(pvalue.CS_comp1), cex=0.5, col=pal.yellow[4], 
     main="CS: Testing (correctly specified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.CS_comp2), cex=0.5, col=pal.blue[8], 
     main="CS: Testing (correctly specified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.CS_comp3), cex=0.5, col=pal.greys[4], 
     main="CS: Testing (correctly specified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS1: Testing (misspecified) reciprocity", 
     xlab="pvalue",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp2), cex=0.5, col=pal.blue[8], 
     main="MS1: Testing (correctly specified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp3), cex=0.5, col=pal.greys[4], 
     main="MS1: Testing (correctly specified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS2: Testing (correctly specified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp2), cex=0.5, col=pal.blue[8], 
     main="MS2: Testing (misspecified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp3), cex=0.5, col=pal.greys[4], 
     main="MS2: Testing (correctly specified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS3: Testing (correctly specified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp2), cex=0.5, col=pal.blue[8], 
     main="MS3: Testing (correctly specified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp3), cex=0.5, col=pal.greys[4], 
     main="MS3: Testing (misspecified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS4: Testing (misspecified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp2), cex=0.5, col=pal.blue[8], 
     main="MS4: Testing (misspecified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp3), cex=0.5, col=pal.greys[4], 
     main="MS4: Testing (correctly specified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS5: Testing (misspecified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp2), cex=0.5, col=pal.blue[8], 
     main="MS5: Testing (correctly specified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp3), cex=0.5, col=pal.greys[4], 
     main="MS5: Testing (misspecified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS6: Testing (correctly specified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp2), cex=0.5, col=pal.blue[8], 
     main="MS6: Testing (misspecified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp3), cex=0.5, col=pal.greys[4], 
     main="MS6: Testing (misspecified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS7: Testing (misspecified) reciprocity", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp2), cex=0.5, col=pal.blue[8], 
     main="MS7: Testing (misspecified) x1", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp3), cex=0.5, col=pal.greys[4], 
     main="MS7: Testing (misspecified) x2", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
dev.off()

#-------------------------------------------------------------------------------
## 3.3. Three-elements-included ####
pal.blue <- brewer.pal(9, "Blues")
pal.yellow <- brewer.pal(9, "YlOrBr")
pal.greys <- brewer.pal(9, "Greys")
load("01-Simulation-Studies/03-Omnibus-Testing/output/global_four_elements.RData")
pdf("03-Supplementary-Materials/pictures/global_power_four_elements_single_A.pdf", width=100, height = 50)

par(mfrow=c(4,8))
plot(ecdf(pvalue.CS_comp1), cex=0.5, col=pal.yellow[4], 
     main="CS: Testing (correctly specified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.CS_comp2), cex=0.5, col=pal.blue[8], 
     main="CS: Testing (correctly specified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.CS_comp3), cex=0.5, col=pal.greys[4], 
     main="CS: Testing (correctly specified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.CS_comp4), cex=0.5, col=pal.greys[4], 
     main="CS: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS1: Testing (misspecified) reciprocity - RE CS", 
     xlab="pvalue",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp2), cex=0.5, col=pal.blue[8], 
     main="MS1: Testing (correctly specified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp3), cex=0.5, col=pal.greys[4], 
     main="MS1: Testing (correctly specified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS1_comp4), cex=0.5, col=pal.greys[4], 
     main="MS1: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS2: Testing (correctly specified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp2), cex=0.5, col=pal.blue[8], 
     main="MS2: Testing (misspecified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp3), cex=0.5, col=pal.greys[4], 
     main="MS2: Testing (correctly specified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS2_comp4), cex=0.5, col=pal.greys[4], 
     main="MS2: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS3: Testing (correctly specified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp2), cex=0.5, col=pal.blue[8], 
     main="MS3: Testing (correctly specified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp3), cex=0.5, col=pal.greys[4], 
     main="MS3: Testing (misspecified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS3_comp4), cex=0.5, col=pal.greys[4], 
     main="MS3: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS4: Testing (misspecified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp2), cex=0.5, col=pal.blue[8], 
     main="MS4: Testing (misspecified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp3), cex=0.5, col=pal.greys[4], 
     main="MS4: Testing (correctly specified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS4_comp4), cex=0.5, col=pal.greys[4], 
     main="MS4: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS5: Testing (misspecified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp2), cex=0.5, col=pal.blue[8], 
     main="MS5: Testing (correctly specified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp3), cex=0.5, col=pal.greys[4], 
     main="MS5: Testing (misspecified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS5_comp4), cex=0.5, col=pal.greys[4], 
     main="MS5: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS6: Testing (correctly specified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp2), cex=0.5, col=pal.blue[8], 
     main="MS6: Testing (misspecified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp3), cex=0.5, col=pal.greys[4], 
     main="MS6: Testing (misspecified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS6_comp4), cex=0.5, col=pal.greys[4], 
     main="MS6: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS7: Testing (misspecified) reciprocity - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp2), cex=0.5, col=pal.blue[8], 
     main="MS7: Testing (misspecified) x1 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp3), cex=0.5, col=pal.greys[4], 
     main="MS7: Testing (misspecified) x2 - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS7_comp4), cex=0.5, col=pal.greys[4], 
     main="MS7: Testing (correctly specified) RE - RE CS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
dev.off()

pdf("03-Supplementary-Materials/pictures/global_power_four_elements_single_B.pdf", width=100, height = 50)
par(mfrow=c(4,8))
plot(ecdf(pvalue.MS8_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS8: Testing (misspecified) reciprocity - RE MS", 
     xlab="pvalue",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS8_comp2), cex=0.5, col=pal.blue[8], 
     main="MS8: Testing (correctly specified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS8_comp3), cex=0.5, col=pal.greys[4], 
     main="MS8: Testing (correctly specified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS8_comp4), cex=0.5, col=pal.greys[4], 
     main="MS8: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS9_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS9: Testing (correctly specified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS9_comp2), cex=0.5, col=pal.blue[8], 
     main="MS9: Testing (misspecified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS9_comp3), cex=0.5, col=pal.greys[4], 
     main="MS9: Testing (correctly specified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS9_comp4), cex=0.5, col=pal.greys[4], 
     main="MS9: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS10_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS10: Testing (correctly specified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS10_comp2), cex=0.5, col=pal.blue[8], 
     main="MS10: Testing (correctly specified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS10_comp3), cex=0.5, col=pal.greys[4], 
     main="MS10: Testing (misspecified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS10_comp4), cex=0.5, col=pal.greys[4], 
     main="MS10: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS11_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS11: Testing (misspecified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS11_comp2), cex=0.5, col=pal.blue[8], 
     main="MS11: Testing (misspecified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS11_comp3), cex=0.5, col=pal.greys[4], 
     main="MS11: Testing (correctly specified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS11_comp4), cex=0.5, col=pal.greys[4], 
     main="MS11: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS12_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS12: Testing (misspecified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS12_comp2), cex=0.5, col=pal.blue[8], 
     main="MS12: Testing (correctly specified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS12_comp3), cex=0.5, col=pal.greys[4], 
     main="MS12: Testing (misspecified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS12_comp4), cex=0.5, col=pal.greys[4], 
     main="MS12: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS13_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS13: Testing (correctly specified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS13_comp2), cex=0.5, col=pal.blue[8], 
     main="MS13: Testing (misspecified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS13_comp3), cex=0.5, col=pal.greys[4], 
     main="MS13: Testing (misspecified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS13_comp4), cex=0.5, col=pal.greys[4], 
     main="MS13: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS14_comp1), cex=0.5, col=pal.yellow[4], 
     main="MS14: Testing (misspecified) reciprocity - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS14_comp2), cex=0.5, col=pal.blue[8], 
     main="MS14: Testing (misspecified) x1 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS14_comp3), cex=0.5, col=pal.greys[4], 
     main="MS14: Testing (misspecified) x2 - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
plot(ecdf(pvalue.MS14_comp4), cex=0.5, col=pal.greys[4], 
     main="MS14: Testing (misspecified) RE - RE MS", 
     xlab="pvalue", 
     ylab="ecdf",
     cex.main=4)
abline(0,1)
dev.off()
