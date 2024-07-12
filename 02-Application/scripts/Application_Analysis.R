#-------------------------------------------------------------------------------
# Script: Application_Analysis.R
# Author: Martina Boschi
# Date: May 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
source("00-Functions/covariate.R")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2. Preliminary analysis on the data
email <- read.csv("02-Application/input/manufacturing.csv", sep = ";")
nrow(email)
range(email$EventDate)

# Do not consider auto-email
sel <- which(email[,1]==email[,2])
email <- email[-sel,]

# Consider only unique email
email <- unique(email)

# Consider only one email per time-point
sel <- unique(email$EventDate)
o <- !duplicated(email$EventDate)
email <- email[o,]

# Numerical format of time 
# (# of days from the first email)
dates <- email$EventDate
tms <- time_length(seconds(ymd_hms(email[,3])))
tms <- tms-min(tms)+1
tms <- tms/(24*60*60)

# Emails as relational events (s,r,t)
email <- cbind(tms, email[,c(1,2)])
colnames(email) <- c("stp","s","r")

# Set of unique senders and receivers
senders <- unique(email$s)
receivers <- unique(email$r)
length(senders)
length(receivers)

# Numerical codification of senders and receivers
actors <- unique(c(email$s, email$r))
email$s <- match(email$s, actors)
email$r <- match(email$r, actors)
p <- length(actors)

# Potential couples of senders and receivers (risk set)
possible.sr <- matrix(1, nrow = p, ncol = p)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3. Case-Control Dataset
dat.gam <- NULL
set.seed(1234)

for (i in 2:nrow(email)){

  if(i%%1000==0){print(i)}

  # Fixed response equal to 1
  y <- 1

  # Event information
  stp <- email[i,1]
  s.ev <- email[i,2]
  r.ev <- email[i,3]
  sr.ev <- (r.ev - 1) * p + s.ev

  # Non-event information
  # sender and receiver are sampled from risk set with equal probility
  sr.nv <- sample(setdiff(which(possible.sr==1), sr.ev), 1)
  s.nv <- (sr.nv-1)%%p+1
  r.nv <- (sr.nv-1)%/%p+1

  # Covariate difference between event and non-event at the time of the event
  # (covariates for specific event and non event are reported as well);
  covariates.ev <- covariate(email, s.ev, r.ev, stp, i)
  covariates.nv <- covariate(email, s.nv, r.nv, stp, i)

  dat.gam <- rbind(dat.gam,
                   c(1, stp,
                     s.ev, r.ev,
                     s.nv, r.nv,
                     covariates.ev[1], covariates.nv[1],
                     covariates.ev[1]- covariates.nv[1],
                     covariates.ev[2], covariates.nv[2],
                     covariates.ev[2]- covariates.nv[2],
                     covariates.ev[3], covariates.nv[3],
                     covariates.ev[3]- covariates.nv[3],
                     covariates.ev[4], covariates.nv[4],
                     covariates.ev[4]- covariates.nv[4],
                     covariates.ev[5], covariates.nv[5],
                     covariates.ev[5]- covariates.nv[5],
                     covariates.ev[6], covariates.nv[6],
                     covariates.ev[6]- covariates.nv[6],
                     covariates.ev[7], covariates.nv[7],
                     covariates.ev[7]- covariates.nv[7],
                     covariates.ev[8], covariates.nv[8],
                     covariates.ev[8]- covariates.nv[8]))
}

dat.gam.complete <- data.frame(
  # Fixed response equal to 1
  y = dat.gam[,1],
  # Event information
  stp = dat.gam[,2],
  s1 = dat.gam[,3],
  r1 = dat.gam[,4],
  # Non-Event information
  s2 = dat.gam[,5],
  r2 = dat.gam[,6],
  # Covariate difference between event and non-event at the time of the event
  # Reciprocity Indicator
  rec1 = dat.gam[,7],
  rec2 = dat.gam[,8],
  rec = dat.gam[,9],
  # Reciprocity Time
  timeRec1 = dat.gam[,10],
  timeRec2 = dat.gam[,11],
  timeRec = dat.gam[,12],
  # Cyclic Closure Indicator
  cyc1 = dat.gam[,13],
  cyc2 = dat.gam[,14],
  cyc = dat.gam[,15],
  # Cyclic Closure Time
  timeCyc1 = dat.gam[,16],
  timeCyc2 = dat.gam[,17],
  timeCyc = dat.gam[,18],
  # Repetition Indicator
  iner1 = dat.gam[,19],
  iner2 = dat.gam[,20],
  iner = dat.gam[,21],
  # Repetition Time
  timeIner1 = dat.gam[,22],
  timeIner2 = dat.gam[,23],
  timeIner = dat.gam[,24],
  # Transitive Closure Indicator
  triad1 = dat.gam[,25],
  triad2 = dat.gam[,26],
  triad = dat.gam[,27],
  # Transitive Closure Time
  timeTriad1 = dat.gam[,28],
  timeTriad2 = dat.gam[,29],
  timeTriad = dat.gam[,30]
)
dat.gam <- dat.gam.complete
# save(dat.gam, file="02-Application/output/dat_gam.RData")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 4. Model selection
# 4.1. Building regression formula 
# Including for each dynamic either the time-based function 
# or the indicator function
covariates_original <- c("rec", "timeRec", "triad", "timeTriad",
                         "iner", "timeIner", "cyc", "timeCyc")

regression_formulas.without.time <- combinations(covariates_original)[[1]]
regression_formulas.with.time <- combinations(covariates_original)[[2]]

# 4.2. Model selection based on the lack of time-related information
gam_fits_aics.without.time <- gam.fitting(regression_formulas.without.time,
                                          dat.gam)
gam.fits.without.time <- gam_fits_aics.without.time[[1]]
aics.without.time <- gam_fits_aics.without.time[[2]]
gam.best.without.time <- gam.fits.without.time[[which.min(aics.without.time)]]

# 4.3. Model selection based on the presence of time-related information
gam_fits_aics.with.time <- gam.fitting(regression_formulas.with.time, dat.gam)
gam.fits.with.time <- gam_fits_aics.with.time[[1]]
aics.with.time <- gam_fits_aics.with.time[[2]]
gam.best.with.time <- gam.fits.with.time[[which.min(aics.with.time)]]
names(aics.with.time) <- regression_formulas.with.time
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 5. Reproducing the asymptotic behaviour of the MR process
BB.stat_9 <- BB.single(dim.k=9)
BB.stat_159 <- BB.single(dim.k=159, n = 2000, n.sim = 10000)

# save(BB.stat_9, BB.stat_159,
#      file="02-Application/output/bb_bridges.RData")
#-------------------------------------------------------------------------------

source("00-Functions/functions.R")
load(file="02-Application/output/dat_gam.RData")
load(file="02-Application/output/fitted_gams_application.RData")
load(file="02-Application/output/bb_bridges.RData")
load(file="02-Application/output/models_and_GOFs.RData")

#-------------------------------------------------------------------------------
# 6. Evaluation of the GOF

# 6.1. Without time - FLE only
GOF_rec.without.time.linear <- GOF_univariate(data = 
                        dat.gam,
                        gam.fit = gam.best.without.time,
                        index = 1)
(GOF_rec.without.time.linear[[1]])
GOF_trs.without.time.linear <- GOF_univariate(data = 
                        dat.gam,
                        gam.fit = gam.best.without.time,
                        index = 2)
(GOF_trs.without.time.linear[[1]])
GOF_rep.without.time.linear <- GOF_univariate(data = 
                        dat.gam,
                        gam.fit = gam.best.without.time,
                        index = 3)
(GOF_rep.without.time.linear[[1]])
GOF_cyc.without.time.linear <- GOF_univariate(data = 
                        dat.gam,
                        gam.fit = gam.best.without.time,
                        index = 4)
(GOF_cyc.without.time.linear[[1]])

# Global

T_g.without.time.linear <- 
  (tan(pi*(0.5-GOF_rec.without.time.linear[[1]])) + 
   tan(pi*(0.5-GOF_trs.without.time.linear[[1]])) + 
   tan(pi*(0.5-GOF_rep.without.time.linear[[1]])) +
   tan(pi*(0.5-GOF_cyc.without.time.linear[[1]])))/4
GOF_global.without.time.linear <- 1/2 - atan(T_g.without.time.linear)/pi

# 6.2. With time - FLE only

GOF_rec.with.time.linear <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.best.with.time,
                            index = 1)
(GOF_rec.with.time.linear[[1]])
GOF_trs.with.time.linear <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.best.with.time,
                            index = 2)
(GOF_trs.with.time.linear[[1]])
GOF_rep.with.time.linear <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.best.with.time,
                            index = 3)
(GOF_rep.with.time.linear[[1]])
GOF_cyc.with.time.linear <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.best.with.time,
                            index = 4)
(GOF_cyc.with.time.linear[[1]])


# Global

T_g.with.time.linear <- 
  (tan(pi*(0.5-GOF_rec.with.time.linear[[1]])) + 
   tan(pi*(0.5-GOF_trs.with.time.linear[[1]])) + 
   tan(pi*(0.5-GOF_rep.with.time.linear[[1]])) +
   tan(pi*(0.5-GOF_cyc.with.time.linear[[1]])))/4
GOF_global.with.time.linear <- 1/2 - atan(T_g.with.time.linear)/pi

#-------------------------------------------------------------------------------
# 7. Enriching the model formulation

# 7.1 Fitting a more complex model formulation
# Including non-linearity for reciprocity and transitivity + random effects
Weight <- Rec_Mat <- cbind(dat.gam$timeRec1,
                           dat.gam$timeRec2)

Trs_Mat <- cbind(dat.gam$timeTriad1,
                 dat.gam$timeTriad2)

Weight[,1] <- 1
Weight[,2] <- -1

s1 <- dat.gam$s1
s2 <- dat.gam$s2
ss <- factor(c(s1,s2))
dim(ss) <- c(length(s1),2)
Ls <- matrix(1,length(s1),2); Ls[,2] <- -1

r1 <- dat.gam$r1
r2 <- dat.gam$r2
rr <- factor(c(r1,r2))
dim(rr) <- c(length(r1),2)
Lr <- matrix(1,length(r1),2); Lr[,2] <- -1

gam.fitNL.re <- gam(y ~ - 1 + s(Rec_Mat, by=Weight) +
                      s(Trs_Mat, by=Weight) + iner + cyc +
                      s(ss, by=Ls, bs="re") + 
                      s(rr, by=Lr, bs="re"),
                    data=dat.gam,
                    family="binomial"(link = 'logit'), 
                    method="REML")

# 7.2. Fitting the relative GOF

coefficients(gam.fitNL.re)[1]
GOF_rep.with.time.nle.re <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 1)
(GOF_rep.with.time.nle.re[[1]])
  
coefficients(gam.fitNL.re)[2]
GOF_cyc.with.time.nle.re <- GOF_univariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 2)
(GOF_cyc.with.time.nle.re[[1]])

coefficients(gam.fitNL.re)[3:11]
GOF_rec.with.time.nle.re <- GOF_multivariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 3:11, 
                            BB.stat = BB.stat_9[[2]])
(GOF_rec.with.time.nle.re[[1]])

coefficients(gam.fitNL.re)[12:20]
GOF_trs.with.time.nle.re <- GOF_multivariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 12:20, 
                            BB.stat = BB.stat_9[[2]])
(GOF_trs.with.time.nle.re[[1]])

GOF_sre.with.time.nle.re <- GOF_multivariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 21:179, 
                            BB.stat = BB.stat_159[[2]])
(GOF_sre.with.time.nle.re[[1]])

GOF_rre.with.time.nle.re <- GOF_multivariate(data = 
                            dat.gam,
                            gam.fit = gam.fitNL.re,
                            index = 180:338, 
                            BB.stat = BB.stat_159[[2]])
(GOF_rre.with.time.nle.re[[1]])

# Global

T_g.with.time.nle.re <- 
  (tan(pi*(0.5-GOF_rec.with.time.nle.re[[1]])) + 
     tan(pi*(0.5-GOF_trs.with.time.nle.re[[1]])) + 
     tan(pi*(0.5-GOF_rep.with.time.nle.re[[1]])) +
     tan(pi*(0.5-GOF_cyc.with.time.nle.re[[1]])))/4
GOF_global.with.time.nle.re <- 1/2 - atan(T_g.with.time.nle.re)/pi

T_g.with.time.nle.re_with_re <- 
  (tan(pi*(0.5-GOF_rec.with.time.nle.re[[1]])) + 
     tan(pi*(0.5-GOF_trs.with.time.nle.re[[1]])) + 
     tan(pi*(0.5-GOF_rep.with.time.nle.re[[1]])) +
     tan(pi*(0.5-GOF_cyc.with.time.nle.re[[1]])) +
     tan(pi*(0.5-GOF_sre.with.time.nle.re[[1]])) +
     tan(pi*(0.5-GOF_rre.with.time.nle.re[[1]])))/6
GOF_global.with.time.nle.re_with_re <- 1/2 - atan(T_g.with.time.nle.re_with_re)/pi

# 7.3 Fitting a even more complex model formulation
# Including non-linearity for cyclic closure and repetition as well

Cyc_Mat <- cbind(dat.gam$timeCyc1,
                 dat.gam$timeCyc2)

Rep_Mat <- cbind(dat.gam$timeIner1,
                 dat.gam$timeIner2)

gam.fitNL_complete.re <- gam(formula = y ~  - 1 + s(Rec_Mat, by=Weight) +
                            s(Trs_Mat, by=Weight) +
                            s(Rep_Mat, by=Weight) +
                            s(Cyc_Mat, by=Weight) +
                            s(ss, by=Ls, bs="re") + 
                            s(rr, by=Lr, bs="re"),
                          family="binomial"(link = 'logit'),
                          method="REML", data=dat.gam)

# 7.4. Fitting the relative GOF

coefficients(gam.fitNL_complete.re)[1:9]
GOF_rec.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 1:9, 
                                     BB.stat = BB.stat_9[[2]])
(GOF_rec.with.time.nle_complete.re[[1]])

coefficients(gam.fitNL_complete.re)[10:18]
GOF_trs.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 10:18, 
                                     BB.stat = BB.stat_9[[2]])
(GOF_trs.with.time.nle_complete.re[[1]])

coefficients(gam.fitNL_complete.re)[19:27]
GOF_rep.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 19:27,
                                     BB.stat = BB.stat_9[[2]])
(GOF_rep.with.time.nle_complete.re[[1]])

coefficients(gam.fitNL_complete.re)[28:36]
GOF_cyc.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 28:36,
                                     BB.stat = BB.stat_9[[2]])
(GOF_cyc.with.time.nle_complete.re[[1]])

coefficients(gam.fitNL_complete.re)[37:195]
GOF_sre.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 37:195, 
                                     BB.stat = BB.stat_159[[2]])
(GOF_sre.with.time.nle_complete.re[[1]])

coefficients(gam.fitNL_complete.re)[196:354]
GOF_rre.with.time.nle_complete.re <- GOF_multivariate(data = 
                                     dat.gam,
                                     gam.fit = gam.fitNL_complete.re,
                                     index = 196:354, 
                                     BB.stat = BB.stat_159[[2]])
(GOF_rre.with.time.nle_complete.re[[1]])

# Global

T_g.with.time.nle_complete.re <- 
  (tan(pi*(0.5-GOF_rec.with.time.nle_complete.re[[1]])) + 
     tan(pi*(0.5-GOF_trs.with.time.nle_complete.re[[1]])) + 
     tan(pi*(0.5-GOF_rep.with.time.nle_complete.re[[1]])) +
     tan(pi*(0.5-GOF_cyc.with.time.nle_complete.re[[1]])))/4
GOF_global.with.time.nle_complete.re <- 1/2 - 
  atan(T_g.with.time.nle_complete.re)/pi

T_g.with.time.nle_complete.re_with_re <- 
  (tan(pi*(0.5-GOF_rec.with.time.nle_complete.re[[1]])) + 
     tan(pi*(0.5-GOF_trs.with.time.nle_complete.re[[1]])) + 
     tan(pi*(0.5-GOF_rep.with.time.nle_complete.re[[1]])) +
     tan(pi*(0.5-GOF_cyc.with.time.nle_complete.re[[1]])) +
     tan(pi*(0.5-GOF_sre.with.time.nle_complete.re[[1]])) +
     tan(pi*(0.5-GOF_rre.with.time.nle_complete.re[[1]])))/6
GOF_global.with.time.nle_complete.re_with_re <- 1/2 - 
  atan(T_g.with.time.nle_complete.re_with_re)/pi

#-------------------------------------------------------------------------------

save(gam.best.without.time,
     gam.best.with.time, 
     gam.fitNL.re, 
     gam.fitNL_complete.re,
     GOF_rec.without.time.linear,
     GOF_rec.with.time.linear,
     GOF_rec.with.time.nle.re,
     GOF_rec.with.time.nle_complete.re,
     GOF_trs.without.time.linear,
     GOF_trs.with.time.linear,
     GOF_trs.with.time.nle.re,
     GOF_trs.with.time.nle_complete.re,
     GOF_rep.without.time.linear,
     GOF_rep.with.time.linear,
     GOF_rep.with.time.nle.re,
     GOF_rep.with.time.nle_complete.re,
     GOF_cyc.without.time.linear,
     GOF_cyc.with.time.linear,
     GOF_cyc.with.time.nle.re,
     GOF_cyc.with.time.nle_complete.re,
     GOF_sre.with.time.nle.re,
     GOF_rre.with.time.nle.re,
     GOF_sre.with.time.nle_complete.re,
     GOF_rre.with.time.nle_complete.re,
     GOF_global.without.time.linear,
     GOF_global.with.time.linear,
     GOF_global.with.time.nle.re,
     GOF_global.with.time.nle.re_with_re,
     GOF_global.with.time.nle_complete.re,
     GOF_global.with.time.nle_complete.re_with_re, 
     file="02-Application/output/models_and_GOFs.RData")

#-------------------------------------------------------------------------------
# 8. Results and Visualization

source("00-Functions/functions.R")
load(file="02-Application/output/dat_gam.RData")
load(file="02-Application/output/fitted_gams_application.RData")
load(file="02-Application/output/bb_bridges.RData")
load(file="02-Application/output/models_and_GOFs.RData")

# 8.1. Towards a correctly specified model
GOF_global.without.time.linear
GOF_global.with.time.linear
GOF_global.with.time.nle.re
GOF_global.with.time.nle_complete.re

# 8.2. According to Likelihood-Based Model Selection
AIC(gam.best.without.time)
AIC(gam.best.with.time)
AIC(gam.fitNL.re)
AIC(gam.fitNL_complete.re)

pdf("02-Application/pictures/application-results.pdf", width = 20, heigh=10)
par(mfrow=c(1,2))

plot(dat.gam$stp,
     apply(GOF_rep.without.time.linear[[2]], 
           1, function(x) crossprod(x,x)),
     lwd=1, lty=1, type="l", 
     col=pal.blue[8], 
     main="Testing GOF of Repetition",
     xlab = "Time", 
     ylab = "Martingale-Residual Type Process")

lines(dat.gam$stp,
      apply(GOF_rep.with.time.linear[[2]], 
            1, function(x) crossprod(x,x)), 
      lwd=1, lty=3, 
      col=pal.blue[8])

lines(dat.gam$stp,
      apply(GOF_rep.with.time.nle.re[[2]], 
            1, function(x) crossprod(x,x)), 
      lwd=1, lty=1, 
      col=pal.yellow[6])

lines(dat.gam$stp,
      apply(GOF_rep.with.time.nle_complete.re[[2]], 
            1, function(x) crossprod(x,x))/9, 
      lwd=1, lty=3, 
      col=pal.yellow[6])

legend("topright", 
       legend=c("id-linear", 
                "time-linear", 
                "NLE for rec&trs", 
                "NLE for all"),
       col=c(pal.blue[8], 
             pal.blue[8],
             pal.yellow[6], 
             pal.yellow[6]), 
       lwd=c(1,1,1,1), 
       lty=c(1,3,1,3),
       cex=0.8)

seq.s <- seq(0,1, length.out=nrow(dat.gam))
plotting.data.rep <- data.frame(seq.s=seq.s, 
                                GOF.curve=
                                  apply(GOF_rep.with.time.nle_complete.re[[2]], 
                                        1, function(x) crossprod(x,x)))

seq.s.theory <- seq(0,1, length.out=2000)

plot(plotting.data.rep,
  type="l", 
  lwd=1,
  col=pal.blue[8], 
  lty=1,
  ylim= c(0,10), 
  main="GOF curves - Model assuming NLE for all the dynamics included", 
  xlab = "Time", 
  ylab = "Martingale-Residual Type Process"
)
for (iter in 1:200){
  lines(
    seq.s.theory,
    BB.stat_9[[1]][,iter], 
    col=pal.greys[3], 
    lwd=0.5, 
    lty=2)
}
lines(
  plotting.data.rep, 
  lwd=1,
  col=pal.blue[8], 
  lty=1,
)
lines(
  seq.s,
  apply(GOF_rec.with.time.nle_complete.re[[2]], 
        1, function(x) crossprod(x,x)),
  lwd=1,
  col=pal.blue[8], 
  lty=3,
)
lines(
  seq.s,
  apply(GOF_trs.with.time.nle_complete.re[[2]], 
        1, function(x) crossprod(x,x)),
  lwd=1,
  col=pal.yellow[6], 
  lty=1,
)
lines(
  seq.s,
  apply(GOF_cyc.with.time.nle_complete.re[[2]], 
        1, function(x) crossprod(x,x)),
  lwd=1,
  col=pal.yellow[6], 
  lty=3,
)

legend("topright", 
       legend=c("Repetition", 
                "Reciprocity", 
                "Transitive Closure", 
                "Cyclic Closure"),
       col=c(pal.blue[8], 
             pal.blue[8],
             pal.yellow[6], 
             pal.yellow[6]), 
       lwd=c(1,1,1,1), 
       lty=c(1,3,1,3), 
       cex=0.8)
dev.off()
#-------------------------------------------------------------------------------
