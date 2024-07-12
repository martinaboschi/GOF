#-------------------------------------------------------------------------------
# Script: Application_Results.R
# Author: Martina Boschi
# Date: July 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
source("00-Functions/covariate.R")

load(file="02-Application/output/dat_gam.RData")
load(file="02-Application/output/fitted_gams_application.RData")
load(file="02-Application/output/bb_bridges.RData")
load(file="02-Application/output/models_and_GOFs.RData")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
GOF_global.without.time.linear
GOF_global.with.time.linear
GOF_global.with.time.nle.re
GOF_global.with.time.nle_complete.re
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
AIC(gam.best.without.time)
AIC(gam.best.with.time)
AIC(gam.fitNL.re)
AIC(gam.fitNL_complete.re)
#-------------------------------------------------------------------------------

pdf("02-Application/pictures/application-finalmodel.pdf", width = 15, heigh=15)

dim.k <- 9
Bridge.data.NLE <- replicate(50,
                            {
                              BB.data <- replicate(dim.k, BBridge(0,0,N=2000-1))
                              process <- apply(BB.data, 1, function(x) crossprod(x,x))
                              return(process)
                            })
Bridge.stat.NLE <- apply(Bridge.data.NLE, 2, function(x) max(abs(x)))
seq.s = seq(0, 1, length.out = nrow(dat.gam))
seq.s_BB = seq(0, 1, length.out = 2000)

plot(
  seq.s,
  apply(GOF_rep.with.time.nle_complete.re[[2]], 1, function(x) crossprod(x,x)),
  type="l", 
  lwd=1.8,
  col=pal.blue[8], 
  lty=1,
  ylim= c(0,10), 
  main="GOF curves - NLE for all the dynamics", 
  xlab = "Time", 
  ylab = "Squared Norm of W",
  cex.main=2.5, 
  cex.lab=2
)
for (iter in 1:50){
  lines(
    seq.s_BB,
    Bridge.data.NLE[,iter], 
    col=pal.greys[3], 
    lwd=0.8, 
    lty=2)
}
lines(
  seq.s,
  apply(GOF_rec.with.time.nle_complete.re[[2]], 1, function(x) crossprod(x,x)),
  lwd=1.8,
  col=pal.blue[8], 
  lty=3,
)
lines(
  seq.s,
  apply(GOF_trs.with.time.nle_complete.re[[2]], 1, function(x) crossprod(x,x)),
  lwd=1.8,
  col=pal.yellow[6], 
  lty=1,
)
lines(
  seq.s,
  apply(GOF_cyc.with.time.nle_complete.re[[2]], 1, function(x) crossprod(x,x)),
  lwd=1.8,
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
       lwd=c(1.8,1.8,1.8,1.8), 
       lty=c(1,3,1,3), 
       cex=1.5)

#-------------------------------------------------------------------------------
dev.off()

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

new.data <- data.frame(Weight=1,
                       Rec_Mat = dat.gam$timeRec1,
                       Cyc_Mat = dat.gam$timeCyc1,
                       Rep_Mat = dat.gam$timeIner1,
                       Trs_Mat = dat.gam$timeTriad1,
                       ss = ss[,1],
                       Ls = 1,
                       rr = rr[,1],
                       Lr = 1)

dat.gam[, c("or.rec1",
            "or.rec2",
            "or.cyclic1",
            "or.cyclic2",
            "or.inert1",
            "or.inert2",
            "or.triadic1",
            "or.triadic2")] <- apply(dat.gam[c(
              "timeRec1",
              "timeRec2",
              "timeCyc1",
              "timeCyc2",
              "timeIner1",
              "timeIner2",
              "timeTriad1",
              "timeTriad2"
            )], 2, function(x) -log(x))
#-------------------------------------------------------------------------------
library(lubridate)
email <- read.csv("02-Application/input/manufacturing.csv", sep = ";")
sel <- which(email[,1]==email[,2])
email <- email[-sel,]
email <- unique(email)
sel <- unique(email$EventDate)
o <- !duplicated(email$EventDate)
email <- email[o,]
hours <- hour(email$EventDate)
hist(hours)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
pred.triadic  <- predict(gam.fitNL_complete.re,
                         exclude = c("s(ss):Ls",
                                     "s(rr):Lr",
                                     "s(Rec_Mat):Weight", 
                                     "s(Rep_Mat):Weight", 
                                     "s(Cyc_Mat):Weight"), 
                         type = "terms", 
                         newdata = new.data)

pred.rec  <- predict(gam.fitNL_complete.re,
                     exclude = c("s(ss):Ls",
                                 "s(rr):Lr",
                                 "s(Trs_Mat):Weight", 
                                 "s(Rep_Mat):Weight", 
                                 "s(Cyc_Mat):Weight"), 
                     type = "terms", 
                     newdata = new.data)

pred.rep  <- predict(gam.fitNL_complete.re,
                     exclude = c("s(ss):Ls",
                                 "s(rr):Lr",
                                 "s(Trs_Mat):Weight", 
                                 "s(Rec_Mat):Weight", 
                                 "s(Cyc_Mat):Weight"), 
                     type = "terms", 
                     newdata = new.data)

pred.cyc  <- predict(gam.fitNL_complete.re,
                     exclude = c("s(ss):Ls",
                                 "s(rr):Lr",
                                 "s(Trs_Mat):Weight", 
                                 "s(Rec_Mat):Weight", 
                                 "s(Rep_Mat):Weight"), 
                     type = "terms", 
                     newdata = new.data)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
xstart = c(0,0.33,0.66,0.916,1.08, 1.67)
xend = c(0.33,0.66,0.916,1.08,1.67,2.66)
labels = c("00-08h",
           "08-16h",
           "16-22h", 
           "22-26h",
           "26-40h",
           "40-64h")
rects = data.frame(xstart = xstart,
                   xend = xend,
                   col = labels)
rec_plot <- ggplot() +
  theme_classic() +
  geom_rect(data = rects, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = -Inf, 
                              ymax = Inf, 
                              fill = col),
            alpha = 0.5) +
  geom_line(data = data.frame(dat.gam$or.rec1[dat.gam$or.rec1<7], 
                              pred.rec[dat.gam$or.rec1<7]), 
            aes(x=dat.gam$or.rec1[dat.gam$or.rec1<7],
                y=pred.rec[dat.gam$or.rec1<7]), 
            size=1)
for(x.line in c(xstart,xend[length(xend)])){
  rec_plot <- rec_plot +
    geom_vline(xintercept = x.line, 
               linetype = "dashed", 
               color = "black")
}

pdf("02-Application/pictures/application-reciprocity.pdf", width =15, height=15)
rec_plot + 
  labs(title = "Non-Linear Effect of Reciprocity",
       x = "Time from the Last Reciprocal Event",
       y = "Effect") +
  scale_fill_discrete(name="Time-Window",
                      labels=labels) + 
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 12, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "grey20",
                                   size = 12, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20",
                                    size = 14, angle = 90,
                                    hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20",
                                    size = 14, angle = 0,
                                    hjust = 1, vjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20",
                                   size = 18,
                                   face = "plain",),
        title = element_text(color = "black",
                             size = 22, angle = 0,
                             hjust = .5, vjust = .5,
                             family="serif", face = "bold"))

dev.off()

