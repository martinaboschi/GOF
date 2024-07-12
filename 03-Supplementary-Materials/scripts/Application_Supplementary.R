#-------------------------------------------------------------------------------
# Script: Application_Supplementary.R
# Author: Martina Boschi
# Date: May 2024
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 1. Importing packages and built-in functions
source("00-Functions/functions.R")
source("00-Functions/covariate.R")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 2. Import the results of the analysis
load(file="02-Application/output/analysis.RData")

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

library(lubridate)
email <- read.csv("02-Application/input/manufacturing.csv", sep = ";")
sel <- which(email[,1]==email[,2])
email <- email[-sel,]
email <- unique(email)
sel <- unique(email$EventDate)
o <- !duplicated(email$EventDate)
email <- email[o,]
hours <- hour(email$EventDate)


hist(hours, freq=FALSE, 
     main="Email hours", 
     xlab="Hour of the day", 
     ylab="Density", 
     col=pal.blue[4])

#-------------------------------------------------------------------------------
# 3. Predict coefficients according to the fitted model

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
# 4. Visualization

## 4.1. Reciprocity

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

pdf("reciprocity.pdf", width =10, height=10)
rec_plot + 
  labs(title = "Non-Linear Effect of Reciprocity",
       x = "Time from the Last Reciprocal Event",
       y = "Effect") +
  scale_fill_discrete(name="Time-Window",
                      labels=labels) + 
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20",
                                    size = 12, angle = 90,
                                    hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20",
                                    size = 12, angle = 0,
                                    hjust = 1, vjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20",
                                   size = 12,
                                   face = "plain",),
        title = element_text(color = "black",
                             size = 14, angle = 0,
                             hjust = .5, vjust = .5,
                             family="serif", face = "bold"))
dev.off()

## 4.2. Repetition

xstart = c(0,0.083,0.167,0.33,0.833,1.67)
xend = c(0.083,0.167,0.33,0.833,1.67,2.5)
labels = c("00-02h",
           "02-04h",
           "04-08h",
           "08-20h",
           "20-40h",
           "40-60h")
rects = data.frame(xstart = xstart,
                   xend = xend,
                   col = labels)

rep_plot <- ggplot() +
  theme_classic() +
  geom_rect(data = rects, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = -Inf, 
                              ymax = Inf, 
                              fill = col),
            alpha = 0.5) +
  geom_line(data = data.frame(dat.gam$or.inert1[dat.gam$or.inert1<7], 
                              pred.rep[dat.gam$or.inert1<7]), 
            aes(x=dat.gam$or.inert1[dat.gam$or.inert1<7],
                y=pred.rep[dat.gam$or.inert1<7]), 
            size=1)
for(x.line in c(xstart,xend[length(xend)])){
  rep_plot <- rep_plot +
    geom_vline(xintercept = x.line, 
               linetype = "dashed", 
               color = "black")
}

pdf("repetition.pdf", width =10, height=10)
rep_plot + 
  labs(title = "Non-Linear Effect of Repetition",
       x = "Time from the Last Event",
       y = "Effect") +
  scale_fill_discrete(name="Time-Window",
                      labels=labels) + 
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20",
                                    size = 12, angle = 90,
                                    hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20",
                                    size = 12, angle = 0,
                                    hjust = 1, vjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20",
                                   size = 12,
                                   face = "plain",),
        title = element_text(color = "black",
                             size = 14, angle = 0,
                             hjust = .5, vjust = .5,
                             family="serif", face = "bold"))
dev.off()

## 4.3. Transitive Closure

xstart = c(0,0.083,0.167,0.33,0.67,0.916,1.08,1.70)
xend = c(0.083,0.167,0.33,0.67,0.916,1.08,1.70,2.5)
labels = c("00-02h",
           "02-04h",
           "04-08h",
           "08-16h",
           "16-22h",
           "22-28h",
           "28-42h",
           "42-56h")
rects = data.frame(xstart = xstart,
                   xend = xend,
                   col = labels)

trs_plot <- ggplot() +
  theme_classic() +
  geom_rect(data = rects, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = -Inf, 
                              ymax = Inf, 
                              fill = col),
            alpha = 0.5) +
  geom_line(data = data.frame(dat.gam$or.triadic1[dat.gam$or.triadic1<7], 
                              pred.triadic[dat.gam$or.triadic1<7]), 
            aes(x=dat.gam$or.triadic1[dat.gam$or.triadic1<7],
                y=pred.triadic[dat.gam$or.triadic1<7]), 
            size=1)
for(x.line in c(xstart,xend[length(xend)])){
  trs_plot <- trs_plot +
    geom_vline(xintercept = x.line, 
               linetype = "dashed", 
               color = "black")
}

pdf("transitive.pdf", width =10, height=10)
trs_plot + 
  labs(title = "Non-Linear Effect of Transtive Closure",
       x = "Time from the Second-Leg Event",
       y = "Effect") +
  scale_fill_discrete(name="Time-Window",
                      labels=labels) + 
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20",
                                    size = 12, angle = 90,
                                    hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20",
                                    size = 12, angle = 0,
                                    hjust = 1, vjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20",
                                   size = 12,
                                   face = "plain",),
        title = element_text(color = "black",
                             size = 14, angle = 0,
                             hjust = .5, vjust = .5,
                             family="serif", face = "bold"))
dev.off()

## 4.4. Cyclic Closure

xstart = c(0,0.5,0.875,1.5)
xend = c(0.5,0.875,1.5,2.5)
labels = c("00-12h",
           "12-21h",
           "21-36h",
           "36-60h")
rects = data.frame(xstart = xstart,
                   xend = xend,
                   col = labels)

cyc_plot <- ggplot() +
  theme_classic() +
  geom_rect(data = rects, aes(xmin = xstart, 
                              xmax = xend, 
                              ymin = -Inf, 
                              ymax = Inf, 
                              fill = col),
            alpha = 0.5) +
  geom_line(data = data.frame(dat.gam$or.cyclic1[dat.gam$or.cyclic1<7], 
                              pred.cyc[dat.gam$or.cyclic1<7]), 
            aes(x=dat.gam$or.cyclic1[dat.gam$or.cyclic1<7],
                y=pred.cyc[dat.gam$or.cyclic1<7]), 
            size=1)
for(x.line in c(xstart,xend[length(xend)])){
  cyc_plot <- cyc_plot +
    geom_vline(xintercept = x.line, 
               linetype = "dashed", 
               color = "black")
}

pdf("cyclic.pdf", width =10, height=10)
cyc_plot + 
  labs(title = "Non-Linear Effect of Cyclic Closure",
       x = "Time from the Second-Leg Event",
       y = "Effect") +
  scale_fill_discrete(name="Time-Window",
                      labels=labels) + 
  theme(axis.text.y = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(color = "grey20",
                                   size = 10, angle = 0,
                                   hjust = 1, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20",
                                    size = 12, angle = 90,
                                    hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(color = "grey20",
                                    size = 12, angle = 0,
                                    hjust = 1, vjust = 0, face = "plain"),
        legend.text = element_text(color = "grey20",
                                   size = 12,
                                   face = "plain",),
        title = element_text(color = "black",
                             size = 14, angle = 0,
                             hjust = .5, vjust = .5,
                             family="serif", face = "bold"))
dev.off()

