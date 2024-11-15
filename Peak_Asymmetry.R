rm(list=ls())
library(MCMCglmm)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(gridExtra)

############################################################
#### Model: Asymmetry in the caterpillar abundance peak ####
############################################################

cater_habitat <- read.csv("~/cater_habitat_data.csv")
load("~/CurveShape23.RData")  
summary(CurveShape23)

#### Checking model fits the data and converged ####

plot(CurveShape23)  

CurveShape23.Sim<-simulate(CurveShape23,nsim=1000) #simulate 1000 times

par(mfcol=c(1,1))
hist(apply(CurveShape23.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(CurveShape23.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data

## Differentiation to calculate peak date
# (-b +/- sqrt(b^2-4*a*c))/(2*a)
cdf <- data.frame(CurveShape23$Sol[,1:4])
cdf$a <- 3*CurveShape23$Sol[,4]
cdf$b <- 2*CurveShape23$Sol[,3]
cdf$c <- CurveShape23$Sol[,2]

cdf$pd <- ((-cdf$b - sqrt(cdf$b^2-4*cdf$a*cdf$c))/(2*cdf$a))

## Calculate peak height
cdf$ph <- CurveShape23$Sol[,1]+CurveShape23$Sol[,2]*cdf$pd+CurveShape23$Sol[,3]*cdf$pd^2+CurveShape23$Sol[,4]*cdf$pd^3

## Calculate width either side at 50% peak height
# Alter 'peak height' to find dates the peak is at a certain proportion of the true peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/2))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.5[i] <- A[1]
  cdf$r2.5[i] <- A[2]
  cdf$r3.5[i] <- A[3]
}

# turn into real numbers
cdf$r1.5 <- Re(cdf$r1.5)[abs(Im(cdf$r1.5)) < 1e-6]
cdf$r2.5 <- Re(cdf$r2.5)[abs(Im(cdf$r2.5)) < 1e-6]
cdf$r3.5 <- Re(cdf$r3.5)[abs(Im(cdf$r3.5)) < 1e-6]

#remove the root that is outside of the data range
cdf$r1.5 <- ifelse(cdf$r1.5<min(cater_habitat$datescaled),NA,cdf$r1.5)
cdf$r2.5 <- ifelse(cdf$r2.5<min(cater_habitat$datescaled),NA,cdf$r2.5)
cdf$r3.5 <- ifelse(cdf$r3.5<min(cater_habitat$datescaled),NA,cdf$r3.5)

#get real dates
roots.5 <- data.frame(r1.5=cdf$r1.5,r2.5=cdf$r2.5,r3.5=cdf$r3.5)
cdf$root1.5 <- ggfortify::unscale((apply(roots.5, 1, min, na.rm=TRUE)), center= 147.4742, scale=14.19027)  
cdf$root2.5 <- ggfortify::unscale((apply(roots.5, 1, max, na.rm=TRUE)), center= 147.4742, scale=14.19027) 

#days either side of peak
cdf$left.5 <- ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027)-cdf$root1.5
cdf$left.5 <- cdf$left.5$var1
cdf$right.5 <- cdf$root2.5-(ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027))
cdf$right.5 <- cdf$right.5$V1

#proportion of width
cdf$width.5 <- cdf$left.5+cdf$right.5
cdf$propleft.5 <- cdf$left.5/cdf$width.5
cdf$propright.5 <- cdf$right.5/cdf$width.5

mean(cdf$propleft.5) # 0.550347  
HPDinterval(mcmc(cdf$propleft.5)) # 0.5391472 0.5606735
mean(cdf$propright.5) # 0.449653
HPDinterval(mcmc(cdf$propright.5)) # 0.4393265 0.4608528
mean(cdf$width.5) # 25.10665
HPDinterval(mcmc(cdf$width.5)) # 23.48493 27.03101

#### Calculate width either side at 25% peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/4))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.25[i] <- A[1]
  cdf$r2.25[i] <- A[2]
  cdf$r3.25[i] <- A[3]
}

cdf$r1.25 <- Re(cdf$r1.25)[abs(Im(cdf$r1.25)) < 1e-6]
cdf$r2.25 <- Re(cdf$r2.25)[abs(Im(cdf$r2.25)) < 1e-6]
cdf$r3.25 <- Re(cdf$r3.25)[abs(Im(cdf$r3.25)) < 1e-6]

cdf$r1.25 <- ifelse(cdf$r1.25<min(cater_habitat$datescaled),NA,cdf$r1.25)
cdf$r2.25 <- ifelse(cdf$r2.25<min(cater_habitat$datescaled),NA,cdf$r2.25)
cdf$r3.25 <- ifelse(cdf$r3.25<min(cater_habitat$datescaled),NA,cdf$r3.25)

roots.25 <- data.frame(r1.25=cdf$r1.25,r2.25=cdf$r2.25,r3.25=cdf$r3.25)
cdf$root1.25 <- ggfortify::unscale((apply(roots.25, 1, min, na.rm=TRUE)), center= 147.4742, scale=14.19027)   
cdf$root2.25 <- ggfortify::unscale((apply(roots.25, 1, max, na.rm=TRUE)), center= 147.4742, scale=14.19027) 

cdf$left.25 <- ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027)-cdf$root1.25
cdf$left.25 <- cdf$left.25$var1
cdf$right.25 <- cdf$root2.25-ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027)
cdf$right.25 <- cdf$right.25$V1
cdf$width.25 <- cdf$left.25+cdf$right.25
cdf$propleft.25 <- cdf$left.25/cdf$width.25
cdf$propright.25 <- cdf$right.25/cdf$width.25

mean(cdf$propleft.25) # 0.5767176   
HPDinterval(mcmc(cdf$propleft.25)) # 0.5575208 0.5964568
mean(cdf$propright.25) # 0.4232824
HPDinterval(mcmc(cdf$propright.25)) # 0.4035432 0.4424792
mean(cdf$width.25) # 36.71111
HPDinterval(mcmc(cdf$width.25)) # 33.99017 39.14467

#### Calculate width either side at 75% peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])*0.75))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.75[i] <- A[1]
  cdf$r2.75[i] <- A[2]
  cdf$r3.75[i] <- A[3]
}

cdf$r1.75 <- Re(cdf$r1.75)[abs(Im(cdf$r1.75)) < 1e-6]
cdf$r2.75 <- Re(cdf$r2.75)[abs(Im(cdf$r2.75)) < 1e-6]
cdf$r3.75 <- Re(cdf$r3.75)[abs(Im(cdf$r3.75)) < 1e-6]

cdf$r1.75 <- ifelse(cdf$r1.75<min(cater_habitat$datescaled),NA,cdf$r1.75)
cdf$r2.75 <- ifelse(cdf$r2.75<min(cater_habitat$datescaled),NA,cdf$r2.75)
cdf$r3.75 <- ifelse(cdf$r3.75<min(cater_habitat$datescaled),NA,cdf$r3.75)

roots.75 <- data.frame(r1.75=cdf$r1.75,r2.75=cdf$r2.75,r3.75=cdf$r3.75)
cdf$root1.75 <- ggfortify::unscale((apply(roots.75, 1, min, na.rm=TRUE)), center= 147.4742, scale=14.19027)   
cdf$root2.75 <- ggfortify::unscale((apply(roots.75, 1, max, na.rm=TRUE)), center= 147.4742, scale=14.19027)  

cdf$left.75 <- ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027)-cdf$root1.75
cdf$left.75 <- cdf$left.75$var1
cdf$right.75 <- cdf$root2.75-ggfortify::unscale(cdf$pd, center= 147.4742, scale=14.19027)
cdf$right.75 <- cdf$right.75$V1
cdf$width.75 <- cdf$left.75+cdf$right.75
cdf$propleft.75 <- cdf$left.75/cdf$width.75
cdf$propright.75 <- cdf$right.75/cdf$width.75

mean(cdf$propleft.75) # 0.531324 
HPDinterval(mcmc(cdf$propleft.75)) # 0.5247186 0.5372034
mean(cdf$propright.75) # 0.468676
HPDinterval(mcmc(cdf$propright.75)) # 0.4627966 0.4752814
mean(cdf$width.75) # 15.92487
HPDinterval(mcmc(cdf$width.75)) # 14.88575 17.18528

#### Plot curves
# colour: 
mycol <- rgb(0, 120, 0, max = 250, alpha = 15, names = "greentrans")

# mean trend
dayscal <- seq(min(cater_habitat$datescaled),max(cater_habitat$datescaled),0.001)
curve <- mean(CurveShape23$Sol[,1])+mean(CurveShape23$Sol[,2])*dayscal+mean(CurveShape23$Sol[,3])*dayscal^2+mean(CurveShape23$Sol[,4])*dayscal^3
days <- unscale(dayscal, center= 147.4742, scale=14.19027)$V1

# metrics for lines 
quart <- data.frame(qd=seq(mean(cdf$root1.25$V1),(mean(cdf$root2.25$V1)),0.1))
quart$qh <- mean(exp(cdf$ph)/4)
half <- data.frame(hd=seq((mean(cdf$root1.5$V1)+0.2),(mean(cdf$root2.5$V1)+0.2),0.1))
half$hh <- mean(exp(cdf$ph)/2)
tquart <- data.frame(tqd=seq((mean(cdf$root1.75$V1)+0.45),(mean(cdf$root2.75$V1)),0.1))
tquart$tqh <- mean(exp(cdf$ph)*0.75)

## Plotting posterior distribution curves with percentage of distribution to either side of the peak
par(mfcol=c(1,1),mar=c(3.9, 3.8, 1, 1), cex=1.4, las=1)
plot(days,exp(curve), type="l", ylim=c(0,0.08), xlab="Ordinal Date", ylab="Abundance", yaxs="i")

for(i in 1:1700){ 
  A <- CurveShape23$Sol[i,1]+CurveShape23$Sol[i,2]*dayscal+CurveShape23$Sol[i,3]*dayscal^2+CurveShape23$Sol[i,4]*dayscal^3
  points(days, exp(A), type="l", col=mycol, lwd=0.5)
}

points(quart$qd, quart$qh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(half$hd, half$hh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(tquart$tqd, tquart$tqh, type="l", lty="dashed", lwd=0.7, col="gray66")
abline(v=mean(unscale(cdf$pd, center= 147.4742, scale=14.19027)$var1), lwd=0.8, lty="dashed", col="gray66")
points(days, exp(curve), type="l")

text(151.8, 0.01, "57.7%", cex=0.9, col="gray40")
text(159.2, 0.01, "42.3%", cex=0.9, col="gray40")
text(151.8, 0.0225, "55.0%", cex=0.9, col="gray40")
text(159.2, 0.0225, "45.0%", cex=0.9, col="gray40")
text(151.8, 0.035, "53.1%", cex=0.9, col="gray40")
text(159.2, 0.035, "46.9%", cex=0.9, col="gray40")

text(124, mean(quart$qh), "0.25", cex=0.9, col="gray66")
text(124.4, mean(half$hh), "0.5", cex=0.9, col="gray66")
text(124, mean(tquart$tqh), "0.75", cex=0.9, col="gray66")

arrows(x0=127.5, y0=mean(half$hh), x1=136, y1=mean(half$hh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(quart$qh), x1=130.25, y1=mean(quart$qh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(tquart$tqh), x1=140, y1=mean(tquart$tqh), length=0.1, col="gray66") #Saved as 8"x8"


############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean

#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(CurveShape23$Sol[,1]),3)," (",
                      round(HPDinterval(CurveShape23$Sol[,1])[1],3)," - ",
                      round(HPDinterval(CurveShape23$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(CurveShape23$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(CurveShape23$Sol[,2]),3)," (",
                          round(HPDinterval(CurveShape23$Sol[,2])[1],3)," - ",
                          round(HPDinterval(CurveShape23$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(CurveShape23$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(CurveShape23$Sol[,3]),3)," (",
                           round(HPDinterval(CurveShape23$Sol[,3])[1],3)," - ",
                           round(HPDinterval(CurveShape23$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(CurveShape23$Sol[,3]))),
  c("Date³ (scaled)",paste(round(mean(CurveShape23$Sol[,4]),3)," (",
                           round(HPDinterval(CurveShape23$Sol[,4])[1],3)," - ",
                           round(HPDinterval(CurveShape23$Sol[,4])[2],3),")",sep=""),
    round(effectiveSize(CurveShape23$Sol[,4]))))

#### random terms
column<-1
siteyear1<-c("SiteYear- Intercept var",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                             round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-2
siteyear2<-c("SiteYear- Intercept:Date slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                          round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-3
siteyear3<-c("SiteYear- Intercept:Date² slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                           round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-4
siteyear4<-c("SiteYear- Intercept:Date³ slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                           round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-6
siteyear5<-c("SiteYear- Date slope var",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                              round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-7
siteyear6<-c("SiteYear- Date slope:Date² slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                            round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-8
siteyear7<-c("SiteYear- Date slope:Date³ slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                            round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-11
siteyear8<-c("SiteYear- Date² slope var",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                               round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-12
siteyear9<-c("SiteYear- Date² slope:Date³ slope covar",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                             round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                             round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape23$VCV[, column])))

column<-16
siteyear10<-c("SiteYear- Date³ slope var",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                                                round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                                                round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(CurveShape23$VCV[, column])))

column<-17
recorder<-c("Recorder",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape23$VCV[, column])))


column<-18
siteday<-c("Site Day",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                            round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                            round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(CurveShape23$VCV[, column])))


column<-19
treeID<-c("Tree ID",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                          round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                          round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(CurveShape23$VCV[, column])))

column<-20
residual<-c("Residual",paste(round(posterior.mode(CurveShape23$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape23$VCV[, column])))

random<-rbind(siteyear1,siteyear2,siteyear3,siteyear4,siteyear5,siteyear6,siteyear7,siteyear8,siteyear9,siteyear10, recorder, siteday, treeID, residual)

write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/TableS7_CurveShape.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)

