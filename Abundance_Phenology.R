rm(list=ls())
library(MCMCglmm)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(gridExtra)

#######################################################################
#### Model: Phenological distribution of abundance among tree taxa ####
#######################################################################

cater_habitat <- read.csv("~/cater_habitat_data.csv")
load("~/ATTC23.RData")  
summary(ATTC23)

#### Checking model fits the data and converged ####

plot(ATTC23) #look at fixed effect and random term trace plots 

ATTC23.Sim<-simulate(ATTC23,nsim=1000) #simulate 1000 times

par(mfcol=c(1,1))
hist(apply(ATTC23.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
#normally not too many rogue values but reduce x axis a bit to see main distribution relative to observed value 
hist(apply(ATTC23.Sim,2,sum), breaks=10000, xlim=c(0,100000))
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(ATTC23.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data


# tree taxon specific coefficients for date^2 (A), date (B) and intercept (C)
AlderC <- (ATTC23$Sol[,1] + ATTC23$Sol[,4])
AshC <- (ATTC23$Sol[,1] + ATTC23$Sol[,5])
BeechC <- (ATTC23$Sol[,1] + ATTC23$Sol[,6])
BirchC <- (ATTC23$Sol[,1] + ATTC23$Sol[,7])
ElmC <- (ATTC23$Sol[,1] + ATTC23$Sol[,8])
HazelC <- (ATTC23$Sol[,1] + ATTC23$Sol[,9])
OakC <- (ATTC23$Sol[,1] + ATTC23$Sol[,10])
RowanC <- (ATTC23$Sol[,1] + ATTC23$Sol[,11])
SycamoreC <- (ATTC23$Sol[,1] + ATTC23$Sol[,12])
WillowC <- (ATTC23$Sol[,1] + ATTC23$Sol[,13])

AlderB <- (ATTC23$Sol[,2] + ATTC23$Sol[,14])
AshB <- (ATTC23$Sol[,2] + ATTC23$Sol[,15])
BeechB <- (ATTC23$Sol[,2] + ATTC23$Sol[,16])
BirchB <- (ATTC23$Sol[,2] + ATTC23$Sol[,17])
ElmB <- (ATTC23$Sol[,2] + ATTC23$Sol[,18])
HazelB <- (ATTC23$Sol[,2] + ATTC23$Sol[,19])
OakB <- (ATTC23$Sol[,2] + ATTC23$Sol[,20])
RowanB <- (ATTC23$Sol[,2] + ATTC23$Sol[,21])
SycamoreB <- (ATTC23$Sol[,2] + ATTC23$Sol[,22])
WillowB <- (ATTC23$Sol[,2] + ATTC23$Sol[,23])

AlderA <- (ATTC23$Sol[,3] + ATTC23$Sol[,24])
AshA <- (ATTC23$Sol[,3] + ATTC23$Sol[,25])
BeechA <- (ATTC23$Sol[,3] + ATTC23$Sol[,26])
BirchA <- (ATTC23$Sol[,3] + ATTC23$Sol[,27])
ElmA <- (ATTC23$Sol[,3] + ATTC23$Sol[,28])
HazelA <- (ATTC23$Sol[,3] + ATTC23$Sol[,29])
OakA <- (ATTC23$Sol[,3] + ATTC23$Sol[,30])
RowanA <- (ATTC23$Sol[,3] + ATTC23$Sol[,31])
SycamoreA <- (ATTC23$Sol[,3] + ATTC23$Sol[,32])
WillowA <- (ATTC23$Sol[,3] + ATTC23$Sol[,33])

#fixed effect coefficients
MeanA <- (ATTC23$Sol[,3])
MeanB <- (ATTC23$Sol[,2])
MeanC <- (ATTC23$Sol[,1])

# unscaled date
dayscal <- seq(min(cater_habitat$datescaled),max(cater_habitat$datescaled),0.001)
days <- ggfortify::unscale(dayscal, center= 147.4742, scale=14.19027)$V1


#mean curve for each tree taxon and the average trend
AlderCurve <- exp(mean(AlderC) +mean(AlderB)*dayscal + mean(AlderA)*dayscal^2)
AshCurve <- exp(mean(AshC) +mean(AshB)*dayscal + mean(AshA)*dayscal^2)
BeechCurve <- exp(mean(BeechC) +mean(BeechB)*dayscal + mean(BeechA)*dayscal^2)
BirchCurve <- exp(mean(BirchC) +mean(BirchB)*dayscal + mean(BirchA)*dayscal^2)
ElmCurve <- exp(mean(ElmC) +mean(ElmB)*dayscal + mean(ElmA)*dayscal^2)
HazelCurve <- exp(mean(HazelC) +mean(HazelB)*dayscal + mean(HazelA)*dayscal^2)
OakCurve <- exp(mean(OakC) +mean(OakB)*dayscal + mean(OakA)*dayscal^2)
RowanCurve <- exp(mean(RowanC) +mean(RowanB)*dayscal + mean(RowanA)*dayscal^2)
SycamoreCurve <- exp(mean(SycamoreC) +mean(SycamoreB)*dayscal + mean(SycamoreA)*dayscal^2)
WillowCurve <- exp(mean(WillowC) +mean(WillowB)*dayscal + mean(WillowA)*dayscal^2)
MeanCurve <- exp(mean(MeanC) +mean(MeanB)*dayscal + mean(MeanA)*dayscal^2)

#### Peak date per taxa ####   -b/2a
MeanPD <- ggfortify::unscale((-MeanB/(2*MeanA)), center= 147.4742, scale=14.19027)$var1
AlderPD <-  ggfortify::unscale((-AlderB/(2*AlderA)), center= 147.4742, scale=14.19027)$var1
AshPD <-  ggfortify::unscale((-AshB/(2*AshA)), center= 147.4742, scale=14.19027)$var1
BeechPD <-  ggfortify::unscale((-BeechB/(2*BeechA)), center= 147.4742, scale=14.19027)$var1
BirchPD <-  ggfortify::unscale((-BirchB/(2*BirchA)), center= 147.4742, scale=14.19027)$var1
ElmPD <-  ggfortify::unscale((-ElmB/(2*ElmA)), center= 147.4742, scale=14.19027)$var1
HazelPD <-  ggfortify::unscale((-HazelB/(2*HazelA)), center= 147.4742, scale=14.19027)$var1
OakPD <-  ggfortify::unscale((-OakB/(2*OakA)), center= 147.4742, scale=14.19027)$var1
RowanPD <-  ggfortify::unscale((-RowanB/(2*RowanA)), center= 147.4742, scale=14.19027)$var1
SycamorePD <-  ggfortify::unscale((-SycamoreB/(2*SycamoreA)), center= 147.4742, scale=14.19027)$var1
WillowPD <-  ggfortify::unscale((-WillowB/(2*WillowA)), center= 147.4742, scale=14.19027)$var1

## Peak date diff to mean
AlderPDdif <- AlderPD-MeanPD
AshPDdif <- AshPD-MeanPD
BeechPDdif <- BeechPD-MeanPD
BirchPDdif <- BirchPD-MeanPD
ElmPDdif <- ElmPD-MeanPD
HazelPDdif <- HazelPD-MeanPD
OakPDdif <- OakPD-MeanPD
RowanPDdif <- RowanPD-MeanPD
SycamorePDdif <- SycamorePD-MeanPD
WillowPDdif <- WillowPD-MeanPD

## Peak date diff to oak
AlderPDOD <- AlderPD-OakPD
AshPDOD <- AshPD-OakPD
BeechPDOD <- BeechPD-OakPD
BirchPDOD <- BirchPD-OakPD
ElmPDOD <- ElmPD-OakPD
HazelPDOD <- HazelPD-OakPD
RowanPDOD <- RowanPD-OakPD
SycamorePDOD <- SycamorePD-OakPD
WillowPDOD <- WillowPD-OakPD

#### Peak height per taxa ####
MeanPH <- exp(MeanC +MeanB*(-MeanB/(2*MeanA)) + MeanA*(-MeanB/(2*MeanA))^2)
AlderPH <- exp(AlderC +AlderB*(-AlderB/(2*AlderA)) + AlderA*(-AlderB/(2*AlderA))^2)
AshPH <- exp(AshC +AshB*(-AshB/(2*AshA)) + AshA*(-AshB/(2*AshA))^2)
BeechPH <- exp(BeechC +BeechB*(-BeechB/(2*BeechA)) + BeechA*(-BeechB/(2*BeechA))^2)
BirchPH <- exp(BirchC +BirchB*(-BirchB/(2*BirchA)) + BirchA*(-BirchB/(2*BirchA))^2)
ElmPH <- exp(ElmC +ElmB*(-ElmB/(2*ElmA)) + ElmA*(-ElmB/(2*ElmA))^2)
HazelPH <- exp(HazelC +HazelB*(-HazelB/(2*HazelA)) + HazelA*(-HazelB/(2*HazelA))^2)
OakPH <- exp(OakC +OakB*(-OakB/(2*OakA)) + OakA*(-OakB/(2*OakA))^2)
RowanPH <- exp(RowanC +RowanB*(-RowanB/(2*RowanA)) + RowanA*(-RowanB/(2*RowanA))^2)
SycamorePH <- exp(SycamoreC +SycamoreB*(-SycamoreB/(2*SycamoreA)) + SycamoreA*(-SycamoreB/(2*SycamoreA))^2)
WillowPH <- exp(WillowC +WillowB*(-WillowB/(2*WillowA)) + WillowA*(-WillowB/(2*WillowA))^2)

## Peak height diff to mean
AlderPHdif <- AlderPH-MeanPH
AshPHdif <- AshPH-MeanPH
BeechPHdif <- BeechPH-MeanPH
BirchPHdif <- BirchPH-MeanPH
ElmPHdif <- ElmPH-MeanPH
HazelPHdif <- HazelPH-MeanPH
OakPHdif <- OakPH-MeanPH
RowanPHdif <- RowanPH-MeanPH
SycamorePHdif <- SycamorePH-MeanPH
WillowPHdif <- WillowPH-MeanPH

## Peak height diff to oak
AlderPHOD <- AlderPH-OakPH
AshPHOD <- AshPH-OakPH
BeechPHOD <- BeechPH-OakPH
BirchPHOD <- BirchPH-OakPH
ElmPHOD <- ElmPH-OakPH
HazelPHOD <- HazelPH-OakPH
RowanPHOD <- RowanPH-OakPH
SycamorePHOD <- SycamorePH-OakPH
WillowPHOD <- WillowPH-OakPH

## Peak height proportion of mean
AlderPHprop <- AlderPH/MeanPH
AshPHprop <- AshPH/MeanPH
BeechPHprop <- BeechPH/MeanPH
BirchPHprop <- BirchPH/MeanPH
ElmPHprop <- ElmPH/MeanPH
HazelPHprop <- HazelPH/MeanPH
OakPHprop <- OakPH/MeanPH
RowanPHprop <- RowanPH/MeanPH
SycamorePHprop <- SycamorePH/MeanPH
WillowPHprop <- WillowPH/MeanPH

## Peak height proportion of oak
AlderPHODprop <- AlderPH/OakPH
AshPHODprop <- AshPH/OakPH
BeechPHODprop <- BeechPH/OakPH
BirchPHODprop <- BirchPH/OakPH
ElmPHODprop <- ElmPH/OakPH
HazelPHODprop <- HazelPH/OakPH 
RowanPHODprop <- RowanPH/OakPH
SycamorePHODprop <- SycamorePH/OakPH
WillowPHODprop <- WillowPH/OakPH

#### Peak width per taxa (set height) ####   x= (-b +/- sqrt(b^2 - 4ac))/2a

# set at 0.01 - roughly half the height of the lowest curve (Alder)

MeanPW1 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.01)))))/(2*MeanA)
MeanPW2 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.01)))))/(2*MeanA)
MeanPW <- ggfortify::unscale((MeanPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((MeanPW1), center= 147.4742, scale=14.19027)
MeanPW <- MeanPW$var1 

AlderPW1 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.01)))))/(2*AlderA)
AlderPW2 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.01)))))/(2*AlderA)
AlderPW <- ggfortify::unscale((AlderPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((AlderPW1), center= 147.4742, scale=14.19027)
AlderPW <- AlderPW$var1 

AshPW1 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(0.01)))))/(2*AshA)
AshPW2 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(0.01)))))/(2*AshA)
AshPW <- ggfortify::unscale((AshPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((AshPW1), center= 147.4742, scale=14.19027)
AshPW <- AshPW$var1 

BeechPW1 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.01)))))/(2*BeechA)
BeechPW2 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.01)))))/(2*BeechA)
BeechPW <- ggfortify::unscale((BeechPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((BeechPW1), center= 147.4742, scale=14.19027)
BeechPW <- BeechPW$var1 

BirchPW1 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.01)))))/(2*BirchA)
BirchPW2 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.01)))))/(2*BirchA)
BirchPW <- ggfortify::unscale((BirchPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((BirchPW1), center= 147.4742, scale=14.19027)
BirchPW <- BirchPW$var1 

ElmPW1 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.01)))))/(2*ElmA)
ElmPW2 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.01)))))/(2*ElmA)
ElmPW <- ggfortify::unscale((ElmPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((ElmPW1), center= 147.4742, scale=14.19027)
ElmPW <- ElmPW$var1 

HazelPW1 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.01)))))/(2*HazelA)
HazelPW2 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.01)))))/(2*HazelA)
HazelPW <- ggfortify::unscale((HazelPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((HazelPW1), center= 147.4742, scale=14.19027)
HazelPW <- HazelPW$var1 

OakPW1 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(0.01)))))/(2*OakA)
OakPW2 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(0.01)))))/(2*OakA)
OakPW <- ggfortify::unscale((OakPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((OakPW1), center= 147.4742, scale=14.19027)
OakPW <- OakPW$var1 

RowanPW1 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.01)))))/(2*RowanA)
RowanPW2 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.01)))))/(2*RowanA)
RowanPW <- ggfortify::unscale((RowanPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((RowanPW1), center= 147.4742, scale=14.19027)
RowanPW <- RowanPW$var1 

SycamorePW1 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.01)))))/(2*SycamoreA)
SycamorePW2 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.01)))))/(2*SycamoreA)
SycamorePW <- ggfortify::unscale((SycamorePW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((SycamorePW1), center= 147.4742, scale=14.19027)
SycamorePW <- SycamorePW$var1 

WillowPW1 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.01)))))/(2*WillowA)
WillowPW2 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.01)))))/(2*WillowA)
WillowPW <- ggfortify::unscale((WillowPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((WillowPW1), center= 147.4742, scale=14.19027)
WillowPW <- WillowPW$var1 

## Peak width diff to mean
AlderPWdif <- AlderPW-MeanPW
AshPWdif <- AshPW-MeanPW
BeechPWdif <- BeechPW-MeanPW
BirchPWdif <- BirchPW-MeanPW
ElmPWdif <- ElmPW-MeanPW
HazelPWdif <- HazelPW-MeanPW
OakPWdif <- OakPW-MeanPW
RowanPWdif <- RowanPW-MeanPW
SycamorePWdif <- SycamorePW-MeanPW
WillowPWdif <- WillowPW-MeanPW

## Peak width diff to oak
AlderPWOD <- AlderPW-OakPW
AshPWOD <- AshPW-OakPW
BeechPWOD <- BeechPW-OakPW
BirchPWOD <- BirchPW-OakPW
ElmPWOD <- ElmPW-OakPW
HazelPWOD <- HazelPW-OakPW
RowanPWOD <- RowanPW-OakPW
SycamorePWOD <- SycamorePW-OakPW
WillowPWOD <- WillowPW-OakPW

#### Peak width per taxa (half each PH) ####   x= (-b +/- sqrt(b^2 - 4ac))/2a

# at half the height of each peak

MeanPW1.5 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW2.5 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW.5 <- ggfortify::unscale((MeanPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((MeanPW1.5), center= 147.4742, scale=14.19027)
MeanPW.5 <- MeanPW.5$var1 

AlderPW1.5 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW2.5 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW.5 <- ggfortify::unscale((AlderPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((AlderPW1.5), center= 147.4742, scale=14.19027)
AlderPW.5 <- AlderPW.5$var1 

AshPW1.5 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW2.5 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW.5 <- ggfortify::unscale((AshPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((AshPW1.5), center= 147.4742, scale=14.19027)
AshPW.5 <- AshPW.5$var1 

BeechPW1.5 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW2.5 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW.5 <- ggfortify::unscale((BeechPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((BeechPW1.5), center= 147.4742, scale=14.19027)
BeechPW.5 <- BeechPW.5$var1 

BirchPW1.5 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW2.5 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW.5 <- ggfortify::unscale((BirchPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((BirchPW1.5), center= 147.4742, scale=14.19027)
BirchPW.5 <- BirchPW.5$var1

ElmPW1.5 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW2.5 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW.5 <- ggfortify::unscale((ElmPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((ElmPW1.5), center= 147.4742, scale=14.19027)
ElmPW.5 <- ElmPW.5$var1 

HazelPW1.5 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW2.5 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW.5 <- ggfortify::unscale((HazelPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((HazelPW1.5), center= 147.4742, scale=14.19027)
HazelPW.5 <- HazelPW.5$var1 

OakPW1.5 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW2.5 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW.5 <- ggfortify::unscale((OakPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((OakPW1.5), center= 147.4742, scale=14.19027)
OakPW.5 <- OakPW.5$var1 

RowanPW1.5 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW2.5 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW.5 <- ggfortify::unscale((RowanPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((RowanPW1.5), center= 147.4742, scale=14.19027)
RowanPW.5 <- RowanPW.5$var1 

SycamorePW1.5 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW2.5 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW.5 <- ggfortify::unscale((SycamorePW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((SycamorePW1.5), center= 147.4742, scale=14.19027)
SycamorePW.5 <- SycamorePW.5$var1 

WillowPW1.5 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW2.5 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW.5 <- ggfortify::unscale((WillowPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((WillowPW1.5), center= 147.4742, scale=14.19027)
WillowPW.5 <- WillowPW.5$var1 

## Peak width diff to mean
AlderPW.5dif <- AlderPW.5-MeanPW.5
AshPW.5dif <- AshPW.5-MeanPW.5
BeechPW.5dif <- BeechPW.5-MeanPW.5
BirchPW.5dif <- BirchPW.5-MeanPW.5
ElmPW.5dif <- ElmPW.5-MeanPW.5
HazelPW.5dif <- HazelPW.5-MeanPW.5
OakPW.5dif <- OakPW.5-MeanPW.5
RowanPW.5dif <- RowanPW.5-MeanPW.5
SycamorePW.5dif <- SycamorePW.5-MeanPW.5
WillowPW.5dif <- WillowPW.5-MeanPW.5

## Peak width diff to oak
AlderPW.5OD <- AlderPW.5-OakPW.5
AshPW.5OD <- AshPW.5-OakPW.5
BeechPW.5OD <- BeechPW.5-OakPW.5
BirchPW.5OD <- BirchPW.5-OakPW.5
ElmPW.5OD <- ElmPW.5-OakPW.5
HazelPW.5OD <- HazelPW.5-OakPW.5
RowanPW.5OD <- RowanPW.5-OakPW.5
SycamorePW.5OD <- SycamorePW.5-OakPW.5
WillowPW.5OD <- WillowPW.5-OakPW.5


#### Data frame of PDHW mean and CIs ####   

AbundCurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                          PD=c(mean(AlderPD),mean(AshPD),mean(BeechPD),mean(BirchPD),mean(ElmPD),mean(HazelPD),mean(OakPD),mean(RowanPD),mean(SycamorePD),mean(WillowPD)),
                          PDLCI=c(HPDinterval(mcmc(AlderPD))[1], HPDinterval(mcmc(AshPD))[1], HPDinterval(mcmc(BeechPD))[1], HPDinterval(mcmc(BirchPD))[1], HPDinterval(mcmc(ElmPD))[1], HPDinterval(mcmc(HazelPD))[1], HPDinterval(mcmc(OakPD))[1], HPDinterval(mcmc(RowanPD))[1], HPDinterval(mcmc(SycamorePD))[1], HPDinterval(mcmc(WillowPD))[1]),
                          PDUCI=c(HPDinterval(mcmc(AlderPD))[2], HPDinterval(mcmc(AshPD))[2], HPDinterval(mcmc(BeechPD))[2], HPDinterval(mcmc(BirchPD))[2], HPDinterval(mcmc(ElmPD))[2], HPDinterval(mcmc(HazelPD))[2], HPDinterval(mcmc(OakPD))[2], HPDinterval(mcmc(RowanPD))[2], HPDinterval(mcmc(SycamorePD))[2], HPDinterval(mcmc(WillowPD))[2]),
                          PDdif=c(mean(AlderPDdif),mean(AshPDdif),mean(BeechPDdif),mean(BirchPDdif),mean(ElmPDdif),mean(HazelPDdif),mean(OakPDdif),mean(RowanPDdif),mean(SycamorePDdif),mean(WillowPDdif)),
                          PDdifLCI=c(HPDinterval(mcmc(AlderPDdif))[1], HPDinterval(mcmc(AshPDdif))[1], HPDinterval(mcmc(BeechPDdif))[1], HPDinterval(mcmc(BirchPDdif))[1], HPDinterval(mcmc(ElmPDdif))[1], HPDinterval(mcmc(HazelPDdif))[1], HPDinterval(mcmc(OakPDdif))[1], HPDinterval(mcmc(RowanPDdif))[1], HPDinterval(mcmc(SycamorePDdif))[1], HPDinterval(mcmc(WillowPDdif))[1]),
                          PDdifUCI=c(HPDinterval(mcmc(AlderPDdif))[2], HPDinterval(mcmc(AshPDdif))[2], HPDinterval(mcmc(BeechPDdif))[2], HPDinterval(mcmc(BirchPDdif))[2], HPDinterval(mcmc(ElmPDdif))[2], HPDinterval(mcmc(HazelPDdif))[2], HPDinterval(mcmc(OakPDdif))[2], HPDinterval(mcmc(RowanPDdif))[2], HPDinterval(mcmc(SycamorePDdif))[2], HPDinterval(mcmc(WillowPDdif))[2]),
                          PH=c(mean(AlderPH),mean(AshPH),mean(BeechPH),mean(BirchPH),mean(ElmPH),mean(HazelPH),mean(OakPH),mean(RowanPH),mean(SycamorePH),mean(WillowPH)),
                          PHLCI=c(HPDinterval(mcmc(AlderPH))[1], HPDinterval(mcmc(AshPH))[1], HPDinterval(mcmc(BeechPH))[1], HPDinterval(mcmc(BirchPH))[1], HPDinterval(mcmc(ElmPH))[1], HPDinterval(mcmc(HazelPH))[1], HPDinterval(mcmc(OakPH))[1], HPDinterval(mcmc(RowanPH))[1], HPDinterval(mcmc(SycamorePH))[1], HPDinterval(mcmc(WillowPH))[1]),
                          PHUCI=c(HPDinterval(mcmc(AlderPH))[2], HPDinterval(mcmc(AshPH))[2], HPDinterval(mcmc(BeechPH))[2], HPDinterval(mcmc(BirchPH))[2], HPDinterval(mcmc(ElmPH))[2], HPDinterval(mcmc(HazelPH))[2], HPDinterval(mcmc(OakPH))[2], HPDinterval(mcmc(RowanPH))[2], HPDinterval(mcmc(SycamorePH))[2], HPDinterval(mcmc(WillowPH))[2]),
                          PHprop=c(mean(AlderPHprop),mean(AshPHprop),mean(BeechPHprop),mean(BirchPHprop),mean(ElmPHprop),mean(HazelPHprop),mean(OakPHprop),mean(RowanPHprop),mean(SycamorePHprop),mean(WillowPHprop)),
                          PHpropLCI=c(HPDinterval(mcmc(AlderPHprop))[1], HPDinterval(mcmc(AshPHprop))[1], HPDinterval(mcmc(BeechPHprop))[1], HPDinterval(mcmc(BirchPHprop))[1], HPDinterval(mcmc(ElmPHprop))[1], HPDinterval(mcmc(HazelPHprop))[1], HPDinterval(mcmc(OakPHprop))[1], HPDinterval(mcmc(RowanPHprop))[1], HPDinterval(mcmc(SycamorePHprop))[1], HPDinterval(mcmc(WillowPHprop))[1]),
                          PHpropUCI=c(HPDinterval(mcmc(AlderPHprop))[2], HPDinterval(mcmc(AshPHprop))[2], HPDinterval(mcmc(BeechPHprop))[2], HPDinterval(mcmc(BirchPHprop))[2], HPDinterval(mcmc(ElmPHprop))[2], HPDinterval(mcmc(HazelPHprop))[2], HPDinterval(mcmc(OakPHprop))[2], HPDinterval(mcmc(RowanPHprop))[2], HPDinterval(mcmc(SycamorePHprop))[2], HPDinterval(mcmc(WillowPHprop))[2]),
                          PW=c(mean(AlderPW, na.rm=TRUE),mean(AshPW, na.rm=TRUE),mean(BeechPW, na.rm=TRUE),mean(BirchPW, na.rm=TRUE),mean(ElmPW, na.rm=TRUE),mean(HazelPW, na.rm=TRUE),mean(OakPW, na.rm=TRUE),mean(RowanPW, na.rm=TRUE),mean(SycamorePW, na.rm=TRUE),mean(WillowPW, na.rm=TRUE)),
                          PWLCI=c(HPDinterval(mcmc(AlderPW))[1], HPDinterval(mcmc(AshPW))[1], HPDinterval(mcmc(BeechPW))[1], HPDinterval(mcmc(BirchPW))[1], HPDinterval(mcmc(ElmPW))[1], HPDinterval(mcmc(HazelPW))[1], HPDinterval(mcmc(OakPW))[1], HPDinterval(mcmc(RowanPW))[1], HPDinterval(mcmc(SycamorePW))[1], HPDinterval(mcmc(WillowPW))[1]),
                          PWUCI=c(HPDinterval(mcmc(AlderPW))[2], HPDinterval(mcmc(AshPW))[2], HPDinterval(mcmc(BeechPW))[2], HPDinterval(mcmc(BirchPW))[2], HPDinterval(mcmc(ElmPW))[2], HPDinterval(mcmc(HazelPW))[2], HPDinterval(mcmc(OakPW))[2], HPDinterval(mcmc(RowanPW))[2], HPDinterval(mcmc(SycamorePW))[2], HPDinterval(mcmc(WillowPW))[2]),
                          PWdif=c(mean(AlderPWdif, na.rm=TRUE),mean(AshPWdif, na.rm=TRUE),mean(BeechPWdif, na.rm=TRUE),mean(BirchPWdif, na.rm=TRUE),mean(ElmPWdif, na.rm=TRUE),mean(HazelPWdif, na.rm=TRUE),mean(OakPWdif, na.rm=TRUE),mean(RowanPWdif, na.rm=TRUE),mean(SycamorePWdif, na.rm=TRUE),mean(WillowPWdif, na.rm=TRUE)),
                          PWdifLCI=c(HPDinterval(mcmc(AlderPWdif))[1], HPDinterval(mcmc(AshPWdif))[1], HPDinterval(mcmc(BeechPWdif))[1], HPDinterval(mcmc(BirchPWdif))[1], HPDinterval(mcmc(ElmPWdif))[1], HPDinterval(mcmc(HazelPWdif))[1], HPDinterval(mcmc(OakPWdif))[1], HPDinterval(mcmc(RowanPWdif))[1], HPDinterval(mcmc(SycamorePWdif))[1], HPDinterval(mcmc(WillowPWdif))[1]),
                          PWdifUCI=c(HPDinterval(mcmc(AlderPWdif))[2], HPDinterval(mcmc(AshPWdif))[2], HPDinterval(mcmc(BeechPWdif))[2], HPDinterval(mcmc(BirchPWdif))[2], HPDinterval(mcmc(ElmPWdif))[2], HPDinterval(mcmc(HazelPWdif))[2], HPDinterval(mcmc(OakPWdif))[2], HPDinterval(mcmc(RowanPWdif))[2], HPDinterval(mcmc(SycamorePWdif))[2], HPDinterval(mcmc(WillowPWdif))[2]),
                          PW.5=c(mean(AlderPW.5, na.rm=TRUE),mean(AshPW.5, na.rm=TRUE),mean(BeechPW.5, na.rm=TRUE),mean(BirchPW.5, na.rm=TRUE),mean(ElmPW.5, na.rm=TRUE),mean(HazelPW.5, na.rm=TRUE),mean(OakPW.5, na.rm=TRUE),mean(RowanPW.5, na.rm=TRUE),mean(SycamorePW.5, na.rm=TRUE),mean(WillowPW.5, na.rm=TRUE)),
                          PW.5LCI=c(HPDinterval(mcmc(AlderPW.5))[1], HPDinterval(mcmc(AshPW.5))[1], HPDinterval(mcmc(BeechPW.5))[1], HPDinterval(mcmc(BirchPW.5))[1], HPDinterval(mcmc(ElmPW.5))[1], HPDinterval(mcmc(HazelPW.5))[1], HPDinterval(mcmc(OakPW.5))[1], HPDinterval(mcmc(RowanPW.5))[1], HPDinterval(mcmc(SycamorePW.5))[1], HPDinterval(mcmc(WillowPW.5))[1]),
                          PW.5UCI=c(HPDinterval(mcmc(AlderPW.5))[2], HPDinterval(mcmc(AshPW.5))[2], HPDinterval(mcmc(BeechPW.5))[2], HPDinterval(mcmc(BirchPW.5))[2], HPDinterval(mcmc(ElmPW.5))[2], HPDinterval(mcmc(HazelPW.5))[2], HPDinterval(mcmc(OakPW.5))[2], HPDinterval(mcmc(RowanPW.5))[2], HPDinterval(mcmc(SycamorePW.5))[2], HPDinterval(mcmc(WillowPW.5))[2]),
                          PW.5dif=c(mean(AlderPW.5dif, na.rm=TRUE),mean(AshPW.5dif, na.rm=TRUE),mean(BeechPW.5dif, na.rm=TRUE),mean(BirchPW.5dif, na.rm=TRUE),mean(ElmPW.5dif, na.rm=TRUE),mean(HazelPW.5dif, na.rm=TRUE),mean(OakPW.5dif, na.rm=TRUE),mean(RowanPW.5dif, na.rm=TRUE),mean(SycamorePW.5dif, na.rm=TRUE),mean(WillowPW.5dif, na.rm=TRUE)),
                          PW.5difLCI=c(HPDinterval(mcmc(AlderPW.5dif))[1], HPDinterval(mcmc(AshPW.5dif))[1], HPDinterval(mcmc(BeechPW.5dif))[1], HPDinterval(mcmc(BirchPW.5dif))[1], HPDinterval(mcmc(ElmPW.5dif))[1], HPDinterval(mcmc(HazelPW.5dif))[1], HPDinterval(mcmc(OakPW.5dif))[1], HPDinterval(mcmc(RowanPW.5dif))[1], HPDinterval(mcmc(SycamorePW.5dif))[1], HPDinterval(mcmc(WillowPW.5dif))[1]),
                          PW.5difUCI=c(HPDinterval(mcmc(AlderPW.5dif))[2], HPDinterval(mcmc(AshPW.5dif))[2], HPDinterval(mcmc(BeechPW.5dif))[2], HPDinterval(mcmc(BirchPW.5dif))[2], HPDinterval(mcmc(ElmPW.5dif))[2], HPDinterval(mcmc(HazelPW.5dif))[2], HPDinterval(mcmc(OakPW.5dif))[2], HPDinterval(mcmc(RowanPW.5dif))[2], HPDinterval(mcmc(SycamorePW.5dif))[2], HPDinterval(mcmc(WillowPW.5dif))[2]))


#### Plotting curves in ggplot ####
Curves <- data.frame(date=days, "Fixed Effects"=MeanCurve, Alder=AlderCurve, Ash=AshCurve, Beech=BeechCurve, Birch=BirchCurve, Elm=ElmCurve, Hazel=HazelCurve, Oak=OakCurve, Rowan=RowanCurve, Sycamore=SycamoreCurve, Willow=WillowCurve)
Curveslong <- gather(Curves, key="TreeTaxa", value="Abundance", 2:12)
AllCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "black","darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

# Curves with legend
ggplot(Curveslong, aes(date, Abundance, col=TreeTaxa))+ #saved as 8"x9"  5x7forBES
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))

#Curves without legend
AbundCurvesplot <- ggplot(Curveslong, aes(date, Abundance, col=TreeTaxa))+ 
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  guides(color = "none", linetype="none")+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))

# Peak Date difference from mean
(PDdifplot <- ggplot(AbundCurves, aes(fct_rev(TT), (PDdif), colour=TT))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymax=(PDdifUCI), ymin=(PDdifLCI), width=0.5))+
    theme_bw()+
    xlab("")+
    ylab("Peak Timing (days)")+
    coord_flip()+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    guides(color = "none")+
    scale_colour_manual(values=AllTaxaCols)) 

# Peak Height prop difference from mean
(PHpropdifplot <- ggplot(AbundCurves, aes(fct_rev(TT), PHprop, colour=TT))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
    theme_bw()+
    xlab("")+
    ylab("Height (proportional)")+
    geom_hline(yintercept=1, linetype="dashed", color = "black")+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    guides(color = "none")+
    coord_flip()+
    scale_colour_manual(values=AllTaxaCols))

# Peak Width.5 difference from mean
PW.5difplot <- ggplot(AbundCurves, aes(fct_rev(TT), PW.5dif, colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=PW.5difUCI, ymin=PW.5difLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Relative Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)


library(gridExtra)

col5 <- grid.arrange(AbundCurvesplot, nrow = 1, widths = 1)
col7 <- grid.arrange(PDdifplot, PHpropdifplot, PW.5difplot, nrow = 3, heights = c(1,1,1))
AbundCurvesFig2 <- grid.arrange(col5, col7,  ncol = 2, widths = c(4,2)) #saved as 8"x10"

#### Supp. mat. figures ####

# Peak Width (set height) difference from mean
ggplot(AbundCurves, aes(fct_rev(TT), PWdif))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWdifUCI, ymin=PWdifLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Duration (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+#+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")# saved as 6"x4"



## Oak Difference PDHW 
AbundODcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                            PD=c(mean(AlderPDOD),mean(AshPDOD),mean(BeechPDOD),mean(BirchPDOD),mean(ElmPDOD),mean(HazelPDOD),mean(RowanPDOD),mean(SycamorePDOD),mean(WillowPDOD)),
                            PDLCI=c(HPDinterval(mcmc(AlderPDOD))[1], HPDinterval(mcmc(AshPDOD))[1], HPDinterval(mcmc(BeechPDOD))[1], HPDinterval(mcmc(BirchPDOD))[1], HPDinterval(mcmc(ElmPDOD))[1], HPDinterval(mcmc(HazelPDOD))[1], HPDinterval(mcmc(RowanPDOD))[1], HPDinterval(mcmc(SycamorePDOD))[1], HPDinterval(mcmc(WillowPDOD))[1]),
                            PDUCI=c(HPDinterval(mcmc(AlderPDOD))[2], HPDinterval(mcmc(AshPDOD))[2], HPDinterval(mcmc(BeechPDOD))[2], HPDinterval(mcmc(BirchPDOD))[2], HPDinterval(mcmc(ElmPDOD))[2], HPDinterval(mcmc(HazelPDOD))[2], HPDinterval(mcmc(RowanPDOD))[2], HPDinterval(mcmc(SycamorePDOD))[2], HPDinterval(mcmc(WillowPDOD))[2]),
                            PHprop=c(mean(AlderPHODprop),mean(AshPHODprop),mean(BeechPHODprop),mean(BirchPHODprop),mean(ElmPHODprop),mean(HazelPHODprop),mean(RowanPHODprop),mean(SycamorePHODprop),mean(WillowPHODprop)),
                            PHpropLCI=c(HPDinterval(mcmc(AlderPHODprop))[1], HPDinterval(mcmc(AshPHODprop))[1], HPDinterval(mcmc(BeechPHODprop))[1], HPDinterval(mcmc(BirchPHODprop))[1], HPDinterval(mcmc(ElmPHODprop))[1], HPDinterval(mcmc(HazelPHODprop))[1], HPDinterval(mcmc(RowanPHODprop))[1], HPDinterval(mcmc(SycamorePHODprop))[1], HPDinterval(mcmc(WillowPHODprop))[1]),
                            PHpropUCI=c(HPDinterval(mcmc(AlderPHODprop))[2], HPDinterval(mcmc(AshPHODprop))[2], HPDinterval(mcmc(BeechPHODprop))[2], HPDinterval(mcmc(BirchPHODprop))[2], HPDinterval(mcmc(ElmPHODprop))[2], HPDinterval(mcmc(HazelPHODprop))[2], HPDinterval(mcmc(RowanPHODprop))[2], HPDinterval(mcmc(SycamorePHODprop))[2], HPDinterval(mcmc(WillowPHODprop))[2]),
                            PW=c(mean(AlderPWOD, na.rm=TRUE),mean(AshPWOD, na.rm=TRUE),mean(BeechPWOD, na.rm=TRUE),mean(BirchPWOD, na.rm=TRUE),mean(ElmPWOD, na.rm=TRUE),mean(HazelPWOD, na.rm=TRUE),mean(RowanPWOD, na.rm=TRUE),mean(SycamorePWOD, na.rm=TRUE),mean(WillowPWOD, na.rm=TRUE)),
                            PWLCI=c(HPDinterval(mcmc(AlderPWOD))[1], HPDinterval(mcmc(AshPWOD))[1], HPDinterval(mcmc(BeechPWOD))[1], HPDinterval(mcmc(BirchPWOD))[1], HPDinterval(mcmc(ElmPWOD))[1], HPDinterval(mcmc(HazelPWOD))[1], HPDinterval(mcmc(RowanPWOD))[1], HPDinterval(mcmc(SycamorePWOD))[1], HPDinterval(mcmc(WillowPWOD))[1]),
                            PWUCI=c(HPDinterval(mcmc(AlderPWOD))[2], HPDinterval(mcmc(AshPWOD))[2], HPDinterval(mcmc(BeechPWOD))[2], HPDinterval(mcmc(BirchPWOD))[2], HPDinterval(mcmc(ElmPWOD))[2], HPDinterval(mcmc(HazelPWOD))[2], HPDinterval(mcmc(RowanPWOD))[2], HPDinterval(mcmc(SycamorePWOD))[2], HPDinterval(mcmc(WillowPWOD))[2]),
                            PW.5=c(mean(AlderPW.5OD, na.rm=TRUE),mean(AshPW.5OD, na.rm=TRUE),mean(BeechPW.5OD, na.rm=TRUE),mean(BirchPW.5OD, na.rm=TRUE),mean(ElmPW.5OD, na.rm=TRUE),mean(HazelPW.5OD, na.rm=TRUE),mean(RowanPW.5OD, na.rm=TRUE),mean(SycamorePW.5OD, na.rm=TRUE),mean(WillowPW.5OD, na.rm=TRUE)),
                            PW.5LCI=c(HPDinterval(mcmc(AlderPW.5OD))[1], HPDinterval(mcmc(AshPW.5OD))[1], HPDinterval(mcmc(BeechPW.5OD))[1], HPDinterval(mcmc(BirchPW.5OD))[1], HPDinterval(mcmc(ElmPW.5OD))[1], HPDinterval(mcmc(HazelPW.5OD))[1], HPDinterval(mcmc(RowanPW.5OD))[1], HPDinterval(mcmc(SycamorePW.5OD))[1], HPDinterval(mcmc(WillowPW.5OD))[1]),
                            PW.5UCI=c(HPDinterval(mcmc(AlderPW.5OD))[2], HPDinterval(mcmc(AshPW.5OD))[2], HPDinterval(mcmc(BeechPW.5OD))[2], HPDinterval(mcmc(BirchPW.5OD))[2], HPDinterval(mcmc(ElmPW.5OD))[2], HPDinterval(mcmc(HazelPW.5OD))[2], HPDinterval(mcmc(RowanPW.5OD))[2], HPDinterval(mcmc(SycamorePW.5OD))[2], HPDinterval(mcmc(WillowPW.5OD))[2]))

NoOakCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "royalblue4", "slateblue2", "orchid")

PDODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PD))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PDUCI, ymin=PDLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Timing (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PHpropODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PHprop))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (proportional)")+
  coord_flip()+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PW))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWUCI, ymin=PWLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Duration (days) ")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWOD.5plot <- ggplot(AbundODcurves, aes(fct_rev(TT), PW.5))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PW.5UCI, ymin=PW.5LCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))


AbundODplots <- grid.arrange(PDODplot, PHpropODplot, PWOD.5plot, PWODplot, ncol = 4, widths=c(1,1,1,1)) #saved 6" by 13"


#### Posterior distribution sample size (without NaNs) ####
AbundCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                                    PD=c(length(which(is.na(AlderPD)==FALSE)),length(which(is.na(AshPD)==FALSE)),length(which(is.na(BeechPD)==FALSE)),length(which(is.na(BirchPD)==FALSE)),length(which(is.na(ElmPD)==FALSE)),length(which(is.na(HazelPD)==FALSE)),length(which(is.na(OakPD)==FALSE)),length(which(is.na(RowanPD)==FALSE)),length(which(is.na(SycamorePD)==FALSE)),length(which(is.na(WillowPD)==FALSE))),
                                    PDdif=c(length(which(is.na(AlderPDdif)==FALSE)),length(which(is.na(AshPDdif)==FALSE)),length(which(is.na(BeechPDdif)==FALSE)),length(which(is.na(BirchPDdif)==FALSE)),length(which(is.na(ElmPDdif)==FALSE)),length(which(is.na(HazelPDdif)==FALSE)),length(which(is.na(OakPDdif)==FALSE)),length(which(is.na(RowanPDdif)==FALSE)),length(which(is.na(SycamorePDdif)==FALSE)),length(which(is.na(WillowPDdif)==FALSE))),
                                    PH=c(length(which(is.na(AlderPH)==FALSE)),length(which(is.na(AshPH)==FALSE)),length(which(is.na(BeechPH)==FALSE)),length(which(is.na(BirchPH)==FALSE)),length(which(is.na(ElmPH)==FALSE)),length(which(is.na(HazelPH)==FALSE)),length(which(is.na(OakPH)==FALSE)),length(which(is.na(RowanPH)==FALSE)),length(which(is.na(SycamorePH)==FALSE)),length(which(is.na(WillowPH)==FALSE))),
                                    PHprop=c(length(which(is.na(AlderPHprop)==FALSE)),length(which(is.na(AshPHprop)==FALSE)),length(which(is.na(BeechPHprop)==FALSE)),length(which(is.na(BirchPHprop)==FALSE)),length(which(is.na(ElmPHprop)==FALSE)),length(which(is.na(HazelPHprop)==FALSE)),length(which(is.na(OakPHprop)==FALSE)),length(which(is.na(RowanPHprop)==FALSE)),length(which(is.na(SycamorePHprop)==FALSE)),length(which(is.na(WillowPHprop)==FALSE))),
                                    PW=c(length(which(is.na(AlderPW)==FALSE)),length(which(is.na(AshPW)==FALSE)),length(which(is.na(BeechPW)==FALSE)),length(which(is.na(BirchPW)==FALSE)),length(which(is.na(ElmPW)==FALSE)),length(which(is.na(HazelPW)==FALSE)),length(which(is.na(OakPW)==FALSE)),length(which(is.na(RowanPW)==FALSE)),length(which(is.na(SycamorePW)==FALSE)),length(which(is.na(WillowPW)==FALSE))),
                                    PWdif=c(length(which(is.na(AlderPWdif)==FALSE)),length(which(is.na(AshPWdif)==FALSE)),length(which(is.na(BeechPWdif)==FALSE)),length(which(is.na(BirchPWdif)==FALSE)),length(which(is.na(ElmPWdif)==FALSE)),length(which(is.na(HazelPWdif)==FALSE)),length(which(is.na(OakPWdif)==FALSE)),length(which(is.na(RowanPWdif)==FALSE)),length(which(is.na(SycamorePWdif)==FALSE)),length(which(is.na(WillowPWdif)==FALSE))),
                                    PW.5=c(length(which(is.na(AlderPW.5)==FALSE)),length(which(is.na(AshPW.5)==FALSE)),length(which(is.na(BeechPW.5)==FALSE)),length(which(is.na(BirchPW.5)==FALSE)),length(which(is.na(ElmPW.5)==FALSE)),length(which(is.na(HazelPW.5)==FALSE)),length(which(is.na(OakPW.5)==FALSE)),length(which(is.na(RowanPW.5)==FALSE)),length(which(is.na(SycamorePW.5)==FALSE)),length(which(is.na(WillowPW.5)==FALSE))),
                                    PW.5dif=c(length(which(is.na(AlderPW.5dif)==FALSE)),length(which(is.na(AshPW.5dif)==FALSE)),length(which(is.na(BeechPW.5dif)==FALSE)),length(which(is.na(BirchPW.5dif)==FALSE)),length(which(is.na(ElmPW.5dif)==FALSE)),length(which(is.na(HazelPW.5dif)==FALSE)),length(which(is.na(OakPW.5dif)==FALSE)),length(which(is.na(RowanPW.5dif)==FALSE)),length(which(is.na(SycamorePW.5dif)==FALSE)),length(which(is.na(WillowPW.5dif)==FALSE))))

AbundODCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                                      PD=c(length(which(is.na(AlderPDOD)==FALSE)),length(which(is.na(AshPDOD)==FALSE)),length(which(is.na(BeechPDOD)==FALSE)),length(which(is.na(BirchPDOD)==FALSE)),length(which(is.na(ElmPDOD)==FALSE)),length(which(is.na(HazelPDOD)==FALSE)),length(which(is.na(RowanPDOD)==FALSE)),length(which(is.na(SycamorePDOD)==FALSE)),length(which(is.na(WillowPDOD)==FALSE))),
                                      PHprop=c(length(which(is.na(AlderPHODprop)==FALSE)),length(which(is.na(AshPHODprop)==FALSE)),length(which(is.na(BeechPHODprop)==FALSE)),length(which(is.na(BirchPHODprop)==FALSE)),length(which(is.na(ElmPHODprop)==FALSE)),length(which(is.na(HazelPHODprop)==FALSE)),length(which(is.na(RowanPHODprop)==FALSE)),length(which(is.na(SycamorePHODprop)==FALSE)),length(which(is.na(WillowPHODprop)==FALSE))),
                                      PW=c(length(which(is.na(AlderPWOD)==FALSE)),length(which(is.na(AshPWOD)==FALSE)),length(which(is.na(BeechPWOD)==FALSE)),length(which(is.na(BirchPWOD)==FALSE)),length(which(is.na(ElmPWOD)==FALSE)),length(which(is.na(HazelPWOD)==FALSE)),length(which(is.na(RowanPWOD)==FALSE)),length(which(is.na(SycamorePWOD)==FALSE)),length(which(is.na(WillowPWOD)==FALSE))),
                                      PW.5=c(length(which(is.na(AlderPW.5OD)==FALSE)),length(which(is.na(AshPW.5OD)==FALSE)),length(which(is.na(BeechPW.5OD)==FALSE)),length(which(is.na(BirchPW.5OD)==FALSE)),length(which(is.na(ElmPW.5OD)==FALSE)),length(which(is.na(HazelPW.5OD)==FALSE)),length(which(is.na(RowanPW.5OD)==FALSE)),length(which(is.na(SycamorePW.5OD)==FALSE)),length(which(is.na(WillowPW.5OD)==FALSE))))

meanmetrics <- data.frame(metric=c("PD", "PH", "PW.5", "PW.01"), 
                          mean=c(mean(MeanPD), mean(MeanPH),mean(MeanPW.5,na.rm=TRUE), mean(MeanPW,na.rm=TRUE)), 
                          lci=c(HPDinterval(mcmc(MeanPD))[1], HPDinterval(mcmc(MeanPH))[1],HPDinterval(mcmc(MeanPW.5))[1], HPDinterval(mcmc(MeanPW))[1]),
                          uci=c(HPDinterval(mcmc(MeanPD))[2], HPDinterval(mcmc(MeanPH))[2],HPDinterval(mcmc(MeanPW.5))[2], HPDinterval(mcmc(MeanPW))[2]),
                          SS=c(length(which(is.na(MeanPD)==FALSE)),length(which(is.na(MeanPH)==FALSE)),length(which(is.na(MeanPW.5)==FALSE)),length(which(is.na(MeanPW)==FALSE))))



############################
#### Model output table ####  
############################


#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(ATTC23$Sol[,1]),3)," (",
                      round(HPDinterval(ATTC23$Sol[,1])[1],3)," - ",
                      round(HPDinterval(ATTC23$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(ATTC23$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(ATTC23$Sol[,2]),3)," (",
                          round(HPDinterval(ATTC23$Sol[,2])[1],3)," - ",
                          round(HPDinterval(ATTC23$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(ATTC23$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(ATTC23$Sol[,3]),3)," (",
                           round(HPDinterval(ATTC23$Sol[,3])[1],3)," - ",
                           round(HPDinterval(ATTC23$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(ATTC23$Sol[,3]))))

#### random terms
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                             round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                          round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-3
treetaxa3<-c("TreeTaxa- Intercept:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                           round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-5
treetaxa5<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                              round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-6
treetaxa6<-c("TreeTaxa- Date slope:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                            round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-9
treetaxa9<-c("TreeTaxa- Date² slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                               round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTC23$VCV[, column])))

column<-10
site10<-c("Site- Intercept var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                      round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-11
site11<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                   round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                   round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-12
site12<-c("Site- Intercept:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                    round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                    round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-14
site14<-c("Site- Date slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                       round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-15
site15<-c("Site- Date slope:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                     round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                     round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-18
site18<-c("Site- Date² slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                        round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                        round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-19
year10<-c("Year- Intercept var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                      round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-20
year11<-c("Year- Intercept:Date slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                   round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                   round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-21
year12<-c("Year- Intercept:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                    round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                    round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-23
year14<-c("Year- Date slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                       round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-24
year15<-c("Year- Date slope:Date² slope covar",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                                     round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                                     round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-27
year18<-c("Year- Date² slope var",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                                        round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                                        round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-31
recorder<-c("Recorder",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                             round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTC23$VCV[, column])))


column<-30
siteday<-c("Site-Day",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                            round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                            round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(ATTC23$VCV[, column])))

column<-28
siteyear<-c("Site-Year",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                              round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                              round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTC23$VCV[, column])))

column<-29
treeID<-c("Tree ID",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                          round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                          round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTC23$VCV[, column])))

column<-32
residual<-c("Residual",paste(round(posterior.mode(ATTC23$VCV[, column]),3)," (",
                             round(HPDinterval(ATTC23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTC23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTC23$VCV[, column])))

random<-rbind(treetaxa1,treetaxa2,treetaxa3,treetaxa5,treetaxa6,treetaxa9,site10,site11,site12,site14,site15,site18,year10,year11,year12,year14,year15,year18, siteyear, recorder, siteday, treeID, residual)

#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/TableS4_ATTC2.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
