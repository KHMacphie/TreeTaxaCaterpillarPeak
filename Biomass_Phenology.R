rm(list=ls())
library(MCMCglmm)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(gridExtra)

##################################################################
#### Model: Phenological distribution of biomass among tree taxa ####
##################################################################

cater_habitat <- read.csv("~/cater_habitat_data.csv")
cater_habitat_1723 <- subset(cater_habitat, year!="2014")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2015")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2016")

load("~/Mass23_1.RData")  
load("~/Mass23_2.RData")
load("~/Mass23_3.RData")

Sol <- rbind(Mass23_1$Sol,Mass23_2$Sol,Mass23_3$Sol)
dim(Sol)

VCV <- rbind(Mass23_1$VCV,Mass23_2$VCV,Mass23_3$VCV)
dim(VCV)


Mass23 <- list(Sol= Sol,VCV=VCV)
rm(Mass23_1,Mass23_2,Mass23_3, Sol,VCV)

load("~/ATTC23.RData")  

#########################################################
#### Combining abundance and mass to predict biomass ####
#########################################################

#### Extract required columns ####
sampleMass <- data.frame(Mass23$Sol[,1:23]) 
sampleATTC <- data.frame(ATTC23$Sol[,1:33]) 
rm(ATTC23,Mass23)

dim(sampleATTC)
dim(sampleMass)
#unscale date
dayscal <- seq(min(cater_habitat$datescaled),max(cater_habitat$datescaled),0.001)
days <- ggfortify::unscale(dayscal, center= 147.4742, scale=14.19027)$V1

#### Combining PD for each model to get biomass date^2, date and intercept coefficients for each tree taxon
AlderC <- (sampleMass[,1] + sampleMass[,4] + sampleATTC[,1] + sampleATTC[,4])
AshC <- (sampleMass[,1] + sampleMass[,5] + sampleATTC[,1] + sampleATTC[,5])
BeechC <- (sampleMass[,1] + sampleMass[,6] + sampleATTC[,1] + sampleATTC[,6])
BirchC <- (sampleMass[,1] + sampleMass[,7] + sampleATTC[,1] + sampleATTC[,7])
ElmC <- (sampleMass[,1] + sampleMass[,8] + sampleATTC[,1] + sampleATTC[,8])
HazelC <- (sampleMass[,1] + sampleMass[,9] + sampleATTC[,1] + sampleATTC[,9])
OakC <- (sampleMass[,1] + sampleMass[,10] + sampleATTC[,1] + sampleATTC[,10])
RowanC <- (sampleMass[,1] + sampleMass[,11] + sampleATTC[,1] + sampleATTC[,11])
SycamoreC <- (sampleMass[,1] + sampleMass[,12] + sampleATTC[,1] + sampleATTC[,12])
WillowC <- (sampleMass[,1] + sampleMass[,13] + sampleATTC[,1] + sampleATTC[,13])

AlderB <- (sampleMass[,2] + sampleMass[,14] + sampleATTC[,2] + sampleATTC[,14])
AshB <- (sampleMass[,2] + sampleMass[,15] + sampleATTC[,2] + sampleATTC[,15])
BeechB <- (sampleMass[,2] + sampleMass[,16] + sampleATTC[,2] + sampleATTC[,16])
BirchB <- (sampleMass[,2] + sampleMass[,17] + sampleATTC[,2] + sampleATTC[,17])
ElmB <- (sampleMass[,2] + sampleMass[,18] + sampleATTC[,2] + sampleATTC[,18])
HazelB <- (sampleMass[,2] + sampleMass[,19] + sampleATTC[,2] + sampleATTC[,19])
OakB <- (sampleMass[,2] + sampleMass[,20] + sampleATTC[,2] + sampleATTC[,20])
RowanB <- (sampleMass[,2] + sampleMass[,21] + sampleATTC[,2] + sampleATTC[,21])
SycamoreB <- (sampleMass[,2] + sampleMass[,22] + sampleATTC[,2] + sampleATTC[,22])
WillowB <- (sampleMass[,2] + sampleMass[,23] + sampleATTC[,2] + sampleATTC[,23])

AlderA <- (sampleATTC[,3] + sampleATTC[,24] + sampleMass[,3])
AshA <- (sampleATTC[,3] + sampleATTC[,25] + sampleMass[,3])
BeechA <- (sampleATTC[,3] + sampleATTC[,26] + sampleMass[,3])
BirchA <- (sampleATTC[,3] + sampleATTC[,27] + sampleMass[,3])
ElmA <- (sampleATTC[,3] + sampleATTC[,28] + sampleMass[,3])
HazelA <- (sampleATTC[,3] + sampleATTC[,29] + sampleMass[,3])
OakA <- (sampleATTC[,3] + sampleATTC[,30] + sampleMass[,3])
RowanA <- (sampleATTC[,3] + sampleATTC[,31] + sampleMass[,3])
SycamoreA <- (sampleATTC[,3] + sampleATTC[,32] + sampleMass[,3])
WillowA <- (sampleATTC[,3] + sampleATTC[,33] + sampleMass[,3])

# Fixed effect coefficients 
MeanA <- (sampleATTC[,3] + sampleMass[,3])
MeanB <- (sampleATTC[,2] + sampleMass[,2])
MeanC <- (sampleATTC[,1] + sampleMass[,1])

# tree taxon curves
AlderBM2 <- exp(mean(AlderC) +mean(AlderB)*dayscal + mean(AlderA)*dayscal^2)
AshBM2 <- exp(mean(AshC) +mean(AshB)*dayscal + mean(AshA)*dayscal^2)
BeechBM2 <- exp(mean(BeechC) +mean(BeechB)*dayscal + mean(BeechA)*dayscal^2)
BirchBM2 <- exp(mean(BirchC) +mean(BirchB)*dayscal + mean(BirchA)*dayscal^2)
ElmBM2 <- exp(mean(ElmC) +mean(ElmB)*dayscal + mean(ElmA)*dayscal^2)
HazelBM2 <- exp(mean(HazelC) +mean(HazelB)*dayscal + mean(HazelA)*dayscal^2)
OakBM2 <- exp(mean(OakC) +mean(OakB)*dayscal + mean(OakA)*dayscal^2)
RowanBM2 <- exp(mean(RowanC) +mean(RowanB)*dayscal + mean(RowanA)*dayscal^2)
SycamoreBM2 <- exp(mean(SycamoreC) +mean(SycamoreB)*dayscal + mean(SycamoreA)*dayscal^2)
WillowBM2 <- exp(mean(WillowC) +mean(WillowB)*dayscal + mean(WillowA)*dayscal^2)
MeanBM2 <- exp(mean(MeanC) +mean(MeanB)*dayscal + mean(MeanA)*dayscal^2)

#### Plotting curves in ggplot ####
Curves <- data.frame(date=days, "Fixed Effects"=MeanBM2, Alder=AlderBM2, Ash=AshBM2, Beech=BeechBM2, Birch=BirchBM2, Elm=ElmBM2, Hazel=HazelBM2, Oak=OakBM2, Rowan=RowanBM2, Sycamore=SycamoreBM2, Willow=WillowBM2)
Curveslong <- gather(Curves, key="TreeTaxa", value="Biomass", 2:12)
AllCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "black","darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

# Curves with legend
ggplot(Curveslong, aes(date, Biomass, col=TreeTaxa))+ #saved as 8"x9"
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Biomass (g)")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  #guides(color = "none")+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))

#Curves without legend
BiomassCurvesplot <- ggplot(Curveslong, aes(date, Biomass, col=TreeTaxa))+ #saved as 8"x9"
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Biomass (g)")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  guides(color = "none", linetype="none")+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))


#### Peak date per taxa ####   -b/2a
MeanPD <- ggfortify::unscale((-MeanB/(2*MeanA)), center= 147.4742, scale=14.19027)$V1
AlderPD <-  ggfortify::unscale((-AlderB/(2*AlderA)), center= 147.4742, scale=14.19027)$V1
AshPD <-  ggfortify::unscale((-AshB/(2*AshA)), center= 147.4742, scale=14.19027)$V1
BeechPD <-  ggfortify::unscale((-BeechB/(2*BeechA)), center= 147.4742, scale=14.19027)$V1
BirchPD <-  ggfortify::unscale((-BirchB/(2*BirchA)), center= 147.4742, scale=14.19027)$V1
ElmPD <-  ggfortify::unscale((-ElmB/(2*ElmA)), center= 147.4742, scale=14.19027)$V1
HazelPD <-  ggfortify::unscale((-HazelB/(2*HazelA)), center= 147.4742, scale=14.19027)$V1
OakPD <-  ggfortify::unscale((-OakB/(2*OakA)), center= 147.4742, scale=14.19027)$V1
RowanPD <-  ggfortify::unscale((-RowanB/(2*RowanA)), center= 147.4742, scale=14.19027)$V1
SycamorePD <-  ggfortify::unscale((-SycamoreB/(2*SycamoreA)), center= 147.4742, scale=14.19027)$V1
WillowPD <-  ggfortify::unscale((-WillowB/(2*WillowA)), center= 147.4742, scale=14.19027)$V1


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

## Peak height prop of mean
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

## Peak height prop of oak
AlderPHpropOD <- AlderPH/OakPH
AshPHpropOD <- AshPH/OakPH
BeechPHpropOD <- BeechPH/OakPH
BirchPHpropOD <- BirchPH/OakPH
ElmPHpropOD <- ElmPH/OakPH
HazelPHpropOD <- HazelPH/OakPH
RowanPHpropOD <- RowanPH/OakPH
SycamorePHpropOD <- SycamorePH/OakPH
WillowPHpropOD <- WillowPH/OakPH

#### Peak width per taxa (set height) ####   x= (-b +/- sqrt(b^2 - 4ac))/2a

# set at 0.00025g - half the height of the lowest curve (Alder)

MeanPW1 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.00025)))))/(2*MeanA)
MeanPW2 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.00025)))))/(2*MeanA)
MeanPW <- ggfortify::unscale((MeanPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((MeanPW1), center= 147.4742, scale=14.19027)
MeanPW <- MeanPW$V1 

AlderPW1 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.00025)))))/(2*AlderA)
AlderPW2 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.00025)))))/(2*AlderA)
AlderPW <- ggfortify::unscale((AlderPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((AlderPW1), center= 147.4742, scale=14.19027)
AlderPW <- AlderPW$V1 

AshPW1 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(0.00025)))))/(2*AshA)
AshPW2 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(0.00025)))))/(2*AshA)
AshPW <- ggfortify::unscale((AshPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((AshPW1), center= 147.4742, scale=14.19027)
AshPW <- AshPW$V1 

BeechPW1 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.00025)))))/(2*BeechA)
BeechPW2 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.00025)))))/(2*BeechA)
BeechPW <- ggfortify::unscale((BeechPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((BeechPW1), center= 147.4742, scale=14.19027)
BeechPW <- BeechPW$V1 

BirchPW1 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.00025)))))/(2*BirchA)
BirchPW2 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.00025)))))/(2*BirchA)
BirchPW <- ggfortify::unscale((BirchPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((BirchPW1), center= 147.4742, scale=14.19027)
BirchPW <- BirchPW$V1 

ElmPW1 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.00025)))))/(2*ElmA)
ElmPW2 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.00025)))))/(2*ElmA)
ElmPW <- ggfortify::unscale((ElmPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((ElmPW1), center= 147.4742, scale=14.19027)
ElmPW <- ElmPW$V1 

HazelPW1 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.00025)))))/(2*HazelA)
HazelPW2 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.00025)))))/(2*HazelA)
HazelPW <- ggfortify::unscale((HazelPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((HazelPW1), center= 147.4742, scale=14.19027)
HazelPW <- HazelPW$V1 

OakPW1 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(0.00025)))))/(2*OakA)
OakPW2 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(0.00025)))))/(2*OakA)
OakPW <- ggfortify::unscale((OakPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((OakPW1), center= 147.4742, scale=14.19027)
OakPW <- OakPW$V1 

RowanPW1 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.00025)))))/(2*RowanA)
RowanPW2 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.00025)))))/(2*RowanA)
RowanPW <- ggfortify::unscale((RowanPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((RowanPW1), center= 147.4742, scale=14.19027)
RowanPW <- RowanPW$V1 

SycamorePW1 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.00025)))))/(2*SycamoreA)
SycamorePW2 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.00025)))))/(2*SycamoreA)
SycamorePW <- ggfortify::unscale((SycamorePW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((SycamorePW1), center= 147.4742, scale=14.19027)
SycamorePW <- SycamorePW$V1 

WillowPW1 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.00025)))))/(2*WillowA)
WillowPW2 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.00025)))))/(2*WillowA)
WillowPW <- ggfortify::unscale((WillowPW2), center= 147.4742, scale=14.19027)-ggfortify::unscale((WillowPW1), center= 147.4742, scale=14.19027)
WillowPW <- WillowPW$V1 

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

MeanPW1.5 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW2.5 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW.5 <- ggfortify::unscale((MeanPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((MeanPW1.5), center= 147.4742, scale=14.19027)
MeanPW.5 <- MeanPW.5$V1 

AlderPW1.5 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW2.5 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW.5 <- ggfortify::unscale((AlderPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((AlderPW1.5), center= 147.4742, scale=14.19027)
AlderPW.5 <- AlderPW.5$V1 

AshPW1.5 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW2.5 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW.5 <- ggfortify::unscale((AshPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((AshPW1.5), center= 147.4742, scale=14.19027)
AshPW.5 <- AshPW.5$V1 

BeechPW1.5 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW2.5 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW.5 <- ggfortify::unscale((BeechPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((BeechPW1.5), center= 147.4742, scale=14.19027)
BeechPW.5 <- BeechPW.5$V1 

BirchPW1.5 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW2.5 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW.5 <- ggfortify::unscale((BirchPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((BirchPW1.5), center= 147.4742, scale=14.19027)
BirchPW.5 <- BirchPW.5$V1 

ElmPW1.5 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW2.5 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW.5 <- ggfortify::unscale((ElmPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((ElmPW1.5), center= 147.4742, scale=14.19027)
ElmPW.5 <- ElmPW.5$V1 

HazelPW1.5 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW2.5 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW.5 <- ggfortify::unscale((HazelPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((HazelPW1.5), center= 147.4742, scale=14.19027)
HazelPW.5 <- HazelPW.5$V1 

OakPW1.5 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW2.5 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW.5 <- ggfortify::unscale((OakPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((OakPW1.5), center= 147.4742, scale=14.19027)
OakPW.5 <- OakPW.5$V1 

RowanPW1.5 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW2.5 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW.5 <- ggfortify::unscale((RowanPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((RowanPW1.5), center= 147.4742, scale=14.19027)
RowanPW.5 <- RowanPW.5$V1 

SycamorePW1.5 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW2.5 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW.5 <- ggfortify::unscale((SycamorePW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((SycamorePW1.5), center= 147.4742, scale=14.19027)
SycamorePW.5 <- SycamorePW.5$V1 

WillowPW1.5 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW2.5 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW.5 <- ggfortify::unscale((WillowPW2.5), center= 147.4742, scale=14.19027)-ggfortify::unscale((WillowPW1.5), center= 147.4742, scale=14.19027)
WillowPW.5 <- WillowPW.5$V1 

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
BMcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
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

# Peak Date difference from mean
PDdifplot <- ggplot(BMcurves, aes(fct_rev(TT), (PDdif), colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=(PDdifUCI), ymin=(PDdifLCI), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak timing (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Height difference from mean
PHpropdifplot <- ggplot(BMcurves, aes(fct_rev(TT), PHprop, colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (proportional)")+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  coord_flip()+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width difference from mean
PW.5difplot <- ggplot(BMcurves, aes(fct_rev(TT), PW.5dif, colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=PW.5difUCI, ymin=PW.5difLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Relative Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

library(gridExtra)

col5 <- grid.arrange(BiomassCurvesplot, nrow = 1, widths = 1)
col7 <- grid.arrange(PDdifplot, PHpropdifplot, PW.5difplot, nrow = 3, heights = c(1,1,1))
BMCurvesFig2 <- grid.arrange(col5, col7,  ncol = 2, widths = c(4,2)) #saved as 8"x10"


#### Supp. mat. ####

# Peak Width (set height) difference from mean
ggplot(BMcurves, aes(fct_rev(TT), PWdif))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWdifUCI, ymin=PWdifLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Duration (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))# saved as 6"x4"

# Difference to oak
BMODcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                         PD=c(mean(AlderPDOD),mean(AshPDOD),mean(BeechPDOD),mean(BirchPDOD),mean(ElmPDOD),mean(HazelPDOD),mean(RowanPDOD),mean(SycamorePDOD),mean(WillowPDOD)),
                         PDLCI=c(HPDinterval(mcmc(AlderPDOD))[1], HPDinterval(mcmc(AshPDOD))[1], HPDinterval(mcmc(BeechPDOD))[1], HPDinterval(mcmc(BirchPDOD))[1], HPDinterval(mcmc(ElmPDOD))[1], HPDinterval(mcmc(HazelPDOD))[1], HPDinterval(mcmc(RowanPDOD))[1], HPDinterval(mcmc(SycamorePDOD))[1], HPDinterval(mcmc(WillowPDOD))[1]),
                         PDUCI=c(HPDinterval(mcmc(AlderPDOD))[2], HPDinterval(mcmc(AshPDOD))[2], HPDinterval(mcmc(BeechPDOD))[2], HPDinterval(mcmc(BirchPDOD))[2], HPDinterval(mcmc(ElmPDOD))[2], HPDinterval(mcmc(HazelPDOD))[2], HPDinterval(mcmc(RowanPDOD))[2], HPDinterval(mcmc(SycamorePDOD))[2], HPDinterval(mcmc(WillowPDOD))[2]),
                         PHprop=c(mean(AlderPHpropOD),mean(AshPHpropOD),mean(BeechPHpropOD),mean(BirchPHpropOD),mean(ElmPHpropOD),mean(HazelPHpropOD),mean(RowanPHpropOD),mean(SycamorePHpropOD),mean(WillowPHpropOD)),
                         PHpropLCI=c(HPDinterval(mcmc(AlderPHpropOD))[1], HPDinterval(mcmc(AshPHpropOD))[1], HPDinterval(mcmc(BeechPHpropOD))[1], HPDinterval(mcmc(BirchPHpropOD))[1], HPDinterval(mcmc(ElmPHpropOD))[1], HPDinterval(mcmc(HazelPHpropOD))[1], HPDinterval(mcmc(RowanPHpropOD))[1], HPDinterval(mcmc(SycamorePHpropOD))[1], HPDinterval(mcmc(WillowPHpropOD))[1]),
                         PHpropUCI=c(HPDinterval(mcmc(AlderPHpropOD))[2], HPDinterval(mcmc(AshPHpropOD))[2], HPDinterval(mcmc(BeechPHpropOD))[2], HPDinterval(mcmc(BirchPHpropOD))[2], HPDinterval(mcmc(ElmPHpropOD))[2], HPDinterval(mcmc(HazelPHpropOD))[2], HPDinterval(mcmc(RowanPHpropOD))[2], HPDinterval(mcmc(SycamorePHpropOD))[2], HPDinterval(mcmc(WillowPHpropOD))[2]),
                         PW=c(mean(AlderPWOD, na.rm=TRUE),mean(AshPWOD, na.rm=TRUE),mean(BeechPWOD, na.rm=TRUE),mean(BirchPWOD, na.rm=TRUE),mean(ElmPWOD, na.rm=TRUE),mean(HazelPWOD, na.rm=TRUE),mean(RowanPWOD, na.rm=TRUE),mean(SycamorePWOD, na.rm=TRUE),mean(WillowPWOD, na.rm=TRUE)),
                         PWLCI=c(HPDinterval(mcmc(AlderPWOD))[1], HPDinterval(mcmc(AshPWOD))[1], HPDinterval(mcmc(BeechPWOD))[1], HPDinterval(mcmc(BirchPWOD))[1], HPDinterval(mcmc(ElmPWOD))[1], HPDinterval(mcmc(HazelPWOD))[1], HPDinterval(mcmc(RowanPWOD))[1], HPDinterval(mcmc(SycamorePWOD))[1], HPDinterval(mcmc(WillowPWOD))[1]),
                         PWUCI=c(HPDinterval(mcmc(AlderPWOD))[2], HPDinterval(mcmc(AshPWOD))[2], HPDinterval(mcmc(BeechPWOD))[2], HPDinterval(mcmc(BirchPWOD))[2], HPDinterval(mcmc(ElmPWOD))[2], HPDinterval(mcmc(HazelPWOD))[2], HPDinterval(mcmc(RowanPWOD))[2], HPDinterval(mcmc(SycamorePWOD))[2], HPDinterval(mcmc(WillowPWOD))[2]),
                         PW.5=c(mean(AlderPW.5OD, na.rm=TRUE),mean(AshPW.5OD, na.rm=TRUE),mean(BeechPW.5OD, na.rm=TRUE),mean(BirchPW.5OD, na.rm=TRUE),mean(ElmPW.5OD, na.rm=TRUE),mean(HazelPW.5OD, na.rm=TRUE),mean(RowanPW.5OD, na.rm=TRUE),mean(SycamorePW.5OD, na.rm=TRUE),mean(WillowPW.5OD, na.rm=TRUE)),
                         PW.5LCI=c(HPDinterval(mcmc(AlderPW.5OD))[1], HPDinterval(mcmc(AshPW.5OD))[1], HPDinterval(mcmc(BeechPW.5OD))[1], HPDinterval(mcmc(BirchPW.5OD))[1], HPDinterval(mcmc(ElmPW.5OD))[1], HPDinterval(mcmc(HazelPW.5OD))[1], HPDinterval(mcmc(RowanPW.5OD))[1], HPDinterval(mcmc(SycamorePW.5OD))[1], HPDinterval(mcmc(WillowPW.5OD))[1]),
                         PW.5UCI=c(HPDinterval(mcmc(AlderPW.5OD))[2], HPDinterval(mcmc(AshPW.5OD))[2], HPDinterval(mcmc(BeechPW.5OD))[2], HPDinterval(mcmc(BirchPW.5OD))[2], HPDinterval(mcmc(ElmPW.5OD))[2], HPDinterval(mcmc(HazelPW.5OD))[2], HPDinterval(mcmc(RowanPW.5OD))[2], HPDinterval(mcmc(SycamorePW.5OD))[2], HPDinterval(mcmc(WillowPW.5OD))[2]))

NoOakCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "royalblue4", "slateblue2", "orchid")

PDODplot <- ggplot(BMODcurves, aes(fct_rev(TT), PD))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PDUCI, ymin=PDLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Timing (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PHpropODplot <- ggplot(BMODcurves, aes(fct_rev(TT), PHprop))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (proportional)")+
  coord_flip()+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWODplot <- ggplot(BMODcurves, aes(fct_rev(TT), PW))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWUCI, ymin=PWLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Relative Duration (days) ")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWOD.5plot <- ggplot(BMODcurves, aes(fct_rev(TT), PW.5))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PW.5UCI, ymin=PW.5LCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Relative Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

ODplots <- grid.arrange(PDODplot, PHpropODplot, PWOD.5plot,PWODplot,  ncol = 4, widths=c(1,1,1,1)) #saved 6"x13"


#### Posterior distribution sample size (without NaNs) ####
BMCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                                 PD=c(length(which(is.na(AlderPD)==FALSE)),length(which(is.na(AshPD)==FALSE)),length(which(is.na(BeechPD)==FALSE)),length(which(is.na(BirchPD)==FALSE)),length(which(is.na(ElmPD)==FALSE)),length(which(is.na(HazelPD)==FALSE)),length(which(is.na(OakPD)==FALSE)),length(which(is.na(RowanPD)==FALSE)),length(which(is.na(SycamorePD)==FALSE)),length(which(is.na(WillowPD)==FALSE))),
                                 PDdif=c(length(which(is.na(AlderPDdif)==FALSE)),length(which(is.na(AshPDdif)==FALSE)),length(which(is.na(BeechPDdif)==FALSE)),length(which(is.na(BirchPDdif)==FALSE)),length(which(is.na(ElmPDdif)==FALSE)),length(which(is.na(HazelPDdif)==FALSE)),length(which(is.na(OakPDdif)==FALSE)),length(which(is.na(RowanPDdif)==FALSE)),length(which(is.na(SycamorePDdif)==FALSE)),length(which(is.na(WillowPDdif)==FALSE))),
                                 PH=c(length(which(is.na(AlderPH)==FALSE)),length(which(is.na(AshPH)==FALSE)),length(which(is.na(BeechPH)==FALSE)),length(which(is.na(BirchPH)==FALSE)),length(which(is.na(ElmPH)==FALSE)),length(which(is.na(HazelPH)==FALSE)),length(which(is.na(OakPH)==FALSE)),length(which(is.na(RowanPH)==FALSE)),length(which(is.na(SycamorePH)==FALSE)),length(which(is.na(WillowPH)==FALSE))),
                                 PHdif=c(length(which(is.na(AlderPHdif)==FALSE)),length(which(is.na(AshPHdif)==FALSE)),length(which(is.na(BeechPHdif)==FALSE)),length(which(is.na(BirchPHdif)==FALSE)),length(which(is.na(ElmPHdif)==FALSE)),length(which(is.na(HazelPHdif)==FALSE)),length(which(is.na(OakPHdif)==FALSE)),length(which(is.na(RowanPHdif)==FALSE)),length(which(is.na(SycamorePHdif)==FALSE)),length(which(is.na(WillowPHdif)==FALSE))),
                                 PHprop=c(length(which(is.na(AlderPHprop)==FALSE)),length(which(is.na(AshPHprop)==FALSE)),length(which(is.na(BeechPHprop)==FALSE)),length(which(is.na(BirchPHprop)==FALSE)),length(which(is.na(ElmPHprop)==FALSE)),length(which(is.na(HazelPHprop)==FALSE)),length(which(is.na(OakPHprop)==FALSE)),length(which(is.na(RowanPHprop)==FALSE)),length(which(is.na(SycamorePHprop)==FALSE)),length(which(is.na(WillowPHprop)==FALSE))),
                                 PW=c(length(which(is.na(AlderPW)==FALSE)),length(which(is.na(AshPW)==FALSE)),length(which(is.na(BeechPW)==FALSE)),length(which(is.na(BirchPW)==FALSE)),length(which(is.na(ElmPW)==FALSE)),length(which(is.na(HazelPW)==FALSE)),length(which(is.na(OakPW)==FALSE)),length(which(is.na(RowanPW)==FALSE)),length(which(is.na(SycamorePW)==FALSE)),length(which(is.na(WillowPW)==FALSE))),
                                 PWdif=c(length(which(is.na(AlderPWdif)==FALSE)),length(which(is.na(AshPWdif)==FALSE)),length(which(is.na(BeechPWdif)==FALSE)),length(which(is.na(BirchPWdif)==FALSE)),length(which(is.na(ElmPWdif)==FALSE)),length(which(is.na(HazelPWdif)==FALSE)),length(which(is.na(OakPWdif)==FALSE)),length(which(is.na(RowanPWdif)==FALSE)),length(which(is.na(SycamorePWdif)==FALSE)),length(which(is.na(WillowPWdif)==FALSE))),
                                 PW.5=c(length(which(is.na(AlderPW.5)==FALSE)),length(which(is.na(AshPW.5)==FALSE)),length(which(is.na(BeechPW.5)==FALSE)),length(which(is.na(BirchPW.5)==FALSE)),length(which(is.na(ElmPW.5)==FALSE)),length(which(is.na(HazelPW.5)==FALSE)),length(which(is.na(OakPW.5)==FALSE)),length(which(is.na(RowanPW.5)==FALSE)),length(which(is.na(SycamorePW.5)==FALSE)),length(which(is.na(WillowPW.5)==FALSE))),
                                 PW.5dif=c(length(which(is.na(AlderPW.5dif)==FALSE)),length(which(is.na(AshPW.5dif)==FALSE)),length(which(is.na(BeechPW.5dif)==FALSE)),length(which(is.na(BirchPW.5dif)==FALSE)),length(which(is.na(ElmPW.5dif)==FALSE)),length(which(is.na(HazelPW.5dif)==FALSE)),length(which(is.na(OakPW.5dif)==FALSE)),length(which(is.na(RowanPW.5dif)==FALSE)),length(which(is.na(SycamorePW.5dif)==FALSE)),length(which(is.na(WillowPW.5dif)==FALSE))))

BMODCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                                   PD=c(length(which(is.na(AlderPDOD)==FALSE)),length(which(is.na(AshPDOD)==FALSE)),length(which(is.na(BeechPDOD)==FALSE)),length(which(is.na(BirchPDOD)==FALSE)),length(which(is.na(ElmPDOD)==FALSE)),length(which(is.na(HazelPDOD)==FALSE)),length(which(is.na(RowanPDOD)==FALSE)),length(which(is.na(SycamorePDOD)==FALSE)),length(which(is.na(WillowPDOD)==FALSE))),
                                   PH=c(length(which(is.na(AlderPHOD)==FALSE)),length(which(is.na(AshPHOD)==FALSE)),length(which(is.na(BeechPHOD)==FALSE)),length(which(is.na(BirchPHOD)==FALSE)),length(which(is.na(ElmPHOD)==FALSE)),length(which(is.na(HazelPHOD)==FALSE)),length(which(is.na(RowanPHOD)==FALSE)),length(which(is.na(SycamorePHOD)==FALSE)),length(which(is.na(WillowPHOD)==FALSE))),
                                   PHprop=c(length(which(is.na(AlderPHpropOD)==FALSE)),length(which(is.na(AshPHpropOD)==FALSE)),length(which(is.na(BeechPHpropOD)==FALSE)),length(which(is.na(BirchPHpropOD)==FALSE)),length(which(is.na(ElmPHpropOD)==FALSE)),length(which(is.na(HazelPHpropOD)==FALSE)),length(which(is.na(RowanPHpropOD)==FALSE)),length(which(is.na(SycamorePHpropOD)==FALSE)),length(which(is.na(WillowPHpropOD)==FALSE))),
                                   PW=c(length(which(is.na(AlderPWOD)==FALSE)),length(which(is.na(AshPWOD)==FALSE)),length(which(is.na(BeechPWOD)==FALSE)),length(which(is.na(BirchPWOD)==FALSE)),length(which(is.na(ElmPWOD)==FALSE)),length(which(is.na(HazelPWOD)==FALSE)),length(which(is.na(RowanPWOD)==FALSE)),length(which(is.na(SycamorePWOD)==FALSE)),length(which(is.na(WillowPWOD)==FALSE))),
                                   PW.5=c(length(which(is.na(AlderPW.5OD)==FALSE)),length(which(is.na(AshPW.5OD)==FALSE)),length(which(is.na(BeechPW.5OD)==FALSE)),length(which(is.na(BirchPW.5OD)==FALSE)),length(which(is.na(ElmPW.5OD)==FALSE)),length(which(is.na(HazelPW.5OD)==FALSE)),length(which(is.na(RowanPW.5OD)==FALSE)),length(which(is.na(SycamorePW.5OD)==FALSE)),length(which(is.na(WillowPW.5OD)==FALSE))))

meanmetrics <- data.frame(metric=c("PD", "PH", "PW.5", "PW.01"), 
                          mean=c(mean(MeanPD, na.rm=T), mean(MeanPH, na.rm=T),mean(MeanPW.5,na.rm=TRUE), mean(MeanPW,na.rm=TRUE)), 
                          lci=c(HPDinterval(mcmc(MeanPD))[1], HPDinterval(mcmc(MeanPH))[1],HPDinterval(mcmc(MeanPW.5))[1], HPDinterval(mcmc(MeanPW))[1]),
                          uci=c(HPDinterval(mcmc(MeanPD))[2], HPDinterval(mcmc(MeanPH))[2],HPDinterval(mcmc(MeanPW.5))[2], HPDinterval(mcmc(MeanPW))[2]),
                          SS=c(length(which(is.na(MeanPD)==FALSE)),length(which(is.na(MeanPH)==FALSE)),length(which(is.na(MeanPW.5)==FALSE)),length(which(is.na(MeanPW)==FALSE))))



#### Credible intervals for difference between abundance and biomass timing and duration ####
BMPD <- MeanPD #using biomass pd
BMPW <- MeanPW.5 #using biomass pw at half height so comparable

MeanA <- (sampleATTC[,3])
MeanB <- (sampleATTC[,2])
MeanC <- (sampleATTC[,1])

APD <- ggfortify::unscale((-MeanB/(2*MeanA)), center= 147.4742, scale=14.19027)$V1

APH <- exp(MeanC +MeanB*(-MeanB/(2*MeanA)) + MeanA*(-MeanB/(2*MeanA))^2)

APW1 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(APH/2)))))/(2*MeanA)
APW2 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(APH/2)))))/(2*MeanA)
APW <- ggfortify::unscale((APW2), center= 147.4742, scale=14.19027)$V1-ggfortify::unscale((APW1), center= 147.4742, scale=14.19027)$V1 

mean(mcmc(APW), na.rm=TRUE)
HPDinterval(mcmc(APW), na.rm=TRUE)

mean(BMPD-APD, na.rm=TRUE)
HPDinterval(mcmc(BMPD-APD), na.rm=TRUE)
mean(BMPW-APW, na.rm=TRUE)
HPDinterval(mcmc(BMPW-APW), na.rm=TRUE)
