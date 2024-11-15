rm(list=ls())
library(MCMCglmm)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(colorspace)
library(gridExtra)

##################################################################
#### Model: Phenological distribution of mass among tree taxa ####
##################################################################

cater_habitat <- read.csv("~/cater_habitat_data.csv")
cater_habitat_1723 <- subset(cater_habitat, year!="2014")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2015")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2016")

load("~/Mass23_1.RData")  
load("~/Mass23_2.RData")
load("~/Mass23_3.RData")

## Check chains have converged to same range of values: (rows should have similar values)
fixed_mean <- data.frame(c1=apply(Mass23_1$Sol[,1:3],2,mean),
                         c2=apply(Mass23_2$Sol[,1:3],2,mean),
                         c3=apply(Mass23_3$Sol[,1:3],2,mean))
fixed_sd <- data.frame(c1=apply(Mass23_1$Sol[,1:3],2,sd),
                       c2=apply(Mass23_2$Sol[,1:3],2,sd),
                       c3=apply(Mass23_3$Sol[,1:3],2,sd))
random_mean <- data.frame(c1=apply(Mass23_1$VCV[,1:15],2,mean),
                          c2=apply(Mass23_2$VCV[,1:15],2,mean),
                          c3=apply(Mass23_3$VCV[,1:15],2,mean))
random_sd <- data.frame(c1=apply(Mass23_1$VCV[,1:15],2,sd),
                        c2=apply(Mass23_2$VCV[,1:15],2,sd),
                        c3=apply(Mass23_3$VCV[,1:15],2,sd))

## Join chains
Sol <- rbind(Mass23_1$Sol,Mass23_2$Sol,Mass23_3$Sol)
dim(Sol)
VCV <- rbind(Mass23_1$VCV,Mass23_2$VCV,Mass23_3$VCV)
dim(VCV)
X <- Mass23_1$X
Mass23 <- list(Sol= Sol,VCV=VCV,X=X)

#### Checking model fits the data and converged ####
plot(Mass23$VCV) 
Mass23.Sim<-simulate(Mass23_1,nsim=100) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(Mass23.Sim,2,mean),50) #histogram of simulation predictions for mean mass (log scale)
abline(v=mean((cater_habitat_1723$logmpc1+cater_habitat_1723$logmpc2)/2, na.rm=TRUE),col=2) # red line for mean mass (from the mean of logmpc1 and logmpc2)

# tree taxon specific coefficients for date^2 (A), date (B) and intercept (C)
AlderC <- (Mass23$Sol[,1] + Mass23$Sol[,4])
AshC <- (Mass23$Sol[,1] + Mass23$Sol[,5])
BeechC <- (Mass23$Sol[,1] + Mass23$Sol[,6])
BirchC <- (Mass23$Sol[,1] + Mass23$Sol[,7])
ElmC <- (Mass23$Sol[,1] + Mass23$Sol[,8])
HazelC <- (Mass23$Sol[,1] + Mass23$Sol[,9])
OakC <- (Mass23$Sol[,1] + Mass23$Sol[,10])
RowanC <- (Mass23$Sol[,1] + Mass23$Sol[,11])
SycamoreC <- (Mass23$Sol[,1] + Mass23$Sol[,12])
WillowC <- (Mass23$Sol[,1] + Mass23$Sol[,13])

AlderB <- (Mass23$Sol[,2] + Mass23$Sol[,14])
AshB <- (Mass23$Sol[,2] + Mass23$Sol[,15])
BeechB <- (Mass23$Sol[,2] + Mass23$Sol[,16])
BirchB <- (Mass23$Sol[,2] + Mass23$Sol[,17])
ElmB <- (Mass23$Sol[,2] + Mass23$Sol[,18])
HazelB <- (Mass23$Sol[,2] + Mass23$Sol[,19])
OakB <- (Mass23$Sol[,2] + Mass23$Sol[,20])
RowanB <- (Mass23$Sol[,2] + Mass23$Sol[,21])
SycamoreB <- (Mass23$Sol[,2] + Mass23$Sol[,22])
WillowB <- (Mass23$Sol[,2] + Mass23$Sol[,23])

AlderA <- (Mass23$Sol[,3])
AshA <- (Mass23$Sol[,3])
BeechA <- (Mass23$Sol[,3])
BirchA <- (Mass23$Sol[,3])
ElmA <- (Mass23$Sol[,3])
HazelA <- (Mass23$Sol[,3])
OakA <- (Mass23$Sol[,3])
RowanA <- (Mass23$Sol[,3])
SycamoreA <- (Mass23$Sol[,3])
WillowA <- (Mass23$Sol[,3])

# Fixed effect coefficients 
MeanA <- (Mass23$Sol[,3])
MeanB <- (Mass23$Sol[,2])
MeanC <- (Mass23$Sol[,1])

#### Plotting growth curves ####

preddayscaled <- seq(min(cater_habitat_1723$datescaled),max(cater_habitat_1723$datescaled),0.001)
predday <- ggfortify::unscale(preddayscaled, center= 147.4742, scale=14.19027)$V1


meanslope <- mean(MeanC)+mean(MeanB)*preddayscaled+mean(MeanA)*preddayscaled^2
alder <- mean(AlderC)+mean(AlderB)*preddayscaled+mean(AlderA)*preddayscaled^2
ash <- mean(AshC)+mean(AshB)*preddayscaled+mean(AshA)*preddayscaled^2
beech <- mean(BeechC)+mean(BeechB)*preddayscaled+mean(BeechA)*preddayscaled^2
birch <- mean(BirchC)+mean(BirchB)*preddayscaled+mean(BirchA)*preddayscaled^2
elm <- mean(ElmC)+mean(ElmB)*preddayscaled+mean(ElmA)*preddayscaled^2
hazel <- mean(HazelC)+mean(HazelB)*preddayscaled+mean(HazelA)*preddayscaled^2
oak <- mean(OakC)+mean(OakB)*preddayscaled+mean(OakA)*preddayscaled^2
rowan <- mean(RowanC)+mean(RowanB)*preddayscaled+mean(RowanA)*preddayscaled^2
sycamore <- mean(SycamoreC)+mean(SycamoreB)*preddayscaled+mean(SycamoreA)*preddayscaled^2
willow <- mean(WillowC)+mean(WillowB)*preddayscaled+mean(WillowA)*preddayscaled^2

mycolblack <- rgb(100, 100, 100, max = 250, alpha = 50, names = "blacktrans")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

fixedslope <- data.frame(date=predday, logmass= meanslope, mass=exp(meanslope))
fixedslope <- fixedslope[fixedslope$date < 172, ]

cater_habitat_1723$mean.mpc <- (cater_habitat_1723$mpc1+cater_habitat_1723$mpc2)/2
cater_habitat_1723$log.mean.mpc <- log(cater_habitat_1723$mean.mpc)

(dataplot <- ggplot(cater_habitat_1723, aes(date, mean.mpc))+
    geom_point(col=mycolblack, size=0.5)+
    geom_line(data=fixedslope, aes(date, mass), col=1, lty=5)+
    xlab("Ordinal Date")+
    ylab("Mass (g)")+
    coord_cartesian(xlim=c(117, 172))+ 
    coord_trans(y="log")+
    scale_y_continuous(breaks=c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1), labels=c("0.01", "0.02", "0.03", "0.04", "0.05", "0.10", "0.20", "0.30", "0.40", "0.50","1.00"), minor_breaks = c(0.06, 0.07, 0.08, 0.09, 0.6, 0.7, 0.8, 0.9))+ #expand = c(0,0.02), 
    scale_x_continuous(breaks=seq(120,170,10))+
    theme_bw()+
    theme(text=element_text(size= 15))+
    annotate(geom="text", x=132, y=0.95, label="- - -  Fixed effect",
             color="black"))

taxoncurves <- data.frame(date = predday, FixedEffect=exp(meanslope), Alder=exp(alder), Ash=exp(ash), Beech=exp(beech), Birch=exp(birch), Elm=exp(elm), Hazel=exp(hazel), Oak=exp(oak), Rowan=exp(rowan),Sycamore=exp(sycamore),Willow=exp(willow)) 
taxoncurveslong <- gather(taxoncurves, Taxon, Mass, 2:12)
taxoncurveslong$Taxon <- as.factor(taxoncurveslong$Taxon)
taxoncurveslong$Taxon <- factor(taxoncurveslong$Taxon , c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak","Rowan","Sycamore","Willow" ,"FixedEffect"))

AllCols <- c( "darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid","black")

(taxonplot <- ggplot(taxoncurveslong, aes(date,Mass, col=Taxon))+
    geom_line(aes(linetype=Taxon))+
    xlab("Ordinal Date")+
    ylab("Mass (g)")+
    scale_colour_manual(values=AllCols)+
    scale_linetype_manual(values=c(1,1,1,1,1,1,1,1,1,1,5))+
    xlim(117, 170)+
    geom_vline(xintercept=168, linetype="dotted", colour="gray45", size=0.5)+
    scale_y_continuous(expand = c(0,0.003))+
    guides(color = "none", linetype="none")+
    theme_bw()+
    theme(text=element_text(size= 15)))

#### Mass on 168 ####
latedate <- scale(168,  center= 147.4742, scale=14.19027)[1,1]
Mean168 <- exp(MeanC + MeanB*latedate + MeanA*latedate^2)
Alder168 <- exp(AlderC + AlderB*latedate + AlderA*latedate^2)
Ash168 <- exp(AshC + AshB*latedate + AshA*latedate^2)
Beech168 <- exp(BeechC + BeechB*latedate + BeechA*latedate^2)
Birch168 <- exp(BirchC + BirchB*latedate + BirchA*latedate^2)
Elm168 <- exp(ElmC + ElmB*latedate + ElmA*latedate^2)
Hazel168 <- exp(HazelC + HazelB*latedate + HazelA*latedate^2)
Oak168 <- exp(OakC + OakB*latedate + OakA*latedate^2)
Rowan168 <- exp(RowanC + RowanB*latedate + RowanA*latedate^2)
Sycamore168 <- exp(SycamoreC + SycamoreB*latedate + SycamoreA*latedate^2)
Willow168 <- exp(WillowC + WillowB*latedate + WillowA*latedate^2)

#### Difference to mean on 168 (0.96) ####
Alder168dif <- Alder168-Mean168
Ash168dif <- Ash168-Mean168
Beech168dif <- Beech168-Mean168
Birch168dif <- Birch168-Mean168
Elm168dif <- Elm168-Mean168
Hazel168dif <- Hazel168-Mean168
Oak168dif <- Oak168-Mean168
Rowan168dif <- Rowan168-Mean168
Sycamore168dif <- Sycamore168-Mean168
Willow168dif <- Willow168-Mean168

#### Difference to mean on 168 (0.96) ####
Alder168prop <- Alder168/Mean168
Ash168prop <- Ash168/Mean168
Beech168prop <- Beech168/Mean168
Birch168prop <- Birch168/Mean168
Elm168prop <- Elm168/Mean168
Hazel168prop <- Hazel168/Mean168
Oak168prop <- Oak168/Mean168
Rowan168prop <- Rowan168/Mean168
Sycamore168prop <- Sycamore168/Mean168
Willow168prop <- Willow168/Mean168

# mean coefficient and CIs for mass of 168 for each tree taxon
Mass168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                      median=c(median(Alder168), median(Ash168), median(Beech168), median(Birch168), median(Elm168), median(Hazel168), median(Oak168), median(Rowan168), median(Sycamore168), median(Willow168)),
                      lowci=c(HPDinterval(mcmc(Alder168))[1],HPDinterval(mcmc(Ash168))[1],HPDinterval(mcmc(Beech168))[1],HPDinterval(mcmc(Birch168))[1],HPDinterval(mcmc(Elm168))[1],HPDinterval(mcmc(Hazel168))[1],HPDinterval(mcmc(Oak168))[1],HPDinterval(mcmc(Rowan168))[1],HPDinterval(mcmc(Sycamore168))[1],HPDinterval(mcmc(Willow168))[1]),
                      upci=c(HPDinterval(mcmc(Alder168))[2],HPDinterval(mcmc(Ash168))[2],HPDinterval(mcmc(Beech168))[2],HPDinterval(mcmc(Birch168))[2],HPDinterval(mcmc(Elm168))[2],HPDinterval(mcmc(Hazel168))[2],HPDinterval(mcmc(Oak168))[2],HPDinterval(mcmc(Rowan168))[2],HPDinterval(mcmc(Sycamore168))[2],HPDinterval(mcmc(Willow168))[2]),
                      meanprop=c(mean(Alder168prop), mean(Ash168prop), mean(Beech168prop), mean(Birch168prop), mean(Elm168prop), mean(Hazel168prop), mean(Oak168prop), mean(Rowan168prop), mean(Sycamore168prop), mean(Willow168prop)),
                      proplowci=c(HPDinterval(mcmc(Alder168prop))[1],HPDinterval(mcmc(Ash168prop))[1],HPDinterval(mcmc(Beech168prop))[1],HPDinterval(mcmc(Birch168prop))[1],HPDinterval(mcmc(Elm168prop))[1],HPDinterval(mcmc(Hazel168prop))[1],HPDinterval(mcmc(Oak168prop))[1],HPDinterval(mcmc(Rowan168prop))[1],HPDinterval(mcmc(Sycamore168prop))[1],HPDinterval(mcmc(Willow168prop))[1]),
                      propupci=c(HPDinterval(mcmc(Alder168prop))[2],HPDinterval(mcmc(Ash168prop))[2],HPDinterval(mcmc(Beech168prop))[2],HPDinterval(mcmc(Birch168prop))[2],HPDinterval(mcmc(Elm168prop))[2],HPDinterval(mcmc(Hazel168prop))[2],HPDinterval(mcmc(Oak168prop))[2],HPDinterval(mcmc(Rowan168prop))[2],HPDinterval(mcmc(Sycamore168prop))[2],HPDinterval(mcmc(Willow168prop))[2]))

# Plot of mass difference to the fixed effect prediction
(Mass168propplot <- ggplot(Mass168, aes(fct_rev(TT), meanprop, col=TT))+
    geom_point(size=2, alpha=0.9)+
    geom_errorbar(aes(ymax=propupci, ymin=proplowci, width=0.5))+
    geom_hline(yintercept=1, linetype="longdash", size=0.3)+  
    coord_flip()+
    xlab("Tree Taxon")+
    ylab("Prop. difference to fixed eff.")+
    theme_bw()+
    theme(text=element_text(size= 15))+
    scale_colour_manual(values=AllTaxaCols)+
    guides(color = "none")) # saved as 6"x4.5"

space <- ggplot() + theme_void()
MassPlots <- grid.arrange(dataplot,space, taxonplot,space, Mass168propplot, ncol = 5, widths = c(1,0.1,1,0.1,1.1)) #saved as 6"x12"


#### For supp mat ####
## Difference to oak on 168 (0.96)
AlderOD <- Alder168/Oak168
AshOD <- Ash168/Oak168
BeechOD <- Beech168/Oak168 
BirchOD <- Birch168/Oak168 
ElmOD <- Elm168/Oak168 
HazelOD <- Hazel168/Oak168 
RowanOD <- Rowan168/Oak168
SycamoreOD <- Sycamore168/Oak168
WillowOD <- Willow168/Oak168

# mean and CIs for mass difference to oak on day 168
OD168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                    mean=c(mean(AlderOD), mean(AshOD), mean(BeechOD), mean(BirchOD), mean(ElmOD), mean(HazelOD), mean(RowanOD), mean(SycamoreOD), mean(WillowOD)),
                    lowci=c(HPDinterval(mcmc(AlderOD))[1],HPDinterval(mcmc(AshOD))[1],HPDinterval(mcmc(BeechOD))[1],HPDinterval(mcmc(BirchOD))[1],HPDinterval(mcmc(ElmOD))[1],HPDinterval(mcmc(HazelOD))[1],HPDinterval(mcmc(RowanOD))[1],HPDinterval(mcmc(SycamoreOD))[1],HPDinterval(mcmc(WillowOD))[1]),
                    upci=c(HPDinterval(mcmc(AlderOD))[2],HPDinterval(mcmc(AshOD))[2],HPDinterval(mcmc(BeechOD))[2],HPDinterval(mcmc(BirchOD))[2],HPDinterval(mcmc(ElmOD))[2],HPDinterval(mcmc(HazelOD))[2],HPDinterval(mcmc(RowanOD))[2],HPDinterval(mcmc(SycamoreOD))[2],HPDinterval(mcmc(WillowOD))[2]))

# Plot of mass difference to oak prediction
ggplot(OD168, aes(fct_rev(TT), mean))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=1, linetype="dashed", colour="black", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Prop. difference to oak")+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 15),axis.text.x = element_text(angle = 45, hjust=1)) # saved as 6"x4"


meanmetrics <- data.frame("168mean"=round(mean(Mean168),3), 
                          lci=round(HPDinterval(mcmc(Mean168))[1],3),
                          uci=round(HPDinterval(mcmc(Mean168))[2],3),
                          SS=c(length(which(is.na(Mean168)==FALSE))))

############################
#### Model output table ####   
############################

#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(Mass23$Sol[,1]),3)," (",
                      round(HPDinterval(mcmc(Mass23$Sol[,1]))[1],3)," - ",
                      round(HPDinterval(mcmc(Mass23$Sol[,1]))[2],3),")",sep=""),
    round(effectiveSize(Mass23$Sol[,1]))),
  c("Date scaled",paste(round(mean(Mass23$Sol[,2]),3)," (",
                        round(HPDinterval(mcmc(Mass23$Sol[,2]))[1],3)," - ",
                        round(HPDinterval(mcmc(Mass23$Sol[,2]))[2],3),")",sep=""),
    round(effectiveSize(Mass23$Sol[,2]))),
  c("Date² scaled",paste(round(mean(Mass23$Sol[,3]),3)," (",
                         round(HPDinterval(mcmc(Mass23$Sol[,3]))[1],3)," - ",
                         round(HPDinterval(mcmc(Mass23$Sol[,3]))[2],3),")",sep=""),
    round(effectiveSize(Mass23$Sol[,3]))))

#### random terms 
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
             round(effectiveSize(Mass23$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                                          round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                                          round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
             round(effectiveSize(Mass23$VCV[, column])))


column<-4
treetaxa4<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                              round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                              round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
             round(effectiveSize(Mass23$VCV[, column])))


column<-5
site5<-c("Site- Intercept var",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                     round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                     round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
         round(effectiveSize(Mass23$VCV[, column])))

column<-6
site6<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                                  round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                                  round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
         round(effectiveSize(Mass23$VCV[, column])))


column<-8
site8<-c("Site- Date slope var",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                                      round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                                      round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
         round(effectiveSize(Mass23$VCV[, column])))

column<-9
year<-c("Year",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                     round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                     round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
        round(effectiveSize(Mass23$VCV[, column])))

column <- 10
siteyear<-c("Site-Year",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                              round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                              round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
            round(effectiveSize(Mass23$VCV[, column])))


column<-13
recorder<-c("Recorder",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
            round(effectiveSize(Mass23$VCV[, column])))


column<-12
siteday<-c("Site-Day",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                            round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                            round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
           round(effectiveSize(Mass23$VCV[, column])))


column<-11
treeID<-c("Tree ID",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                          round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                          round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
          round(effectiveSize(Mass23$VCV[, column])))

column<-14
weight<-c("Weighting",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                            round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                            round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
          round(effectiveSize(Mass23$VCV[, column])))

column<-15
residual<-c("Residual",paste(round(posterior.mode(mcmc(Mass23$VCV[, column])),3)," (",
                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[1],3)," - ",
                             round(HPDinterval(mcmc(Mass23$VCV[, column]))[2],3),")",sep=""),
            round(effectiveSize(Mass23$VCV[, column])))




random<-rbind(treetaxa1,treetaxa2,treetaxa4,site5,site6,site8,year,siteyear, recorder, siteday, treeID, weight, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/TableS5_Mass.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
