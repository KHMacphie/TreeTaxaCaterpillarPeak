rm(list=ls())
library(MCMCglmm)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library("colorspace")
library(gridExtra)

#####################################################
#### Model: Habitat and TreeTaxa Abundance model ####
#####################################################

cater_habitat <- read.csv("~/cater_habitat_data.csv")
load("~/TTHA23.RData")  
Habitat_Site <- read.csv("~/Habitat_Site.csv")
summary(TTHA23)

#### Checking model fits the data and converged ####

plot(TTHA23) 

TTHA23.Sim<-simulate(TTHA23,nsim=1000) #simulate 1000 times

par(mfcol=c(1,1))
hist(apply(TTHA23.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(TTHA23.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data


# dataframe for meeaans and CIs for beaten tree taxa effects
TTHA <- TTHA23$Sol[,3:12] # crop to just the columns wanted
TTHA.df <- data.frame(treetaxa=c(colnames(TTHA))) #dataframe with column for beaten tree taxa
TTHA.df$coeff <- apply(TTHA,2, mean) 
TTHA.df$prop_coeff <- apply(exp(TTHA),2, mean) 
for(i in 1:length(TTHA.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA[,i])
  TTHA.df$lowci[i] <- A["var1","lower"] 
  TTHA.df$upci[i] <- A["var1","upper"] 
  B <- HPDinterval(exp(TTHA[,i]))
  TTHA.df$prop_lowci[i] <- B["var1","lower"] 
  TTHA.df$prop_upci[i] <- B["var1","upper"] 
} 

TTHA.df$treetaxa <- gsub("tree.species.","", TTHA.df$treetaxa) # adjust name
TTHA.dfs <- TTHA.df
TTHA.dfs[,2:7] <- round(TTHA.dfs[,2:7],2)

#dataframe for means and CIs for habitat FS effects
TTHA2 <- TTHA23$Sol[,13:24] # crop to just the columns wanted
TTHA2.df <- data.frame(treetaxa=c(colnames(TTHA2))) #dataframe with column for FS tree taxa
TTHA2.df$coeff <- apply(TTHA2,2, mean) 
TTHA2.df$exp_coeff <- apply(exp(TTHA2),2, mean) 
for(i in 1:length(TTHA2.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA2[,i])
  TTHA2.df$lowci[i] <- A["var1","lower"] 
  TTHA2.df$upci[i] <- A["var1","upper"] 
  B <- HPDinterval(exp(TTHA2[,i]))
  TTHA2.df$exp_lowci[i] <- B["var1","lower"] 
  TTHA2.df$exp_upci[i] <- B["var1","upper"] 
} 
TTHA2.df$treetaxa <- gsub("_cent.NA.1","", TTHA2.df$treetaxa) # adjust name
TTHA2.df$treetaxa <- gsub("OthDecid","Other", TTHA2.df$treetaxa) # adjust name
TTHA2.dfs <- TTHA2.df
TTHA2.dfs[,2:6] <- round(TTHA2.dfs[,2:6],4)

TTHA.df[11:12,] <- TTHA2.df[11:12,]
TTHA.df[11:12,2:4] <- NA

#plot beaten tree taxa effects
Colours <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid", "gray57", "gray35")

(plot1 <- ggplot(TTHA.df, aes(fct_rev(fct_inorder(treetaxa)), coeff, col=fct_inorder(treetaxa)))+
    geom_point(size=2, alpha=0.9)+
    geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
    theme_bw()+
    coord_flip()+
    theme(text = element_text(size=15))+
    scale_colour_manual(values=Colours)+
    geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+
    xlab("Host tree taxon")+
    ylab("Coefficient (log scale)")+
    guides(linetype="none", colour="none"))

#plot for habitat FS's effects
(plot2 <- ggplot(TTHA2.df, aes(fct_rev(fct_inorder(treetaxa)), coeff, col=fct_inorder(treetaxa)))+
    geom_point(size=2, alpha=0.9)+
    geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
    theme_bw()+
    theme(text = element_text(size=15))+
    scale_colour_manual(values=Colours)+
    geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+
    coord_flip()+
    xlab("Tree stand composition")+
    ylab("Coefficient (log scale)")+
    guides(linetype="none", colour="none"))

#Proportional difference to oak 
alderOD <- (TTHA23$Sol[,3]-TTHA23$Sol[,9])
ashOD <- (TTHA23$Sol[,4]-TTHA23$Sol[,9])
beechOD <- (TTHA23$Sol[,5]-TTHA23$Sol[,9])
birchOD <- (TTHA23$Sol[,6]-TTHA23$Sol[,9])
elmOD <- (TTHA23$Sol[,7]-TTHA23$Sol[,9])
hazelOD <- (TTHA23$Sol[,8]-TTHA23$Sol[,9])
rowanOD <- (TTHA23$Sol[,10]-TTHA23$Sol[,9])
sycamoreOD <- (TTHA23$Sol[,11]-TTHA23$Sol[,9])
willowOD <- (TTHA23$Sol[,12]-TTHA23$Sol[,9])

# extract means and CIs
TTOD <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"),
                   meanOD=c(mean(exp(alderOD)), mean(exp(ashOD)), mean(exp(beechOD)), mean(exp(birchOD)), mean(exp(elmOD)), mean(exp(hazelOD)), mean(exp(rowanOD)), mean(exp(sycamoreOD)), mean(exp(willowOD))),
                   lowci=c(HPDinterval(exp(alderOD))[1], HPDinterval(exp(ashOD))[1],HPDinterval(exp(beechOD))[1],HPDinterval(exp(birchOD))[1],HPDinterval(exp(elmOD))[1],HPDinterval(exp(hazelOD))[1],HPDinterval(exp(rowanOD))[1],HPDinterval(exp(sycamoreOD))[1],HPDinterval(exp(willowOD))[1]),
                   upci=c(HPDinterval(exp(alderOD))[2], HPDinterval(exp(ashOD))[2],HPDinterval(exp(beechOD))[2],HPDinterval(exp(birchOD))[2],HPDinterval(exp(elmOD))[2],HPDinterval(exp(hazelOD))[2],HPDinterval(exp(rowanOD))[2],HPDinterval(exp(sycamoreOD))[2],HPDinterval(exp(willowOD))[2]))



# plot of proportional differences to oak
ggplot(TTOD, aes(fct_rev(TT), meanOD))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=(upci), ymin=(lowci), width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Difference to oak (proportional)")+
  coord_flip()+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1)) #saved 6"x4"

#### Plot with slope for effect on abundance for FS of each tree taxa
#Plotted for range of foliage scores of each taxon present at sites
summary(rowSums(Habitat_Site[2:13])) #real FS total scale
summary(Habitat_Site$Total_cent) #used in model
summary(as.matrix(Habitat_Site[14:25])) #taxaFS's in model

totcentfs <- seq(-44.561, 45.029, 0.001) # for slope calcs 

# range of FS's for each taxa
aldercent <- seq(min(Habitat_Site$Alder_cent), max(Habitat_Site$Alder_cent), 0.01)
ashcent <- seq(min(Habitat_Site$Ash_cent), max(Habitat_Site$Ash_cent), 0.01)
beechcent <- seq(min(Habitat_Site$Beech_cent), max(Habitat_Site$Beech_cent), 0.01)
birchcent <- seq(min(Habitat_Site$Birch_cent), max(Habitat_Site$Birch_cent), 0.01)
elmcent <- seq(min(Habitat_Site$Elm_cent), max(Habitat_Site$Elm_cent), 0.01)
hazelcent <- seq(min(Habitat_Site$Hazel_cent), max(Habitat_Site$Hazel_cent), 0.01)
oakcent <- seq(min(Habitat_Site$Oak_cent), max(Habitat_Site$Oak_cent), 0.01)
rowancent <- seq(min(Habitat_Site$Rowan_cent), max(Habitat_Site$Rowan_cent), 0.01)
sycamorecent <- seq(min(Habitat_Site$Sycamore_cent), max(Habitat_Site$Sycamore_cent), 0.01)
willowcent <- seq(min(Habitat_Site$Willow_cent), max(Habitat_Site$Willow_cent), 0.01)
conifercent <- seq(min(Habitat_Site$Conifer_cent), max(Habitat_Site$Conifer_cent), 0.01)
othercent <- seq(min(Habitat_Site$OthDecid_cent), max(Habitat_Site$OthDecid_cent), 0.01)

Totalslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2])*totcentfs
# each slope coefficient is the deviation from the total slope from fixed effects
Alderslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,13])*aldercent
Ashslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,14])*ashcent
Beechslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,15])*beechcent
Birchslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,16])*birchcent
Elmslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,17])*elmcent
Hazelslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,18])*hazelcent
Oakslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,19])*oakcent
Rowanslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,20])*rowancent
Sycamoreslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,21])*sycamorecent
Willowslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,22])*willowcent
Coniferslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,23])*conifercent
Otherslope <- mean(TTHA23$Sol[,1])+mean(TTHA23$Sol[,2]+TTHA23$Sol[,24])*othercent


## Plot in ggplot
TTFSlong <- data.frame(FS=c(aldercent,ashcent,beechcent,birchcent,elmcent,hazelcent,oakcent,rowancent,sycamorecent,willowcent,conifercent,othercent),
                       logabund=c(Alderslope,Ashslope,Beechslope,Birchslope,Elmslope,Hazelslope,Oakslope,Rowanslope,Sycamoreslope,Willowslope,Coniferslope,Otherslope),
                       TT=c(rep("Alder",length(aldercent)),rep("Ash",length(ashcent)),rep("Beech",length(beechcent)),rep("Birch",length(birchcent)),rep("Elm",length(elmcent)),
                            rep("Hazel",length(hazelcent)),rep("Oak",length(oakcent)),rep("Rowan",length(rowancent)),rep("Sycamore",length(sycamorecent)),rep("Willow",length(willowcent)),
                            rep("Conifer",length(conifercent)),rep("OtherDecid",length(othercent))))

Colours <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid", "gray57", "gray35")

(plot3 <- ggplot(TTFSlong, aes(FS, exp(logabund), col=fct_inorder(TT)))+
    geom_line(size=0.6,aes(linetype=fct_inorder(TT)))+
    theme_bw()+
    theme(text = element_text(size=15))+
    xlab("Deviation in FS from mean")+
    ylab("Caterpillar Abundance")+
    scale_colour_manual(values=Colours)+
    scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","dashed","solid","dashed","dashed","dashed","dashed","dashed"))+
    guides(color = "none", linetype="none"))
gap <- ggplot()+theme_void()
row1 <- grid.arrange(plot1,gap, plot2,gap,plot3, ncol=5, widths=c(1,0.1,1,0.1,1)) #saved as 6"x12"


### looking at credible intervals on oak slope
oakhabitat <- data.frame(FS = seq(min(Habitat_Site$Oak_cent), max(Habitat_Site$Oak_cent), 0.1))
for(i in 1:length(oakhabitat$FS)){
  oakhabitat$mean[i] <- mean(TTHA23$Sol[,1]+((TTHA23$Sol[,2]+TTHA23$Sol[,19])*oakhabitat$FS[i]))
  oakhabitat$LCI[i] <- HPDinterval(TTHA23$Sol[,1]+((TTHA23$Sol[,2]+TTHA23$Sol[,19])*oakhabitat$FS[i]))[1]
  oakhabitat$UCI[i] <- HPDinterval(TTHA23$Sol[,1]+((TTHA23$Sol[,2]+TTHA23$Sol[,19])*oakhabitat$FS[i]))[2]
}

plot(oakhabitat$FS, exp(oakhabitat$mean), type="l", ylim=c(0,0.3), lwd=2)
points(oakhabitat$FS, exp(oakhabitat$LCI), type="l")
points(oakhabitat$FS, exp(oakhabitat$UCI), type="l")

maxoak <- TTHA23$Sol[,1]+((TTHA23$Sol[,2]+TTHA23$Sol[,19])*max(Habitat_Site$Oak_cent))
nooak <- TTHA23$Sol[,1]+((TTHA23$Sol[,2]+TTHA23$Sol[,19])*min(Habitat_Site$Oak_cent))
mean(exp(maxoak-nooak)) 
HPDinterval(mcmc(exp(maxoak-nooak))) 


############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean

####fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(TTHA23$Sol[,1]),3)," (",
                      round(HPDinterval(TTHA23$Sol[,1])[1],3)," - ",
                      round(HPDinterval(TTHA23$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(TTHA23$Sol[,1]))),
  c("Total Foliage Score",paste(round(mean(TTHA23$Sol[,2]),3)," (",
                                round(HPDinterval(TTHA23$Sol[,2])[1],3)," - ",
                                round(HPDinterval(TTHA23$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(TTHA23$Sol[,2]))))

####random terms 
column<-1
treetaxa<-c("Sampled Tree Taxa",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                                      round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA23$VCV[, column])))

column<-2
habitat<-c("Habitat Composition",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                                       round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(TTHA23$VCV[, column])))

column<-3
site<-c("Site",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA23$VCV[, column])))

column<-4
year<-c("Year",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA23$VCV[, column])))

column<-5
siteyear<-c("Site Year",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                              round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                              round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA23$VCV[, column])))

column<-6
treeID<-c("Tree ID",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                          round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                          round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(TTHA23$VCV[, column])))

column<-7
siteday<-c("Site Day",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                            round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                            round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(TTHA23$VCV[, column])))

column<-8
recorder<-c("Recorder",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA23$VCV[, column])))

column<-9
residual<-c("Residual",paste(round(posterior.mode(TTHA23$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA23$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA23$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA23$VCV[, column])))


random<-rbind(treetaxa, habitat, site, year, siteyear, treeID, siteday,recorder,residual)


#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/TableS3_TTHA.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)

