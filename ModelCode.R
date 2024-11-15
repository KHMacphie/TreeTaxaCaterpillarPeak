rm(list=ls())
setwd("/Users/kmacphie/Dropbox/Kirsty's/Chapter1/Inc2023/For Publication/")
library(MCMCglmm)


cater_habitat <- read.csv("cater_habitat_data.csv")

k <- 10000

################
#### Models ####
################


#### Model: Variance decomposition model ####

prior1<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

AbundVar23<- MCMCglmm(caterpillars~datescaled+I(datescaled^2),
                      random=~recorder+year+siteyear+siteday+site+treeID+tree.species+yearday,
                      family="poisson", data=cater_habitat, prior=prior1, nitt=4000000, burnin=100000, thin=1500)
save(AbundVar23, file = "~/AbundVar23.RData") 
plot(AbundVar23)


#### Model: Habitat and TreeTaxa Abundance model ####

prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


TTHA23<- MCMCglmm(caterpillars~Total_cent,
                  random=~tree.species+idv(~Alder_cent+Ash_cent+Beech_cent+Birch_cent+Elm_cent+Hazel_cent+Oak_cent+Rowan_cent+Sycamore_cent+Willow_cent+Conifer_cent+OthDecid_cent)
                    +site+year+siteyear+treeID+siteday+recorder,
                  family="poisson", data=cater_habitat, prior=prior2, nitt=9000,  pr=TRUE)
save(TTHA23, file = "~/TTHA23.RData")
plot(TTHA23)

#### Model: Asymmetry in the caterpillar abundance peak ####

prior3<-list(R=list(V=diag(1), nu=0.002),
            G=list(G1=list(V=diag(4), nu=4, alpha.mu=c(0,0,0,0), alpha.V=diag(4)*k),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*k),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*k),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*k)))

CurveShape23<- MCMCglmm(caterpillars~datescaled+I(datescaled^2)+I(datescaled^3),
                        random=~us(1+datescaled+I(datescaled^2)+I(datescaled^3)):siteyear+recorder+siteday+treeID,
                        family="poisson", data=cater_habitat, prior=prior3, nitt=3500000, burnin=100000, thin=2000 ,pr=TRUE)
save(CurveShape23, file = "~/CurveShape23.RData") 
plot(CurveShape23)


#### Model: Phenological distribution of abundance among tree taxa ####

prior4<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

ATTC23<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2),
                    random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + siteyear + treeID + siteday + recorder,
                    family="poisson", data=cater_habitat, prior=prior4, nitt=3300000, burnin=300000, pr=TRUE, thin=1000)
save(ATTC23, file = "~/ATTC23.RData") 
plot(ATTC23) 

#### Model: Phenological distribution of mass among tree taxa ####

# Mass only available from 2017
cater_habitat_1723 <- subset(cater_habitat, year!="2014")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2015")
cater_habitat_1723 <- subset(cater_habitat_1723, year!="2016")
 
prior5<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

## Ran as three chains for efficiency

Mass23_1<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2),
                random=~us(1+datescaled):tree.species + us(1+datescaled):site + year + siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units,
                  family="cengaussian", data=cater_habitat_1723, prior=prior5, nitt=2500000, burnin=500000, pr=TRUE, thin=2000)
save(Mass23_1, file = "~/Mass23_1.RData")
plot(Mass23_1)

Mass23_2<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2),
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + year + siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units,
                  family="cengaussian", data=cater_habitat_1723, prior=prior5, nitt=2500000, burnin=500000, pr=TRUE, thin=2000)
save(Mass23_2, file = "~/Mass23_2.RData")
plot(Mass23_2)

Mass23_3<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2),
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + year + siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units,
                  family="cengaussian", data=cater_habitat_1723, prior=prior5, nitt=2500000, burnin=500000, pr=TRUE, thin=2000)
save(Mass23_3, file = "~/Mass23_3.RData")
plot(Mass23_3)
