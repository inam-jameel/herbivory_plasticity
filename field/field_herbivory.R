#libraries
library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
packageVersion("lme4")
library(tidyr)
library(dplyr)
library(broom)
library(car)
library(MASS)
library(emmeans)
library(visreg)
library(fitdistrplus)
library(gamlss)
library(rjags)
library(gamlss)
library(effects)

setwd("~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/")

#### merging data ####
#schofield <- read.csv("Field2021_schofield.csv")
#s_data <- read.csv("schofield_plant_data.csv")

#schofield_a<- merge(schofield,s_data,by="plantID")

#write.csv(schofield_a,"~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/Field2021_schofield_merge.csv")

#g_2020 <- read.csv("Field2021_gothic2020.csv")
#g_data <- read.csv("gothic_2020_plant_data.csv")

#gothic_a<- merge(g_2020,g_data,by="plantID")

#write.csv(gothic_a,"~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/Field2021_gothic2020_merge.csv")


#### read in data ####
data<- read.csv("Field2021_merged.csv")


data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)

data$genotype<-as.factor(data$genotype)
data$treatment<-as.factor(data$Treatment)
data$block<-as.factor(data$block)
data$cohort<-as.factor(data$cohort)

# creating a df for LAR
LAR <- drop_na(data,LAR) #this removes the rows without LAR values

##convert LAR to proportion
LAR$LAR_prop<- LAR$LAR/100
##Check that the conversion worked
hist(LAR$LAR_prop)

LAR<-na.omit(LAR)
 
gamlss.model_1a<- gamlss (formula=LAR_prop~S_elev*treatment*garden+ random(block)+random(genotype), family=BEZI, data= LAR)
summary(gamlss.model_1a)
plot(gamlss.model_1a)

Anova(gamlss.model_1a,type="III") #error regarding aliased coefficients

LSmeans_LAR <- emmeans(gamlss.model_1a,~ S_elev:treatment:garden,pbkrtest.limit = 3612)
LSmeans_LAR
write.csv(LSmeans_germ,file="LSmeans_germ_21apr21.csv")

# use predictor effect to plot
#LAR_fig <-predictorEffect("S_elev",  partial.residuals=TRUE, gamlss.model_1a)
#plot(LAR_fig, lwd=2,xlab="Elevation of origin", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),      partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))
#does not work

#looking at LAR v elevation ggplot
plot(LAR$elevation, LAR$LAR, xlab="elev", ylab="LAR")

LAR = ggplot(LAR, aes(x= elevation,y=LAR))+geom_point(size=5) +scale_x_continuous("LAR")+ scale_y_continuous("elevation")
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6)
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(), legend.position = "bottom")+ geom_smooth(method="glm",size=1.6
                                                                                                      , formula=y~poly(x,2)) 


# zoib --------------------------------------------------------------------

#trying zero inflated model package "zoib"
#this was not working, cant specify random effects
#zoib(LAR ~ Treatment*elevation, x2021_cohort_061820, random = population) 
#zoib(LAR ~ Treatment*elevation |1,
#  data = x2021_cohort_061820, random = 1, EUID= x2021_cohort_061820$PlantID,
#  zero.inflation = TRUE, one.inflation = FALSE, joint = FALSE)



#zoib example
#eg3 <- zoib(Percentage ~ Grade*Gender+MedDays|1|Grade*Gender+MedDays|1,
#            data = AlcoholUse, random = 1, EUID= AlcoholUse$County,
#            zero.inflation = TRUE, one.inflation = FALSE, joint = FALSE,
#            n.iter=5000, n.thin=20, n.burn=1000)

