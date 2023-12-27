######## PROJECT: Greenhouse experiment: Fitness and phenotypes in response to herbivory
#### PURPOSE:Examine fitness and traits in response to herbivory treatment, removed the data for the grandparental generation
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 18 Dec 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(tidyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for plotting
require(glmmTMB) # for running survival model, have to load twice
require(gamlss) # for running phenology model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(vioplot) #for violin plots
library(RColorBrewer) #for 3d visreg plots



# set working directory

setwd("~/OneDrive - University of Georgia/Inam_experiments/Herbivory_data/greenhouse/")  #this is where you specify the folder where you have the data on your computer

#read in data 
gh2 <- read.csv("gen2_GHSummary2022_no_parent.csv")

gh2 <- filter(gh2, Include == "yes" )

gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)



# convert LAR  to proportion

#maternal leaf damage

gh2$mat_avgLAR<- gh2$mat_avgLAR/100
hist(gh2$mat_avgLAR) 

#summary
gh2$LAR_avg_prop<- gh2$LAR_avg/100
hist(gh2$LAR_avg_prop) 

gh2$LAR_max_prop<- gh2$LAR_max/100 
hist(gh2$LAR_max_prop) 

gh2$LAR_avg_S1_prop<- gh2$LAR_avg_S1/100
hist(gh2$LAR_avg_S1_prop) 

gh2$LAR_max_S1_prop<- gh2$LAR_max_S1/100
hist(gh2$LAR_max_S1_prop)

gh2$LAR_avg_S2_prop<- gh2$LAR_avg_S2/100
hist(gh2$LAR_avg_S2_prop) 

gh2$LAR_max_S2_prop<- gh2$LAR_max_S2/100 
hist(gh2$LAR_max_S2_prop) 



#season 2
gh2$S2_LAR1_prop<- gh2$S2_LAR1/100
hist(gh2$S2_LAR1_prop) 

gh2$S2_LAR2_prop<- gh2$S2_LAR2/100 #most recent
hist(gh2$S2_LAR2_prop) 


#season 1

gh2$S1_LAR1_prop<- gh2$S1_LAR1/100 
hist(gh2$S1_LAR1_prop)

gh2$S1_LAR3_prop<- gh2$S1_LAR3/100 #most damage
hist(gh2$S1_LAR3_prop) 

gh2$S1_LAR5_prop<- gh2$S1_LAR5/100 
hist(gh2$S1_LAR5_prop) 

gh2$S1_LAR13_prop<- gh2$S1_LAR13/100 #last in S1
hist(gh2$S1_LAR13_prop) 

##variables from characters to factors
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$exp_ID<-as.factor(gh2$exp_ID)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line
gh2 $treat_mat<-interaction(gh2 $treatment, gh2 $mat_treat,sep = "_")
gh2 $S_LAR_avg_prop<-scale(gh2$LAR_avg_prop,center=TRUE, scale=TRUE)

#gh2$trichomes <- gh2$avg_trichomes +0.001 #added to do gamma glmer


#*******************************************************************************
#### 1.SLA for greenhouse experiment #####
#*******************************************************************************

#### SLA individual for full dataset ####

SLA_S1 =ggplot(gh2, aes(x= elevation,y= S1_SLA, color = treat_mat))+ geom_point(size=3.5)+ ggtitle("SLA season 1 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S1

SLA_S2 =ggplot(gh2, aes(x= elevation,y= S2_SLA, color = treat_mat))+ geom_point(size=3.5)+ ggtitle("SLA season 2 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S2

#converting wide to long for repeated measures

SLA_data<- gh2 %>% pivot_longer(cols=c('S2_SLA_include', 'S1_SLA_include'),
                 names_to='season',
                 values_to='SLA')

SLA_data <- dplyr::select(SLA_data, SLA, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, mat_avgLAR)
SLA_data$season[SLA_data$season == "S1_SLA_include"] <- "1"
SLA_data$season[SLA_data$season == "S2_SLA_include"] <- "2"
SLA_data$season<-as.numeric(SLA_data$season)

ggplot(SLA_data, aes(x= SLA))+ geom_histogram(color="black", fill="white")+ facet_grid(treatment ~  .)

#with treatment/mat_treat
SLA_RM_a <- lmer(SLA ~ treatment*mat_treat*elev+season + (1|block)+(1|genotype)+(1|exp_ID),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)

plot(SLA_RM_a)

Anova(SLA_RM_a, type = "III")

visreg(SLA_RM_a,"elev", by="treatment", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

#with avg damage/mat avg damage
SLA_RM_b <- lmer(SLA ~ LAR_avg_prop*mat_avgLAR*elev+season + (1|block)+(1|genotype)+(1|exp_ID),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)

plot(SLA_RM_b)

Anova(SLA_RM_b, type = "III")

summary(SLA_RM_b)

visreg2d(SLA_RM_b,"elev","LAR_avg_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")



#Testing the random effect of exp_ID

SLAnogeno <- lm(SLA ~ treatment*mat_treat*elev+year 
                 # +(1|exp_ID)
                     ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)


anova(SLA_RM, SLAnogeno)


#### selection, repeated measures ####


#have to run after datafiles are made for SLA and the fitness components
SLA_test <- dplyr::select(SLA_data, SLA, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, mat_avgLAR)
repro_test <- dplyr::select(repro_data, reproduced, exp_ID,season)
fruit_test<- dplyr::select(fruit_data, fruit_length, exp_ID,season)

fit_test <- merge(repro_test, fruit_test, by=c("exp_ID","season"))
selection <- merge(SLA_test, fit_test, by=c("exp_ID","season"))

selection $s_SLA<-scale(selection$SLA,center=TRUE, scale=TRUE)


hurdle_Model_repeated_SLAa <- glmmTMB(fruit_length ~treatment*s_SLA+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=selection, zi=~treatment*s_SLA+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows significiant source elevation.
Anova(hurdle_Model_repeated_SLAa,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_repeated_SLAa,type="III", component="cond")

require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model_repeated_SLAa)

plotQQunif(hurdle_Model_repeated_SLAa)
plotResiduals(hurdle_Model_repeated_SLAa)

summary(hurdle_Model_repeated_SLAa)


library(ggeffects)

plot(ggpredict(hurdle_Model_repeated_SLAa, terms = "elev"))

result <- ggpredict(hurdle_Model_repeated_SLAa, c("s_SLA","treatment"))
plot(result)

#damage

hurdle_Model_repeated_SLAa <- glmmTMB(fruit_length ~LAR_avg_prop*s_SLA+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=selection, zi=~LAR_avg_prop*s_SLA+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows significiant source elevation.
Anova(hurdle_Model_repeated_SLAa,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_repeated_SLAa,type="III", component="cond")

require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model_repeated_SLAa)

plotQQunif(hurdle_Model_repeated_SLAa)
plotResiduals(hurdle_Model_repeated_SLAa)

summary(hurdle_Model_repeated_SLAa)


library(ggeffects)

plot(ggpredict(hurdle_Model_repeated_SLAa, terms = "elev"))

result <- ggpredict(hurdle_Model_repeated_SLAa, c("s_SLA","LAR_avg_prop"))
plot(result)




#### individual years ####


SLA_data_S1 $s_LAR_avg_S1_prop<-scale(SLA_data_S1$LAR_avg_S1_prop,center=TRUE, scale=TRUE)
SLA_data_S1 $s_S1_SLA_include<-scale(SLA_data_S1$S1_SLA_include,center=TRUE, scale=TRUE)


#looking at season 1
#SLA for first season is S1_SLA
#also have to filter for the indiv that we have appropriate values for SLA 
SLA_data_S1=filter(gh2,S1_include_SLA == 1)


repro_SLA_S1_a <-glmmTMB(S1_reproduced~treatment*S1_SLA_include*mat_treat +(1|block)+(1|genotype),data=SLA_data_S1,family=binomial(link="logit"))

Anova(repro_SLA_S1_a,type="III") 

visreg(repro_SLA_S1_a,"S1_SLA_include",  scale = "response", xlab="SLA", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))


repro_SLA_S1_b <-glmmTMB(S1_reproduced~LAR_avg_S1_prop*S1_SLA_include*mat_avgLAR +(1|block)+(1|genotype),data=SLA_data_S1,family=binomial(link="logit"))

Anova(repro_SLA_S1_b,type="III") 

visreg(repro_SLA_S1,"S1_SLA_include", by="LAR_avg_S1_prop", overlay=FALSE,  scale = "response", xlab="SLA", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

#looking at season 2
#SLA for second season is S2_SLA
#also have to filter for the indiv that we have appropriate values for SLA 
SLA_data_S2=filter(gh2,S2_include_SLA == 1)


SLA_data_S2 $s_LAR_avg_S2_prop<-scale(SLA_data_S2$LAR_avg_S2_prop,center=TRUE, scale=TRUE)
SLA_data_S2 $s_S2_SLA_include<-scale(SLA_data_S2$S2_SLA_include,center=TRUE, scale=TRUE)

repro_SLA_S2_a <-glmmTMB(S2_reproduced~treatment*mat_treat*S2_SLA_include +(1|block)+(1|genotype),data=SLA_data_S2,family=binomial(link="logit"))

Anova(repro_SLA_S2_a,type="III") 

visreg(repro_SLA_S2_a,"S2_SLA_include",  scale = "response", xlab="SLA", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))


repro_SLA_S2_b <-glmmTMB(S2_reproduced~LAR_avg_S2_prop*mat_avgLAR*S2_SLA_include +(1|block)+(1|genotype),data=SLA_data_S2,family=binomial(link="logit"))

Anova(repro_SLA_S2_b,type="III") 

visreg(repro_SLA_S2_b,"S2_SLA_include",  scale = "response", xlab="SLA", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))


#### linear regression of fecundity as a fitness component ####

repro_dataS1 <- filter(SLA_data_S1, S1_Mature_length_siliques > 0 )
repro_dataS2 <- filter(SLA_data_S2, S2_Mature_length_siliques > 0 )



#looking at season 1

repro_dataS1 $S_S1_SLA_include<-scale(repro_dataS1$S1_SLA_include,center=TRUE, scale=TRUE)

fecundity_SLA_S1a <-glmmTMB(S1_Mature_length_siliques~treatment*S_LAR_avg_S1_prop*mat_treat+(1|block)+(1|genotype),data=repro_dataS1,family=ziGamma(link="log"))

Anova(fecundity_SLA_S1a,type="III")

summary(fecundity_SLA_S1a)

fecundSLA <-predictorEffect("S_LAR_avg_S1_prop",  partial.residuals=TRUE, fecundity_SLA_S1a)
plot(fecundSLA, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"),
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,2000))


repro_dataS1 $mat_avgLAR<-scale(repro_dataS1$mat_avgLAR,center=TRUE, scale=TRUE)
repro_dataS1 $S_LAR_avg_S1_prop<-scale(repro_dataS1$LAR_avg_S1_prop,center=TRUE, scale=TRUE)

fecundity_SLA_S1b <-glmmTMB(S1_Mature_length_siliques~LAR_avg_S1_prop*S_S1_SLA_include*mat_avgLAR +(1|block)+(1|genotype),data=repro_dataS1,family=ziGamma(link="log"))

Anova(fecundity_SLA_S1b,type="III")

fecundSLA <-predictorEffect("mat_avgLAR_pzero",  partial.residuals=TRUE, fecundity_SLA_S1b)
plot(fecundSLA, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(200,300))



#looking at season 2
repro_dataS2 $S_S2_SLA_include<-scale(repro_dataS2$S2_SLA_include,center=TRUE, scale=TRUE)

fecundity_SLA_S2 <-glmmTMB(S2_Mature_length_siliques~treatment*S_S2_SLA_include*mat_treat+(1|block)+(1|genotype),data=repro_dataS2,family=ziGamma(link="log"))

Anova(fecundity_SLA_S2,type="III")

summary(fecundity_SLA_S2)

fecundSLA_S2 <-predictorEffect("S_S2_SLA_include",  partial.residuals=TRUE, fecundity_SLA_S2)
plot(fecundSLA_S2, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1000))



fecundity_SLA_S2 <-glmmTMB(S2_Mature_length_siliques~LAR_avg_S2_prop*S_S2_SLA_include*mat_avgLAR+(1|block)+(1|genotype),data=repro_dataS2,family=ziGamma(link="log"))

Anova(fecundity_SLA_S2,type="III")

summary(fecundity_SLA_S2)

fecundSLA_S2 <-predictorEffect("S2_SLA_include",  partial.residuals=TRUE, fecundity_SLA_S2)
plot(fecundSLA_S2, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(200,1000))



#*******************************************************************************
#### 2.LAR for greenhouse experiment #####
#*******************************************************************************


herb <- subset(gh2, treatment == "Herbivorized")

#repeated measures 

LAR_data<- herb %>% pivot_longer(cols=c('S1_LAR1_prop', 'S1_LAR3_prop', 'S1_LAR13_prop','S2_LAR1_prop', 'S2_LAR2_prop'),
                                 names_to='exposure',
                                 values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, exposure,mat_avgLAR)

LAR_data$census<-LAR_data$exposure
LAR_data$census[LAR_data$census == "S1_LAR1_prop"] <- "1"
LAR_data$census[LAR_data$census == "S1_LAR3_prop"] <- "2"
LAR_data$census[LAR_data$census == "S1_LAR13_prop"] <-"3"
LAR_data$census[LAR_data$census == "S2_LAR1_prop"] <- "4"
LAR_data$census[LAR_data$census == "S2_LAR2_prop"] <- "5"

ggplot(LAR_data, aes(x= LAR))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

LAR_data <- drop_na(LAR_data,LAR) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

LAR_data $y_beta<- (LAR_data $LAR*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(LAR_data $y_beta)

max(LAR_data $y_beta)



#Then, the analysis for treatment/mat_treat

library(betareg)
beta_modela<- betareg ( y_beta ~elev*mat_treat+census, data=LAR_data)
Anova(beta_modela,type="III")
visreg(beta_modela, "elev", scale="response")

visreg(beta_modela, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 



summary(beta_modela)

plot(beta_modela)
gamlss_moda<- gamlss (formula= y_beta ~elev*mat_treat+census + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_moda)
plot(gamlss_moda)
Anova(gamlss_moda)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(gamlss_moda, re.form = NA),
         
         fit.c = predict(gamlss_moda, re.form = NULL), #all random effects

         resid = residuals(gamlss_moda))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 



LAR_fig =ggplot(newdf2,aes(x= elev,y= LAR,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  geom_point(aes(shape= mat_treat),size=4)+scale_shape_manual(values = c(0,17,20)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))

LAR_fig




##And if you wanted to plot the partial residuals (instead of raw data) – again, you could customize this to your heart’s delight:


predictedProb_fig =ggplot(newdf2,aes(x= elev,y= fit.m_trans + resid_trans,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
 geom_smooth(method = "lm", se = TRUE,formula=y~x) + 
  
  theme_bw()  

predictedProb_fig

##### Then, the analysis for avg damage/mat avg damage #####

beta_modelb<- betareg ( y_beta ~elev*mat_avgLAR+census, data=LAR_data)
Anova(beta_modelb,type="III")
visreg(beta_modelb, "elev", scale="response")

visreg(beta_modelb, 'elev', by= "mat_avgLAR", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 



plot(beta_modelb)
gamlss_modb<- gamlss (formula= y_beta ~elev*mat_avgLAR+census + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_modb)
plot(gamlss_mod)



#pull the fitted values out and plot them in ggplot

newdf3 <- LAR_data %>% 
  
  mutate(fit.m = predict(gamlss_modb, re.form = NA),
         
         fit.c = predict(gamlss_modb, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_modb))

##Convert fit.m and fit.c back to the proportional scale

newdf3 $fit.m_trans<-1/(1+exp(-(newdf3 $fit.m))) 

newdf3 $fit.c_trans<-1/(1+exp(-(newdf3 $fit.m))) 

newdf3 $resid_trans<-1/(1+exp(-(newdf3 $resid))) 



LAR_fig =ggplot(newdf3,aes(x= elev,y= LAR,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  geom_point(aes(shape= mat_treat),size=4)+scale_shape_manual(values = c(0,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))

LAR_fig




##And if you wanted to plot the partial residuals (instead of raw data) – again, you could customize this to your heart’s delight:


predictedProb_fig =ggplot(newdf3,aes(x= elev,y= fit.m_trans + resid_trans,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  geom_smooth(method = "lm", se = TRUE,formula=y~x) + 
  
  theme_bw()  

predictedProb_fig





#*******************************************************************************
#### 3.Fitness models: hurdle ~ mat treatment x treatment x elev #####
#*******************************************************************************


#repeated measures
fruit_data<- gh2 %>% pivot_longer(cols=c('S1_Mature_length_siliques', 'S2_Mature_length_siliques'),
                                  names_to='season',
                                  values_to='fruit_length')

fruit_data <- dplyr::select(fruit_data, fruit_length, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, LAR_max_prop,mat_avgLAR)
fruit_data$season[fruit_data$season == "S1_Mature_length_siliques"] <- "1"
fruit_data$season[fruit_data$season == "S2_Mature_length_siliques"] <- "2"
fruit_data$season<-as.numeric(fruit_data$season)



#### treatment ####
hist(gh2$Overall_Mature_length_siliques)

hurdle_Model<- glmmTMB(round(Overall_Mature_length_siliques) ~treatment*elev*mat_treat + (1|block)+ (1| genotype), zi=~treatment*elev*mat_treat + (1|block)+ (1| genotype),data = gh2 ,family=truncated_nbinom2)

Anova(hurdle_Model,type="III", component="zi")
Anova(hurdle_Model,type="III", component="cond")



hurdle_Model_repeateda <- glmmTMB(fruit_length ~treatment*elev*mat_treat+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=fruit_data, zi=~treatment*elev*mat_treat+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))


##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows significiant source elevation.
Anova(hurdle_Model_repeateda,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_repeateda,type="III", component="cond")

require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model_repeateda)

check_zeroinflation(hurdle_Model_repeateda)

plotQQunif(hurdle_Model_repeateda)
plotResiduals(hurdle_Model_repeateda)

summary(hurdle_Model_repeateda)

# predictor effects is not working with zero inflated due to gamma?
#zeroinflated <-predictorEffect("elev",  partial.residuals=FALSE, hurdle_Model_repeated,component="zi")
#plot(zeroinflated, lwd=2,xlab="Source elevation (km)", ylab="Fecundity", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1500))

library(ggeffects)

plot(ggpredict(hurdle_Model_repeateda, terms = "elev"))

result <- ggpredict(hurdle_Model_repeateda, c("elev","mat_avgLAR_pzero","LAR_avg_prop"))
plot(result)

# Extract the predicted probabilities of zero counts
predicted_probs_zero <- predict(hurdle_Model_repeateda, type = "response", re.form = NA)

# Create a data frame for plotting
plot_data <- data.frame(
  elev = fruit_data$elev,
  treatment = fruit_data$treatment,
  mat_treat = fruit_data$mat_treat,
  season = fruit_data$season,
  predicted_probs_zero = predicted_probs_zero
)

# Scatterplot with regression line
plot(plot_data$elev, plot_data$predicted_probs_zero, 
     xlab = "Covariate 1", ylab = "Predicted Probability of Zero Count",
     main = "Zero-Inflated Component Regression Plot")
abline(lm(predicted_probs_zero ~ elev+treatment+mat_treat+season, data = plot_data), col = "red")



# Manually subset to exclude cases with zero counts
positive_counts_indices <- fruit_data$fruit_length > 0
predicted_counts_positive <- predicted_probs_zero[positive_counts_indices]


# Create a data frame for plotting
plot_data_positive <- data.frame(
  elev = fruit_data$elev[positive_counts_indices],
  treatment = fruit_data$treatment[positive_counts_indices],
  mat_treat = fruit_data$mat_treat[positive_counts_indices],
  season = fruit_data$season[positive_counts_indices],
  predicted_counts_positive = predicted_counts_positive
)

# Scatterplot with regression line
plot(plot_data_positive$elev, plot_data_positive$predicted_counts_positive, 
     xlab = "Covariate 1", ylab = "Predicted Counts for Positive Counts",
     main = "Count Data Component Regression Plot")
abline(lm(predicted_counts_positive ~ elev*treatment*mat_treat+season, data = plot_data_positive), col = "blue")

##### avg LAR ####

hurdle_Modelb<- glmmTMB(round(Overall_Mature_length_siliques) ~LAR_avg_prop*mat_avgLAR + (1|block)+ (1| genotype), zi=~LAR_avg_prop*elev*mat_avgLAR_pzero + (1|block)+ (1| genotype),data = gh2 ,family=truncated_nbinom2)

hurdle_Model_repeatedb <- glmmTMB(fruit_length ~LAR_avg_prop*elev+mat_avgLAR+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=fruit_data, zi=~LAR_avg_prop*elev+mat_avgLAR+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))


##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows a signifciant source elevation
Anova(hurdle_Model_repeatedb,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced).
Anova(hurdle_Model_repeatedb,type="III", component="cond")

diagnose(hurdle_Model_repeatedb)

testDispersion(hurdle_Model_repeatedb)
plotQQunif(hurdle_Model_repeatedb)
plotResiduals(hurdle_Model_repeatedb)

summary(hurdle_Model_repeatedb)


library(ggeffects)

plot(ggpredict(hurdle_Model_repeated, terms = "mat_avgLAR_pzero"))

result <- ggpredict(hurdle_Model_repeatedb, c("elev"))
plot(result)


#*******************************************************************************
#### 3.Fitness models: Probability of Reproduction ~ mat treatment x treatment x elev #####
#*******************************************************************************

# Probability of reproduction 


gh_pf= ggplot(gh2, aes(x= elev,y= Overall_flowered, group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Probability of flowering")+ scale_x_continuous("Source Elevation")  
gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))


#repeated measures
repro_data<- gh2 %>% pivot_longer(cols=c('S1_reproduced', 'S2_reproduced'),
                                names_to='season',
                                values_to='reproduced')

repro_data <- dplyr::select(repro_data, reproduced, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, LAR_max_prop,mat_avgLAR)
repro_data$season[repro_data$season == "S1_reproduced"] <- "1"
repro_data$season[repro_data$season == "S2_reproduced"] <- "2"
repro_data$season<-as.numeric(repro_data$season)


####treatment####

mod_repro_REavgLARa<-glmmTMB(reproduced~elev*treatment*mat_treat+season+(1|block)+(1|genotype)+(1|exp_ID),data=repro_data,family=binomial(link="logit"))

Anova(mod_repro_REavgLARa, type="III") #sig interaction btw mat_treat and elev

visreg(mod_repro_REavgLARa,"elev", by="mat_treat", overlay=FALSE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg(mod_repro_REavgLAR, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Probability of flowering", partial=FALSE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255"))) 

#### Avg LAR####

repro_data $s_LAR_avg_prop<-scale(repro_data$LAR_avg_prop,center=TRUE, scale=TRUE) #doesnt need to be scaled, but i thought it would be good to be consistant
repro_data $s_mat_avgLAR<-scale(repro_data$mat_avgLAR,center=TRUE, scale=TRUE)
#have to scale the maternal damage

mod_repro_REavgLARb<-glmmTMB(reproduced~elev*s_mat_avgLAR*s_LAR_avg_prop+season+(1|block)+(1|genotype)+(1|exp_ID),data=repro_data,family=binomial(link="logit"))

Anova(mod_repro_REavgLARb,type="III") #sig interaction btw mat_treat and elev

visreg(mod_repro_REavgLAR,"elev", by="mat_avgLAR_pzero", overlay=FALSE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(mod_repro_REavgLARb,"elev","s_mat_avgLAR", scale = "response", xlab="Source elevation (Km)", ylab="maternal damage (scaled)",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,1, by=.25))




#*******************************************************************************
#### 4.Fitness models: Fecundity ~ garden x treatment x elev #####
#*******************************************************************************

repro <- filter(gh2, Overall_Mature_length_siliques > 0 )

gh_pf= ggplot(repro, aes(x= elev,y= round(Overall_Mature_length_siliques), group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Fecundity(summed silique length)")+ scale_x_continuous("Source Elevation")  

gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))

vioplot(round(Overall_Mature_length_siliques) ~ treatment, data= repro, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255","grey"), col=c("#6699cc","#882255","grey"), ylab="Fecundity(summed silique length)", xlab="Herbivore treatment")+stripchart(round(Overall_Mature_length_siliques) ~ mat_treat, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)


only_repro <- filter(fruit_data, fruit_length > 0 )

hist(only_repro$fruit_length)

#### treatment####
mod_fecundity_avgLARa <- glmmTMB (fruit_length ~ treatment*elev*mat_treat+season  + (1| block)+(1|genotype)+(1|exp_ID), family=Gamma(link="log"), data = only_repro)


Anova(mod_fecundity_avgLARa,type="III")

summary(mod_fecundity_avgLARa)

#### avg damage ####
mod_fecundity_avgLARb <- glmmTMB (fruit_length ~ LAR_avg_prop*elev*mat_avgLAR+season  + (1| block)+(1|genotype)+(1|exp_ID), family=Gamma(link="log"), data = only_repro)

Anova(mod_fecundity_avgLARb,type="III")

summary(mod_fecundity_avgLARb)

fecundmoda <-predictorEffect("elev",  partial.residuals=TRUE, mod_fecundity_avgLARb)
plot(fecundmoda, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(200,400))

library(ggeffects)

plot(ggpredict(mod_fecundity_avgLARb, terms = "elev"))

result <- ggpredict(mod_fecundity_avgLARb, c("elev","LAR_avg_prop","mat_avgLAR"))
plot(result)










