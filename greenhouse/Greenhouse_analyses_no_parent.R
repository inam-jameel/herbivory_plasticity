######## PROJECT: Greenhouse experiment: Fitness and phenotypes in response to herbivory
#### PURPOSE:Examine fitness and traits in response to herbivory treatment, removed the data for the grandparental generation
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 5 apr 24

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
library(ggeffects) #for plotting
library(DHARMa)

# set working directory

setwd("~/OneDrive - University of Georgia/Inam_experiments/Herbivory_data/greenhouse/data")  #this is where you specify the folder where you have the data on your computer

#read in data 
gh2 <- read.csv("greenhouse_noparent_long_04Apr24.csv")

gh2 <- filter(gh2, Include == "yes" )

gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev_km<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)



# convert LAR  to proportion

#maternal leaf damage

gh2$mat_avgLAR<- gh2$mat_avgLAR/100
hist(gh2$mat_avgLAR) 

#summary
gh2$LAR_avg_prop<- gh2$avg_LAR/100
hist(gh2$LAR_avg_prop) 

gh2$LAR_max_prop<- gh2$max_LAR/100 
hist(gh2$LAR_max_prop) 


gh2$LAR_1_prop<- gh2$LAR_1/100 
hist(gh2$LAR_1_prop)

gh2$LAR_2_prop<- gh2$LAR_2/100 
hist(gh2$LAR_2_prop) 

gh2$LAR_3_prop<- gh2$LAR_3/100 
hist(gh2$LAR_3_prop) 


##variables from characters to factors

gh2$season<-as.factor(gh2$season)
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$exp_ID<-as.factor(gh2$exp_ID)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line

##Change the baseline for maternal treatment
gh2 $mat_treat <-factor(gh2 $mat_treat, levels = c("herb", "cont"))

##Change the baseline for offspring treatment
gh2 $treatment <-factor(gh2 $treatment, levels = c("Herbivorized", "Control"))


gh2 $treat_mat<-interaction(gh2 $treatment, gh2 $mat_treat,sep = "_")

gh2 $S_LAR_avg_prop<-scale(gh2$LAR_avg_prop,center=TRUE, scale=TRUE) #scale avg damage

#gh2$trichomes <- gh2$avg_trichomes +0.001 #added to do gamma glmer

cols=c("darkred","#56B4E9")

#*******************************************************************************
#### 1.SLA for greenhouse experiment #####
#*******************************************************************************

#### SLA  ####

foliar<-subset(gh2,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

#Box_plot
SLA_box <-ggplot(foliar, aes(x = treatment, y = SLA, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Offspring treatment")+ scale_y_continuous("Specific leaf area") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Herbivorized", "Control")) +  scale_fill_manual(values = c("darkred","#56B4E9"), name = "offspring treatment", labels = c("Included","excluded"))+facet_wrap(~mat_treat)


SLA_plot =ggplot(foliar, aes(x= elevation,y= SLA, color = treatment))+ geom_point(size=1)+ ggtitle("SLA by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_classic()+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)+facet_wrap(~mat_treat)
SLA_plot

##

sla_model <- glmmTMB(SLA ~ treatment*mat_treat*S_elev+season+(1|exp_ID)+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") #treatment is significant

#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


SLA_means<-emmeans(sla_model, ~ treatment:mat_treat, type="response", adjust = "sidak")
cld(SLA_means, details=TRUE)



library(ggeffects)
pred_sla <- ggpredict(sla_model, terms = c("S_elev[all]","treatment","mat_treat"), type = "re", interval="confidence")

sla_cline <-plot(pred_sla, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = c("darkred","#56B4E9"), facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Specific leaf area")


sla_cline


#with avg damage/mat avg damage

sla_model <- glmmTMB(SLA ~ LAR_max_prop*mat_treat*S_elev+season+(1|exp_ID)+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") 

visreg2d(sla_model,"S_elev","LAR_avg_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")


sla_model_b <- glmmTMB(SLA ~ LAR_avg_prop*mat_avgLAR*S_elev+season+(1|exp_ID)+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))

Anova(sla_model_b, type = "III") 

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


plot(sla_model_b)

summary(SLA_RM_b)

visreg2d(SLA_RM_b,"elev","LAR_avg_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")


####


visreg(sla_model,"S_elev", by="treatment", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))




#### selection, repeated measures ####


#have to run after datafiles are made for SLA and the fitness components
SLA_test <- dplyr::select(foliar, SLA, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, mat_avgLAR)
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

LAR_data<- gh2 %>% pivot_longer(cols=c('LAR_1_prop', 'LAR_2_prop', 'LAR_3_prop'),
                                 names_to='census',
                                 values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, S_elev,elevation, genotype, treatment, mat_treat, exp_ID, S_init_size, block,season, census,mat_avgLAR)


LAR_data$census[LAR_data$census == "LAR_1_prop"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2_prop"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3_prop"] <-"3"

ggplot(LAR_data, aes(x= LAR))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

LAR_data$census<-as.factor(LAR_data$census)

LAR_data <- drop_na(LAR_data,LAR) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

LAR_data $y_beta<- (LAR_data $LAR*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(LAR_data $y_beta)

max(LAR_data $y_beta)

########################
#Box_plot
LAR_box <-ggplot(LAR_data, aes(x = treatment, y = y_beta, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivory Treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

LAR_box + theme_classic() + theme(text = element_text(size=14),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Herbivorized", "Control")) +  scale_fill_manual(values = c(cols), name = "Herbivory treatment", labels = c("Herbivorized","Control"))+facet_wrap(~mat_treat)


LAR_plot =ggplot(LAR_data, aes(x= elevation,y= y_beta, color = treatment))+ geom_point(size=1)+ ggtitle("SLA by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("LAR") +theme_classic()+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)+facet_wrap(~mat_treat+season)
LAR_plot



## model
LAR_Model<- gamlss (formula= y_beta ~S_elev*season+S_elev*treatment*mat_treat +random(census) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>%  
  mutate(fit.m = predict(LAR_Model, se.fit=FALSE),              
         resid = residuals(LAR_Model))

##Convert coefficients to probabilities
newdf2 $predicted<-1/(1+exp(-(newdf2 $fit.m))) 
newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


LAR_fig= ggplot(newdf2, aes(x= S_elev,y= predicted, group= treatment, colour= treatment))+geom_point(size=2,aes(x= S_elev,y= LAR)) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")  +scale_colour_manual(values = c(cols), name = "Herbivore treatment", labels = c("Herbivorized","Control")) +geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) + facet_wrap(~ mat_treat:season)




#Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*treatment*mat_treat+S_elev*treatment*season+S_elev* mat_treat*season+ treatment* mat_treat*season+ random(census) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)


LAR_Model_two<- gamlss (formula= y_beta ~S_elev*treatment+S_elev*mat_treat+ treatment* mat_treat+season*treatment+season*mat_treat+S_elev*season+ random(census) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)

LAR_Model_one<- gamlss (formula= y_beta ~S_elev+treatment+mat_treat+ season+ random(census) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)




# not great
#damage_model <- glmmTMB(y_beta ~S_elev*mat_treat*treatment+season + (1| census)+ (1|exp_ID)+ (1| block)+ (1| genotype), data=LAR_data, family=beta_family())

#simulationOutput <- simulateResiduals(fittedModel= damage_model, plot = T, re.form = NULL,allow.new.levels =T)

#plot(LAR_Model)
#summary(damage_model)







###########################

#Then, the analysis for treatment/mat_treat

library(betareg)
beta_modela<- betareg ( y_beta ~S_elev*mat_treat*treatment+census, data=LAR_data)
Anova(beta_modela,type="III")
visreg(beta_modela, "S_elev", scale="response")

visreg(beta_modela, 'S_elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 


summary(beta_modela)

plot(beta_modela)
gamlss_moda<- gamlss (formula= y_beta ~S_elev*mat_treat+census + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

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



#### treatment ####
hist(gh2$mature_length_siliques)

hurdle_Model<- glmmTMB(round(mature_length_siliques) ~treatment*S_elev*mat_treat + (1|block)+ (1| genotype), zi=~treatment*S_elev*mat_treat + (1|block)+ (1| genotype),data = gh2 ,family=truncated_nbinom2)

Anova(hurdle_Model,type="III", component="zi")
Anova(hurdle_Model,type="III", component="cond")


hurdle_Model_repeateda <- glmmTMB(mature_length_siliques ~treatment*S_elev*mat_treat+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=gh2, zi=~treatment*S_elev*mat_treat+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))


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


gh_pf= ggplot(gh2, aes(x= elevation,y= reproduced, group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Source Elevation")  
gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))

Prob_violin <-ggplot(gh2, aes(x = treatment, y = reproduced, fill = treatment)) +
  geom_violin() +xlab("Herbivory treatment")+ 
  stat_summary(fun.data = "mean_cl_normal")+
  scale_y_continuous("Probability of reproduction") +
  geom_point(pch = 21, position = position_jitterdodge())

repro_treatment<-Prob_violin + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~mat_treat)



####treatment####

mod_repro<-glmmTMB(reproduced~S_elev*treatment*mat_treat+season+(1|block)+(1|genotype)+(1|exp_ID),data=gh2,family=binomial(link="logit"))

Anova(mod_repro, type="III") #sig interaction btw mat_treat and elev

simulationOutput <- simulateResiduals(fittedModel= mod_repro, plot = T, re.form = NULL,allow.new.levels =T)





#### Avg LAR####

repro_data $s_LAR_avg_prop<-scale(repro_data$LAR_avg_prop,center=TRUE, scale=TRUE) #doesnt need to be scaled, but i thought it would be good to be consistant
repro_data $s_mat_avgLAR<-scale(repro_data$mat_avgLAR,center=TRUE, scale=TRUE)
#have to scale the maternal damage

mod_repro_REavgLARb<-glmmTMB(reproduced~S_elev*mat_avgLAR*LAR_avg_prop+season+(1|block)+(1|genotype)+(1|exp_ID),data=gh2,family=binomial(link="logit"))

Anova(mod_repro_REavgLARb,type="III") #sig interaction btw mat_treat and elev

pred_repro <- ggpredict(mod_repro_REavgLARb, terms = c("S_elev[all]","mat_avgLAR"), type = "re", interval="confidence")


test <-plot(pred_repro, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,1)


visreg2d(mod_repro_REavgLARb,"S_elev","mat_avgLAR", scale = "response", xlab="Source elevation (Km)", ylab="maternal damage (scaled)",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,1, by=.25))




#*******************************************************************************
#### 4.Fitness models: Fecundity ~ garden x treatment x elev #####
#*******************************************************************************

repro <- filter(gh2, mature_length_siliques > 0 )

gh_pf= ggplot(repro, aes(x= elevation,y= round(mature_length_siliques), group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Fecundity(summed silique length)")+ scale_x_continuous("Source Elevation")  

gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))



only_repro <- filter(fruit_data, fruit_length > 0 )

hist(only_repro$fruit_length)

#### treatment####
mod_fecundity_avgLARa <- glmmTMB (mature_length_siliques ~ treatment+S_elev+mat_treat+season  + (1| block)+(1|genotype)+(1|exp_ID), family=Gamma(link="log"), data = repro)

simulationOutput <- simulateResiduals(fittedModel= mod_fecundity_avgLARa, plot = T, re.form = NULL,allow.new.levels =T)


Anova(mod_fecundity_avgLARa,type="III")

summary(mod_fecundity_avgLARa)

#### avg damage ####
mod_fecundity_avgLARb <- glmmTMB (mature_length_siliques ~ LAR_avg_prop*S_elev*mat_avgLAR+season  + (1| block)+(1|genotype)+(1|exp_ID), family=Gamma(link="log"), data = repro)

simulationOutput <- simulateResiduals(fittedModel= mod_fecundity_avgLARb, plot = T, re.form = NULL,allow.new.levels =T)


Anova(mod_fecundity_avgLARb,type="III")

summary(mod_fecundity_avgLARb)

fecundmoda <-predictorEffect("elev",  partial.residuals=TRUE, mod_fecundity_avgLARb)
plot(fecundmoda, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(200,400))

library(ggeffects)

plot(ggpredict(mod_fecundity_avgLARb, terms = "elev"))

result <- ggpredict(mod_fecundity_avgLARb, c("elev","LAR_avg_prop","mat_avgLAR"))
plot(result)










