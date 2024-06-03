######## PROJECT: greenhouse experiment: variation in herbivore damage due to treatment
#### PURPOSE:Examine fitness, traits in response to herbivory.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 06 May 24



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
library(DHARMa)
library(ggeffects)
library(ggpubr)

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/greenhouse/data")

#this is where you specify the folder where you have the data on your computer

#read in data 
greenhouse <- read.csv("2gen_all_season_long_050324.csv",stringsAsFactors=T)
greenhouse <- filter(greenhouse, Include == "yes" )
sapply(greenhouse,class)
##Some  variables are being read as characters not factors. Let's fix that
greenhouse$block<-as.factor(greenhouse$block)
greenhouse$Season<-as.factor(greenhouse$Season)
greenhouse$exp_ID<-as.factor(greenhouse$exp_ID)



##Change the baseline for offspring
greenhouse $treatment <-factor(greenhouse $treatment, levels = c("Herbivory", "Naïve"))

##Change the baseline for maternal treatment
greenhouse $mat_treat <-factor(greenhouse $mat_treat, levels = c("Herbivory", "Naïve","Parent"))

#if want to remove the parent maternal environment
greenhouse <- filter(greenhouse, mat_treat != "Parent")


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
greenhouse $S_elev<-scale(greenhouse $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
greenhouse $S_initdiam<-scale(greenhouse $ini_size,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
greenhouse$elev_km<-greenhouse $elevation/1000



#for vegetative traits
foliar<-subset(greenhouse,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

#set colors
cols=c("#882255","#56B4E9")
#cols=c("#CC79A7","lightblue","gray")
#cols2=c("#882255","#77DDAA")




#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************

#reformat datafile


LAR_data<- greenhouse %>% pivot_longer(cols=c("LAR_1","LAR_2"
                                              ,"LAR_3"
                                              ),
                                        names_to='census',
                                        values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elevation, genotype, mat_treat, treatment, block, exp_ID, ini_size, S_initdiam, elev_km, S_elev,census, Season)

LAR_data$census[LAR_data$census == "LAR_1"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3"] <- "3"


LAR_data $census <-as.factor(LAR_data $census)

LAR_data $Season <-as.factor(LAR_data $Season)

LAR_data$censusSeason <- interaction(LAR_data$census, LAR_data$Season)


LAR_data$LAR_prop<-LAR_data $LAR/100
hist(LAR_data$LAR_prop)
LAR_data <- drop_na(LAR_data,LAR_prop) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

LAR_data $y_beta<- (LAR_data $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

min(LAR_data $y_beta)

max(LAR_data $y_beta)



#Box_plot
LAR_box <-ggplot(LAR_data, aes(x = treatment, y = LAR_prop, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Offsping herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Herbivory", "Naïve")) +  scale_fill_manual(values = cols, name = "Offspring herbivory treatment", labels = c("Herbivory","Naïve"))+facet_grid(~Season+mat_treat)




## model
LAR_Model<- gamlss (formula= y_beta ~S_elev*treatment*mat_treat*Season+ random(censusSeason) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))


LAR <-emmeans(LAR_Model, ~ treatment*mat_treat*Season, type="response", adjust = "sidak")


plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)


#glmmtmb model is not a good fit
#damage_model <- glmmTMB(y_beta ~S_elev*treatment*mat_treat*Season + (1| censusSeason)+ (1| exp_ID)+ (1| block)+ (1| genotype), data=LAR_data, family=beta_family())

#Anova(damage_model)
#simulationOutput <- simulateResiduals(fittedModel= damage_model, plot = T, re.form = NULL,allow.new.levels =TRUE) #not 



#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>%  
  mutate(fit.m = predict(LAR_Model, se.fit=FALSE),              
         resid = residuals(LAR_Model))

##Convert coefficients to probabilities
newdf2 $predicted<-1/(1+exp(-(newdf2 $fit.m))) 
newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 

##Box plot
LAR_box <-ggplot(LAR_data, aes(x = treatment, y = LAR_prop, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Offspring herbivore treatment")+ 
  scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

damage_treatment<-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~ Season+mat_treat)



#cline, plotting predicted data
damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted, group= treatment, 
                                 colour= treatment))+geom_point(pch = 20) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + facet_grid(mat_treat~ Season) + theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Water availability", labels = c("Restricted","Supplemental"))



#cline, plotting predicted data
#damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted))+geom_point(pch = 20) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + facet_grid(~ Season) + theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Water availability", labels = c("Restricted","Supplemental"))




#cline, plotting predicted data
#damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted))+geom_point(pch = 20) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + facet_grid(mat_treat~ Season) + theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Water availability", labels = c("Restricted","Supplemental"))




##Conversation of standardized elevation back to raw units
-1*sd(greenhouse $elevation,na.rm=TRUE)+mean(greenhouse $elevation,na.rm=TRUE)
0*sd(greenhouse $elevation,na.rm=TRUE)+mean(greenhouse $elevation,na.rm=TRUE)
1*sd(greenhouse $elevation,na.rm=TRUE)+mean(greenhouse $elevation,na.rm=TRUE)
2*sd(greenhouse $elevation,na.rm=TRUE)+mean(greenhouse $elevation,na.rm=TRUE)

#Figure2

library(ggpubr)
figure2 <- ggarrange(
  damage_cline,
  damage_treatment,
  
  #LAR_fecundity,
  
  labels = c("A", "B"),
  ncol = 1, nrow = 2)
figure2



#env concatenates water and herbivore, needed for slope calculations

LAR_data $env<-interaction(LAR_data $treatment, LAR_data $mat_treat)

LAR_data $env_season<-interaction(LAR_data $env, LAR_data $Season)

LAR_ModelE<- gamlss (formula= y_beta ~S_elev* env_season-1+ random(censusSeason) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

#just season x s_elev
LAR_ModelE<- gamlss (formula= y_beta ~S_elev* Season-1+ random(censusSeason) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(LAR_ModelE)

plot(LAR_Model)


#Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*treatment*mat_treat+S_elev*treatment*Season+S_elev* mat_treat*Season+ treatment* mat_treat*Season+ random(censusSeason) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)

#Subsets of models for drop 1
LAR_Model_two<- gamlss (formula= y_beta ~S_elev*treatment+S_elev*mat_treat+ treatment* mat_treat+Season*treatment+Season*mat_treat+S_elev*Season+ random(censusSeason) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)

LAR_Model_one<- gamlss (formula= y_beta ~S_elev+treatment+mat_treat+ Season+ random(census) + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)


#*******************************************************************************
#### Specific Leaf Area #####
#*******************************************************************************
#*


foliar<-subset(greenhouse,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

hist(foliar$SLA)


sla_model <- glmmTMB(SLA ~ treatment*mat_treat*S_elev+Season+(1|exp_ID)+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


SLA <-emmeans(sla_model, ~ treatment*mat_treat, type="response", adjust = "sidak")
cld(SLA, details=TRUE)
Anova(sla_model, type = "III") 


##Box plot
SLA_box <-ggplot(foliar, aes(x = treatment, y = SLA, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ 
  scale_y_continuous("Specific leaf area") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

SLA_treatment<-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~mat_treat)


#significance of random effects

sla_model_no_plantID <- glmmTMB(SLA ~ treatment*mat_treat*S_elev+Season+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))

sla_model_no_genotype <- glmmTMB(SLA ~ treatment*mat_treat*S_elev+Season+(1|exp_ID)+(1|block), data = foliar, family= lognormal(link="log"))

sla_model_no_block <- glmmTMB(SLA ~ treatment*mat_treat*S_elev+Season+(1|exp_ID)+(1|genotype), data = foliar, family= lognormal(link="log"))


anova(sla_model,sla_model_no_plantID)
anova(sla_model,sla_model_no_genotype)
anova(sla_model,sla_model_no_block)




#*******************************************************************************
#### Flowering phenology #####
#*******************************************************************************


##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(greenhouse, Date_flowering!="NA")

hist(flowering$Date_flowering_exp)


FT_model<-glmmTMB(Date_flowering_exp ~ S_elev*treatment*mat_treat*Season+(1|genotype)+(1|block)+(1| exp_ID),data= flowering,family=lognormal(link="log"))
Anova(FT_model, type = "III") # 

#Use the DHARMa package to examine the residuals

simulationOutput <- simulateResiduals(fittedModel= FT_model, plot = T, re.form = NULL,allow.new.levels =T)


#FT_model_1<-glmmTMB(Date_flowering_exp ~ S_elev*treatment*mat_treat*Season+(1|genotype)+(1|block)+(1| exp_ID),data= flowering,family=Gamma(link="sqrt"))
#Anova(FT_model_1, type = "III") # 

#Use the DHARMa package to examine the residuals, which are reasonable

#simulationOutput <- simulateResiduals(fittedModel= FT_model_1, plot = T, re.form = NULL,allow.new.levels =T)


Ft <-emmeans(FT_model, ~ mat_treat:treatment:Season, type="response", adjust = "sidak")
cld(Ft, details=TRUE)

##Box plot

##Box plot
ft_box <-ggplot(flowering, aes(x = treatment, y = Date_flowering_exp, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ 
  scale_y_continuous("Day of flowering") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

ft_treatment<-ft_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~Season)

#significance of random effects

FT_model_no_plantID <- glmmTMB(Date_flowering_exp ~ S_elev*treatment*mat_treat*Season+(1|genotype)+(1|block),data= flowering,family=lognormal(link="log"))

FT_model_no_genotype <- glmmTMB(Date_flowering_exp ~ S_elev*treatment*mat_treat*Season+(1|block)+(1| exp_ID),data= flowering,family=lognormal(link="log"))

FT_model_no_block <- glmmTMB(Date_flowering_exp ~ S_elev*treatment*mat_treat*Season+(1|genotype)+(1| exp_ID),data= flowering,family=lognormal(link="log"))


anova(FT_model,FT_model_no_plantID)
anova(FT_model,FT_model_no_genotype)
anova(FT_model,FT_model_no_block)



#*******************************************************************************
#### Tallest bolt at flowering #####
#*******************************************************************************

#This calculates flowering duration
#flowering$summed_height_flowering<-(flowering $Height1_flowering + flowering $Height2_flowering)


height_model<-glmmTMB(Height1_flowering ~ S_elev*treatment*mat_treat+Season+(1|genotype)+(1|block)+(1| exp_ID),data= flowering,family=lognormal(link="log"))

hist(flowering$Height1_flowering)

Anova(height_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= height_model, plot = T, re.form = NULL)


height <-emmeans(height_model, ~ treatment, type="response", adjust = "sidak")
cld(height, details=TRUE)

max_height_nogeno<-glmmTMB(Height1_flowering ~ S_elev*treatment*mat_treat+Season+(1|block)+(1| exp_ID),data= flowering,family=lognormal(link="log"))

max_height_noblock<-glmmTMB(Height1_flowering ~ S_elev*treatment*mat_treat+Season+(1|genotype)+(1| exp_ID),data= flowering,family=lognormal(link="log"))

max_height_nopid<-glmmTMB(Height1_flowering ~ S_elev*treatment*mat_treat+Season+(1|genotype)+(1|block),data= flowering,family=lognormal(link="log"))

anova(height_model,max_height_nogeno)
anova(height_model,max_height_noblock)
anova(height_model,max_height_nopid)



##Box plot

height_box <-ggplot(flowering, aes(x = treatment, y = Height1_flowering, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ 
  scale_y_continuous("Height of tallest stem at flowering") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

height_treatment<-height_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_grid(~Season)



#Figure3

library(ggpubr)
figure3 <- ggarrange(
  SLA_treatment,
  ft_treatment,
  height_treatment,
  
  #LAR_fecundity,
  
  labels = c("A", "B","C"),
  ncol = 1, nrow = 3)
figure3



#*******************************************************************************
####  Probability of reproduction #####
#*******************************************************************************
greenhouse2 <- filter(greenhouse, Season != "1")


prob_repro <- glmmTMB(Reproduction ~S_initdiam+ treatment * S_elev *mat_treat +Season+(1|exp_ID) + (1|block)+(1|genotype), data= greenhouse2,family=binomial(link="logit"))


Anova(prob_repro, type="III")
simulationOutput <- simulateResiduals(fittedModel= prob_repro, plot = T, re.form = NULL,allow.new.levels =TRUE)

testOutliers(simulationOutput, type = c(
  "bootstrap"), nBoot = 100, plot = T)

prob_repro_nogeno <-glmmTMB(Reproduction ~S_initdiam+ treatment * S_elev *mat_treat +Season+(1|exp_ID) + (1|block), data= greenhouse2,family=binomial(link="logit"))

prob_repro_nopid <- glmmTMB(Reproduction ~S_initdiam+ treatment * S_elev *mat_treat +Season+ (1|block)+(1|genotype), data= greenhouse2,family=binomial(link="logit"))

prob_repro_block <- glmmTMB(Reproduction ~S_initdiam+ treatment * S_elev *mat_treat +Season+(1|exp_ID)+(1|genotype), data= greenhouse2,family=binomial(link="logit"))


anova(prob_repro,prob_repro_nogeno)
anova(prob_repro,prob_repro_nopid)
anova(prob_repro,prob_repro_block)


repro<-emmeans(prob_repro, ~ treatment:mat_treat, type="response", adjust = "sidak")
cld(repro, details=TRUE)



##Box plot
Prob_violin <-ggplot(greenhouse, aes(x = treatment, y = Reproduction, fill = treatment)) +
  geom_violin() +xlab("Herbivore treatment")+ 
  stat_summary(fun.data = "mean_cl_normal")+
  scale_y_continuous("Probability of reproduction") +
  geom_point(pch = 21, position = position_jitterdodge())

repro_treatment<-Prob_violin + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~mat_treat)


#*******************************************************************************
####  Fecundity #####
#*******************************************************************************
#*
repro<-subset(greenhouse, Reproduction=="1")



#repro $Water1 <-factor(repro $Water, levels = c("Supplemental", "Restricted"))
#repro $Water <-factor(repro $Water, levels = c("Restricted","Supplemental"))


fecundity <- glmmTMB(Mature_length_siliques ~S_initdiam+treatment*mat_treat* S_elev+Season+(1|exp_ID) + (1|block)+(1|genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))



Anova(fecundity,type="III")

simulationOutput <- simulateResiduals(fittedModel= fecundity, plot = T, re.form = NULL,allow.new.levels =TRUE)

summary(fecundity)



Fecmeans<-emmeans(fecundity, ~ mat_treat:treatment , type="response", adjust = "sidak")
cld(Fecmeans, details=TRUE)

plot(Fecmeans)


#pull the fitted values out and plot them in ggplot

newdf2 <- repro %>%  
  mutate(fit.m = predict(fecundity, se.fit=FALSE),              
         resid = residuals(fecundity))

##Convert coefficients to probabilities
newdf2 $predicted<-exp((newdf2 $fit.m)) 
newdf2 $resid_trans<-exp((newdf2 $resid)) 


fecundity_nogeno <-glmmTMB(Mature_length_siliques ~S_initdiam+treatment*mat_treat* S_elev+Season+(1|exp_ID) + (1|block), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

fecundity_nopid <- glmmTMB(Mature_length_siliques ~S_initdiam+treatment*mat_treat* S_elev+Season + (1|block)+(1|genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

fecundity_block <- glmmTMB(Mature_length_siliques ~S_initdiam+treatment*mat_treat* S_elev+Season+(1|exp_ID) +(1|genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))


anova(fecundity,fecundity_nogeno)
anova(fecundity,fecundity_nopid)
anova(fecundity,fecundity_block)




##Box plot
fec_box <-ggplot(repro, aes(x = treatment, y = Mature_length_siliques, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ 
  scale_y_continuous("Fecundity (mature fruit length)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

fecundity_treatment<-fec_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~mat_treat)



#Figure4

library(ggpubr)
figure4 <- ggarrange(
  repro_treatment,
  fecundity_treatment,
  
  labels = c("A", "B"),
  ncol = 1, nrow = 2)
figure4


#################################################################################

#*******************************************************************************
#### Succulence #####
#*******************************************************************************

foliar$succulence_mg<-foliar$succulence*1000
foliar$Suc_mg<- foliar$succulence_mg+0.01

succulence_data <- dplyr::select(foliar, succulence, elevation, genotype, treatment, mat_treat, exp_ID, S_initdiam, block, elev_km, S_elev,  Season, Suc_mg)

succulence_data <- drop_na(succulence_data,succulence) 

succ_model <- glmmTMB(Suc_mg ~ treatment*mat_treat*S_elev+Season+(1|exp_ID)+(1|genotype)+(1|block), data = succulence_data, family=lognormal(link="log"))
Anova(succ_model, type = "III") # 

#Use the DHARMa package to examine the residuals
simulationOutput <- simulateResiduals(fittedModel= succ_model, plot = T, re.form = NULL,allow.new.levels =T)




succulence <-emmeans(succ_model, ~ treatment:mat_treat, type="response", adjust = "sidak")
cld(succulence, details=TRUE)



##Box plot

Suc_box <-ggplot(succulence_data, aes(x = treatment, y = Suc_mg, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("herbivore treatment")+ 
  #scale_y_continuous("Leaf succulence (mg)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

Suc_treatment<-Suc_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "herbivore treament", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~mat_treat)


#*******************************************************************************
####  Leaf number #####
#*******************************************************************************


foliar<-subset(greenhouse,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)


leafnumber_model <- glmmTMB(avg_leafnumber ~ treatment*mat_treat*S_elev+Season+(1|exp_ID)+(1|genotype)+(1|block), data = foliar, family= lognormal(link="log"))
Anova(leafnumber_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable


emtrends(leafnumber_model, c("mat_treat"), var = "S_elev")


simulationOutput <- simulateResiduals(fittedModel= leafnumber_model, plot = T, re.form = NULL,allow.new.levels =T)


pred_leaf <- ggpredict(leafnumber_model, terms = c("S_elev[all]","mat_treat"), type = "re", interval="confidence")


leaf_cline <-plot(pred_leaf, show_data=TRUE, show_title =FALSE, show_legend=FALSE, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Leaf number")

leaf_cline 




leaf <-emmeans(leafnumber_model, ~ treatment*mat_treat, type="response", adjust = "sidak")
cld(leaf, details=TRUE)

##Box plot

leaf_box <-ggplot(foliar, aes(x = treatment, y = avg_leafnumber, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Leaf succulence (mg)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

leaf_treatment<-leaf_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water availability", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~mat_treat)



