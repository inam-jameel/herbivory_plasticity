######## PROJECT: grasshopper experiment: variation in herbivore damage due to water availability 
#### PURPOSE:Examine fitness, traits in response to water availability and herbivory .
#### AUTHORS: Inam Jameel & Jill Anderson
#### DATE LAST MODIFIED: 14 feb 24

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

##setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")
setwd("~/Documents/personnel/Jameel/grasshopper")

 #this is where you specify the folder where you have the data on your computer


#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long.csv",stringsAsFactors=T)

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper$year<-as.factor(grasshopper$year)

##Change the baseline for Water Water
grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))



##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

#Elevational distance in km
grasshopper$elev_dist_km<-grasshopper $elev_dist/1000

#Adjust flowering time based on final_model in flowering_time_adjustment.R
grasshopper$FT_Adj<-round((grasshopper$Ordinal_Date_flowering - (grasshopper$Silique_length_flowering /2.5901)),0)

plot(grasshopper$FT_Adj~grasshopper$Ordinal_Date_flowering)

grasshopper$Snowmelt_FT_Adj<-grasshopper$FT_Adj-grasshopper$Day_of_snowmelt
plot(grasshopper$Snowmelt_FT_Adj ~grasshopper$Snowmelt_Date_flowering)

hist(grasshopper$Ordinal_Date_flowering)
hist(grasshopper$FT_Adj)

#Let's concatenate herbivore and watering Waters, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Water,sep = "_")

grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))

##Let's correlate rosette and bolt leaf data to see if we can come up with composite figures
plot(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
plot(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
plot(grasshopper$rosette_lwc~grasshopper$bolt_lwc)

mod1<-lm(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
summary(mod1)

LAR_Mod<-lm(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
summary(LAR_Mod)

mod3<-lm(grasshopper$rosette_lwc~grasshopper$bolt_lwc)
summary(mod3)

##Create composite leaf SLA, succuclence and lwc variables based on the regressions above. This gives us foliar trait data for the 67 plants for which we have bolt but not rosette collections
grasshopper$SLA <- ifelse(is.na(grasshopper$rosette_SLA), (71.8177 + 0.7319*grasshopper$bolt_SLA), grasshopper$rosette_SLA)
grasshopper$succulence <- ifelse(is.na(grasshopper$rosette_succulence), (0.0023948 + 0.5609208*grasshopper$bolt_succulence), grasshopper$rosette_succulence)
grasshopper$lwc <- ifelse(is.na(grasshopper$rosette_lwc), (0.19571 + 0.67021*grasshopper$bolt_lwc), grasshopper$rosette_lwc)

#This calculates flowering durationgrasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $FT_Adj)


##Remove plants inadvertently killed by experimentors or gophers
grasshopper_exclude <- filter(grasshopper, Exclude == "Include")


#*******************************************************************************
#### Treatment Water effect #####
#*******************************************************************************

#just looking at cage/block. no plants

VWC_gg<-ggplot(grasshopper, aes(x = Water, y = avg_vwc, fill = Water,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Watering Water")+ scale_y_continuous("Volumetric Water Content") +
  geom_point(pch = 21, position = position_jitterdodge())
VWC_gg + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                               panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Ample water","Water-restricted")) +  scale_fill_manual(values = c( "lightblue","#CC79A7"))


VWC_repeated <- lmer(avg_vwc~ Water*Herbivore+year + (1|Cage_Block), data=grasshopper )

visreg(VWC_repeated, 'Water', type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Volumetric Water Content", partial=TRUE) 


Anova(VWC_repeated, type = "III")


#pull the fitted values out and plot them in ggplot

vwc_df <- grasshopper %>% 
  
  mutate(fit.m = predict(VWC_repeated, re.form = NA),
         
         fit.c = predict(VWC_repeated, re.form = NULL), #all random effects
         
         resid = residuals(VWC_repeated))

vioplot(fit.c ~ Water, data= vwc_df, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,ylim=c(0,20),colMed = "black",colMed2 = c("lightblue","#CC79A7"), col=c("lightblue","#CC79A7"), ylab="Volumetric water content", xlab="Watering Water") +stripchart(fit.c ~ Water, data= vwc_df, col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)

ggplot(vwc_df, aes(x = Water, y = fit.c, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("#CC79A7", "lightblue")) +
  labs(y = "Volumetric Water Content") +
  labs(x = "Watering Water") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 


#*******************************************************************************
####  Survival #####
#*******************************************************************************
##survA <- glmmTMB(Season_survival ~ S_initdiam +Water*Herbivore * S_elev +Water*Herbivore * I(S_elev^2) +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data=grasshopper,family=binomial(link="logit"))
survA <- glmmTMB(Season_survival ~ S_initdiam +Water*Herbivore * S_elev +I(S_elev^2) +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data=grasshopper_exclude,family=binomial(link="logit"))
Anova(survA, type="III")

surv <- glmmTMB(Season_survival ~ S_initdiam+ Water*Herbivore * S_elev *year+Water*Herbivore * I(S_elev^2) +(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopper_exclude,family=binomial(link="logit"))
Anova(surv, type="III")

library(effects)
cols=c("#CC79A7","blue")
pred_FT <- ggpredict(surv, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Survival")
phen_plot


#*******************************************************************************
####  Probability of reproduction #####
#*******************************************************************************
grasshopperFT <- filter(grasshopper, year != "2021")

##reproA <- glmmTMB(Reproduced ~ S_initdiam+Water*Herbivore * S_elev +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
#reproA <- glmmTMB(Reproduced ~ S_initdiam+Water*Herbivore * S_elev+Water*Herbivore * I(S_elev^2) +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
reproA <- glmmTMB(Reproduced ~ S_initdiam+Water*Herbivore * S_elev+ I(S_elev^2) +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
Anova(reproA, type="III")

##repro <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
#repro <- glmmTMB(Reproduced ~ S_initdiam+Water*Herbivore * S_elev *year+Water*Herbivore * I(S_elev^2)+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
repro <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))
Anova(repro, type="III")

emmip(repro, Herbivore ~ year*Water, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Probability of reproduction")+ aes(shape = Water)+theme(legend.position = c(0.8, 0.8))+scale_color_manual(values=cols)+facet_grid(~year)


emmip(repro, Herbivore ~ year, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Probability of reproduction")+ theme(legend.position = c(0.8, 0.8))+scale_color_manual(values=cols)




#*******************************************************************************
####  Fecundity #####
#*******************************************************************************
reproduced<-subset(grasshopperFT, Reproduced=="1")

#fec <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*S_elev +I(S_elev ^2)*Water*Herbivore +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= reproduced,family=Gamma(link="log"))
fec <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*S_elev +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= reproduced,family=Gamma(link="log"))
Anova(fec,type="III")
simulationOutput <- simulateResiduals(fittedModel= fec, plot = T, re.form = NULL,allow.new.levels =TRUE)

emmip(fec, Herbivore ~ Water, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Fecundity")+ aes(shape = Water)+theme(legend.position = c(0.5, 0.8))+scale_color_manual(values=cols)

simulationOutput <- simulateResiduals(fittedModel= fec, plot = T, re.form = NULL,allow.new.levels =TRUE)



#fec2 <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*S_elev*year +I(S_elev ^2)*Water*Herbivore +(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= reproduced,family=Gamma(link="log"))
fec2 <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*S_elev*year  +(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= reproduced,family=Gamma(link="log"))
Anova(fec2,type="III")
simulationOutput <- simulateResiduals(fittedModel= fec2, plot = T, re.form = NULL,allow.new.levels =TRUE)




#*******************************************************************************
####  Hurdle #####
#*******************************************************************************

##hurdle_ModelLA <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore * S_elev +year+Water*Herbivore * I(S_elev^2)+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT, zi=~Water*Herbivore * S_elev+year+Water*Herbivore * I(S_elev^2)+(1|PlantID) + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_ModelLA <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore * S_elev +year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data=grasshopperFT, zi=~Water*Herbivore * S_elev +year +(1|PlantID) + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
summary(hurdle_ModelLA)
simulationOutput <- simulateResiduals(fittedModel= hurdle_ModelLA, plot = T, re.form = NULL,allow.new.levels =TRUE)

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_ModelLA,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_ModelLA,type="III", component="cond")


##hurdle_ModelLA_full <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore * S_elev*year+Water*Herbivore * I(S_elev^2)+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT, zi=~Water*Herbivore * S_elev*year+Water*Herbivore * I(S_elev^2)+(1|PlantID) + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_ModelLA_full <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore * S_elev*year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT, zi=~Water*Herbivore * S_elev*year +(1|PlantID) + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
summary(hurdle_ModelLA_full)

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_ModelLA_full,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_ModelLA_full,type="III", component="cond")
simulationOutput <- simulateResiduals(fittedModel= hurdle_ModelLA_full, plot = T, re.form = NULL,allow.new.levels =TRUE)


#*******************************************************************************
#### Lifetime Fitness #####
#*******************************************************************************

#read in data 
grasshopperLF <- read.csv("Grasshopper_lifetimefitness.csv", stringsAsFactors=T)

sapply(grasshopperLF,class)
##Change the baseline for Water
grasshopperLF $Water <-factor(grasshopperLF $Water, levels = c("Control", "watered"))
grasshopperLF$Cage <-as.factor(grasshopperLF$Cage)
grasshopperLF$Block <-as.factor(grasshopperLF$Block)


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopperLF $S_elev<-scale(grasshopperLF $elevation,center=TRUE, scale=TRUE)

#This rescales source elevation from meters to km
grasshopperLF$elev_km<-grasshopperLF $elevation/1000
#This calculates elevational distance
grasshopperLF$elev_dist<-grasshopperLF$elev_km - 2.890
hist(grasshopperLF$elev_dist)
grasshopperLF $S_initdiam<-scale(grasshopperLF $init.diam,center=TRUE, scale=TRUE)

#Let's concatenate herbivore and watering Waters, which is helpful for some models.
grasshopperLF $treat<-interaction(grasshopperLF$Herbivore, grasshopperLF$Water,sep = "_")
##Remove plants inadvertently killed by experimentors or gophers
grasshopperLF_exclude <- filter(grasshopperLF, Exclude == "Include")



##### Survival ####

#mod_surv_elevdist <-glmmTMB(Season_survival ~ S_initdiam + elev_dist*Water*Herbivore+ I(elev_dist^2)*Water+(1|Cage_Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))
mod_surv_elevdist <-glmmTMB(Season_survival ~ S_initdiam + elev_dist*Water*Herbivore+ I(elev_dist^2)*Water+(1|Cage_Block)+(1|Genotype),data=grasshopperLF,family=binomial(link="logit"))
Anova(mod_surv_elevdist,type="III")
simulationOutput <- simulateResiduals(fittedModel= mod_surv_elevdist, plot = T, re.form = NULL,allow.new.levels =TRUE)

pred_surv_elevdist <- ggpredict(mod_surv_elevdist, terms = c("elev_dist[all]", "Herbivore","Water"), type = "re", interval ="confidence")

plot(pred_surv_elevdist, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevational distance")+ scale_y_continuous("Probability of survival")

#mod_surv <-glmmTMB(Season_survival ~ S_initdiam +S_elev*Water*Herbivore+I(S_elev^2)*Water+(1|Cage:Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))
mod_surv <-glmmTMB(Season_survival ~ S_initdiam +S_elev*Water*Herbivore+I(S_elev^2)*Water+(1|Cage_Block)+(1|Genotype),data=grasshopperLF,family=binomial(link="logit"))
simulationOutput <- simulateResiduals(fittedModel= mod_surv, plot = T, re.form = NULL,allow.new.levels =TRUE)

Anova(mod_surv,type="III")

pred_surv <- ggpredict(mod_surv, terms = c("S_elev[all]", "Herbivore","Water"), type = "re", interval ="confidence")

plot(pred_surv, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevational distance")+ scale_y_continuous("Probability of survival")


##### Reproduction ####

#mod_repro_elevdist <-glmmTMB(reproduced ~ S_initdiam + elev_dist*Water*Herbivore+ I(elev_dist^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))

##mod_repro_elevdist <-glmmTMB(reproduced ~ S_initdiam + elev_dist*Water*Herbivore+ I(elev_dist^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF,family=binomial(link="logit"))

##mod_repro_elevdist <-glmmTMB(reproduced ~ S_initdiam + elev_dist*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))

mod_repro_elevdist <-glmmTMB(reproduced ~ S_initdiam + elev_dist*Water*Herbivore+ (1|Cage_Block)+(1|Genotype),data= grasshopperLF,family=binomial(link="logit"))
Anova(mod_repro_elevdist,type="III")
simulationOutput <- simulateResiduals(fittedModel= mod_repro_elevdist, plot = T, re.form = NULL,allow.new.levels =TRUE)

#mod_repro <-glmmTMB(reproduced ~ S_initdiam + S_elev*Water*Herbivore+ I(S_elev ^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))

##mod_repro <-glmmTMB(reproduced ~ S_initdiam + S_elev*Water*Herbivore+ I(S_elev ^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF,family=binomial(link="logit"))

mod_repro <-glmmTMB(reproduced ~ S_initdiam + S_elev*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= grasshopperLF_exclude,family=binomial(link="logit"))

mod_repro <-glmmTMB(reproduced ~ S_initdiam + S_elev*Water*Herbivore+  (1|Cage_Block)+(1|Genotype),data= grasshopperLF,family=binomial(link="logit"))


##### Fecundity, amongst individuals that reproduced ####
repro <-subset(grasshopperLF_exclude, reproduced=="1")

##fecund_elevdist <-glmmTMB(total_silique_length ~ S_initdiam + elev_dist*Water*Herbivore+ I(elev_dist^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= repro,family=Gamma(link="log"))

fecund_elevdist <-glmmTMB(total_silique_length ~ S_initdiam + elev_dist*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= repro,family=Gamma(link="log"))
Anova(fecund_elevdist,type="III")
simulationOutput <- simulateResiduals(fittedModel= fecund_elevdist, plot = T, re.form = NULL,allow.new.levels =TRUE)


emmip(fecund_elevdist, Water ~ Herbivore, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Fecundity")+ aes(shape = Herbivore)+theme(legend.position = c(0.5, 0.8))+scale_color_manual(values=cols)


##fecundity <-glmmTMB(total_silique_length ~ S_initdiam + S_elev*Water*Herbivore+ I(S_elev ^2)*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= repro,family=Gamma(link="log"))

fecundity <-glmmTMB(total_silique_length ~ S_initdiam + S_elev*Water*Herbivore+(1|Cage_Block)+(1|Genotype),data= repro,family=Gamma(link="log"))
Anova(fecundity,type="III")
simulationOutput <- simulateResiduals(fittedModel= fecundity, plot = T, re.form = NULL,allow.new.levels =TRUE)

emmip(fecundity, Water ~ Herbivore, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Fecundity")+ aes(shape = Herbivore)+theme(legend.position = c(0.5, 0.8))+scale_color_manual(values=cols)




##### Hurdle model####
hurdle_Model <- glmmTMB(total_silique_length ~ S_initdiam +Water*S_elev*Herbivore + (1|Genotype)+ (1|Cage_Block), data=grasshopperLF_exclude, zi=~ S_initdiam +Water*S_elev*Herbivore + (1| Genotype)+ (1| Cage_Block),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model,type="III", component="cond")
simulationOutput <- simulateResiduals(fittedModel= hurdle_Model, plot = T, re.form = NULL,allow.new.levels =TRUE)



####HERE 
cols=c("#56B4E9","#D55E00")

model1_num <- brm(bf(total_silique_number ~ S_initdiam + S_elev*Water*Herbivore + (1|Genotype)+ (1|Cage_Block),hu ~ S_initdiam + S_elev*Water*Herbivore+ (1|Genotype)+ (1|Cage_Block)), data = grasshopperLF, family = hurdle_gamma())

summary(model1_num)
pred_model1_num <- ggpredict(model1_num, terms = c("S_elev[all]", "Water","Herbivore"))

plot(pred_model1_num, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Number of siliques")
result_num <- ggpredict(model1_num, c("Herbivore","Water"))
plot(result_num)


model1_length <- brm(bf(total_silique_length ~ S_initdiam + S_elev*Water*Herbivore + (1|Genotype)+ (1|Cage_Block),hu ~ S_initdiam + S_elev*Water*Herbivore+ (1|Genotype)+ (1|Cage_Block)), data = grasshopperLF, family = hurdle_gamma())

summary(model1_length)
pred_model1_length <- ggpredict(model1_length, terms = c("S_elev[all]", "Water","Herbivore"))

plot(pred_model1_length, show_residuals=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity (silique length)")
result_bayes <- ggpredict(model1_length, c("Herbivore","Water"))
plot(result_bayes)






#*******************************************************************************
#### 2.leaf damage across censuses #####
#*******************************************************************************
#####repeated measures with all damage data####


#reformat datafile

LAR_data<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3","LAR_4","LAR_5","LAR_6"),
                                        names_to='census',
                                        values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev, treat,census, year)

LAR_data$census[LAR_data$census == "LAR_1"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3"] <-"3"
LAR_data$census[LAR_data$census == "LAR_4"] <- "4"
LAR_data$census[LAR_data$census == "LAR_5"] <- "5"
LAR_data$census[LAR_data$census == "LAR_6"] <- "6"

LAR_data $census <-as.factor(LAR_data $census)

LAR_data $year <-as.factor(LAR_data $year)


LAR_data$LAR_prop<-LAR_data $LAR/100
hist(LAR_data$LAR_prop)

ggplot(LAR_data, aes(x= LAR_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

LAR_data <- drop_na(LAR_data,LAR_prop) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

LAR_data $y_beta<- (LAR_data $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(LAR_data $y_beta)

max(LAR_data $y_beta)



LAR_Mod<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(LAR_Mod)
drop1(LAR_Mod)


LAR_Model<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(LAR_Model)
drop1(LAR_Model)

##Full models with year crossed with other factors is slightly better
AIC(LAR_Mod, LAR_Model, c=True)


#Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+S_elev*Water*year+S_elev* Herbivore*year+ Water* Herbivore*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)

#Subsets of models for drop 1
LAR_Model_two<- gamlss (formula= y_beta ~S_elev*Water+S_elev*Herbivore+ Water* Herbivore+year*Water+year*Herbivore+S_elev*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)

LAR_Model_one<- gamlss (formula= y_beta ~S_elev+Water+Herbivore+ year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(LAR_Mod, re.form = NA),
         
         fit.c = predict(LAR_Mod, re.form = NULL), #all random effects
         
         resid = residuals(LAR_Mod))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


LAR_box <-ggplot(newdf2, aes(x = Herbivore, y = fit.c_trans, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Watering Water", labels = c("Water-restricted","Supplemental watering"))

LAR_fig= ggplot(newdf2, aes(x= elev_km,y= fit.c_trans, group= Water, 
                                   colour= Water))+geom_point(size=2) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#CC79A7","lightblue"), name = "Water Water", labels = c("Restricted","Supplemental"))



##With year


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(LAR_Model, re.form = NA),
         
         fit.c = predict(LAR_Model, re.form = NULL), #all random effects
         
         resid = residuals(LAR_Model))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


LAR_box <-ggplot(newdf2, aes(x = Herbivore, y = fit.c_trans, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Watering Water", labels = c("Water-restricted","Supplemental watering"))

LAR_fig= ggplot(newdf2, aes(x= elev_km,y= fit.c_trans, group= Water, 
                                   colour= Water))+geom_point(size=2) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#CC79A7","lightblue"), name = "Water Water", labels = c("Restricted","Supplemental"))
                            
                            
                            
##Try glmmTMB

lar_model <- glmmTMB( y_beta ~S_elev*Water*Herbivore+year + (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data, family=beta_family)
Anova(lar_model,type="III")
simulationOutput <- simulateResiduals(fittedModel= lar_model, plot = T, re.form = NULL,allow.new.levels =TRUE)

library(ggeffects)
cols=c("#CC79A7","blue")
pred_lar <- ggpredict(lar_model, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
pred_lar <-plot(pred_dam, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores (proportion)")
pred_lar

damage_model <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year + (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data, family=beta_family())
Anova(damage_model,type="III")
simulationOutput <- simulateResiduals(fittedModel= damage_model, plot = T, re.form = NULL,allow.new.levels =TRUE)

library(ggeffects)
cols=c("#CC79A7","blue")
pred_dam <- ggpredict(damage_model, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
dam_plot <-plot(pred_dam, show_residuals=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores (proportion)")
dam_plot


##REturn here
library(brms)


mode=get_prior(y_beta ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), LAR_data, family = beta())


zoib_model <- brm(
bf(LAR ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype),
  phi ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype),
  zoi ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype),
  coi ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), 
  family = zero_one_inflated_beta(), LAR_data, prior=prior, chains=4, iter=110000,warmup=10000,control = list(max_treedepth = 100)))
  
  
  model_beta_bayes <- brm(
  bf(y_beta ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype),
     phi ~  S_elev*Water*Herbivore*year+ (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype)),
  data = LAR_data,
  family = Beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234, 
  file = "model_beta_bayes"  # Save this so it doesn't have to always rerun
)


summary(model_beta_bayes)
library(ggeffects)
cols=c("#CC79A7","blue")
pred_dam <- ggpredict(model_beta_bayes, terms = c("S_elev[all]", "Water","Herbivore"), interval="confidence")
dam_plot <-plot(pred_dam, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores (proportion)")
dam_plot




#*******************************************************************************
#### Flowering time #####
#*******************************************************************************

#filter out the two plants that flowered in the first season. This step is likely not necessary because those plants didn't have flowering times.
grasshopperFT <- filter(grasshopper, year != "2021")

##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(grasshopperFT, Ordinal_Date_flowering!="NA")
hist(flowering$Ordinal_Date_flowering)
hist(flowering$Snowmelt_Date_flowering)
min(flowering$Snowmelt_Date_flowering)
# actual model

flowering_time<-glmmTMB(Snowmelt_Date_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

flowering_time<-glmmTMB(Snowmelt_FT_Adj ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

flowering_time<-glmmTMB(Snowmelt_FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

plot(flowering_time)

Anova(flowering_time, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel= flowering_time, plot = T, re.form = NULL,allow.new.levels =TRUE)

library(ggeffects)
cols=c("#CC79A7","blue")
pred_FT <- ggpredict(flowering_time, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_residuals=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (day since snowmelt)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")
phen_plot




##### ordinal floweringtime ####
flowering_time_ordinal<-glmmTMB(Ordinal_Date_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(flowering_time_ordinal, type = "III") # 
simulationOutput <- simulateResiduals(fittedModel= flowering_time_ordinal, plot = T, re.form = NULL)


FT_ordinal<-glmmTMB(Ordinal_Date_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

Anova(FT_ordinal, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= FT_ordinal, plot = T, re.form = NULL)


cols=c("#CC79A7","blue")
pred_FT <- ggpredict(FT_ordinal, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")


# Interactions with year
FT_interact<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(FT_interact, type = "III") # 
library(effects)
cols=c("#CC79A7","blue")
pred_FT <- ggpredict(FT_interact, terms = c("S_elev[all]", "Water","Herbivore", "year"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")
phen_plot


##### peak floweringtime ####
flowering_time_peak<-glmmTMB(Date_peak_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(flowering_time_peak, type = "III") # 
simulationOutput <- simulateResiduals(fittedModel= flowering_time_peak, plot = T, re.form = NULL)

cols=c("#CC79A7","blue")
pred_FT <- ggpredict(flowering_time_peak, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Peak flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")
phen_plot



#*******************************************************************************
#### Height at flowering #####
#*******************************************************************************

hist(flowering$Sum_height_flowering)

height<-glmmTMB(Sum_height_flowering ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)


height<-glmmTMB(Sum_height_flowering ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=Gamma(link="log"))

height<-glmmTMB(Sum_height_flowering ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(height, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= height, plot = T, re.form = NULL)

cols=c("#CC79A7","blue")
pred_height <- ggpredict(height, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
height_plot <-plot(pred_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Height at flowering (cm)")
height_plot


##Max_height_flowering
max_height<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))


max_height<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(max_height, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= max_height, plot = T, re.form = NULL)

cols=c("#CC79A7","blue")
pred_max_height <- ggpredict(max_height, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
max_height_plot <-plot(pred_max_height, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Maximum at flowering (cm)")
max_height_plot

#*******************************************************************************
#### Flowering duration #####
#*******************************************************************************
hist(flowering$flowering_duration)

flowering_duration<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)

Anova(flowering_duration, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)

library(ggeffects)
cols=c("#CC79A7","blue")
pred_FT <- ggpredict(flowering_duration, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_plot <-plot(pred_FT, show_residuals=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering duration (elapsed days)")
phen_plot


emmip(flowering_duration, Water ~ Herbivore, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Flowering duration (elapsed days)")+ aes(shape = Water)+theme(legend.position = c(0.2, 0.2))+scale_color_manual(values=cols)



##### All interactions ####
duration<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)
Anova(duration, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= duration, plot = T, re.form = NULL)


emmip(duration, Water ~ Herbivore, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ ylab("Flowering duration (elapsed days)")+ aes(shape = Water)+theme(legend.position = c(0.2, 0.2))+scale_color_manual(values=cols)

ggplot(flowering, aes(x = Herbivore, y = flowering_duration, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Water) +
  theme_light() +
  scale_fill_manual(values = c("ivory", "#117733")) +
  labs(y = "Duration of flowering") +
  labs(x = "Herbivore treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

ggplot(flowering, aes(x = Water, y = flowering_duration, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("ivory", "#117733")) +
  labs(y = "Duration of flowering") +
  labs(x = "Water treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 



#*******************************************************************************
#### Specific leaf area #####
#*******************************************************************************

###### SLA  ####
foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)



SLA_RM_box<-ggplot(grasshopper, aes(x = Herbivore, y = rosette_SLA, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Specific Leaf Area") +
  geom_point(pch = 21, position = position_jitterdodge())

SLA_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "#CC79A7","lightblue"), name = "Watering Water", labels = c("Water-restricted","Ample water"))


hist(foliar$SLA)
sla <- glmmTMB(SLA ~ Water*Herbivore*S_elev+year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family=lognormal(link="log"))

Anova(sla, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= sla, plot = T, re.form = NULL,allow.new.levels =T)

cols=c("#CC79A7","blue")
pred_sla <- ggpredict(sla, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
sla_plot <-plot(pred_sla, show_residuals=TRUE, show_title =TRUE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Specific leaf area")
sla_plot

sla_full <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))
Anova(sla_full, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= sla_full, plot = T, re.form = NULL,allow.new.levels =T)
cols=c("#CC79A7","blue")
pred_sla <- ggpredict(sla_full, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
sla_plot <-plot(pred_sla, show_residuals=TRUE, show_title =TRUE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Specific leaf area")
sla_plot

rosette_SLA <- glmmTMB(rosette_SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))
Anova(rosette_SLA, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= rosette_SLA, plot = T, re.form = NULL,allow.new.levels =T)
cols=c("#CC79A7","blue")
pred_sla <- ggpredict(rosette_SLA, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
sla_plot <-plot(pred_sla, show_residuals=TRUE, show_title =TRUE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Specific leaf area (rosette leaves only)")
sla_plot


#*******************************************************************************
#### Leaf water content #####
#*******************************************************************************
lwc_data <- dplyr::select(foliar, lwc, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,  year)
hist(lwc_data$lwc)

lwc_data <- drop_na(lwc_data,lwc) 

n<-nrow(lwc_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

lwc_data $y_beta<- (lwc_data $lwc*(n-1) + 0.5)/n

hist(lwc_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(lwc_data $y_beta)

max(lwc_data $y_beta)


LWC_Mod<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=lwc_data,control = gamlss.control(n.cyc = 500))
summary(LWC_Mod)
drop1(LWC_Mod)

cols=c("#CC79A7","blue")
pred_lwc <- ggpredict(LWC_Mod, terms = c("S_elev[all]", "Herbivore","Water"), type = "re", interval="confidence")
lwc_plot <-plot(pred_lwc, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Leaf water content (cm)")
lwc_plot

ggplot(lwc_data, aes(x = Water, y = lwc, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#CC79A7","lightblue")) +
  labs(y = "Leaf water content") +
  labs(x = "Water treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 


#pull the fitted values out and plot them in ggplot

newdf2 <- lwc_data %>% 
  
  mutate(fit.m = predict(LWC_Mod, re.form = NA),
         
         fit.c = predict(LWC_Mod, re.form = NULL), #all random effects
         
         resid = residuals(LWC_Mod))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


lwc_box <-ggplot(newdf2, aes(x = Herbivore, y = fit.c_trans, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Leaf water content") +
  geom_point(pch = 21, position = position_jitterdodge())

lwc_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Watering Water", labels = c("Water-restricted","Supplemental watering"))

lwc_fig= ggplot(newdf2, aes(x= elev_km,y= fit.c_trans, group= Water, 
                                   colour= Water))+geom_point(size=2) + scale_y_continuous("Leaf water content")+ scale_x_continuous("Source Elevation")  
lwc_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#CC79A7","lightblue"), name = "Water Water", labels = c("Restricted","Supplemental"))
                            
                            
                            
LWC_Mod_full<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=lwc_data,control = gamlss.control(n.cyc = 500))
summary(LWC_Mod_full)
drop1(LWC_Mod_full)

#*******************************************************************************
#### Succulence (fresh weight - dried weight)/leaf area. We calculated it in g/cm^2. I've multiplied by 1000 to have mg/cm^2 #####
#*******************************************************************************

##Check distribution first
hist(foliar$succulence)
min(foliar$succulence, na.rm=TRUE)

max(foliar$succulence, na.rm=TRUE)
foliar$succulence_mg<-foliar$succulence*1000
hist(foliar$succulence)
hist(foliar$succulence_mg)

succulence_data <- dplyr::select(foliar, succulence, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,  year, succulence_mg)


succulence_data <- drop_na(succulence_data,succulence) 

n<-nrow(succulence_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

succulence_data $y_beta<- (succulence_data $succulence*(n-1) + 0.5)/n

hist(succulence_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(succulence_data $y_beta)

max(succulence_data $y_beta)



Model_succulence<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=succulence_data,control = gamlss.control(n.cyc = 500))
summary(Model_succulence)
drop1(Model_succulence)
plot(Model_succulence)

Model_succulence_two<- gamlss (formula= y_beta ~Water*Herbivore+S_elev*Herbivore+S_elev*Water+year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=succulence_data,control = gamlss.control(n.cyc = 500))
drop1(Model_succulence_two)

Model_succulence_one<- gamlss (formula= y_beta ~Water+Herbivore+S_elev+year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=succulence_data,control = gamlss.control(n.cyc = 500))
drop1(Model_succulence_one)



cols=c("#CC79A7","blue")
pred_succulence <- ggpredict(Model_succulence, terms = c("S_elev[all]", "Herbivore","Water"), type = "re", interval="confidence")

pred_succulence <- ggpredict(Model_succulence, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
succulence_plot <-plot(pred_succulence, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Succulence")
succulence_plot

#pull the fitted values out and plot them in ggplot

newdf2 <- succulence_data %>% 
  
  mutate(fit.m = predict(Model_succulence, re.form = NA),
         
         fit.c = predict(Model_succulence, re.form = NULL), #all random effects
         
         resid = residuals(Model_succulence))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


succulence_box <-ggplot(newdf2, aes(x = Herbivore, y = fit.c_trans, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Succulence") +
  geom_point(pch = 21, position = position_jitterdodge())

succulence_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Watering Water", labels = c("Water-restricted","Supplemental watering"))

succulence_fig= ggplot(newdf2, aes(x= elev_km,y= fit.c_trans, group= Water, 
                                   colour= Water))+geom_point(size=2) + scale_y_continuous("Succulence")+ scale_x_continuous("Source Elevation")  
succulence_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#CC79A7","lightblue"), name = "Water Water", labels = c("Restricted","Supplemental"))
                            
    ggplot(succulence_data, aes(x = Water, y = succulence, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#CC79A7","lightblue")) +
  labs(y = "Leaf succulence") +
  labs(x = "Water treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 
  
  
 Model_succulence_full<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+  random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=succulence_data,control = gamlss.control(n.cyc = 500))
summary(Model_succulence_full)
drop1(Model_succulence_full)
 
  
##The residuals don't look great if I use a normal distribution. We can't use values of 0 with a lognormal distribution,so I've added a small number to succulence. 
  succulence_data$Suc_mg<- succulence_data$succulence_mg+0.01
succ <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev+year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))
                      
 Anova(succ, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= succ, plot = T, re.form = NULL,allow.new.levels =T)

cols=c("#CC79A7","blue")
pred_succ <- ggpredict(succ, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot_succ <-plot(pred_succ, show_residuals=TRUE, show_title =TRUE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Leaf succulence")+ylim(0,20)
plot_succ
            
succ_full <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))
Anova(succ_full, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= succ_full, plot = T, re.form = NULL,allow.new.levels =T)



    ggplot(succulence_data, aes(x = Herbivore, y = succulence, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("#F05039","#1F449C")) +
  labs(y = "Leaf succulence") +
  labs(x = "Herbivore treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

