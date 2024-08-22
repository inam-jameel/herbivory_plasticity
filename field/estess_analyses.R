######## PROJECT: Field experiment: variation in herbivore damage due to treatment in Estess
#### PURPOSE:Examine fitness, traits in response to herbivory.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 8 July 24

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
require(DHARMa)

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field")
#setwd("~/Documents/personnel/Jameel/chapter3")

#this is where you specify the folder where you have the data on your computer


#read in data 
estess <- read.csv("estess2021_wide.csv", stringsAsFactors = TRUE)

sapply(estess,class)
##Some  variables are being read as characters not factors. Let's fix that
estess $Block <-factor(estess $Block)



estess <- filter(estess, Exclude_2021 == "include")



##Change the baseline for pesticide treatment
estess $Treatment <-factor(estess $Treatment, levels = c("Control", "Pesticide"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
estess $S_elev<-scale(estess $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
estess $S_initdiam<-scale(estess $initial_size,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
estess$elev_km<-estess $elevation/1000


#Subsetting data for plants that flowered
repro <- filter(estess, Flowered_2021 == 1) #only plants that flowered


cols=c("#882255","#6699cc")




#*******************************************************************************
####  fitness: Standardize elevation #####
#*******************************************************************************

#survival

#reproduction ####

field_R= ggplot(estess, aes(x= elev_km,y= Season_survival_2021, group= Treatment, 
                           colour= Treatment))+geom_point(size=5) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Source Elevation")  
field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x) +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2)) +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

#####survival

##mod_survR_LA <-glmmTMB(Season_survival~ S_initdiam+ S_elev*Treatment*garden+I(S_elev^2)*Treatment*garden+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

mod_survR_LA <-glmmTMB(Season_survival_2021~ S_initdiam+S_elev*Treatment+(1|PlantID)+(1|Block)+(1|Genotype),data=estess,family=binomial(link="logit"))
Anova(mod_survR_LA,type="III")

simulationOutput <- simulateResiduals(fittedModel= mod_survR_LA, plot = T, re.form = NULL,allow.new.levels =T)


surv_means<-emmeans(mod_survR_LA, ~ Treatment, type="response", adjust = "sidak")
cld(surv_means, details=TRUE)



surv_prob <- ggpredict(mod_survR_LA, terms = c("S_elev[all]","Treatment"), type = "re", interval="confidence")
surv_plot <-plot(surv_prob, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Probability of survival")+scale_x_continuous("Elevation")
surv_plot


##Box plot
##Box plot
surv_violin <-ggplot(estess, aes(x = Treatment, y = Season_survival_2021, fill = Treatment)) +
  geom_violin() +xlab("Herbivore treatment")+ 
  stat_summary(fun.data = "mean_cl_normal")+
  scale_y_continuous("Probability of survival") +
  geom_point(pch = 21, position = position_jitterdodge())

surv_treatment<-surv_violin + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols))





# Hurdle model #### 
##hurdle_Model <- glmmTMB(Mature_length_siliques ~ S_initdiam + Treatment*S_elev*garden+Treatment*I(S_elev^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~ S_initdiam +Treatment * S_elev*garden+Treatment*I(S_elev^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))
hurdle_Model <- glmmTMB(Mature_length_siliques ~ Year+cohort+S_initdiam + Treatment*S_elev*garden+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~ Year+cohort+S_initdiam +Treatment * S_elev*garden+I(S_elev^2)+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))


hurdle_Model <- glmmTMB(Mature_length_siliques ~ Year+cohort+S_initdiam + Treatment*S_elev*garden+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~ S_initdiam +Treatment * S_elev*garden+I(S_elev^2)+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))


##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model,type="III", component="cond")

simulationOutput <- simulateResiduals(fittedModel= hurdle_Model, plot = T, re.form = NULL,allow.new.levels =T)


coefficients_fec <- emtrends(hurdle_Model, specs = c("Treatment","garden"), var = "S_elev",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('Treatment',"garden",'S_elev.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table

cld(coefficients_fec, details=TRUE)

library(ggeffects)
pred_fec <- ggpredict(hurdle_Model, terms = c("S_elev[all]","Treatment","garden"), type = "zero_inflated_random", interval="confidence")
Local_adaptation <-plot(pred_fec, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,1000)
Local_adaptation


#reproduction ####

field_R= ggplot(field, aes(x= elev_km,y= Reproduction, group= Treatment, 
                           colour= Treatment))+geom_point(size=5) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Source Elevation")  
field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

#prob of repro 

#mod_repro <-glmmTMB(Reproduction~ S_initdiam + Treatment*S_elev*garden+Treatment*I(S_elev^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
mod_repro <-glmmTMB(Reproduction~ S_initdiam +Year+cohort+ Treatment*S_elev*garden+I(S_elev^2)+(1|PlantID) + (1| Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
Anova(mod_repro, type="III") 
simulationOutput <- simulateResiduals(fittedModel= mod_repro, plot = T, re.form = NULL,allow.new.levels =T)

pred_repro <- ggpredict(mod_repro, terms = c("S_elev[all]","Treatment","garden"), type = "re", interval="confidence")
repro_plot <-plot(pred_repro, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Probability of reproduction")+scale_x_continuous("Source elevation")
repro_plot


#fecundity ####

repro <- filter(field, Reproduction == 1) #only plants that reproduced

#length siliques
##mod_repro_f <- glmmTMB(Mature_length_siliques ~ S_initdiam +S_elev*Treatment*garden+I(S_elev^2)*Treatment*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))
mod_repro_f <- glmmTMB(Mature_length_siliques ~ S_initdiam+Year+cohort + S_elev*Treatment*garden +(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))
Anova(mod_repro_f, type="III") 
simulationOutput <- simulateResiduals(fittedModel= mod_repro_f, plot = T, re.form = NULL,allow.new.levels =T)

fec_pred <- ggpredict(mod_repro_f, terms = c("S_elev[all]","Treatment","garden"), type = "re", interval="confidence")
fec_plot <-plot(fec_pred, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Fecundity")+scale_x_continuous("Elevation")
fec_plot

fec_pred2 <- ggpredict(mod_repro_f, terms = c("S_elev[all]","Treatment"), type = "re", interval="confidence")
fec_plot2 <-plot(fec_pred2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Fecundity")+scale_x_continuous("Elevation")
fec_plot2


coefficients_fec <- emtrends(mod_repro_f, specs = c("Treatment","garden"), var = "S_elev",type="response")


result <- ggpredict(mod_repro_f, c("S_elev","Treatment","garden"))
plot(result, show_residuals=TRUE,)





#*******************************************************************************
####  fitness:Elevational difference #####
#*******************************************************************************
# Hurdle model #### 
##hurdle_Model <- glmmTMB(Mature_length_siliques ~ S_initdiam + Treatment*elev_dist_km*garden+Treatment*I(elev_dist_km^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~ S_initdiam +Treatment * elev_dist_km*garden+Treatment*I(elev_dist_km^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))
hurdle_Model <- glmmTMB(Mature_length_siliques ~ S_initdiam + Treatment*elev_dist_km*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~ S_initdiam +Treatment * elev_dist_km*garden+I(elev_dist_km^2)+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model,type="III", component="cond")

simulationOutput <- simulateResiduals(fittedModel= hurdle_Model, plot = T, re.form = NULL,allow.new.levels =T)


coefficients_fec <- emtrends(hurdle_Model, specs = c("Treatment","garden"), var = "elev_dist_km",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('Treatment',"garden",'elev_dist_km.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(elev_dist_km.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table

Fecmeans<-emmeans(hurdle_Model, ~ Water:Herbivore, type="response")
cld(coefficients_fec, details=TRUE)

library(ggeffects)
pred_fec <- ggpredict(hurdle_Model, terms = c("elev_dist_km[all]","Treatment","garden"), type = "zero_inflated_random", interval="confidence")
Local_adaptation <-plot(pred_fec, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,1000)
Local_adaptation


#reproduction ####

field_R= ggplot(field, aes(x= elev_km,y= Reproduction, group= Treatment, 
                           colour= Treatment))+geom_point(size=5) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Source Elevation")  
field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

field_R + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

#prob of repro 

#mod_repro <-glmmTMB(Reproduction~ S_initdiam + Treatment*elev_dist_km*garden+Treatment*I(elev_dist_km^2)*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
mod_repro <-glmmTMB(Reproduction~ S_initdiam + Treatment*elev_dist_km*garden+I(elev_dist_km^2)+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

Anova(mod_repro, type="III") 

pred_repro <- ggpredict(mod_repro, terms = c("elev_dist_km[all]","Treatment","garden"), type = "re", interval="confidence")
repro_plot <-plot(pred_repro, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme_classic(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Probability of reproduction")
repro_plot



#fecundity ####

repro <- filter(field, Reproduction == 1) #only plants that reproduced

field_f= ggplot(repro, aes(x= elev_km,y= Mature_length_siliques, group= Treatment, 
                           colour= Treatment))+geom_point(size=5) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
field_f + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


field_f + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


field_R_box<-ggplot(repro, aes(x = garden, y = Mature_length_siliques, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Mature length siliques") +
  geom_point(pch = 21, position = position_jitterdodge())


field_R_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


#length siliques
##mod_repro_f <- glmmTMB(Mature_length_siliques ~ S_initdiam +elev_dist_km*Treatment*garden+I(elev_dist_km^2)*Treatment*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))
mod_repro_f <- glmmTMB(Mature_length_siliques ~ S_initdiam +elev_dist_km*Treatment*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))

Anova(mod_repro_f, type="III") 
simulationOutput <- simulateResiduals(fittedModel= mod_repro_f, plot = T, re.form = NULL,allow.new.levels =T)

fec_pred <- ggpredict(mod_repro_f, terms = c("elev_dist_km[all]","Treatment","garden"), type = "re", interval="confidence")
fec_plot <-plot(fec_pred, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Fecundity")+scale_x_continuous("Elevational distance")
fec_plot

fec_pred2 <- ggpredict(mod_repro_f, terms = c("elev_dist_km[all]","Treatment"), type = "re", interval="confidence")
fec_plot2 <-plot(fec_pred2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Fecundity")+scale_x_continuous("Elevational distance")
fec_plot2


plot(predictorEffects(mod_repro_f, ~ elev_dist_km), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Fecundity",line=list(multiline=TRUE, lty=1:2,col=c(cols)))


coefficients_fec <- emtrends(mod_repro_f, specs = c("garden"), var = "elev_dist_km",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('garden','elev_dist_km.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(elev_dist_km.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

result <- ggpredict(mod_repro_f, c("elev_dist_km","garden"))
plot(result, show_residuals=TRUE,)



mod_reproF_LA_2 <- glmmTMB(Mature_length_siliques ~ Treatment * S_elev  *garden + Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))
Anova(mod_reproF_LA_2,type="III")

fec_pred3 <- ggpredict(mod_repro_f, terms = c("S_elev[all]","Treatment"), type = "re", interval="confidence")
fec_plot3 <-plot(fec_pred2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Fecundity")+scale_x_continuous("Standardized elevation")
fec_plot3


Local_adaptation <-plot(fec_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = c("#882255","#6699cc"), facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,1700)

Local_adaptation


#survival

##mod_survR_LA <-glmmTMB(Season_survival~ S_initdiam+ elev_dist_km*Treatment*garden+I(elev_dist_km^2)*Treatment*garden+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
mod_survR_LA <-glmmTMB(Season_survival~ S_initdiam+ elev_dist_km*Treatment*garden+I(elev_dist_km^2)+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

Anova(mod_survR_LA)

mod_surv <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_survR_LA)

plot(mod_surv, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=TRUE, pch=19, col="black"))


survived <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_survR_LA)
plot(survived, lwd=2,xlab="Distance from Source Elevation", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#*******************************************************************************
#### 2.leaf damage across censuses #####
#*******************************************************************************
#####repeated measures with all damage data####
LARdat<-estess[,c(1:10,49:58)]


#reformat datafile

LAR_data<- LARdat %>% pivot_longer(cols=c("LAR_1_2021","LAR_2_2021","LAR_3_2021","LAR_4_2021","LAR_5_2021","LAR_6_2021","LAR_7_2021","LAR_8_2021","LAR_9_2021","LAR_10_2021"),
                                   names_to='census',
                                   values_to='LAR')

LAR_data$census[LAR_data$census == "LAR_1_2021"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2_2021"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3_2021"] <-"3"
LAR_data$census[LAR_data$census == "LAR_4_2021"] <- "4"
LAR_data$census[LAR_data$census == "LAR_5_2021"] <- "5"
LAR_data$census[LAR_data$census == "LAR_6_2021"] <- "6"
LAR_data$census[LAR_data$census == "LAR_7_2021"] <- "7"
LAR_data$census[LAR_data$census == "LAR_8_2021"] <- "8"
LAR_data$census[LAR_data$census == "LAR_9_2021"] <- "9"
LAR_data$census[LAR_data$census == "LAR_10_2021"] <- "10"

LAR_data $census <-as.factor(LAR_data $census)

#LAR_data $Year <-as.factor(LAR_data $Year)


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

LAR_data $S_elev<-scale(LAR_data $elevation,center=TRUE, scale=TRUE)
#LAR_data $Garden_Block<-interaction(LAR_data $garden, LAR_data $block,sep = "_")


mod_LAR<- gamlss (formula= y_beta ~S_elev*Treatment+ random(census) + random(PlantID)+ random(Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(mod_LAR)
# Save  model to .rda file 
save(mod_LAR, file='mod_LAR.rda')   

plot(mod_LAR)
drop1(mod_LAR)


pred_lar <- ggpredict(mod_LAR, terms = c("S_elev[all]", "Treatment"), type = "re", interval="confidence")
lar_plot <-plot(pred_lar, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("LAR")
lar_plot


mod_LAR_no_interact<- gamlss (formula= y_beta ~Treatment+S_elev+ random(census) + random(PlantID)+ random(Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(mod_LAR_no_interact)

drop1(mod_LAR_no_interact)


LAR_means<-emmeans(mod_LAR_no_interact, ~ Treatment, type="response", adjust = "sidak")
cld(surv_means, details=TRUE)

#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>%  
  mutate(fit.m = predict(mod_LAR, se.fit=FALSE),              
         resid = residuals(mod_LAR))

##Convert coefficients to probabilities
newdf2 $predicted<-1/(1+exp(-(newdf2 $fit.m))) 
newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


##Box plot
LAR_box <-ggplot(LAR_data, aes(x = Treatment, y = LAR_prop, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ 
  scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

damage_treatment<-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))



#cline, plotting predicted data
damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted, group= Treatment, 
                                 colour= Treatment))+geom_point(pch = 20, size = 1) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + #facet_grid(~ garden) + 
  theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Herbivore availability", labels = c("Pesticide","Control"))

damage_cline #fig 2

#env concatenates water and herbivore. This model allows us to extract means and slopes for each combination of year, treatment and garden
LAR_data $env<-interaction(LAR_data $Treatment, LAR_data $garden)
LAR_data $env_year<-interaction(LAR_data $env, LAR_data $Year)

LAR_ModelE<- gamlss (formula= y_beta ~S_elev*env_year+cohort+ random(census) + random(PlantID)+ random(Garden_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(LAR_ModelE)


##Conversation of standardized elevation back to raw units
-2*sd(newdf2 $elevation,na.rm=TRUE)+mean(newdf2 $elevation,na.rm=TRUE)
-1*sd(newdf2 $elevation,na.rm=TRUE)+mean(newdf2 $elevation,na.rm=TRUE)
0*sd(newdf2 $elevation,na.rm=TRUE)+mean(newdf2 $elevation,na.rm=TRUE)
1*sd(newdf2 $elevation,na.rm=TRUE)+mean(newdf2 $elevation,na.rm=TRUE)
2*sd(newdf2 $elevation,na.rm=TRUE)+mean(newdf2 $elevation,na.rm=TRUE)

#####
