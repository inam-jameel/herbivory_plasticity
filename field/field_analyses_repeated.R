######## PROJECT: Field experiment: variation in herbivore damage due to treatment
#### PURPOSE:Examine fitness, traits in response to herbivory.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 5 Apr 24

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
setwd("~/Documents/personnel/Jameel/chapter3")

 #this is where you specify the folder where you have the data on your computer


#read in data 
field <- read.csv("field_long.csv", stringsAsFactors = TRUE)

sapply(field,class)
##Some  variables are being read as characters not factors. Let's fix that
field $cohort <-as.factor(field $cohort)
field $Year <-as.factor(field $Year)
field <- filter(field, Exclude == "Include")

#Let's concatenate garden and block
field $Garden_Block<-interaction(field$garden, field$block,sep = "_")

#Let's concatenate garden and treatment
field $treat<-interaction(field$garden, field$Treatment,sep = "_")

#convert elevation distance (source population elevation - garden elevation) to km
field$elev_dist_km<-field$elev_dist/1000
hist(field$elev_dist_km)


#This standardizes elevation distance to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
field $S_elev_dist<-scale(field $elev_dist,center=TRUE, scale=TRUE)


##Change the baseline for pesticide treatment
field $Treatment <-factor(field $Treatment, levels = c("Control", "Pesticide"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
field $S_elev<-scale(field $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
field $S_initdiam<-scale(field $initial_size,center=TRUE, scale=TRUE)

#This makes number of year as factor
field $Year <-as.factor(field $Year)

#This rescales source elevation from meters to km
field$elev_km<-field $elevation/1000


#Subsetting data for plants that reproduced (had siliques)
repro <- filter(field, Reproduction == 1) #only plants that reproduced


cols=c("#CC79A7","lightblue")


#*******************************************************************************
####  fitness: Standardize elevation #####
#*******************************************************************************
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

Fecmeans<-emmeans(hurdle_Model, ~ Water:Herbivore, type="response")
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

pred_repro <- ggpredict(mod_repro_f, terms = c("S_elev[all]","Treatment","garden"), type = "re", interval="confidence")
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


coefficients_fec <- emtrends(mod_repro_f, specs = c("Treatment"), var = "S_elev",type="response")


result <- ggpredict(mod_repro_f, c("S_elev","garden"))
plot(result, show_residuals=TRUE,)


#survival

##mod_survR_LA <-glmmTMB(Season_survival~ S_initdiam+ S_elev*Treatment*garden+I(S_elev^2)*Treatment*garden+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
mod_survR_LA <-glmmTMB(Season_survival~ S_initdiam+Year+cohort+ S_elev*Treatment*garden+I(S_elev^2)+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))
Anova(mod_survR_LA,type="III")

surv_prob <- ggpredict(mod_survR_LA, terms = c("S_elev[all]","Treatment","garden"), type = "re", interval="confidence")
surv_plot <-plot(surv_prob, show_data=TRUE, show_title =FALSE, show_legend=TRUE, facet=FALSE,jitter=0.025)+theme(text = element_text(size=10),axis.ticks = element_line(colour = "black"),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position="right")+ scale_y_continuous("Probability of survival")+scale_x_continuous("Elevation")
surv_plot



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
LARdat<-field[,c(1:14,51:59)]


#reformat datafile

LAR_data<- LARdat %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3","LAR_4","LAR_5","LAR_6","LAR_7","LAR_8","LAR_9"),
                                        names_to='census',
                                        values_to='LAR')

LAR_data$census[LAR_data$census == "LAR_1"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3"] <-"3"
LAR_data$census[LAR_data$census == "LAR_4"] <- "4"
LAR_data$census[LAR_data$census == "LAR_5"] <- "5"
LAR_data$census[LAR_data$census == "LAR_6"] <- "6"
LAR_data$census[LAR_data$census == "LAR_7"] <- "7"
LAR_data$census[LAR_data$census == "LAR_8"] <- "8"
LAR_data$census[LAR_data$census == "LAR_9"] <- "8"

LAR_data $census <-as.factor(LAR_data $census)

LAR_data $Year <-as.factor(LAR_data $Year)


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
LAR_data $Garden_Block<-interaction(LAR_data $garden, LAR_data $block,sep = "_")


mod_LAR<- gamlss (formula= y_beta ~S_elev*Treatment*garden+Year+cohort+ random(census) + random(PlantID)+ random(Garden_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(mod_LAR)
# Save  model to .rda file 
save(mod_LAR, file='mod_LAR.rda')   


drop1(mod_LAR)


pred_lar <- ggpredict(mod_LAR, terms = c("S_elev[all]", "Treatment","garden"), type = "re", interval="confidence")
lar_plot <-plot(pred_lar, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("LAR")
lar_plot


mod_LAR_two_way<- gamlss (formula= y_beta ~Treatment*garden+S_elev*garden+S_elev*Treatment+Year+cohort+ random(census) + random(PlantID)+ random(Garden_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(mod_LAR_two_way)

mod_LAR_no_interact<- gamlss (formula= y_beta ~Treatment+garden+S_elev+Year+cohort+ random(census) + random(PlantID)+ random(Garden_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(mod_LAR_no_interact)



#*******************************************************************************
#### 3. traits #####
#*******************************************************************************

##### Flowering time  ####

##### ordinal flowering time ####
##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(field, Date_flowering!="NA")

# actual model: Neither of these models have great residuals and the spread of the data seems odd. I'm not sure I trust these data.

##flowering_time_ordinal<-glmmTMB(Date_flowering ~ S_initdiam+Year+cohort+ S_elev*Treatment*garden+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data= flowering,family=lognormal(link="log"))

flowering_time_ordinal<-lmer(Date_flowering ~ Year+cohort+ S_elev*Treatment*garden+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data= flowering)
Anova(flowering_time_ordinal, type = "III") 
simulationOutput <- simulateResiduals(fittedModel= flowering_time_ordinal, plot = T, re.form = NULL)
cols=c("#CC79A7","blue")
pred_ordinal<- ggpredict(flowering_time_ordinal, terms = c("S_elev[all]", "Treatment","garden"), type = "re", interval="confidence")
ordinal_fig <-plot(pred_ordinal, show_residuals=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")
ordinal_fig


##### Duration of flowering####
flowering$duration<-flowering$Date_silique-flowering$Date_flowering+0.01
hist(flowering$duration)


flowering_duration<-glmer(duration ~ Year+cohort+ S_elev*Treatment*garden+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data= flowering, family=Gamma(link="log"))
Anova(flowering_duration, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)


##### Height at flowering####
flowering$Max_height_flowering <-pmax(flowering$Height1_flowering, flowering$Height2_flowering, flowering$Height3_flowering, na.rm=TRUE)
hist(flowering$Max_height_flowering )

height<-lmer(Max_height_flowering ~ Year+cohort+ S_elev*Treatment*garden+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data= flowering)
Anova(height, type = "III") # elevation is significant, not year

#*******************************************************************************
#### Leaf traits #####
#*******************************************************************************
##Let's correlate rosette and bolt leaf data to see if we can come up with composite figures

hist(field$rosette_succulence)
hist(field$bolt_succulence)

plot(field$rosette_succulence~field$bolt_succulence)
plot(field$rosette_SLA~field$bolt_SLA)
plot(field$rosette_LWC~field$bolt_LWC)

##mod1<-lm(field$rosette_succulence~field$bolt_succulence)
##summary(mod1)
##Anova(mod1,type="III")

modSLA<-lm(field$rosette_SLA~field$bolt_SLA)
summary(modSLA)

mod3<-lm(field$rosette_LWC~field$bolt_LWC)
summary(mod3)

##Create composite leaf SLA, succuclence and LWC variables based on the regressions above. This gives us foliar trait data for the 67 plants for which we have bolt but not rosette collections
field$SLA <- ifelse(is.na(field$rosette_SLA), (32.29571  + 0.49145*field$bolt_SLA), field$rosette_SLA)
field$LWC <- ifelse(is.na(field$rosette_LWC), (0.05897 + 0.66460*field$bolt_LWC), field$rosette_LWC)

hist(field$SLA)
hist(field$LWC)
###### SLA  ####
##To enable simulate residuals to work, we have to exclude plants that did not flower
leaf<-subset(field, SLA!="NA")
SLA_RM <- lmer(SLA ~ Year+cohort+Treatment*garden*S_elev+(1|PlantID) +(1|Genotype)+(1|Garden_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = leaf)
simulationOutput <- simulateResiduals(fittedModel= SLA_RM, plot = T, re.form = NULL)
Anova(SLA_RM, type = "III") # garden and year are significant

visreg(SLA_RM,"S_elev", by="garden", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))


###### LWC  ####
leaf2<-subset(field, LWC!="NA")
LWC_RM <- lmer(LWC ~ Year+cohort+Treatment*garden*S_elev+(1|PlantID) +(1|Genotype)+(1|Garden_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = leaf2)
Anova(LWC_RM, type = "III") # cohort, year is significant
simulationOutput <- simulateResiduals(fittedModel= LWC_RM, plot = T, re.form = NULL)


