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

 #this is where you specify the folder where you have the data on your computer


#read in data 
field <- read.csv("field_summary_long.csv", stringsAsFactors = TRUE)

sapply(field,class)
##Some  variables are being read as characters not factors. Let's fix that
#field$Genotype<-as.factor(field$Genotype)
#field$garden<-as.factor(field$garden)
#field$Treatment<-as.factor(field$Treatment)
#field$elevation<-as.numeric(field$elevation)
#field$block<-as.factor(field$block)
#field $PlantID <-as.factor(field $PlantID)

field $cohort <-as.factor(field $cohort)

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

field <- filter(field, Exclude == "Include")

#Subsetting data for plants that reproduced (had siliques)
repro <- filter(field, Reproduction == 1) #only plants that reproduced


cols=c("#CC79A7","lightblue")

#*******************************************************************************
####  fitness #####
#*******************************************************************************


# Hurdle model #### 

#include intital size?

hurdle_Model <- glmmTMB(Mature_length_siliques ~Treatment*elev_dist_km*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=field, zi=~Treatment * elev_dist_km*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))


simulationOutput <- simulateResiduals(fittedModel= hurdle_Model, plot = T, re.form = NULL,allow.new.levels =T)


##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model,type="III", component="cond")

coefficients_fec <- emtrends(hurdle_Model, specs = c("Treatment","garden"), var = "elev_dist_km",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('Treatment',"garden",'elev_dist_km.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(elev_dist_km.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table

Fecmeans<-emmeans(fecundity, ~ Water:Herbivore, type="response")
cld(coefficients_fec, details=TRUE)

pred_fec <- ggpredict(hurdle_Model, terms = c("elev_dist_km[all]","Treatment"), type = "re", interval="confidence")


Local_adaptation <-plot(pred_fec, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,1000)

library(ggeffects)

result <- ggpredict(hurdle_Model, c("elev_dist_km","Treatment"))
plot(result, show_residuals=TRUE,)

result <- ggpredict(hurdle_Model, c("elev_dist_km","Treatment","garden"))
plot(result,show_residuals=TRUE,jitter = FALSE)

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

mod_repro <-glmmTMB(Reproduction~ elev_km*Treatment*garden+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

Anova(mod_repro, type="III") 

visreg(mod_repro, "Year", type="conditional", scale = "response", 
       xlab="Year", ylab="Probability of Reproduction", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))



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
mod_repro_f <- glmmTMB(Mature_length_siliques ~elev_dist_km*Treatment*garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=ziGamma(link="log"))

Anova(mod_repro_f, type="III") 

plot(predictorEffects(mod_repro_f, ~ elev_dist_km), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Fecundity",line=list(multiline=TRUE, lty=1:2,col=c(cols)))


coefficients_fec <- emtrends(mod_repro_f, specs = c("garden"), var = "elev_dist_km",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('garden','elev_dist_km.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(elev_dist_km.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

result <- ggpredict(mod_repro_f, c("elev_dist_km","garden"))
plot(result, show_residuals=TRUE,)


###### local adaptation models: fitness ~ T x elevational distance (source – transplant) ######
#survival, reproduction, fecundity 

#hurdle_ModelLA = glmmTMB (Mature_length_siliques ~ Treatment * elev_dist_km*garden+I(elev_dist_km^2)*Treatment*garden  +Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), zi=~ Treatment * elev_dist_km*garden+I(elev_dist_km^2)*Treatment*garden +Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype),data = repro ,family=Gamma(link=log))

#Warning message:In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) : Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')

#Anova(hurdle_ModelLA,type="III") 

#summary(hurdle_ModelLA)

#fruit <-predictorEffect("elev_dist_km",  partial.residuals=FALSE,hurdle_ModelLA)
#plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2, col="black"), 
#     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1000))

#fecundity

mod_reproF_LA <- glmmTMB(Mature_length_siliques ~ Treatment * S_elev_dist  *garden + I(S_elev_dist^2) * Treatment *garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))

Anova(mod_reproF_LA,type="III")

fruitA <-predictorEffect("S_elev_dist",  partial.residuals=TRUE,mod_reproF_LA)
plot(fruitA, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=TRUE, lty=1:2, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,800))


plot(fruitA, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=1:2, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1500))



mod_reproF_LA_2 <- glmmTMB(Mature_length_siliques ~ Treatment * S_elev  *garden + I(S_elev^2) * Treatment *garden+Year+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=repro,family=Gamma(link="log"))

fruit <-predictorEffect("S_elev",  partial.residuals=TRUE,mod_reproF_LA_2)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,400))
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,400))

Anova(mod_reproF_LA,type="III")


fruit <-predictorEffect("elev_dist_km",  partial.residuals=TRUE,mod_reproF_LA)




##Then, we examine local adaptation using the probability of reproduction as the fitness component. Here, values of 0 include individuals that failed to flower along with those that died


mod_reproR_LA <-glmmTMB(Reproduction~ elev_dist_km*Treatment*garden+I(elev_dist_km^2)*Treatment*garden+cohort+Year+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

Anova(mod_reproR_LA)

pred_repro <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_reproR_LA)

plot(pred_repro, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", pch=19, type="response",lines=list(multiline=TRUE, lty=1:2, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))


plot(pred_repro, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", pch=19, type="response",lines=list(multiline=FALSE, lty=1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))



#survival

mod_survR_LA <-glmmTMB(Season_survival~ elev_dist_km*Treatment*garden+I(elev_dist_km^2)*Treatment*garden+Year+cohort+(1|PlantID)+(1|Garden_Block)+(1|Genotype),data=field,family=binomial(link="logit"))

Anova(mod_survR_LA)

mod_surv <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_survR_LA)

plot(mod_surv, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=TRUE, pch=19, col="black"))


survived <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_survR_LA)
plot(survived, lwd=2,xlab="Distance from Source Elevation", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))




#*******************************************************************************
#### Lifetime Fitness #####
#*******************************************************************************


#read in data 
fieldLF <- read.csv("field_lifetime_fitness.csv",stringsAsFactors = TRUE)

sapply(fieldLF,class)
##Some  variables are being read as characters not factors. Let's fix that
fieldLF$Genotype<-as.factor(fieldLF$Genotype)
fieldLF$Treatment<-as.factor(fieldLF$Treatment)
fieldLF $PlantID <-as.factor(fieldLF $PlantID)
fieldLF $cohort <-as.factor(fieldLF $cohort)


#Let's concatenate garden and block
fieldLF $Garden_Block<-interaction(fieldLF$garden, fieldLF$block,sep = "_")

#Let's concatenate garden and treatment
fieldLF $treat<-interaction(fieldLF$garden, fieldLF$Treatment,sep = "_")

#convert elevation distance (source population elevation - garden elevation) to km
fieldLF$elev_dist_km<-fieldLF$elev_dist/1000
hist(fieldLF$elev_dist_km)

##Change the baseline for pesticide treatment
fieldLF $Treatment <-factor(fieldLF $Treatment, levels = c("Control", "Pesticide"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
fieldLF $S_elev<-scale(fieldLF $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
fieldLF $S_initdiam<-scale(fieldLF $initial_size,center=TRUE, scale=TRUE)

#This rescales source elevation from meters to km
fieldLF$elev_km<-fieldLF $elevation/1000


#fieldLF $Garden_Block <-as.factor(fieldLF $Garden_Block)

reproLF <- filter(fieldLF, Seasons_reproduced != 0) #only plants that reproduced




#multinomial regression

library(nnet)

ggplot(grasshopperLF %>% filter(years_reproduced == 0 | 
                                  years_reproduced == 1),
       aes(x = as.numeric(elev_km), y = as.numeric(years_reproduced))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic()


ggplot(grasshopperLF %>% filter(years_reproduced == "0" | 
                                  years_reproduced == "2"),
       aes(x = as.numeric(elev_km), y = as.numeric(years_reproduced))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic()

fit_basic <- multinom(years_reproduced ~ Treatment*S_elev*Herbivore, data = grasshopperLF)

summary(fit_basic)

library(broom)
tidy(fit_basic, conf.int = TRUE) # nothing significant

##### Hurdle model ####

hurdle_Model_LF <- glmmTMB(Overall_silique_length ~Treatment*elev_dist_km*garden+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=fieldLF, zi=~Treatment * elev_dist_km*garden+cohort +(1|PlantID) + (1| Garden_Block)+(1|Genotype),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model_LF,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_LF,type="III", component="cond")

library(ggeffects)


result <- ggpredict(hurdle_Model_LF, c("elev_dist_km","Treatment"))
plot(result, show_residuals=TRUE,)

result <- ggpredict(hurdle_Model_LF, c("elev_dist_km","Treatment","garden"))
plot(result,show_residuals=TRUE,jitter = FALSE)



require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model)

plotQQunif(hurdle_Model)
plotResiduals(hurdle_Model)

summary(hurdle_Model)



##### Probability of reproduction ####

field_f= ggplot(reproLF, aes(x= elev_km,y= Mature_length_siliques, group= Treatment, 
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

mod_reproF_LF <- glmmTMB(Overall_silique_length ~ Treatment * S_elev  *garden +cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=reproLF,family=Gamma(link="log"))

Anova(mod_reproF_LF, type="III") 

plot(predictorEffects(mod_repro_f, ~ elev_km), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Fecundity",line=list(multiline=TRUE, lty=1:4,col=c("#6699cc","#882255")))

visreg(mod_repro_f, overlay = TRUE, "elev_km", by="Treatment", type="conditional", scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


visreg(mod_repro_f, overlay = FALSE, "elev_km", by="Treatment", type="conditional", scale = "response", 
       xlab="Source elevation (km)", ylab="silique length")


########treatment####

mod_repro <-glmmTMB(reproduced~S_elev*Treatment*Herbivore+(1|Cage_Block)+(1|Genotype),data=grasshopperLF,family=binomial(link="logit"))

Anova(mod_repro, type="III") 

repro_df <- field %>% 
  
  mutate(fit.m = predict(mod_repro, re.form = NA),
         
         fit.c = predict(mod_repro, re.form = NULL), #all random effects
         
         resid = residuals(mod_repro))

##Convert fit.m and fit.c back to the proportional scale

repro_df $fit.m_trans<-1/(1+exp(-(repro_df $fit.m))) 

repro_df $fit.c_trans<-1/(1+exp(-(repro_df $fit.m))) 

repro_df $resid_trans<-1/(1+exp(-(repro_df $resid))) 



vioplot(fit.c_trans ~ Herbivore*year, data= repro_df, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,ylim=c(0,0.2),colMed = "black",colMed2 = c("lightblue","#CC79A7"), col=c("lightblue","#CC79A7"), ylab="Volumetric water content", xlab="Watering treatment") +stripchart(fit.c_trans ~ Treatment*year, data= repro_df, col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)

ggplot(repro_df, aes(x = Herbivore, y = fit.c_trans, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("ivory", "#117733")) +
  labs(y = "Probability of Reproduction") +
  labs(x = "Herbivore Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

#### fecundity ####


field_LF= ggplot(reproLF, aes(x= elev_km,y= Overall_silique_length, group= Treatment, 
                           colour= Treatment))+geom_point(size=5) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
field_LF + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


field_LF + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


field_LF_box<-ggplot(reproLF, aes(x = garden, y = Overall_silique_length, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Mature length siliques") +
  geom_point(pch = 21, position = position_jitterdodge())


field_LF_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Herbivore treatment", labels = c("Control","Pesticide"))


#length siliques
mod_repro_lf <- glmmTMB(Overall_silique_length ~elev_km*Treatment+garden+cohort+(1|PlantID) + (1| Garden_Block)+(1|Genotype), data=reproLF,family=ziGamma(link="log"))

Anova(mod_repro_lf, type="III") 

plot(predictorEffects(mod_repro_lf, ~ elev_km), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Fecundity",line=list(multiline=TRUE, lty=1:4,col=c("#6699cc","#882255")))

visreg(mod_repro_lf, overlay = TRUE, "elev_km", by="Treatment", type="conditional", scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


visreg(mod_repro_f, overlay = FALSE, "elev_km", by="Treatment", type="conditional", scale = "response", 
       xlab="Source elevation (km)", ylab="silique length")

####### treatment####
mod_fecundity <- glmmTMB (total_silique_length ~ S_elev*Treatment*Herbivore  + (1|Cage_Block)+(1|Genotype), family=Gamma(link="log"), data = repro)


Anova(mod_fecundity,type="III") #signficiant watered x herbivore treatment

summary(mod_fecundity)

fecundity_df <- repro %>% 
  
  mutate(fit.m = predict(mod_fecundity, re.form = NA),
         
         fit.c = predict(mod_fecundity, re.form = NULL), #all random effects
         
         resid = residuals(mod_fecundity))


ggplot(fecundity_df, aes(x = Treatment, y = fit.c, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Mature Length Siliques") +
  labs(x = "Watering Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 


mod_fecundity_num <- glmmTMB (total_silique_number ~ S_elev*Treatment*Herbivore  + (1|Cage_Block)+(1|Genotype), family=Gamma(link="log"), data = repro)


Anova(mod_fecundity_num,type="III") #signficiant watered x herbivore treatment

summary(mod_fecundity_num)

fecundity_num_df <- repro %>% 
  
  mutate(fit.m = predict(mod_fecundity_num, re.form = NA),
         
         fit.c = predict(mod_fecundity_num, re.form = NULL), #all random effects
         
         resid = residuals(mod_fecundity_num))


ggplot(fecundity_num_df, aes(x = Treatment, y = fit.c, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Mature Number Siliques") +
  labs(x = "Watering Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

#*******************************************************************************
#### 2.Repeated measures Fitness #####
#*******************************************************************************

##### Hurdle model ####

#filter out the first year fitness, two plants that flowered in the first season
grasshopperRM <- filter(field, year_num != 1)

hurdle_ModelRM <- glmmTMB(Mature_length_siliques ~Treatment*S_elev*Herbivore+year_num + (1|PlantID)+ (1| Genotype)+ (1| Cage_Block), data=field, zi=~Treatment*S_elev*Herbivore+year_num + (1|PlantID)+ (1| Genotype)+ (1| Cage_Block),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_ModelRM,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_ModelRM,type="III", component="cond")



library(ggeffects)


result <- ggpredict(hurdle_ModelRM, c("Herbivore","Treatment"))
plot(result)



require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model)

plotQQunif(hurdle_Model)
plotResiduals(hurdle_Model)

summary(hurdle_Model)


hurdle_ModelRM <- glmmTMB(Mature_length_siliques ~avg_vwc*S_elev*avg_LAR+year_num + (1|PlantID)+ (1| Genotype)+ (1| Cage_Block), data=field, zi=~avg_vwc*S_elev*avg_LAR+year_num + (1|PlantID)+ (1| Genotype)+ (1| Cage_Block),family=ziGamma(link="log"))

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_ModelRM,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_ModelRM,type="III", component="cond")



##### Probability of reproduction ####


grasshopper_PR= ggplot(grasshopperRM, aes(x= elev_km,y= Reproduced, group= Herbivore, 
                                          colour= Herbivore))+geom_point(size=5) + scale_y_continuous("Probability of Reproduction")+ scale_x_continuous("Source Elevation")  
grasshopper_PR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                    axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                                    panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Treatment, scales="free_x") +scale_colour_manual(values = c( "gray","darkgreen"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))


########treatment####

mod_repro_RM <-glmmTMB(Reproduced~S_elev*Treatment*Herbivore+year_num+(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=grasshopperRM,family=binomial(link="logit"))

Anova(mod_repro_RM, type="III") 

plot(predictorEffects(mod_repro_RM, ~ year_num), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Probability of seed germination",line=list(multiline=TRUE, lty=1:4,col=c("#882255","darkorange","#6699cc","#117733")))

visreg(mod_repro_RM, overlay = TRUE, "year_num", by="Herbivore", type="conditional", scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

ggplot(grasshopperRM, aes(x = Herbivore, y = Reproduced, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  #geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Probability of Reproduction") +
  labs(x = "Herbivory Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "left") 


########avgdam/vwc####

mod_repro_RM <-glmmTMB(Reproduced~S_elev*avg_vwc*avg_LAR+year_num+(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=grasshopperRM,family=binomial(link="logit"))

Anova(mod_repro_RM, type="III") # significance of herbivore*year

#### fecundity ####

reproRM <- filter(grasshopperRM, Reproduced == 1 )
hist(reproRM$Mature_length_siliques)

grasshopper_fRM= ggplot(reproRM, aes(x= elev_km,y= Mature_length_siliques, group= Treatment, 
                                 colour= Treatment))+geom_point(size=5) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
grasshopper_fRM + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                   axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                                   panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Watering treatment", labels = c("Ample Water","Water Restricted"))


ggplot(reproRM, aes(x = Treatment, y = Mature_length_siliques, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Volumetric Water Content") +
  labs(x = "") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "left") 

####### treatment####
mod_fecundityRM <-glmmTMB (Mature_length_siliques ~ S_elev*Water*Herbivore+year_num  + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = reproRM)


Anova(mod_fecundityRM,type="III") #signficiant watered x herbivore treatment


fecundity_dfRM <- reproRM %>% 
  
  mutate(fit.m = predict(mod_fecundityRM, re.form = NA),
         
         fit.c = predict(mod_fecundityRM, re.form = NULL), #all random effects
         
         resid = residuals(mod_fecundityRM))

ggplot(fecundity_dfRM, aes(x = Treatment, y = fit.c, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Mature length Siliques") +
  labs(x = "Watering Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 



ggplot(fecundity_dfRM, aes(x = treat, y = fit.c, group = treat)) +
  geom_violin(aes(fill = treat), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("darkorange","#882255","#6699cc","#117733")) +
  labs(y = "Mature Length Siliques") +
  labs(x = "Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 


dfRM_box <-ggplot(fecundity_dfRM, aes(x = Herbivore, y = fit.c, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Mature Length Siliques") +
  geom_point(pch = 21, position = position_jitterdodge())

dfRM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                             panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("field removal", "field addition")) +  scale_fill_manual(values = c("darkred","lightblue"), name = "Watering treatment", labels = c("Water-restricted","Supplemental watering"))

df_fig= ggplot(fecundity_dfRM, aes(x= elev_km,y= Mature_length_siliques, group= Herbivore, 
                            colour= Herbivore))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
df_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Water, scales="free_x") +scale_colour_manual(values = c( "darkgreen","gray"), name = "Herbivore treatment", labels = c("Addition","Removal"))

df_fig= ggplot(fecundity_dfRM, aes(x= elev_km,y= Mature_length_siliques, group= Water, 
                                   colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
df_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

####### avgdam/vwc####
mod_fecundityRM <-glmmTMB (Mature_length_siliques ~ S_elev*avg_LAR*avg_vwc+year_num  + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = reproRM)


Anova(mod_fecundityRM,type="III") #signficiant watered x herbivore treatment




mod_fecundity_numRM<-glmmTMB (Mature_silique_number ~ S_elev*Treatment*Herbivore+year_num  + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = reproRM)


Anova(mod_fecundity_numRM,type="III") #signficiant watered x herbivore treatment

summary(mod_fecundity_numRM)

fecundity_num_df_RM <- reproRM %>% 
  
  mutate(fit.m = predict(mod_fecundity_numRM, re.form = NA),
         
         fit.c = predict(mod_fecundity_numRM, re.form = NULL), #all random effects
         
         resid = residuals(mod_fecundity_numRM))


ggplot(fecundity_num_df_RM, aes(x = Treatment, y = fit.c, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Mature Number Siliques") +
  labs(x = "Watering Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 





#*******************************************************************************
#### 2.leaf damage across treatments #####
#*******************************************************************************

#repeated measures with average damage data across years 

field1 <- drop_na(field,avg_LAR) 


##convert LAR (measured as % of leaf area removed by herbivores) to proportion
grasshopper1 $LAR_prop<- grasshopper1 $avg_LAR/100
##Check that the conversion worked
hist(grasshopper1 $LAR_prop)
head(grasshopper1)  
sapply(grasshopper1,class)    


#First, we will use ggplot to look at the data:


LAR_RM= ggplot(grasshopper1, aes(x= elevation,y= LAR_prop, group= Treatment, 
                                 colour= Treatment))+geom_point(size=5) + scale_y_continuous("Leaf area removed by herbivores (%)")+ scale_x_continuous("Source Elevation")  
LAR_RM + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                   axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                                   panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "#882255","#6699cc"), name = "Watering treatment", labels = c("Ample Water","Water Restricted"))



LAR_RM_box<-ggplot(grasshopper1, aes(x = Herbivore, y = LAR_prop, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("field exclosure", "field enclosure")) +  scale_fill_manual(values = c("darkred","lightblue"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))

ggplot(grasshopper1, aes(x = Treatment, y = LAR_prop, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = c("#882255","#6699cc")) +
  labs(y = "Leaf area removed by herbivores (%)") +
  labs(x = "Watering Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

#Gamlss models fail if there are any NA values in the entire dataset. So, I exclude NAs here

grasshopper2 <- dplyr::select(grasshopper1, LAR_prop, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, elev_km, S_elev,block, treat,Cage_Block,avg_vwc,S_year,year)



##We can proceed to the zero inflated beta regression.

gamlss.model<- gamlss (formula=LAR_prop~Treatment* S_elev*Herbivore+S_year+ random(Cage_Block)+random(population), 
                       sigma.formula=LAR_prop~Treatment* S_elev*Herbivore+ random(Cage_Block)+random(Genotype), nu.formula=LAR_prop~Treatment* S_elev*Herbivore+S_year+ random(Cage_Block)+random(population),family=BEZI, data= grasshopper2)
summary(gamlss.model)

#has not converged



visreg(gamlss.model, overlay = TRUE, "S_elev", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Elevation (KM)", ylab="Leaf Area Herbivorized (Scaled)", partial=TRUE,
       fill=list(col=grey
                 (c(0.99), alpha=0)
       ), band = TRUE, gg = TRUE,
       
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

dropterm(gamlss.model)

mod2<-stepGAIC(gamlss.model) # suggests LAR_prop ~ Treatment + S_elev + Herbivore + random(Cage_Block) +    Treatment:S_elev 

mod2$anova

summary(mod2)


# beta transformation

ggplot(grasshopper1, aes(x= LAR_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(year ~  .)

n<-nrow(grasshopper1)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

grasshopper1 $y_beta<- (grasshopper1 $LAR_prop*(n-1) + 0.5)/n

hist(grasshopper1 $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(grasshopper1 $y_beta)

max(grasshopper1 $y_beta)


#Then, the analysis for Treatment/Herbivore

library(betareg)
beta_modela<- betareg ( y_beta ~S_elev*Treatment*Herbivore+S_year, data=grasshopper1)
Anova(beta_modela,type="III")
visreg(beta_modela, "elev_km", scale="response")

visreg(beta_modela, 'Treatment', by= "Herbivore", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=TRUE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=.2,col=c("#6699cc","#882255"))) 

summary(beta_modela)

plot(beta_modela)

grasshopper3 <- dplyr::select(grasshopper1, LAR_prop, y_beta, S_year,year, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, elev_km, S_elev,block, treat,Cage_Block,avg_vwc)


gamlss_moda<- gamlss (formula= y_beta ~S_elev*Treatment*Herbivore+S_year + random(PlantID)+ random(Cage_Block)+random(population),family=BE(mu.link = "logit"), data=grasshopper3,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_moda)
plot(gamlss_moda)
Anova(gamlss_moda)


#pull the fitted values out and plot them in ggplot

newdf1 <- grasshopper3 %>% 
  
  mutate(fit.m = predict(gamlss_moda, re.form = NA),
         
         fit.c = predict(gamlss_moda, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_moda))

##Convert fit.m and fit.c back to the proportional scale

newdf1 $fit.m_trans<-1/(1+exp(-(newdf1 $fit.m))) 

newdf1 $fit.c_trans<-1/(1+exp(-(newdf1 $fit.m))) 

newdf1 $resid_trans<-1/(1+exp(-(newdf1 $resid))) 



LAR_fig =ggplot(newdf1,aes(x= elev_km,y= fit.c_trans,shape= treat_herb, linetype= treat_herb,color= treat_herb, group= treat_herb)) + 
  
  geom_point(aes(shape= treat_herb),size=4)+scale_shape_manual(values = c(0,2,15,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c("dotdash", "dotdash","solid","solid"))+
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","lightblue","darkred"))

LAR_fig


LAR_RM_box<-ggplot(newdf1, aes(x = Herbivore, y = fit.c_trans, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("field exclosure", "field enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))






#####repeated measures with all damage data####


#reformat datafile

LAR_data<- field %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3","LAR_4","LAR_5","LAR_6","LAR_7","LAR_8","LAR_9"),
                                        names_to='census',
                                        values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR,census, elevation, Genotype, PlantID, initial_size, S_initdiam, Garden_Block, elev_km, S_elev, Year,Treatment,garden)

LAR_data$census[LAR_data$census == "LAR_1"] <- "1"
LAR_data$census[LAR_data$census == "LAR_2"] <- "2"
LAR_data$census[LAR_data$census == "LAR_3"] <-"3"
LAR_data$census[LAR_data$census == "LAR_4"] <- "4"
LAR_data$census[LAR_data$census == "LAR_5"] <- "5"
LAR_data$census[LAR_data$census == "LAR_6"] <- "6"
LAR_data$census[LAR_data$census == "LAR_7"] <- "7"
LAR_data$census[LAR_data$census == "LAR_8"] <- "8"
LAR_data$census[LAR_data$census == "LAR_9"] <- "9"

LAR_data $census <-as.factor(LAR_data $census)


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

#plotting LAR

LAR_box <-ggplot(LAR_data, aes(x = garden, y = y_beta, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                             panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c("darkred","lightblue"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

LAR_fig= ggplot(LAR_data, aes(x= elev_km,y= y_beta, group= Treatment, 
                            colour= Treatment))+geom_point(size=2) + scale_y_continuous("Leaf area removed by herbivores (%)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Herbivore treatment", labels = c("Control","Pesticide"))



##### Final(?) LAR Model ####


Mod2<- gamlss (formula= y_beta ~S_elev*Treatment+garden+Year+ random(census) + random(PlantID)+ random(Garden_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(Mod2)
drop1(Mod2)




#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(Mod2, re.form = NA),
         
         fit.c = predict(Mod2, re.form = NULL), #all random effects
         
         resid = residuals(Mod2))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


LAR_box <-ggplot(newdf2, aes(x = Treatment, y = fit.c_trans, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Control", "Pesticide")) +  scale_fill_manual(values = c("darkred","lightblue"))#, name = "Watering treatment", labels = c("Water-restricted","Supplemental watering"))

LAR_fig= ggplot(newdf2, aes(x= elev_km,y= fit.c_trans, group= Treatment, 
                                   colour= Treatment))+geom_point(size=2) + scale_y_continuous("Leaf area removed by herbivores (%)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Herbivore treatment", labels = c("Control","Pesticide"))



#*******************************************************************************
#### 3. traits #####
#*******************************************************************************

##### Flowering time  ####

##### ordinal floweringtime ####


FT_fig =ggplot(field,aes(x= elev_km,y= Date_flowering_2023,shape= treat, linetype= treat,color= treat, group= treat)) + 
  
  geom_point(aes(shape= treat),size=4)+scale_shape_manual(values = c(0,2,15,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Day of Flowering (Ordinal)") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c("dotdash", "dotdash","solid","solid"))+
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","lightblue","darkred"))

FT_fig


FT_box<-ggplot(field, aes(x = garden, y = Date_flowering_2023, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Day of Flowering (Ordinal)") +
  geom_point(pch = 21, position = position_jitterdodge())

FT_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                             panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Herbivore treatment", labels = c("Control","Pesticide"))

# actual model

FT_RM <- lmer(Ordinal_Date_flowering ~ Treatment*Herbivore*S_elev*year_num+(1|PlantID)+(1|Genotype)+(1|Cage_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = field)

plot(FT_RM)

Anova(FT_RM, type = "III") # 

library(ggeffects)

FT_df <- ggpredict(FT_RM, c("S_elev","Treatment","Herbivore"))

pFT_fig =ggplot(FT_df,aes(x= x,y= predicted,shape= group, linetype= group,color= group, group= group)) + 
  
  scale_x_continuous("Elevation (km)")+ scale_y_continuous("Day of Flowering (Ordinal)") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c(1:2))+
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("darkred","lightblue"))+ facet_grid( ~ facet)

pFT_fig


##### Leaf traits ####

###### SLA  ####



SLA_fig =ggplot(field,aes(x= elev_km,y= rosette_SLA_2023,shape= treat, linetype= treat,color= treat, group= treat)) + 
  
  geom_point(aes(shape= treat),size=4)+scale_shape_manual(values = c(0,2,15,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c("dotdash", "dotdash","solid","solid"))+
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","lightblue","darkred"))

SLA_fig


SLA_RM_box<-ggplot(field, aes(x = garden, y = rosette_SLA_2023, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Specific Leaf Area") +
  geom_point(pch = 21, position = position_jitterdodge())


SLA_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Herbivore treatment", labels = c("Control","Pesticide"))





# damage and VWC 
SLA_RM <- lmer(rosette_SLA ~ avg_vwc*LAR_1*S_elev+year_num+(1|PlantID)+(1|population)+(1|Cage_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = field)

plot(SLA_RM)

Anova(SLA_RM, type = "III") # elevation is significant, not year

visreg(SLA_RM,"avg_vwc", by="LAR_1", overlay=TRUE,   scale = "response", xlab="volumetric water content (%)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

visreg(SLA_RM,"LAR_1", overlay=TRUE,   scale = "response", xlab="volumetric water content (%)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))


# treatment and herbivore
SLA_RMa <- lmer(rosette_SLA ~ Water*Herbivore*S_elev+year_num+(1|PlantID) +(1|Genotype)+(1|Cage_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = field)

plot(SLA_RMa)

Anova(SLA_RMa, type = "III") # elevation, year is significant

visreg(SLA_RMa,"S_elev", by="Herbivore", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

SLA_fig= ggplot(field, aes(x= elev_km,y= rosette_SLA, group= treat, 
                            colour= treat))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  
SLA_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                             panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ year_num, scales="free_x") #+scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))


ggplot(field, aes(x = Herbivore, y = rosette_SLA, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("ivory", "#117733")) +
  labs(y = "Probability of Reproduction") +
  labs(x = "Herbivore Treatment") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 


###### LWC  ####



lwc_fig =ggplot(field,aes(x= elev_km,y= rosette_lwc,shape= treat, linetype= treat,color= treat, group= treat)) + 
  
  geom_point(aes(shape= treat),size=4)+scale_shape_manual(values = c(0,2,15,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c("dotdash", "dotdash","solid","solid"))+
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","lightblue","darkred"))

lwc_fig


lwc_RM_box<-ggplot(field, aes(x = Herbivore, y = rosette_lwc, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Specific Leaf Area") +
  geom_point(pch = 21, position = position_jitterdodge())

lwc_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("field exclosure", "field enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))


# damage and VWC 
SLA_RM <- lmer(rosette_lwc ~ avg_vwc*LAR_1*S_elev+year_num+(1|PlantID)+(1|population)+(1|Cage_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = field)

plot(SLA_RM)

Anova(SLA_RM, type = "III") # elevation is significant, not year

visreg(SLA_RM,"avg_vwc", by="LAR_1", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

visreg(SLA_RM,"avg_vwc", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))


# treatment and herbivore
lwc_RMa <- lmer(rosette_lwc ~ Treatment+Herbivore*S_elev*year_num+Herbivore*I(S_elev^2)*year_num+(1|PlantID) +(1|population)+(1|Cage_Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = field)

plot(lwc_RMa)

Anova(lwc_RMa, type = "III") # elevation, year is significant

visreg(lwc_RMa,"S_elev",by="Herbivore", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))



