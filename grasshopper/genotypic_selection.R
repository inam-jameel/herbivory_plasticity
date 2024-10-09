######## PROJECT: Common garden experiment: Examining clines and plasticity in response to water availability and grasshopper abundance
#### PURPOSE: calculate genotypic means for genotypic selection analyses
#### DATE LAST MODIFIED: 6 Oct 2024



## remove objects and clear workspace
rm(list = ls(all=TRUE))


##require packages
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(tidyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for pairwise comparisons
require(glmmTMB) # for running models, have to load twice occasionally
require(gamlss) # for running herbivor damage model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(multcompView) #for pairwise comparisons
require(vioplot) #for violin plots
library(DHARMa) #for assessing model fits
library(ggeffects) #for plotting
library(ggpubr) #for arranging panels into figures



##this is where you specify the folder where you have the data on your computer

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/Grasshopper_manuscript_files/Grasshopper_manuscript_Submission/scripts_data/")

##setwd("~/Documents/personnel/Jameel/grasshopper")

##read in data 
LSmeans <- read.csv("LSmeans_foliar2.csv",stringsAsFactors=T)

##Change the baseline for Water availability
LSmeans $Water <-factor(LSmeans $Water, levels = c("Restricted", "Supplemental"))

##Change the baseline for Herbivore manipulation
LSmeans $Herbivore <-factor(LSmeans $Herbivore, levels = c("Addition", "Removal"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
LSmeans $S_elev<-scale(LSmeans $elevation,center=TRUE, scale=TRUE)


#### probability of reproduction


LSmeans $sLAR<-scale(LSmeans $LAR,center=TRUE,scale=TRUE)
LSmeans $sSLA<-scale(LSmeans $SLA,center=TRUE,scale=TRUE)
LSmeans $sSUC<-scale(LSmeans $succulence,center=TRUE,scale=TRUE)

selection2 <-glmmTMB(cbind(reproduced , overwintered -reproduced )~ year + Water*Herbivore
                    + sLAR*Water*Herbivore #+ I(sLAR^2)
                    + sSLA*Water*Herbivore #+ I(sSLA^2)
                    + sSUC*Water*Herbivore + I(sSUC^2)
                    + (1|Genotype), data= LSmeans, family=binomial(link="logit"))

simulationOutput <- simulateResiduals(fittedModel= selection2, plot = T, re.form = NULL)

Anova(selection2,type="III") 

##selection on leaf succulence

suc_repro = visreg(selection2,"sSUC", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
suc_repro_plot<-ggplot(suc_repro $fit, aes(sSUC, visregFit)) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= LSmeans, aes(sSUC, survived/overwintered),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Leaf succulence (mg/cm2)")+ scale_y_continuous("Probability of reproduction")
suc_repro_plot

coefficients_SUC <- emtrends(selection, specs = c("sSUC"), var = "sSUC",type="response")
SUC_table<- as.data.frame(summary(coefficients_SUC))[c('sSUC.trend', 'SE')]
SUC_table <- SUC_table%>% mutate(
  slopes = exp(sSUC.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))


##To extract coefficients, we need to create the quadratic effects outside of the model
LSmeans$sSLA2<-LSmeans$sSLA*LSmeans$sSLA
LSmeans$sLAR2<-LSmeans$sLAR*LSmeans$sLAR
LSmeans$sSUC2<-LSmeans$sSUC*LSmeans$sSUC


repro_model_quad <-glmmTMB(cbind(reproduced , overwintered -reproduced )~Water*Herbivore+year+
                             Water*Herbivore*sLAR+ 
                             Water*Herbivore*sSUC+ sSUC2 +
                             Water*Herbivore*sSLA#+ sSLA2 
                           +(1|Genotype),data=LSmeans,family=binomial(link="logit"))
Anova(repro_model_quad,type="III") 

##Slope of selection surfaces 

emtrends(repro_model_quad, specs = c("sSLA"), var = "sSLA")


emtrends(repro_model_quad, specs = c("sSUC"), var = "sSUC")


emtrends(repro_model_quad, specs = c("sSLA"), var = "sSLA")

coefficients_succ <- emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSUC",type="response")
Succ_table<- as.data.frame(summary(coefficients_succ))[c('sSUC.trend', 'SE')]
Succ_table <- Succ_table%>% mutate(
  slopes = exp(sSUC.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

coefficients_succ <- emtrends(repro_model_quad, specs = c("sSUC2"), var = "sSUC2",type="response")
Succ_table<- as.data.frame(summary(coefficients_succ))[c('sSUC2.trend', 'SE')]
Succ_table <- Succ_table%>% mutate(
  slopes = exp(sSUC2.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))



#### probability of survival


LSmeans $sLAR<-scale(LSmeans $LAR,center=TRUE,scale=TRUE)
LSmeans $sSLA<-scale(LSmeans $SLA,center=TRUE,scale=TRUE)
LSmeans $sSUC<-scale(LSmeans $succulence,center=TRUE,scale=TRUE)


selection <-glmmTMB(cbind(survived , overwintered -survived )~ year + Water*Herbivore
                    + sLAR*Water*Herbivore #+ I(sLAR^2)
                    + sSLA*Water*Herbivore #+ I(sSLA^2)
                    + sSUC*Water*Herbivore #+ I(sSUC^2)
                    + (1|Genotype), data= LSmeans, family=binomial(link="logit"))

simulationOutput <- simulateResiduals(fittedModel= selection, plot = T, re.form = NULL)

Anova(surv_model,type="III") 


#Succulence 
emtrends(selection, specs = c("sSUC"), var = "sSUC") # succulence is not significant?


#LAR 
emtrends(selection, specs = c("Water"), var = "sLAR") # LAR * water is not significant?

##selection on leaf succulence

suc_repro = visreg(selection,"sSUC", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
suc_repro_plot<-ggplot(suc_repro $fit, aes(sSUC, visregFit)) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= LSmeans, aes(sSUC, survived/overwintered),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Leaf succulence (mg/cm2)")+ scale_y_continuous("Probability of survival")
suc_repro_plot

coefficients_SUC <- emtrends(selection, specs = c("sSUC"), var = "sSUC",type="response")
SUC_table<- as.data.frame(summary(coefficients_SUC))[c('sSUC.trend', 'SE')]
SUC_table <- SUC_table%>% mutate(
  slopes = exp(sSUC.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))



##selection on Leaf area removed

coefficients_LAR <- emtrends(selection, specs = c("Water"), var = "sLAR",type="response")
LAR_table<- as.data.frame(summary(coefficients_LAR))[c('sLAR.trend', 'Water', 'SE')]
LAR_table <- LAR_table%>% mutate(
  slopes = exp(sLAR.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))


LAR_repro<-visregList(visreg(selection,"sLAR", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                visreg(selection,"sSUC", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


LAR_repro<-visreg(selection,"sLAR", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
                      
LAR_repro_plot<-ggplot(LAR_repro $fit, aes(sLAR, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= LSmeans, aes(sLAR, survived/overwintered, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf area removed (percent)")+ scale_y_continuous("Probability of survival",limits=c(0,1))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))#+facet_wrap(~Herbivore)


LAR_repro = visreg(selection,"sLAR", by = "Water", overlay = TRUE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

LAR_repro_plot<-ggplot(LAR_repro $fit, aes(sLAR, visregFit)) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= LSmeans, aes(sLAR, survived/overwintered),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of survival")
LAR_repro_plot
## univariate models


repro_model_SLA <-glmmTMB(cbind(survived , overwintered -survived )~ year+Water*Herbivore +
                                      sLAR*Water*Herbivore + #Water*Herbivore* I(sLAR^2)
                                   +(1|Genotype), data= LSmeans, family=binomial(link="logit"))
Anova(repro_model_SLA,type="III") 

#sla_pred_uni <- ggpredict(repro_model_SLA, terms = c("sSLA[all]"), type = "re", interval="confidence")
#SLA_reproduction_uni <- plot(sla_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=FALSE,  facet=FALSE,dot_alpha=0.75)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")
#SLA_reproduction_uni

repro_model_SUC <-glmmTMB(cbind(survived , overwintered -survived )~Water*Herbivore+year+
                            Water*Herbivore*sSUC+  I(sSUC^2)+
                            (1|Genotype),data=LSmeans,family=binomial(link="logit")) 

#repro_model_SUC <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
#                             Water*Herbivore*sSUC+I(sSUC^2)+ 
#                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model_SUC,type="III") 
#sSUC_pred_uni <- ggpredict(repro_model_SUC, terms = c("sSUC[all]","Herbivore"), type = "re", interval="confidence")

#Succulence_reproduction_herbivore_uni <-plot(sSUC_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols2,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")
#Succulence_reproduction_herbivore_uni


repro_model_LAR <-glmmTMB(cbind(survived , overwintered -survived )~Water*Herbivore+year+
                             Water*Herbivore+sLAR+ I(sLAR ^2)+
                            (1|Genotype),data=LSmeans,family=binomial(link="logit"))                            

#Anova(repro_model_LAR,type="III") 



#*******************************************************************************
#### Calculating LSmeans #######
#*******************************************************************************


##read in data 
grasshopper <- read.csv("Common_garden_experiment_data.csv",stringsAsFactors=T)

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper$year<-as.factor(grasshopper$year)

##Change the baseline for Water availability
grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))

##Change the baseline for Herbivore manipulation
grasshopper $Herbivore <-factor(grasshopper $Herbivore, levels = c("Addition", "Removal"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)


##This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

##Elevational distance in km
grasshopper$elev_dist_km<-grasshopper $elev_dist/1000

##Adjust flowering time based on final_model in flowering_time_adjustment.R
grasshopper$FT_Adj<-round((grasshopper$Ordinal_Date_flowering - (grasshopper$Silique_length_flowering /2.5901)),0)

plot(grasshopper$FT_Adj~grasshopper$Ordinal_Date_flowering)

hist(grasshopper$Ordinal_Date_flowering)
hist(grasshopper$FT_Adj)

##Correlate rosette and stem leaf data to come up with composite figures
#plot(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
#plot(grasshopper$rosette_SLA~grasshopper$bolt_SLA)

#mod1<-lm(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
#summary(mod1)
#Anova(mod1,type="III")

#SLA_Mod<-lm(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
#summary(SLA_Mod)


##Create composite leaf SLA and succulence based on the regressions above. This gives us foliar trait data for the 67 plants for which we have stem but not rosette collections
grasshopper$SLA <- ifelse(is.na(grasshopper$rosette_SLA), (71.8177 + 0.7319*grasshopper$bolt_SLA), grasshopper$rosette_SLA)
grasshopper$succulence <- ifelse(is.na(grasshopper$rosette_succulence), (0.0023948 + 0.5609208*grasshopper$bolt_succulence), grasshopper$rosette_succulence)


##This calculates flowering duration
grasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $FT_Adj)

##create a dataframe for the vegetative traits
foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)
foliar $env<-interaction(foliar $Water, foliar $Herbivore)


##Create a dataframe for reproductive phenology traits which filters out the two plants that flowered without vernalization in the first season.
grasshopperFT <- filter(grasshopper, year != "2021")
##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(grasshopperFT, Ordinal_Date_flowering!="NA")

##set colors
cols=c("#CC79A7","lightblue") #for water treatment
cols2=c("#882255","#77DDAA") #for the grasshopper treatment


# aggregate number over wintered, survived, and reproduced

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")

overwintered <- aggregate(grasshopper_no2021$Overwinter_survival_2022, list(Genotype = grasshopper_no2021$Genotype, year = grasshopper_no2021$year,  Water = grasshopper_no2021$Water,  Herbivore = grasshopper_no2021$Herbivore), FUN=sum)

survived <- aggregate(grasshopper_no2021$Season_survival, list(Genotype = grasshopper_no2021$Genotype, year = grasshopper_no2021$year,  Water = grasshopper_no2021$Water,  Herbivore = grasshopper_no2021$Herbivore), FUN=sum)

reproduced <- aggregate(grasshopper_no2021$Reproduced, list(Genotype = grasshopper_no2021$Genotype, year = grasshopper_no2021$year,  Water = grasshopper_no2021$Water,  Herbivore = grasshopper_no2021$Herbivore), FUN=sum) 

overwintered_surv <- merge(overwintered, survived, by.x=c('Genotype', 'year', 'Water', 'Herbivore'), 
      by.y=c('Genotype', 'year', 'Water', 'Herbivore')) 

ps_repro <- merge(plant_surv, reproduced, by.x=c('Genotype', 'year', 'Water', 'Herbivore'), 
                    by.y=c('Genotype', 'year', 'Water', 'Herbivore')) 

write.csv(ps_repro,file="summary_planted.csv")


#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")

grasshopper_no2021$avg_LAR_prop<-grasshopper_no2021 $avg_LAR/100 #convert to proportion

grasshopper_no2021 <- drop_na(grasshopper_no2021,avg_LAR_prop) 

##this is the beta transformation, which transforms all values of 0 to a small value.
n<-nrow(grasshopper_no2021)
grasshopper_no2021 $y_beta<- (grasshopper_no2021 $avg_LAR_prop*(n-1) + 0.5)/n

hist(grasshopper_no2021 $y_beta)

min(grasshopper_no2021 $y_beta)

max(grasshopper_no2021 $y_beta)

LAR_model <- glmmTMB(y_beta ~ Water*Herbivore*Genotype*year
                     +(1|Cage_Block), data = grasshopper_no2021, 
                     family=beta_family())

lsmeans_LAR <-emmeans (LAR_model, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

LAR <- as.data.frame(lsmeans_LAR)[c('Genotype', 'year', 'Water', 'Herbivore','response', 'SE')]


test <- merge(ps_repro, LAR, by.x=c('Genotype', 'year', 'Water', 'Herbivore'), 
                  by.y=c('Genotype', 'year', 'Water', 'Herbivore')) 

#write.csv(lsmeans_LAR,file="LSmeans_LAR.csv")


#*******************************************************************************
#### SLA #######
#*******************************************************************************
#*

sla_model_2 <- glmmTMB(SLA ~ Water*Herbivore*Genotype*year+(1|Cage_Block), data = foliar, family= lognormal(link="log"),control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))

lsmeans_SLA <-emmeans (sla_model_2, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_SLA,file="LSmeans_SLA.csv")


SLA <- as.data.frame(lsmeans_SLA)[c('Genotype', 'year', 'Water', 'Herbivore','response', 'SE')]


test2 <- merge(test, SLA, by.x=c('Genotype', 'year', 'Water', 'Herbivore'), 
              by.y=c('Genotype', 'year', 'Water', 'Herbivore')) 


#*******************************************************************************
#### Succulence #######
#*******************************************************************************


foliar$succulence_mg<-foliar$succulence*1000 #transform succulence values for model convergence

succulence_data <- dplyr::select(foliar, succulence, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, Cage_Block, elev_km, S_elev,  year, elev_dist_km) #create separate data frame with relevant data columns

succulence_data <- drop_na(succulence_data,succulence) #remove missing data

hist(succulence_data$succulence)

#beta transformation
n<-nrow(succulence_data)
succulence_data $y_beta<- (succulence_data $succulence*(n-1) + 0.5)/n

succ_model <- glmmTMB(y_beta ~ Water*Herbivore*Genotype*year
                      #+(1|PlantID)
                      #+(1|Genotype)
                      +(1|Cage_Block), data = succulence_data, 
                      family=beta_family(),
                      control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))

lsmeans_succulence <-emmeans (succ_model, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_succulence,file="LSmeans_succulence.csv")

succ <- as.data.frame(lsmeans_succulence)[c('Genotype', 'year', 'Water', 'Herbivore','response', 'SE')]


test3 <- merge(test2, succ, by.x=c('Genotype', 'year', 'Water', 'Herbivore'), 
               by.y=c('Genotype', 'year', 'Water', 'Herbivore')) 

write.csv(test3,file="LSmeans_foliar2.csv")






