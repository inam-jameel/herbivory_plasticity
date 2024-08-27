### Purpose: Phenotypic selection on traits
### Author: Jill Anderson

rm(list = ls(all=TRUE))


library(lme4)
library(ggplot2)
library(car)
library(visreg)
library(effects)
library(glmmTMB)
library(bbmle)
library(dplyr)
library(tidyr)
library(ggeffects)
library(scales)
library(DHARMa)
library(MuMIn)

cols=c("#D55E00", "blue")

##### read in trait data #####
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")
#setwd("~/Documents/personnel/Jameel/grasshopper")

#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_updated_8March24.csv", stringsAsFactors = TRUE)

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper$year<-as.factor(grasshopper$year)

##Change the baseline for treatment
grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))

#Let's concatenate herbivore and watering treatments, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Water,sep = "_")



#This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000
##Let's correlate rosette and bolt leaf data to see if we can come up with composite figures
plot(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
plot(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
plot(grasshopper$rosette_lwc~grasshopper$bolt_lwc)

mod1<-lm(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
summary(mod1)

mod2<-lm(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
summary(mod2)

mod3<-lm(grasshopper$rosette_lwc~grasshopper$bolt_lwc)
summary(mod3)

##Create composite leaf SLA, succuclence and lwc variables based on the regressions above. This gives us foliar trait data for the 67 plants for which we have bolt but not rosette collections
grasshopper$SLA <- ifelse(is.na(grasshopper$rosette_SLA), (71.8177 + 0.7319*grasshopper$bolt_SLA), grasshopper$rosette_SLA)
grasshopper$succulence <- ifelse(is.na(grasshopper$rosette_succulence), (0.0023948 + 0.5609208*grasshopper$bolt_succulence), grasshopper$rosette_succulence)
grasshopper$lwc <- ifelse(is.na(grasshopper$rosette_lwc), (0.19571 + 0.67021*grasshopper$bolt_lwc), grasshopper$rosette_lwc)

#Adjust flowering time based on final_model in flowering_time_adjustment.R
grasshopper$FT_Adj<-round((grasshopper$Ordinal_Date_flowering - (grasshopper$Silique_length_flowering /2.5901)),0)

plot(grasshopper$FT_Adj~grasshopper$Ordinal_Date_flowering)

grasshopper$Snowmelt_FT_Adj<-grasshopper$FT_Adj-grasshopper$Day_of_snowmelt
plot(grasshopper$Snowmelt_FT_Adj ~grasshopper$Snowmelt_Date_flowering)


#This calculates flowering duration
grasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $FT_Adj)
grasshopper <- filter(grasshopper, Exclude == "Include")

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced_updated
grasshopper_no2021<-subset(grasshopper, year!="2021")
# retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper,FT_Adj,Max_height_flowering, Max_height_peak ,Season_survival, avg_leafnumber, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block, elevation,elev_km,  year,Stem_number_flowering, rosette_succulence,Ordinal_Date_flowering,Sum_height_flowering,Date_peak_flowering,Stem_number_peak,Mature_length_siliques,Mature_length_siliques_updated,Reproduced_updated,Reproduced,flowering_duration,days_until_mortality, SLA, lwc, succulence, rosette_SLA, rosette_lwc,FT_Adj,Snowmelt_FT_Adj, treat,max_LAR) #added max LAR observed in a year to the dataset

##only plants that Reproduced_updated
Reproduced_updated<-subset(traitdat,Reproduced_updated=="1")
ggplot(Reproduced_updated, aes(x= Mature_length_siliques_updated))+ geom_histogram()+ facet_grid(treat ~ .)
ggplot(Reproduced_updated, aes(x= Mature_length_siliques_updated))+ geom_histogram(color="black", fill="white")


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
traitdat $Selev<-scale(traitdat $elevation,center=TRUE, scale=TRUE)
##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
traitdat $S_initdiam<-scale(traitdat $init.diam,center=TRUE, scale=TRUE)

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)
traitdat $sFP<-scale(traitdat $Ordinal_Date_flowering,center=TRUE,scale=TRUE)
traitdat $sduration<-scale(traitdat $flowering_duration,center=TRUE,scale=TRUE)
traitdat $sheight<-scale(traitdat $Sum_height_flowering,center=TRUE,scale=TRUE)
traitdat $s_max_height<-scale(traitdat $Max_height_flowering,center=TRUE,scale=TRUE)
traitdat $s_max_height_peak<-scale(traitdat $Max_height_peak,center=TRUE,scale=TRUE)
traitdat $speak<-scale(traitdat $Date_peak_flowering,center=TRUE,scale=TRUE)
traitdat $sleaf<-scale(traitdat $avg_leafnumber,center=TRUE,scale=TRUE)

traitdat $sLAR<-scale(traitdat $avg_LAR,center=TRUE,scale=TRUE)
traitdat $sSucc<-scale(traitdat $rosette_succulence,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sLWC<-scale(traitdat $lwc,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)
traitdat $sFT<-scale(traitdat $FT_Adj,center=TRUE,scale=TRUE)

head(traitdat)


#Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Removal","Addition"))

#*******************************************************************************
####  Probability of reproduction for vegetative traits only  #####
#*******************************************************************************


#repro_model_full_SUC <-glmmTMB(Reproduced~S_initdiam
#                               +Water*Herbivore*sLAR*year + Water*Herbivore* I(sLAR^2)
#                               +Water*Herbivore*sSUC + I(sSUC^2)
#                               +Water*Herbivore*sSLA + I(sSLA^2)   
#                               +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

#repro_model_full_SUC_a <-glmmTMB(Reproduced_updated~S_initdiam+year+
#                                   Water*Herbivore*sLAR+I(sLAR^2)+
#                                   Water*Herbivore*sSUC+I(sSUC^2)+
#                                   Water*Herbivore*sSLA+I(sSLA^2) 
#                                 +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))



repro_model_full_SUC_b <-glmmTMB(Reproduced_updated~S_initdiam+
                                   Water*Herbivore*sLAR*year+Herbivore*I(sLAR^2)+
                                   Water*Herbivore*sSUC+I(sSUC^2)+
                                   Water*Herbivore*sSLA+I(sSLA^2) 
                                 +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


model.sel(repro_model_full_SUC_a,repro_model_full_SUC_b) #repro_model_full_SUC_b has better AIC value


Anova(repro_model_full_SUC_b,type="III")


simulationOutput <- simulateResiduals(fittedModel= repro_model_full_SUC_b, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")





full_sSUC <- ggpredict(repro_model_full_SUC_b, terms = c("sSUC[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(full_sSUC, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols,facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Probability of reproduction")

full_sla <- ggpredict(repro_model_full_SUC_b, terms = c("sSLA[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(full_sla, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")

full_LAR <- ggpredict(repro_model_full_SUC_b, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(full_LAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area herbivorized")+ scale_y_continuous("Probability of reproduction")


### Univariate models ####

#succulence

#repro_model_SUC_1 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSUC + Water*Herbivore+(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

#model.sel(repro_model_SUC_1,repro_model_SUC_2) 


#repro_model_SUC_2 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSUC*year  +Water*Herbivore*year +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

#repro_model_SUC_3b <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSUC*year +Water*Herbivore*year+
#                             +I(sSUC^2)
#                             +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


repro_model_SUC_3a <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSUC*year #+Water*Herbivore*year+
                           +I(sSUC^2)
                          +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


simulationOutput <- simulateResiduals(fittedModel= repro_model_SUC, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")

Anova(repro_model_SUC_3,type="III") 

#model.sel(repro_model_SUC_3a,repro_model_SUC_3b) #repro_model_SUC_3 has a higher AIC, including  #+Water*Herbivore*year does not effect aic



#standard_succ1 <- ggpredict(repro_model_SUC_3a, terms = c("sSUC[all]", "Water","Herbivore"), type = "re", interval="confidence")
#plot(standard_succ1, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Probability of reproduction")

#standard_succ2 <- ggpredict(repro_model_SUC_3a, terms = c("sSUC[all]", "Water","Herbivore","year"), type = "re", interval="confidence")
#plot(standard_succ2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Probability of reproduction")



#SLA


#repro_model_SLA_2 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSLA + #Herbivore*year
#                          +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


#repro_model_SLA_1 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSLA + Herbivore*year
#                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

#repro_model_SLA_3a <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSLA + Herbivore*year
#                           +Water*Herbivore*I(sSLA^2)
#                          +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


repro_model_SLA_3b <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSLA + Herbivore*year
                            +I(sSLA^2)
                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

#repro_model_SLA_3b2 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sSLA
#                           +I(sSLA^2)
#                          +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


#Anova(repro_model_SLA_3b2,type="III") 
Anova(repro_model_SLA_3b,type="III") 



model.sel(repro_model_SLA_3b,repro_model_SLA_3b2) #repro_model_SLA_3b has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= repro_model_SLA_3b, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")


standard_sla <- ggpredict(repro_model_SLA_3b, terms = c("sSLA[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_sla, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")

#LAR

##Here,we have to exclude 2021 because there were only 2 plants who flowered in that year
traitdat_no2021<-subset(traitdat,year!="2021") 

#repro_model_LAR_1 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sLAR*year  +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat_no2021,family=binomial(link="logit"))


#repro_model_LAR_2 <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sLAR+Herbivore*year +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat_no2021,family=binomial(link="logit"))


repro_model_LAR_3a <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sLAR+Herbivore*year 
                            +I(sLAR^2) 
                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat_no2021,family=binomial(link="logit"))


#repro_model_LAR_3b <-glmmTMB(Reproduced_updated~S_initdiam+Water*Herbivore*sLAR+year+ 
#                               I(sLAR^2) 
#                             +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat_no2021,family=binomial(link="logit"))


Anova(repro_model_LAR_3a,type="III") 
#Anova(repro_model_LAR_3b,type="III") 

model.sel(repro_model_LAR_3a,repro_model_LAR_3b) #repro_model_SLA_3a has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= repro_model_LAR_3a, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")



standard_lar <- ggpredict(repro_model_LAR_3a, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_lar, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of reproduction")


standard_lar2 <- ggpredict(repro_model_LAR_3a, terms = c("sLAR[all]", "Water","Herbivore","year"), type = "re", interval="confidence")
plot(standard_lar2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of reproduction")




#*******************************************************************************
####  Fecundity  #####
#*******************************************************************************


traitdatRepro <- filter(traitdat, Reproduced_updated == 1 ) 


fecund_modeltraits_a <-glmmTMB(Mature_length_siliques_updated~ S_initdiam+year+
                              Water*Herbivore* sFT  
                             +Water*Herbivore* s_max_height 
                             +Water*Herbivore* sduration +Water*I(sduration^2)
                             +Water*Herbivore* sSLA + Water*Herbivore* I(sSLA^2)
                             +Water*Herbivore* sSUC 
                             +Water*Herbivore* sLAR 
                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#fecund_modeltraits_b <-glmmTMB(Mature_length_siliques_updated~ S_initdiam+year+
#                               Water*Herbivore* sFT  
#                             +Water*Herbivore* s_max_height 
#                             +Water*Herbivore* sduration #+Water*I(sduration^2)
#                             +Water*Herbivore* sSLA #+ Water*Herbivore* I(sSLA^2)
#                             +Water*Herbivore* sSUC #+ I(sSUC^2)
#                             +Water*Herbivore* sLAR #+ I(sLAR^2)
#                             +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))


Anova(fecund_modeltraits_a,type="III") 
#Anova(fecund_modeltraits_b,type="III") 

#model.sel(fecund_modeltraits_a,fecund_modeltraits_b) #fecund_modeltraits_a has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits_a, plot = T, re.form = NULL) # residuals have a significant quantile test
testOutliers(simulationOutput, type="bootstrap")




standard_FT <- ggpredict(fecund_modeltraits_a, terms = c("sFT[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering time")+ scale_y_continuous("Fecundity")

standard_duration<- ggpredict(fecund_modeltraits_a, terms = c("sduration[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_duration, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering duration")+ scale_y_continuous("Fecundity")+ylim(0,2000)

standard_sSLA <- ggpredict(fecund_modeltraits_a, terms = c("sSLA[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_sSLA, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Fecundity")

standard_sSUC <- ggpredict(fecund_modeltraits_a, terms = c("sSUC[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_sSUC, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Fecundity")


standard_sLAR <- ggpredict(fecund_modeltraits_a, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_sLAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("LAR")+ scale_y_continuous("Fecundity")



#### Univariate ####

#tallest bolt height at flowering

#fecund_model_maxheight_1 <-glmmTMB(Mature_length_siliques_updated~ S_initdiam+
#                                     Water*Herbivore*s_max_height
#                                   +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

fecund_model_maxheight_2 <-glmmTMB(Mature_length_siliques_updated~ S_initdiam+year+
                                   Water*Herbivore*s_max_height
                                 +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#Anova(fecund_model_maxheight_1,type="III") 
Anova(fecund_model_maxheight_2,type="III") 

#model.sel(fecund_model_maxheight_1,fecund_model_maxheight_2) #repro_model_SLA_3a has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= fecund_model_maxheight_2, plot = T, re.form = NULL) # residuals have a significant quantile test
testOutliers(simulationOutput, type="bootstrap")

standard_height <- ggpredict(fecund_model_maxheight_2, terms = c("s_max_height[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Height at flowering")+ scale_y_continuous("Probability of reproduction")+ylim(0,1500)


#flowering duration

fecund_model_dur_a <-glmmTMB(Mature_length_siliques_updated~S_initdiam+
                             Water*Herbivore*sduration+year 
                             +Water*I(sduration^2)
                             +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#fecund_model_dur_b <-glmmTMB(Mature_length_siliques_updated~S_initdiam+
#                            Water*Herbivore*sduration+Water*year 
#                            +Water*I(sduration^2)
#                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))


Anova(fecund_model_dur_a,type="III") 
#Anova(fecund_model_dur_b,type="III") 

#model.sel(fecund_model_dur_a,fecund_model_dur_b) #fecund_model_dur_a has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= fecund_model_dur_a, plot = T, re.form = NULL) # residuals have a significant quantile test
testOutliers(simulationOutput, type="bootstrap")


standard_dur <- ggpredict(fecund_model_dur_a, terms = c("sduration[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_dur, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+ theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering duration")+ scale_y_continuous("Fecundity (silique length)")



fecund_model_FT_a <-glmmTMB(Mature_length_siliques_updated~S_initdiam+
                              Water*Herbivore*sFT*year
                          +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#fecund_model_FT_b <-glmmTMB(Mature_length_siliques_updated~S_initdiam
#                              Water*Herbivore*sFT*year
#                           +I(sFT^2)
#                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))


Anova(fecund_model_FT_a,type="III") 
#Anova(fecund_model_FT_b,type="III") 

#model.sel(fecund_model_FT_a,fecund_model_FT_b) #fecund_model_dur_a has a higher AIC

simulationOutput <- simulateResiduals(fittedModel= fecund_model_FT_b, plot = T, re.form = NULL) # residuals have a significant quantile test
testOutliers(simulationOutput, type="bootstrap")

standard_FTa <- ggpredict(fecund_model_FT_a, terms = c("sFT[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_FTa, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+ theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering phenology")+ scale_y_continuous("Fecundity (silique length)")

#standard_FTb <- ggpredict(fecund_model_FT_b, terms = c("sFT[all]", "Water","Herbivore"), type = "re", interval="confidence")
#plot(standard_FTb, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+ theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering phenology")+ scale_y_continuous("Fecundity (silique length)")

