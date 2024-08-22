######## PROJECT: greenhouse experiment: variation in herbivore damage due to treatment
#### PURPOSE:Examine fitness, traits in response to herbivory.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 28 Jun 24


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


##Change the baselifor offspring
greenhouse $treatment <-factor(greenhouse $treatment, levels = c("Herbivory", "Naïve"))

##Change the baseline for maternal treatment
greenhouse $mat_treat <-factor(greenhouse $mat_treat, levels = c("Herbivory", "Naïve","Parent"))

greenhouse <- filter(greenhouse, mat_treat != "Parent")

##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
greenhouse $S_elev<-scale(greenhouse $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
greenhouse $S_initdiam<-scale(greenhouse $ini_size,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
greenhouse$elev_km<-greenhouse $elevation/1000


#Adjust flowering time based on final_model in flowering_time_adjustment.R
#greenhouse$FT_Adj<-round((greenhouse$Date_flowering_exp - (greenhouse$Silique_length_flowering /0.54955)),0)

#plot(greenhouse$FT_Adj~greenhouse$Date_flowering_exp)

#grasshopper$Snowmelt_FT_Adj<-grasshopper$FT_Adj-grasshopper$Day_of_snowmelt
#plot(grasshopper$Snowmelt_FT_Adj ~grasshopper$Snowmelt_Date_flowering)

#hist(grasshopper$Ordinal_Date_flowering)
#hist(grasshopper$FT_Adj)

#This calculates flowering duration
#grasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $FT_Adj)

#for vegetative traits
foliar<-subset(greenhouse,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

#set colors
cols=c("#882255","#56B4E9")
#


## selection ######

##Exclude season 1 because we only have LAR data from that year and no plants reproduced
greenhouse_no1<-subset(greenhouse, Season!="1")
# retain only those traits to be included in the models;
colnames(greenhouse);

traitdat <- dplyr::select(greenhouse_no1,Season_survival, avg_leafnumber, avg_LAR_2, genotype, treatment, mat_treat, exp_ID, ini_size, block, elevation,elev_km,  Season,Date_flowering_exp, Height1_flowering,Mature_length_siliques,Reproduction, SLA, succulence,max_LAR, Mature_silique_number)


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
traitdat $Selev<-scale(traitdat $elevation,center=TRUE, scale=TRUE)
##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
traitdat $S_initdiam<-scale(traitdat $ini_size,center=TRUE, scale=TRUE)

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)

#traitdat $sleaf<-scale(traitdat $avg_leafnumber,center=TRUE,scale=TRUE)
#traitdat $sLAR_max<-scale(traitdat $max_LAR,center=TRUE,scale=TRUE)
traitdat $sLAR<-scale(traitdat $avg_LAR_2,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)
traitdat $sFT<-scale(traitdat $Date_flowering_exp,center=TRUE,scale=TRUE)
traitdat $sheight<-scale(traitdat $Height1_flowering,center=TRUE,scale=TRUE)


head(traitdat)


#Change baseline for plotting purposes



#*******************************************************************************
####  Probability of reproduction for vegetative traits only  #####
#*******************************************************************************


repro_model_full <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season+
                             #treatment*mat_treat*sLAR+
                             treatment*mat_treat*sSUC+#I(sSUC^2)*treatment+
                             treatment*mat_treat*sSLA+
                             #treatment*mat_treat*sleaf#+I(sleaf^2) 
                           +(1|exp_ID)+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
Anova(repro_model_full,type="III") 

simulationOutput <- simulateResiduals(fittedModel= repro_model_full, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")


repro_model_full_nogeno <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season+ treatment*mat_treat*sLAR+ treatment*mat_treat*sSLA+(1|exp_ID)+(1|block),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

repro_model_full_nopid <- glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season+ treatment*mat_treat*sLAR+ treatment*mat_treat*sSLA+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

repro_model_full_block <- glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season+ treatment*mat_treat*sLAR+ treatment*mat_treat*sSLA+(1|exp_ID)+(1|genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))


anova(repro_model_full,repro_model_full_nogeno)
anova(repro_model_full,repro_model_full_nopid)
anova(repro_model_full,repro_model_full_block)




#SLA, seelction for increased SLA across trt

sSLA_pred <- ggpredict(repro_model_full, terms = c("sSLA[all]"#, "treatment","mat_treat"
                                                   ), type = "re", interval="confidence")

SLA_reproduction_treatment <-plot(sSLA_pred, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")+ ylim(0,1)


emtrends(repro_model_full, c("treatment"), var = "sSLA")

repro <- emtrends(repro_model_full, specs = c("mat_treat","treatment"), var = "sSLA",type="response")
repro_table<- as.data.frame(summary(repro))[c("mat_treat",'treatment','sSLA.trend', 'SE')]
repro_table <- repro_table%>% mutate(
  oddsratio = exp(sSLA.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE))



#LAR treatment effect, slope varies

sLAR_pred <- ggpredict(repro_model_full, terms = c("sLAR[all]","treatment","mat_treat"), type = "re", interval="confidence")


LAR_reproduction_treatment <-plot(sLAR_pred, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Leaf area removed (proportion)")+ scale_y_continuous("Probability of reproduction")+ ylim(0,1)


emtrends(repro_model_full, specs = c("treatment","mat_treat"), var = "sLAR")

repro <- emtrends(repro_model_full, specs = c("treatment","mat_treat"), var = "sLAR",type="response")
repro_table<- as.data.frame(summary(repro))[c('treatment',"mat_treat",'sLAR.trend', 'SE')]
repro_table <- repro_table%>% mutate(
  oddsratio = exp(sLAR.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE))

#Figure5

library(ggpubr)
figure5 <- ggarrange(
  SLA_reproduction_treatment,
  LAR_reproduction_treatment,
  
  labels = c("A", "B"),
  ncol = 2, nrow = 1)
figure5

##Conversation of standardized elevation back to raw units
sd(traitdat $SLA,na.rm=TRUE)
mean(traitdat$SLA,na.rm=TRUE)

sd(traitdat$avg_LAR,na.rm=TRUE)
mean(traitdat$avg_LAR,na.rm=TRUE)


-1*sd(traitdat $sLAR,na.rm=TRUE)+mean(sLAR $elevation,na.rm=TRUE)
0*sd(traitdat $sLAR,na.rm=TRUE)+mean(sLAR $elevation,na.rm=TRUE)
1*sd(traitdat $sLAR,na.rm=TRUE)+mean(sLAR $elevation,na.rm=TRUE)
2*sd(traitdat $sLAR,na.rm=TRUE)+mean(sLAR $elevation,na.rm=TRUE)



#succulence, context dependent, selection for reduced succulence under control

sSUC_pred <- ggpredict(repro_model_full, terms = c("sSUC[all]","treatment"), type = "re", interval="confidence")

Succulence_reproduction <-plot(sSUC_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")

emtrends(repro_model_full, specs = c("treatment"), var = "sSUC")

#repro <- emtrends(repro_model_full, specs = c("treatment"), var = "sSUC",type="response")
#repro_table<- as.data.frame(summary(repro))[c('treatment','sSUC.trend', 'SE')]
#repro_table <- repro_table%>% mutate(
#  oddsratio = exp(sSUC.trend),
#  Lower95 = oddsratio * exp(-1.96*SE),
#  Upper95 = oddsratio * exp(1.96*SE)) #reduced succulence 


#leaf number divergent selection

#sleaf_pred <- ggpredict(repro_model_full, terms = c("sleaf[all]","treatment"), type = "re", interval="confidence")

#leaf_reproduction <-plot(sleaf_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf number")+ scale_y_continuous("Probability of reproduction")

#emtrends(repro_model_full, specs = c("treatment"), var = "sleaf")

#repro <- emtrends(repro_model_full, specs = c("treatment"), var = "sleaf",type="response")
#repro_table<- as.data.frame(summary(repro))[c('treatment','sleaf.trend', 'SE')]
#repro_table <- repro_table%>% mutate(
#  oddsratio = exp(sleaf.trend),
#  Lower95 = oddsratio * exp(-1.96*SE),
#  Upper95 = oddsratio * exp(1.96*SE))




########### Univariate ################
###
###
###

#succulence
                             

repro_model_SUC <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season
                          +sSUC*treatment*mat_treat
                          #+I(sSUC^2)*treatment
                          +(1|exp_ID)+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"))

Anova(repro_model_SUC,type="III") 


simulationOutput <- simulateResiduals(fittedModel= repro_model_SUC, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_succ <- ggpredict(repro_model_SUC, terms = c("sSUC[all]","treatment"), type = "re", interval="confidence")
plot(standard_succ, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Probability of reproduction")


#SLA

repro_model_SLA <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season
                          +sSLA*treatment*mat_treat
                          #+I(sSLA^2)
                          +(1|exp_ID)+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"))

Anova(repro_model_SLA,type="III") 

simulationOutput <- simulateResiduals(fittedModel= repro_model_SLA, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_SLA <- ggpredict(repro_model_SLA, terms = c("sSLA[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_SLA, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")

simulationOutput <- simulateResiduals(fittedModel= repro_model_SLA, plot = T, re.form = NULL,allow.new.levels =TRUE)




#LAR

repro_model_LAR <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season
                          +sLAR+treatment*mat_treat
                          #+I(sLAR^2)+
                          +(1|exp_ID)+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"))

Anova(repro_model_LAR,type="III") 

simulationOutput <- simulateResiduals(fittedModel= repro_model_LAR, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_LAR <- ggpredict(repro_model_LAR, terms = c("sLAR[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_LAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of reproduction")


#Leaf number

repro_model_Leaf <-glmmTMB(Reproduction~S_initdiam+treatment*mat_treat*Season
                          +sleaf+treatment*mat_treat
                          #+I(sleaf^2)+
                          +(1|exp_ID)+(1|block)+(1|genotype),data=traitdat,family=binomial(link="logit"))

simulationOutput <- simulateResiduals(fittedModel= repro_model_Leaf, plot = T, re.form = NULL,allow.new.levels =TRUE)

Anova(repro_model_Leaf,type="III") 

standard_Leaf <- ggpredict(repro_model_Leaf, terms = c("sleaf[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_Leaf, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("average Leaf number")+ scale_y_continuous("Probability of reproduction")





#*******************************************************************************
####  Fecundity  #####
#*******************************************************************************



traitdatRepro <- filter(traitdat, Reproduction == 1 ) 
##Rescale once we have removed the non-reproductive plants

#traitdatRepro $sLAR_max<-scale(traitdatRepro $max_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sLAR<-scale(traitdatRepro $avg_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sSLA<-scale(traitdatRepro $SLA,center=TRUE,scale=TRUE)
traitdatRepro $sSUC<-scale(traitdatRepro $succulence,center=TRUE,scale=TRUE)
traitdatRepro $sFT<-scale(traitdatRepro $Date_flowering_exp,center=TRUE,scale=TRUE)
traitdatRepro $sheight<-scale(traitdatRepro $Height1_flowering,center=TRUE,scale=TRUE)

fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                               treatment*mat_treat* sSLA+
                              # treatment*mat_treat* sLAR+
                               treatment*mat_treat*sSUC+
                               #sFT*treatment*mat_treat+
                               #sheight*treatment*mat_treat
                             
                               +(1|exp_ID)+(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))


#fecund_modeltraits <-glmmTMB(Mature_silique_number~ S_initdiam+treatment*mat_treat*Season +
#                               treatment*mat_treat* sSLA+
                               #treatment*mat_treat* sSUC+
#                               treatment*mat_treat* sLAR+
#                               treatment*mat_treat* sleaf+
#                               sFT*treatment*mat_treat+
#                               sheight*treatment*mat_treat+
#                              (1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))


Anova(fecund_modeltraits,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)




fecund_modeltraits_nogeno <-fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                                                           treatment*mat_treat* sSLA+
                                                           treatment*mat_treat* sLAR+
                                                           sFT*treatment*mat_treat+
                                                           sheight*treatment*mat_treat+ (1|exp_ID)+
                                                           (1|block),data=traitdatRepro,family=Gamma(link="log"))


fecund_modeltraits_nopid <- fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                                                           treatment*mat_treat* sSLA+
                                                           treatment*mat_treat* sLAR+
                                                           sFT*treatment*mat_treat+
                                                           sheight*treatment*mat_treat+
                                                           (1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))


fecund_modeltraits_block <- fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                                                           treatment*mat_treat* sSLA+
                                                           treatment*mat_treat* sLAR+
                                                           sFT*treatment*mat_treat+
                                                           sheight*treatment*mat_treat+ (1|exp_ID)+
                                                           (1|genotype),data=traitdatRepro,family=Gamma(link="log"))


anova(fecund_modeltraits,fecund_modeltraits_nogeno)
anova(fecund_modeltraits,fecund_modeltraits_nopid)
anova(fecund_modeltraits,fecund_modeltraits_block)



##Conversation of standardized elevation back to raw units
sd(traitdatRepro $SLA,na.rm=TRUE)
mean(traitdatRepro$SLA,na.rm=TRUE)

sd(traitdatRepro$avg_LAR,na.rm=TRUE)
mean(traitdatRepro$avg_LAR,na.rm=TRUE)

sd(traitdatRepro$Date_flowering_exp,na.rm=TRUE)
mean(traitdatRepro$Date_flowering_exp,na.rm=TRUE)

sd(traitdatRepro$Height1_flowering,na.rm=TRUE)
mean(traitdatRepro$Height1_flowering,na.rm=TRUE)



#SLA
emtrends(fecund_modeltraits, specs = c("treatment"), var = "sSLA")

#Context dependent selection
standard_sSLA <- ggpredict(fecund_modeltraits, terms = c("sSLA[all]","treatment"), type = "re", interval="confidence")

SLA_fecundity <-plot(standard_sSLA, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme()+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)

emtrends(fecund_modeltraits, specs = c("treatment"), var = "sSLA")

fecund <- emtrends(fecund_modeltraits, specs = c("treatment"), var = "sSLA",type="response")
fecund_table<- as.data.frame(summary(fecund))[c('treatment','sSLA.trend', 'SE')]
fecund_table <- fecund_table%>% mutate(
  oddsratio = exp(sSLA.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE))


#LAR 
emtrends(fecund_modeltraits, specs = c("treatment"), var = "sLAR")

# no selection
standard_sLAR <- ggpredict(fecund_modeltraits, terms = c("sLAR[all]","treatment"), type = "re", interval="confidence")

LAR_fecundity <-plot(standard_sLAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Leaf area removed (proportion)")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)




#Context dependent selection
standard_FT <- ggpredict(fecund_modeltraits, terms = c("sFT[all]","treatment","mat_treat"), type = "re", interval="confidence")

FT_fecundity <-plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("flowering time")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)

emtrends(fecund_modeltraits, specs = c("treatment","mat_treat"), var = "sFT")

fecund <- emtrends(fecund_modeltraits, specs = c("treatment","mat_treat"), var = "sFT",type="response")
fecund_table<- as.data.frame(summary(fecund))[c('treatment',"mat_treat",'sFT.trend', 'SE')]
fecund_table <- fecund_table%>% mutate(
  oddsratio = exp(sFT.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE))


#height
standard_height <- ggpredict(fecund_modeltraits, terms = c("sheight[all]","treatment","mat_treat"), type = "re", interval="confidence")

height_fecundity <-plot(standard_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="bottom")+scale_x_continuous("Tallest stem at flowering")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)

emtrends(fecund_modeltraits, var = "sheight")

fecund <- emtrends(fecund_modeltraits, specs = c("treatment","mat_treat"), var = "sheight",type="response")
fecund_table<- as.data.frame(summary(fecund))[c('treatment',"mat_treat",'sheight.trend', 'SE')]
fecund_table <- fecund_table%>% mutate(
  oddsratio = exp(sheight.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE))


#Figure6

library(ggpubr)
figure6 <- ggarrange(
  SLA_fecundity,
  LAR_fecundity,
  FT_fecundity,
  height_fecundity,
  labels = c("A", "B","C","D"),
  ncol = 2, nrow = 2)
figure6




#Context dependent selection
#standard_leaf <- ggpredict(fecund_modeltraits, terms = c("sleaf[all]","treatment","mat_treat"), type = "re", interval="confidence")

#leaf_fecundity <-plot(standard_leaf, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("flowering time")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)

#emtrends(fecund_modeltraits, specs = c("treatment","mat_treat"), var = "sFT")

#fecund <- emtrends(fecund_modeltraits, specs = c("treatment","mat_treat"), var = "sFT",type="response")
#fecund_table<- as.data.frame(summary(fecund))[c('treatment',"mat_treat",'sFT.trend', 'SE')]
#fecund_table <- fecund_table%>% mutate(
#  oddsratio = exp(sFT.trend),
#  Lower95 = oddsratio * exp(-1.96*SE),
#  Upper95 = oddsratio * exp(1.96*SE))




##No selection
standard_sLAR <- ggpredict(fecund_modeltraits, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")

LAR_fecundity <-plot(standard_sLAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="none")+ ylim(0,1200)+scale_y_continuous("Fecundity (Summed length of siliques)")+scale_x_continuous("Leaf area removed by herbivores (Proportion)")


########### Univariate ################
###
###
###

#succulence no selection

fecund_model_SUC <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                               +sSUC*treatment*mat_treat+
                               #+I(sSUC^2)*treatment*mat_treat+
                               +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_SUC,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_model_SUC, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_succ <- ggpredict(fecund_model_SUC, terms = c("sSUC[all]","treatment"), type = "re", interval="confidence")
plot(standard_succ, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("mature siliques")


#SLA

fecund_model_SLA <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                             +sSLA*treatment*mat_treat+
                             +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_SLA,type="III") 


standard_SLA <- ggpredict(fecund_model_SLA, terms = c("sSLA[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_SLA, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("mature siliques")


#LAR no evidence of selection

fecund_model_LAR <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                             +sLAR*treatment*mat_treat+
                             #+I(sLAR^2)*treatment*mat_treat+
                             +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_LAR,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_model_LAR, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_LAR <- ggpredict(fecund_model_LAR, terms = c("sLAR[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_LAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("mature siliques")


#Leaf number #no evidence of selection

fecund_model_leaf <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                             +sleaf*treatment*mat_treat+
                             
                             +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_leaf,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_model_leaf, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_Leaf <- ggpredict(fecund_model_leaf, terms = c("sleaf[all]","treatment","mat_treat"), type = "re", interval="confidence")
plot(standard_Leaf, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("average Leaf number")+ scale_y_continuous("Probability of reproduction")


#FT #no evidence of selection

fecund_model_FT <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                              +sFT*treatment*mat_treat+
                              +I(sFT^2)*treatment
                              +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_FT,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_model_FT, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_FT <- ggpredict(fecund_model_FT, terms = c("sFT[all]","treatment"), type = "re", interval="confidence")
plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("average Leaf number")+ scale_y_continuous("Probability of reproduction")

standard_FT <- ggpredict(fecund_model_FT, terms = c("sFT[all]","mat_treat"), type = "re", interval="confidence")
plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("average Leaf number")+ scale_y_continuous("Probability of reproduction")

#height 

fecund_model_height <-glmmTMB(Mature_length_siliques~ S_initdiam+treatment*mat_treat*Season +
                            +sheight*treatment*mat_treat+
#                            +I(sheight^2)
                            +(1|block)+(1|genotype),data=traitdatRepro,family=Gamma(link="log"))

Anova(fecund_model_height,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_model_height, plot = T, re.form = NULL,allow.new.levels =TRUE)


standard_height <- ggpredict(fecund_model_height, terms = c("sheight[all]"), type = "re", interval="confidence")
plot(standard_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("average Leaf number")+ scale_y_continuous("Probability of reproduction")

