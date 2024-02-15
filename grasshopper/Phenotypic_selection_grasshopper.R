### Purpose: Phenotypic selection on traits
### Author: Inam Jameel
### uses code originally written by Jill Anderson


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


##### read in trait data #####
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")

#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_18Jan24.csv")

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Genotype<-as.factor(grasshopper$Genotype)
grasshopper$population<-as.factor(grasshopper$population)
grasshopper$Water<-as.factor(grasshopper$Water)
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Herbivore<-as.factor(grasshopper$Herbivore)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper $PlantID <-as.factor(grasshopper $PlantID)
grasshopper $Cage_Block <-as.factor(grasshopper $Cage_Block)

##Change the baseline for treatment
grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)

#This makes number of year as factor
grasshopper $year_num<-as.factor(grasshopper $year_num)

#This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

#This calculates flowering duration
grasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $Ordinal_Date_flowering)

#Let's concatenate herbivore and watering treatments, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Water,sep = "_")


grasshopper <- filter(grasshopper, Exclude == "Include")



# retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper, avg_leafnumber, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev, treat, year_num, rosette_succulence,Ordinal_Date_flowering,Sum_height_flowering,Date_peak_flowering,Stem_number_peak,Mature_length_siliques,Reproduced,flowering_duration,days_until_mortality)



##Let's plot histograms of each trait

ggplot(traitdat, aes(x= Sum_height_flowering))+ geom_histogram(color="black", fill="white")+ facet_grid(treat ~  .)

ggplot(traitdat, aes(x= rosette_succulence))+ geom_histogram(color="black", fill="white")+ facet_grid(treat ~  .)


ggplot(traitdat, aes(x= Stem_number_peak))+ geom_histogram()+ facet_grid(treat ~ .)

ggplot(traitdat, aes(x= Ordinal_Date_flowering,))+ geom_histogram()+ facet_grid(treat ~ .)
ggplot(traitdat, aes(x= Ordinal_Date_flowering))+ geom_histogram(color="black", fill="white")

ggplot(traitdat, aes(x= flowering_duration))+ geom_histogram()+ facet_grid(treat ~ .)
ggplot(traitdat, aes(x= flowering_duration))+ geom_histogram(color="black", fill="white")

ggplot(traitdat, aes(x= Mature_length_siliques))+ geom_histogram()+ facet_grid(treat ~ .)
ggplot(traitdat, aes(x= Mature_length_siliques))+ geom_histogram(color="black", fill="white")


## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)


traitdat $sFP<-scale(traitdat $Ordinal_Date_flowering,center=TRUE,scale=TRUE)
traitdat $sduration<-scale(traitdat $flowering_duration,center=TRUE,scale=TRUE)
traitdat $sheight<-scale(traitdat $Sum_height_flowering,center=TRUE,scale=TRUE)
traitdat $speak<-scale(traitdat $Date_peak_flowering,center=TRUE,scale=TRUE)
traitdat $sleaf<-scale(traitdat $avg_leafnumber,center=TRUE,scale=TRUE)

traitdat $sLAR<-scale(traitdat $avg_LAR,center=TRUE,scale=TRUE)
traitdat $sSucc<-scale(traitdat $rosette_succulence,center=TRUE,scale=TRUE)



head(traitdat)


#Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Removal","Addition"))

traitdat $treat <-factor(traitdat $treat, levels = c("Removal_Restricted","Removal_Supplemental","Addition_Restricted", "Addition_Supplemental"))

traitdatRepro <- filter(traitdat, Reproduced == 1 ) #for fecundity models


########## leaf succulence  ###############

#prob repro
succ_fig_r= ggplot(traitdat, aes(x= rosette_succulence,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("probability of reproduction")+ scale_x_continuous("Rosette Succulence")  

succ_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

succ_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#over time
succ_fig= ggplot(traitdat, aes(x= rosette_succulence,y= Reproduced, group= treat, colour= treat))+geom_point(size=2) + scale_y_continuous("probability of reproduction")+ scale_x_continuous("Source Elevation")  
succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ year_num, scales="free_x") +scale_colour_manual(values = c( "darkorange","lightblue","darkred","#56B4E9"), name = "Water treatment", labels = c("Restricted, Removal","Supplemental, Removal","Restricted, addition","Supplemental, addition"))


succ_fig= ggplot(traitdat, aes(x= rosette_succulence,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Rosette Succulence")  

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))


#over time
succ_fig= ggplot(traitdat, aes(x= rosette_succulence,y= Mature_length_siliques, group= treat, colour= treat))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),  panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ year_num, scales="free_x") +scale_colour_manual(values = c( "darkorange","lightblue","darkred","#56B4E9"), name = "Water treatment", labels = c("Restricted, Removal","Supplemental, Removal","Restricted, addition","Supplemental, addition"))


#only reproduced
succ_fig= ggplot(traitdatRepro, aes(x= rosette_succulence,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Rosette Succulence")

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#over time
succ_fig= ggplot(traitdatRepro, aes(x= rosette_succulence,y= Mature_length_siliques, group= treat, colour= treat))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Source Elevation")  

succ_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),  panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ year_num, scales="free_x") +scale_colour_manual(values = c( "darkorange","lightblue","darkred","#56B4E9"), name = "Water treatment", labels = c("Restricted, Removal","Supplemental, Removal","Restricted, addition","Supplemental, addition"))


repro_succ <-glmmTMB(Reproduced~Water*Herbivore*sSucc
                              +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

Anova(repro_succ,type="III") # significant effect water, herbivore,  but not succulence

#no zeros
fecund_succ <-glmmTMB (Mature_length_siliques ~S_elev+year_num
                                +Water*Herbivore*succ 
                                #+Water*Herbivore*I(sSucc^2) #quad didnt work
                                + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

Anova(repro_succ,type="III") # significant effect water, herbivore,  but not succulence


########## leaf area removed  ###############

LAR_fig_r= ggplot(traitdat, aes(x= avg_LAR,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Leaf Area Removed (%)")  

LAR_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ year, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

LAR_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))




LAR_fig= ggplot(traitdat, aes(x= avg_LAR,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Leaf Area Removed (%)")  

LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#only reproduced

LAR_fig= ggplot(traitdatRepro, aes(x= avg_LAR,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Leaf Area Removed (%)")  

LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

repro_LAR <-glmmTMB(Reproduced~#S_elev+year_num+
                     Water*Herbivore*sLAR
                     +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

Anova(repro_LAR,type="III") # significant effect water,  but not LAR

fecund_LAR <-glmmTMB (Mature_length_siliques ~#S_elev+year_num+
                       Water*Herbivore*sLAR 
                       #+Water*Herbivore*I(sLAR^2) #quad didnt work
                       + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

Anova(fecund_LAR,type="III") #significant LAR

visreg(fecund_LAR,"sLAR", overlay=FALSE, scale = "response", xlab="Leaf Area Removed (standardized)", ylab="silique length", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))

########## first flowering  ###############

flower_fig_r = ggplot(traitdat, aes(x= Ordinal_Date_flowering,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("First day of flowering (ordinal)") 

flower_fig_r + theme_bw() + ylim (0,1)+ theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

flower_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



flower_fig= ggplot(traitdat, aes(x= Ordinal_Date_flowering,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("First day of flowering (ordinal)")  

flower_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

flower_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#only reproduced

flower_fig= ggplot(traitdatRepro, aes(x= Ordinal_Date_flowering,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("First day of flowering (ordinal)")  

flower_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

flower_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))


fecund_flowering <-glmmTMB (Mature_length_siliques ~#S_elev+year_num+
                      Water*Herbivore*sFP 
                      #+Water*Herbivore*I(speak^2) #quad not significant
                      + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

Anova(fecund_flowering,type="III") #significant flowering

visreg(fecund_flowering,"sFP", overlay=FALSE, scale = "response", xlab="Day of Flowering (standardized)", ylab="silique length", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))

########## peak flowering  ###############

peak_fig_r= ggplot(traitdat, aes(x= Date_peak_flowering,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Peak day of flowering (ordinal)")  

peak_fig_r + theme_bw() + ylim (0,1)+ theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

peak_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



peak_fig= ggplot(traitdat, aes(x= Date_peak_flowering,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Peak day of flowering (ordinal)")  

peak_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

peak_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#only reproduced

peak_fig= ggplot(traitdatRepro, aes(x= Date_peak_flowering,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Peak day of flowering (ordinal)")  

peak_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

peak_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

fecund_peak <-glmmTMB (Mature_length_siliques ~#S_elev+year_num+
                            Water*Herbivore*speak 
                            #+Water*Herbivore*I(speak^2) #quad not significant
                            + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

#Warning message:
#  In fitTMB(TMBStruc) :
#  Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')

Anova(fecund_peak,type="III")

#visreg(fecund_peak,"speak", overlay=FALSE, scale = "response", xlab="Day of Flowering (standardized)", ylab="silique length", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))



##########  flowering duration  ###############

duration_fig_r= ggplot(traitdat, aes(x= flowering_duration,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Duration of flowering (days)")  

duration_fig_r + theme_bw() + ylim(0,1) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

duration_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



duration_fig= ggplot(traitdat, aes(x= flowering_duration,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Duration of flowering (days)")  

duration_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

duration_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#only reproduced

duration_fig= ggplot(traitdatRepro, aes(x= flowering_duration,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Duration of flowering (days)")  

duration_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

duration_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



repro_duration <-glmmTMB (Mature_length_siliques ~#S_elev+year_num+
                            Water*Herbivore*sduration 
                            #+Water*Herbivore*I(sduration^2) #quad single term is almost significant
                            + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

Anova(repro_duration,type="III") #water and herbivory significant 


##########  height at flowering  ###############

height_fig_r= ggplot(traitdat, aes(x= Sum_height_flowering,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Summed hieght at flowering (ordinal)")  

height_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+ylim (0,1) + geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

height_fig_r + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



height_fig= ggplot(traitdat, aes(x=Sum_height_flowering ,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Summed height at flowering")  

height_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

height_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))


# only reproduced
height_fig= ggplot(traitdatRepro, aes(x=Sum_height_flowering ,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Summed height at flowering")  

height_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

height_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))



repro_height <-glmmTMB (Mature_length_siliques ~#S_elev+year_num+
                          Water*Herbivore*sheight 
                          #+Water*Herbivore*I(sheight^2) #quad didnt work
                          + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)

Anova(repro_height,type="III") #water and herbivory significant 


########## Multivariate   ###############

####logistic regression####
#using all traits are not related to flowering

mod_repro_selection <-glmmTMB(Reproduced~year_num
                              +Water*Herbivore*sLAR
                              +Water*Herbivore*sSucc
                              +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))


Anova(mod_repro_selection,type="III") # significant effect of water, but nothing else


#####Linear regression####
traitdatRepro <- filter(traitdat, Reproduced == 1 )


mod_fecundityRM <-glmmTMB (Mature_length_siliques ~year_num
                           +Water*Herbivore*sLAR
                           #+Water*Herbivore*sSucc #model does not run with succulence
                           +Water*Herbivore*sFP
                           +Water*Herbivore*sduration
                           +Water*Herbivore*sheight
                           + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)
Anova(mod_fecundityRM,type="III")

#Warning messages:
#1: In .checkRankX(TMBStruc, control$rank_check) :
#  fixed effects in conditional model are rank deficient
#2: In fitTMB(TMBStruc) :
#  Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
#3: In fitTMB(TMBStruc) :
#  Model convergence problem; false convergence (8). See vignette('troubleshooting')

mod_fecund_selection <-glmmTMB (Mature_length_siliques ~S_elev+year_num
                           +Water*Herbivore*sLAR 
                           #+Water*Herbivore*sSLA
                           #+Water*Herbivore*sLWC 
                           #+Water*Herbivore*sFP
                           +Water*Herbivore*sduration
                           +Water*Herbivore*sheight  
                           +Water*Herbivore*speak 
                           + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)
Anova(mod_fecund_selection,type="III")



mod_fecunditytreat <-glmmTMB (Mature_length_siliques ~S_elev+year_num
                           +treat*sleaf
                           +treat*sLAR
                           #+Water*Herbivore*sSLA
                           #+Water*Herbivore*sLWC 
                           +treat*sFP
                           +treat*sduration
                           +treat*sheight
                           +treat*speak 
                           + (1|Cage_Block)+(1|Genotype)+(1|PlantID), family=Gamma(link="log"), data = traitdatRepro)
Anova(mod_fecunditytreat,type="III")
library(DHARMa)
simoutput<- simulateResiduals(fittedmodel= mod_fecund_selection, plot = T )

simulationOutput <- simulateResiduals( mod_fecund_selection, plot = T) 

visreg(mod_fecund_selection,"sduration", by="Herbivore",cond=list(Water="Supplemental"), overlay=FALSE, scale = "response", xlab="flowering time (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))
visreg(mod_fecund_selection,"sduration", by="Herbivore",cond=list(Water="Restricted"), overlay=FALSE, scale = "response", xlab="flowering time (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))


visreg(mod_fecund_selection,"speak", by="Herbivore",cond=list(Water="Supplemental"), overlay=FALSE, scale = "response", xlab="flowering time (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))
visreg(mod_fecund_selection,"speak", by="Herbivore",cond=list(Water="Restricted"), overlay=FALSE, scale = "response", xlab="flowering time (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))


# flowering duration
duration_figa= ggplot(traitdat, aes(x= flowering_duration,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Duration of Flowering (days)")  

duration_figa + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Herbivore, scales="free_x") +scale_colour_manual(values = c( "darkred","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental"))

#flowering time, herbivore by flowering time
flowering_fig= ggplot(traitdat, aes(x= Ordinal_Date_flowering,y= Mature_length_siliques, group= Herbivore, colour= Herbivore))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Day of flowering (ordinal)")  

flowering_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+scale_colour_manual(values = c( "darkgreen","black"), name = "Herbivore treatment", labels = c("Addition","Removal"))



# peak flowering
# herbivore by flowering time
Hpeakflowering_fig= ggplot(traitdat, aes(x= Date_peak_flowering,y= Mature_length_siliques, group= Herbivore, colour= Herbivore))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Day of peak flowering (ordinal)")  

Hpeakflowering_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+scale_colour_manual(values = c( "darkgreen","black"), name = "Herbivore treatment", labels = c("Addition","Removal"))

# water by flowering time
Wpeakflowering_fig= ggplot(traitdat, aes(x= Date_peak_flowering,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Day of peak flowering (ordinal)")  

Wpeakflowering_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+scale_colour_manual(values = c( "darkred","lightblue"), name = "Watering treatment", labels = c("Restricted","Supplemental"))

Wpeakflowering_fig= ggplot(traitdat, aes(x= Date_peak_flowering,y= Reproduced, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Probability of reproduction")+ scale_x_continuous("Day of peak flowering (ordinal)")  

Wpeakflowering_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",method.args=list(family = "binomial"), size=1.6, formula=y~x)+scale_colour_manual(values = c( "darkred","lightblue"), name = "Watering treatment", labels = c("Restricted","Supplemental"))


#LAR

Wpeakflowering_fig= ggplot(traitdat, aes(x= avg_LAR,y= Mature_length_siliques, group= Water, colour= Water))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Day of peak flowering (ordinal)")  

Wpeakflowering_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+scale_colour_manual(values = c( "darkred","lightblue"), name = "Watering treatment", labels = c("Restricted","Supplemental"))


lar_fig= ggplot(traitdat, aes(x= avg_LAR,y= Mature_length_siliques, group= Herbivore, 
                                        colour= Herbivore))+geom_point(size=2) + scale_y_continuous("Mature Length Siliques")+ scale_x_continuous("Leaf area removed by Herbivores (%)")  
lar_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                       axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                                       panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Water, scales="free_x") +scale_colour_manual(values = c( "darkgreen","black"), name = "Herbivore treatment", labels = c("Addition","Removal"))


lar_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                       axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                                       panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Water, scales="free_x") +scale_colour_manual(values = c( "darkgreen","gray"), name = "Herbivore treatment", labels = c("Addition","Removal"))



 
##### Hurdle model and older code from Jill and Rachel's paper ####

 hurdle_Model = glmmTMB (Mature_length_siliques ~ elev +Treatment*Herbivore*sFP  +Treatment*Herbivore*sduration + Treatment*Herbivore*sheight +Treatment*Herbivore*sleaf  + (1|Cage_Block)+ (1| Genotype), zi=~elev +Treatment*Herbivore*sFP + Treatment*Herbivore*sduration +Treatment*Herbivore*sheight + Treatment*Herbivore*sleaf  + (1|Cage_Block)+ (1| Genotype),data = traitdat ,family=ziGamma(link="log"))
 
 hurdle_Model_ft = glmmTMB (Mature_length_siliques ~ elev +Treatment*Herbivore*sFP + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID), zi=~elev +Treatment*Herbivore*sFP + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID),data = traitdat ,family=ziGamma(link="log"))
 
 visreg(hurdle_Model,"sFP", scale = "response")
 
 diagnose(hurdle_Model)
 check_zeroinflation(hurdle_Model)
check_overdispersion(hurdle_Model)
Anova(hurdle_Model,type="III", component="cond")
summary(hurdle_Model)
car::Anova(hurdle_Model,component="cond")
car::Anova(hurdle_Model,component="zi")


hurdle_Model_fd = glmmTMB (Mature_length_siliques ~ elev +Treatment*Herbivore*sduration + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID), zi=~elev +Treatment*Herbivore*sduration + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID),data = traitdat ,family=ziGamma(link="log"))

p_repro_fd <-glmmTMB(Reproduced~ Treatment*Herbivore*sduration + year_num+(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

visreg(hurdle_Model_fd,"sduration", by="Treatment", overlay=FALSE, scale = "response")

Anova(p_repro_fd,type="III", component="cond")
summary(hurdle_Model)
car::Anova(hurdle_Model_fd,component="cond")
car::Anova(hurdle_Model_fd,component="zi")


hurdle_Model_fpeak = glmmTMB (Mature_length_siliques ~ elev +Treatment*Herbivore*Date_peak_flowering + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID), zi=~elev +Treatment*Herbivore*Date_peak_flowering + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID),data = traitdat ,family=ziGamma(link="log"))

p_repro_fd <-glmmTMB(Reproduced~ Treatment*Herbivore*sduration + year_num+(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))

visreg(hurdle_Model_fpeak,"Date_peak_flowering", by="Treatment", overlay=FALSE, scale = "response")

Anova(hurdle_Model_fpeak,type="III", component="cond")
summary(hurdle_Model)
car::Anova(hurdle_Model_fpeak,component="zi")



hurdle_Model_SLA = glmmTMB (Mature_length_siliques ~ elev +Treatment*Herbivore*rosette_lwc + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID), zi=~elev +Treatment*Herbivore*rosette_lwc + year_num+(1|Cage_Block)+ (1| Genotype)+ (1| PlantID),data = traitdat ,family=ziGamma(link="log"))



Anova(hurdle_Model_SLA,type="III", component="zi")
summary(hurdle_Model)
car::Anova(hurdle_Model_SLA,component="zi")




rosette_SLA




hurdle_model = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf + Drought*Nutrient*I(sleaf^2)  + (1|Block)+ (1| Population), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1|Block)+ (1| Population),data = traitdat ,family=truncated_nbinom2)
car::Anova(hurdle_model,component="cond")
car::Anova(hurdle_model,component="zi")

HM = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf + I(sleaf^2)*Nutrient + (1|Block)+ (1| Population), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1|Block)+ (1| Population),data = traitdat ,family=truncated_nbinom2)
car::Anova(HM,component="cond")
car::Anova(HM,component="zi")




##Without block and population in zi
noblock_zi = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf  + (1|Block)+ (1| Population), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1| Population),data = traitdat ,family=truncated_nbinom2)
nopop_zi = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf  + (1|Block)+ (1| Population), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1|Block),data = traitdat ,family=truncated_nbinom2)

anova(hurdle_Model,noblock_zi)
anova(hurdle_Model,nopop_zi)

##Without block and population in non-zero-inflated
noblock = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf + (1| Population), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1|Block)+ (1| Population),data = traitdat ,family=truncated_nbinom2)
nopop = glmmTMB (SeedCount_total ~ elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration + Drought*Nutrient*sheight+ Drought*Nutrient*sleaf  + (1|Block), zi=~elev +Drought*Nutrient*sFP + Drought*Nutrient*sduration +Drought*Nutrient*sheight + Drought*Nutrient*sleaf  + (1|Block)+ (1| Population),data = traitdat ,family=truncated_nbinom2)

anova(hurdle_Model,noblock)
anova(hurdle_Model,nopop)

##Visualizing needs to happen outside of glmmTMB
plotrep <- glmer(Reproduction~elev +treatB*sFP + treatB*sduration + treatB*sheight+ treatB*sleaf + (1| Population)+(1|Block),data = traitdat, family=binomial(link=logit), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

visreg(plotrep,"sFP", by="treatB", overlay=FALSE, scale = "response", xlab="Timing of first flowering (standardized)", ylab="Probability of reproduction", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))
visreg(plotrep,"sduration", by="treatB", overlay=FALSE, scale = "response", xlab="Duration of flowering  (standardized)", ylab="Probability of reproduction", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"), ylim=c(-0.1,1.1))

fecund<-subset(traitdat, Reproduction=="1")

fecund $sFP<-scale(fecund $Flowering_time,center=TRUE,scale=TRUE)
fecund $sduration<-scale(fecund $Flowering_duration,center=TRUE,scale=TRUE)
fecund $sheight<-scale(fecund $Height,center=TRUE,scale=TRUE)
fecund $elev<-scale(fecund $Elevation,center=TRUE,scale=TRUE)
fecund $sleaf<-scale(fecund $Leaf_number_first_flowering,center=TRUE,scale=TRUE)

Hard_select <- glmer.nb(SeedCount~elev + treatB*sFP +  treatB*sduration + treatB*sheight+ treatB*sleaf + (1| Population)+(1|Block),data = fecund, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Hard_select,type="III")

visreg(Hard_select,"sheight", by="treatB", overlay=FALSE, scale = "response", xlab="Height at flowering (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))
visreg(Hard_select,"sleaf", by="treatB", overlay=FALSE, scale = "response", xlab="Leaf number at flowering  (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial= TRUE,points=list(cex=1.2, col="black"))


Hard_selectB <- glmer.nb(SeedCount~elev + treatB*sFP +  treatB*sduration + treatB*sheight+ treatB*sleaf+ treatB*I(sleaf^2) + (1| Population)+(1|Block),data = fecund, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(Hard_selectB,type="III")

visreg(Hard_selectB,"sheight", by="treatB", overlay=FALSE, scale = "response", xlab="Height at flowering (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))
visreg(Hard_selectB,"sleaf", by="treatB", overlay=FALSE, scale = "response", xlab="Leaf number at flowering  (standardized)", ylab="Seed count, among individuals that set seeds", line=list(col="black"),partial= TRUE,points=list(cex=1.2, col="black"))


