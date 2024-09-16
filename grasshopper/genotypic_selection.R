######## PROJECT: Common garden experiment: Examining clines and plasticity in response to water availability and grasshopper abundance
#### PURPOSE: calculate genotypic means for genotypic selection analyses
#### DATE LAST MODIFIED: 5 Sept 2024



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
LSmeans <- read.csv("LS_means.csv",stringsAsFactors=T)

##Change the baseline for Water availability
LSmeans $Water <-factor(LSmeans $Water, levels = c("Restricted", "Supplemental"))

##Change the baseline for Herbivore manipulation
LSmeans $Herbivore <-factor(LSmeans $Herbivore, levels = c("Addition", "Removal"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
LSmeans $S_elev<-scale(LSmeans $elevation,center=TRUE, scale=TRUE)


#### probability of reproduction

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)

LSmeans $sLAR<-scale(LSmeans $LAR,center=TRUE,scale=TRUE)
LSmeans $sSLA<-scale(LSmeans $SLA,center=TRUE,scale=TRUE)
LSmeans $sSUC<-scale(LSmeans $succulence,center=TRUE,scale=TRUE)

repro_model <-glmmTMB(prob_repro~Water*Herbivore+
                        Water*Herbivore*sLAR+ 
                        Water*Herbivore*sSUC+ #I(sSUC^2)+
                        Water*Herbivore*sSLA+ #I(sSLA^2) 
                      +(1|Genotype),data=LSmeans,family=gaussian(link="identity"))

hist(LSmeans$prob_repro)

simulationOutput <- simulateResiduals(fittedModel= repro_model, plot = T, re.form = NULL)

Anova(repro_model,type="III") 

#### fecundity

##Create datafile with only reproductive plants
LSmeansRepro <- filter(LSmeans, mature_length_siliques != 0 )

##Rescale after removing the non-reproductive plants
LSmeansRepro $sduration<-scale(LSmeansRepro $duration,center=TRUE,scale=TRUE)
LSmeansRepro $s_max_height<-scale(LSmeansRepro $max_height,center=TRUE,scale=TRUE)
LSmeansRepro $sLAR<-scale(LSmeansRepro $LAR,center=TRUE,scale=TRUE)
LSmeansRepro $sSLA<-scale(LSmeansRepro $SLA,center=TRUE,scale=TRUE)
LSmeansRepro $sSUC<-scale(LSmeansRepro $succulence,center=TRUE,scale=TRUE)
LSmeansRepro $sFT<-scale(LSmeansRepro $phenology,center=TRUE,scale=TRUE)
#LSmeansRepro $S_initdiam <-scale(LSmeansRepro $init.diam,center=TRUE,scale=TRUE)



fecund_modeltraits <-glmmTMB(mature_length_siliques~ Water*Herbivore +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +#Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ #Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+#Water*Herbivore*I(sLAR^2)
                             +(1|Genotype),data=LSmeansRepro)#,family=gaussian(link="identity"))

hist(LSmeansRepro$mature_length_siliques)


Anova(fecund_modeltraits,type="III") 

##To check residuals
simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)


### Univariate models####


###Flowering time###


fecund_modeltraits_FT <-glmmTMB(mature_length_siliques~ Water*Herbivore +
                                  Water*Herbivore*sFT
                               +I(sFT^2)+
                                +(1|Genotype),data=LSmeansRepro,family=gaussian(link="identity"))

Anova(fecund_modeltraits_FT,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits_FT, plot = T, re.form = NULL,allow.new.levels =TRUE)



##### calculating LS means #######


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


#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************

##* Censuses 1 and 2 occurred for all years. Census 3 was for 2021 and 2022 only

##reformat datafile

LAR_data_long_form<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2",
#                                                         "LAR_3"
                                                         ),
                                                  names_to='census',
                                                  values_to='LAR')

LAR_data_long_form <- dplyr::select(LAR_data_long_form, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,census, year)

LAR_data_long_form$census[LAR_data_long_form$census == "LAR_1"] <- "1"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_2"] <- "2"
#LAR_data_long_form$census[LAR_data_long_form$census == "LAR_3"] <-"3"

LAR_data_long_form $census <-as.factor(LAR_data_long_form $census)

LAR_data_long_form $year <-as.factor(LAR_data_long_form $year)
##Let's concatenate census and year
LAR_data_long_form $census_year<-interaction(LAR_data_long_form$census, LAR_data_long_form$year,sep = "_")


LAR_data_long_form$LAR_prop<-LAR_data_long_form $LAR/100
hist(LAR_data_long_form$LAR_prop)
LAR_data_long_form <- drop_na(LAR_data_long_form,LAR_prop) 

n<-nrow(LAR_data_long_form)

##this is the beta transformation, which transforms all values of 0 to a small value.

LAR_data_long_form $y_beta<- (LAR_data_long_form $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data_long_form $y_beta)

min(LAR_data_long_form $y_beta)

max(LAR_data_long_form $y_beta)


LAR_Model_three<- gamlss (formula= y_beta ~Genotype*Water*Herbivore+year+ random(census_year)+ random(Cage_Block),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

LAR_model <- glmmTMB(y_beta ~ Water*Herbivore*Genotype+year
                      #+(1|census_year)
                      +(1|Cage_Block), data = LAR_data_long_form, 
                      family=beta_family())


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

write.csv(lsmeans_LAR,file="LSmeans_LAR.csv")


#*******************************************************************************
#### SLA #######
#*******************************************************************************
#*

sla_model_2 <- glmmTMB(SLA ~ Water*Herbivore*Genotype*year+(1|Cage_Block), data = foliar, family= lognormal(link="log"),control=glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3)))

lsmeans_SLA <-emmeans (sla_model_2, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_SLA,file="LSmeans_SLA.csv")


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

#*******************************************************************************
#### Flowering phenology #######
#*******************************************************************************

FT_model<-glmmTMB(FT_Adj ~ Genotype*Water*Herbivore*year+(1|Cage_Block),data= flowering)#,family=lognormal(link="log"))

lsmeans_FT <-emmeans (FT_model, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_FT,file="LSmeans_phenology.csv")

#*******************************************************************************
#### Flowering duration #######
#*******************************************************************************

flowering_duration<-glmmTMB(flowering_duration ~ Genotype*Water*Herbivore*year+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

lsmeans_duration <-emmeans (flowering_duration, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_duration,file="LSmeans_duration.csv")

#*******************************************************************************
#### Height at flowering #######
#*******************************************************************************

max_height<-glmmTMB(Max_height_flowering ~ Genotype+Water*Herbivore*year+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

lsmeans_height <-emmeans (max_height, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_height,file="LSmeans_height.csv")


#*******************************************************************************
#### probability of reproduction #######
#*******************************************************************************

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")

prob_repro <-glmmTMB(Reproduced~Genotype+Water*Herbivore*year+
                      (1|Cage_Block),data=grasshopper_no2021,family=binomial(link="logit"))

lsmeans_prob_repro <-emmeans (prob_repro, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_prob_repro,file="LSmeans_PR.csv")


#*******************************************************************************
#### fecundity #######
#*******************************************************************************

##Create datafile with only reproductive plants
traitdatRepro <- filter(grasshopper_no2021, Reproduced == 1 )

fecundity <-glmmTMB(Mature_length_siliques~Genotype*Water*Herbivore*year+
                       (1|Cage_Block),data=traitdatRepro,family=Gamma(link="log"))

lsmeans_fecundity <-emmeans (fecundity, ~Genotype*Water*Herbivore*year, type="response",adjust = "sidak")

write.csv(lsmeans_fecundity,file="LSmeans_fecundity.csv")



