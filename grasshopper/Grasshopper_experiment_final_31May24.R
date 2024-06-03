######## PROJECT: Common garden experiment: Examining clines and plasticity in response to water availability and grasshopper abundance
#### PURPOSE:Examine clines, trait values and fitness in response towater availability and grasshopper abundance in the field.
#### AUTHOR: Jill Anderson and Inam Jameel
#### DATE LAST MODIFIED: 31 May 2024


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

#setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/Grasshopper_manuscript_files/Grasshopper_data/")
#setwd("~/Documents/personnel/Jameel/grasshopper")

#this is where you specify the folder where you have the data on your computer


#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_updated_1May24.csv",stringsAsFactors=T)

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

#Let's concatenate herbivore and watering , which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Water,sep = "_")

grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))

##Let's correlate rosette and stem leaf data to see if we can come up with composite figures
plot(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
plot(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
plot(grasshopper$rosette_lwc~grasshopper$bolt_lwc)

mod1<-lm(grasshopper$rosette_succulence~grasshopper$bolt_succulence)
summary(mod1)
Anova(mod1,type="III")

LAR_Mod<-lm(grasshopper$rosette_SLA~grasshopper$bolt_SLA)
summary(LAR_Mod)

mod3<-lm(grasshopper$rosette_lwc~grasshopper$bolt_lwc)
summary(mod3)

##Create composite leaf SLA, succuclence and lwc variables based on the regressions above. This gives us foliar trait data for the 67 plants for which we have stem but not rosette collections
grasshopper$SLA <- ifelse(is.na(grasshopper$rosette_SLA), (71.8177 + 0.7319*grasshopper$bolt_SLA), grasshopper$rosette_SLA)
grasshopper$succulence <- ifelse(is.na(grasshopper$rosette_succulence), (0.0023948 + 0.5609208*grasshopper$bolt_succulence), grasshopper$rosette_succulence)
grasshopper$lwc <- ifelse(is.na(grasshopper$rosette_lwc), (0.19571 + 0.67021*grasshopper$bolt_lwc), grasshopper$rosette_lwc)

#This calculates flowering duration
grasshopper$flowering_duration<-(grasshopper $Date_silique - grasshopper $FT_Adj)

#for vegetative traits
foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

#set colors
cols=c("#CC79A7","lightblue")
cols2=c("#882255","#77DDAA")

#*******************************************************************************
#### Treatment Water effect ####
#*******************************************************************************

#read in data 
VWC <- read.csv("Grasshopper_VWC.csv",stringsAsFactors=T)
VWC $Year <-as.factor(VWC $Year)

#model with multiple censuses across years

#reformat datafile

VWC_data<- VWC %>% pivot_longer(cols=c("VWC_1","VWC_2","VWC_3","VWC_4","VWC_5","VWC_6","VWC_7","VWC_8"),
                                names_to='census',
                                values_to='VWC')

VWC_data$census[VWC_data$census == "VWC_1"] <- "1"
VWC_data$census[VWC_data$census == "VWC_2"] <- "2"
VWC_data$census[VWC_data$census == "VWC_3"] <-"3"
VWC_data$census[VWC_data$census == "VWC_4"] <- "4"
VWC_data$census[VWC_data$census == "VWC_5"] <- "5"
VWC_data$census[VWC_data$census == "VWC_6"] <- "6"
VWC_data$census[VWC_data$census == "VWC_7"] <- "7"
VWC_data$census[VWC_data$census == "VWC_8"] <- "8"

VWC_data $census <-as.factor(VWC_data $census)

VWC_data $Year <-as.factor(VWC_data $Year)

VWC_data$yearWater <- interaction(VWC_data$Year, VWC_data$Water)
VWC_data$yearHerbivore <- interaction(VWC_data$Year, VWC_data$Herbivore)



VWC_mod <- glmmTMB(VWC ~Water*Year*Herbivore+(1|census) + (1|Cage_Block), data= VWC_data,family=Gamma(link="log"))

simulationOutput <- simulateResiduals(fittedModel= VWC_mod, plot = T, re.form = NULL,allow.new.levels =TRUE)

Anova(VWC_mod,type="III")

VWC_means<-emmeans(VWC_mod, ~ Water:Year, type="response", adjust = "sidak")
cld(VWC_means, details=TRUE)


VWC_mod_no_census <- glmmTMB(VWC ~Water*Year*Herbivore+(1|census), data= VWC_data,family=Gamma(link="log"))

VWC_mod_no_cageblock <- glmmTMB(VWC ~Water*Year*Herbivore + (1|Cage_Block), data= VWC_data,family=Gamma(link="log"))


anova(VWC_mod,VWC_mod_no_cageblock)
anova(VWC_mod,VWC_mod_no_census)

##Box plot
VWC_box<-ggplot(VWC_data, aes(x = Water, y = VWC, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ 
  scale_y_continuous("Volumetric water content") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

VWC_plot<-VWC_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none") +  scale_fill_manual(values = c(cols), name = "Water treatment")+facet_wrap(~Year)


VWC_plot # Fig S3



### Hypothesis one: Context dependency of clines and plasticity #####

#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************



#* Censuses 1 and 2 occurred for all years. Census 3 was for 2021 and 2022 only

#reformat datafile

LAR_data_long_form<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3"),
                                                names_to='census',
                                                values_to='LAR')

LAR_data_long_form <- dplyr::select(LAR_data_long_form, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev, treat,census, year)

LAR_data_long_form$census[LAR_data_long_form$census == "LAR_1"] <- "1"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_2"] <- "2"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_3"] <-"3"

LAR_data_long_form $census <-as.factor(LAR_data_long_form $census)

LAR_data_long_form $year <-as.factor(LAR_data_long_form $year)
#Let's concatenate census and year
LAR_data_long_form $census_year<-interaction(LAR_data_long_form$census, LAR_data_long_form$year,sep = "_")


LAR_data_long_form$LAR_prop<-LAR_data_long_form $LAR/100
hist(LAR_data_long_form$LAR_prop)
LAR_data_long_form <- drop_na(LAR_data_long_form,LAR_prop) 

n<-nrow(LAR_data_long_form)

#this is the beta transformation, which transforms all values of 0 to a small value.

LAR_data_long_form $y_beta<- (LAR_data_long_form $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data_long_form $y_beta)

min(LAR_data_long_form $y_beta)

max(LAR_data_long_form $y_beta)

##Conversation of standardized elevation back to raw units
-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)


## model
LAR_Model<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)
save(LAR_Model, file='LAR_Model.rda')   


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data_long_form %>%  
  mutate(fit.m = predict(LAR_Model, se.fit=FALSE),              
         resid = residuals(LAR_Model))

##Convert coefficients to probabilities
newdf2 $predicted<-1/(1+exp(-(newdf2 $fit.m))) 
newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


	##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
(2600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
(3100-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
(3600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)


damage_cline<-visregList(visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2021"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2021"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2022"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2022"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2023"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2023"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)
head(damage_cline$fit)		
damfit<-data.frame(damage_cline$fit)
		head(damfit)
damfit$visregFITexp<-exp(damfit$visregFit)
		
damage_clinal_variation<-ggplot(damfit, aes(S_elev, visregFITexp,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= LAR_data_long_form, aes(S_elev, LAR_prop, color= Water, shape=Water), alpha=0.50, color="black")+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.557947,-0.1085093, 1.340929))+ scale_y_continuous("Leaf area removed by herbivores (proportion)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+ geom_line(aes(group=Water),size=1)+scale_linetype_manual(values=c("solid", "dotted"))+facet_grid(Herbivore~year)
 damage_clinal_variation



#Box_plot. Fig. 2B
LAR_box <-ggplot(LAR_data_long_form, aes(x = Herbivore, y = LAR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = .75,position = position_jitterdodge(0.3))

LAR_box <-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
   panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~year)
LAR_box


#Box_plot. Fig. 2B
LAR_box <-ggplot(LAR_data_long_form, aes(x = Water, y = LAR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

LAR_box <-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(Herbivore~year)
LAR_box

Figure2 <- ggarrange(
  damage_clinal_variation,
  LAR_box,
   labels = c("A", "B"),
  ncol = 2, nrow = 1, common.legend = TRUE, legend="right")
Figure2

Fig2 <- ggarrange(
  damage_clinal_variation,
  LAR_box,
   labels = c("A", "B"),
  ncol = 1, nrow = 2, common.legend = TRUE, legend="right")
Fig2



##Conversation of standardized elevation back to raw units
-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)



#env concatenates water and herbivore. This model allows us to extract means and slopes for each combination of year, treatment and garden
LAR_data_long_form $env<-interaction(LAR_data_long_form $Water, LAR_data_long_form $Herbivore)
LAR_data_long_form $env_year<-interaction(LAR_data_long_form $env, LAR_data_long_form $year)
LAR_ModelE<- gamlss (formula= y_beta ~S_elev* env_year-1+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
summary(LAR_ModelE)
save(LAR_ModelE, file='LAR_ModelE.rda')   


#Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+S_elev*Water*year+S_elev* Herbivore*year+ Water* Herbivore*year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)
save(LAR_Model_three, file='LAR_Model_three.rda')   

#Subsets of models for drop 1
LAR_Model_two<- gamlss (formula= y_beta ~S_elev*Water+S_elev*Herbivore+ Water* Herbivore+year*Water+year*Herbivore+S_elev*year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)
save(LAR_Model_two, file='LAR_Model_two.rda')   

LAR_Model_one<- gamlss (formula= y_beta ~S_elev+Water+Herbivore+ year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)
save(LAR_Model_one, file='LAR_Model_one.rda')




#*******************************************************************************
#### Specific Leaf Area #####
#*******************************************************************************

foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)


sla_model <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


SLA <-emmeans(sla_model, ~ Water:Herbivore, type="response", adjust = "sidak")
cld(SLA, details=TRUE)

##Slope of the relationship between source elevation and SLA. The beta values need to be exponentiated
emtrends(sla_model, specs = c("year"), var = "S_elev")

#To test the random effeccts
sla_model_no_plantID <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))

sla_model_no_genotype <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))

sla_model_no_cageblock <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype), data = foliar, family= lognormal(link="log"))


anova(sla_model,sla_model_no_plantID)
anova(sla_model,sla_model_no_genotype)
anova(sla_model,sla_model_no_cageblock)

	##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
(2600-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)
(3100-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)
(3600-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)


SLA_clinal<-visreg(sla_model,"S_elev", by="year", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
SLA_clinal_variation<-ggplot(SLA_clinal $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= foliar, aes(S_elev, SLA), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("SLA (cm2/g)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))+facet_wrap(~year)
SLA_clinal_variation


SLA_clinal<-visreg(sla_model,"S_elev", by="year", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
SLA_clinal_variation<-ggplot(SLA_clinal $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= foliar, aes(S_elev, SLA), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("SLA (cm2/g)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))+facet_wrap(~year)
SLA_clinal_variation



SLA_box <-ggplot(foliar, aes(x = Water, y = SLA, fill = Water,shape=Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Specific leaf area") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

SLA_treatment<-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "Top")+ scale_x_discrete()+ scale_fill_manual(values = c(cols), name = "Water availability", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_wrap(~Herbivore)



#Box_plot
SLA_box <-ggplot(foliar, aes(x = Water, y = SLA, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

SLA_treatment <-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_grid(Herbivore~year)
SLA_treatment



#*******************************************************************************
#### Succulence #####
#*******************************************************************************

foliar$succulence_mg<-foliar$succulence*1000
foliar$Suc_mg<- foliar$succulence_mg+0.01

succulence_data <- dplyr::select(foliar, succulence, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, Cage_Block, elev_km, S_elev,  year, Suc_mg, elev_dist_km)

succulence_data <- drop_na(succulence_data,succulence) 

succ_model <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))
Anova(succ_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= succ_model, plot = T, re.form = NULL,allow.new.levels =T)

succ_model_no_plantID <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))

succ_model_no_genotype <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))

succ_model_no_cageblock <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype), data = succulence_data, family=lognormal(link="log"))

anova(succ_model,succ_model_no_plantID)
anova(succ_model,succ_model_no_genotype)
anova(succ_model,succ_model_no_cageblock)


succulence <-emmeans(succ_model, ~ Herbivore:Water, type="response", adjust = "sidak")
cld(succulence, details=TRUE)



Succulence_cline <-plot(pred_succ, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Leaf succulence (mg)")



succ_cline<-visregList(visreg(succ_model,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(succ_model,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
succ_clinal_variation<-ggplot(succ_cline $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= succulence_data, aes(S_elev, Suc_mg, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))+ scale_y_continuous("Leaf succulence (mg)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
succ_clinal_variation


##Box plot

Suc_box <-ggplot(succulence_data, aes(x = Water, y = Suc_mg, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Leaf succulence (mg)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

Suc_treatment<-Suc_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water availability", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)



Suc_box_herb <-ggplot(succulence_data, aes(x = Herbivore, y = Suc_mg, fill = Herbivore)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  scale_y_continuous("Leaf succulence (mg)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

FigS4_Suc_treatment_herb <-Suc_box_herb + theme_classic() + theme(text = element_text(size=14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Grasshopper addition", "Grasshopper removal")) +  scale_fill_manual(values = c(cols2), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~year)

FigS4_Suc_treatment_herb #supplement figure


#*******************************************************************************
#### Flowering phenology #####
#*******************************************************************************

#filter out the two plants that flowered in the first season. This step is likely not necessary because those plants didn't have flowering times.
grasshopperFT <- filter(grasshopper, year != "2021")

##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(grasshopperFT, Ordinal_Date_flowering!="NA")


FT_model<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(FT_model, type = "III") # 

 emtrends(FT_model, specs = c("Water","Herbivore"), var = "S_elev")





FT_model_no_geno<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
FT_model_no_cageblock<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1| PlantID),data= flowering,family=lognormal(link="log"))
FT_model_noplantID<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

anova(FT_model,FT_model_no_geno)
anova(FT_model,FT_model_no_cageblock)
anova(FT_model,FT_model_noplantID)

	##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
(2600-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)
(3100-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)
(3600-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)



flowering_time_cline<-visregList(visreg(FT_model,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(FT_model,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
phen_cline <-ggplot(flowering_time_cline $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowering, aes(S_elev, FT_Adj, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.630394,-0.1998366,  1.23072))+ scale_y_continuous("Flowering phenology (ordinal day of year)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
phen_cline


FT_treatment <-FT_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)



#Box_plot
FT_box <-ggplot(flowering, aes(x = Water, y = FT_Adj, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous(element_blank()) +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

FT_treatment <-FT_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)
FT_treatment


pred_FT2 <- ggpredict(FT_model, terms = c("S_elev[all]", "Water","year"), type = "re", interval="confidence")
phen_cline_year <-plot(pred_FT2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")

####


phen_cline_year #supplemental figure, S5



#*******************************************************************************
#### Flowering duration #####
#*******************************************************************************


flowering_duration<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)

Anova(flowering_duration, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)


flowering_duration_nogeno<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Cage_Block)+(1| PlantID),data= flowering)

flowering_duration_nocageblock<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+(1|Genotype)+(1| PlantID),data= flowering)

flowering_duration_noplantid<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block),data= flowering)


anova(flowering_duration,flowering_duration_nogeno)
anova(flowering_duration,flowering_duration_nocageblock)
anova(flowering_duration,flowering_duration_noplantid)




floweringduration_cline<-visregList(visreg(flowering_duration,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(flowering_duration,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
Fduration_clinal_variation<-ggplot(floweringduration_cline $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowering, aes(S_elev, flowering_duration, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.630394,-0.1998366,  1.23072))+ scale_y_continuous("Flowering duration (days)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
Fduration_clinal_variation

Duration_cline #no significant interaction


Duration_box <-ggplot(flowering, aes(x = Water, y = flowering_duration, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

Duration_treatment <-Duration_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)

#Box_plot
Duration_box <-ggplot(flowering, aes(x = Water, y = flowering_duration, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous(element_blank()) +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

Duration_treatment <-Duration_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                  panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)
Duration_treatment



#*******************************************************************************
#### Tallest stem at flowering #####
#*******************************************************************************

max_height<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

Anova(max_height, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= max_height, plot = T, re.form = NULL)

max_height_nogeno<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

max_height_nocageblock<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1| PlantID),data= flowering,family=lognormal(link="log"))

max_height_nopid<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

anova(max_height,max_height_nogeno)
anova(max_height,max_height_nocageblock)
anova(max_height,max_height_nopid)

height_means<-emmeans(max_height, ~year, type="response", adjust = "sidak")
height_means
pairs(height_means)



height_cline<-visregList(visreg(max_height,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(max_height,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
		
height_clinal_variation<-ggplot(height_cline $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowering, aes(S_elev, flowering_duration, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.630394,-0.1998366,  1.23072))+ scale_y_continuous("Height of tallest stem at flowering")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
height_clinal_variation


height_box <-ggplot(flowering, aes(x = Water, y = Max_height_flowering, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

height_treatment <-height_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)





#*******************************************************************************
####   Hypothesis Two: Selection via two fitness components on traits   #####
#*******************************************************************************

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced_updated
grasshopper_no2021<-subset(grasshopper, year!="2021")
# retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper_no2021,FT_Adj,Max_height_flowering, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block,  year, Mature_length_siliques,Reproduced,flowering_duration, SLA, lwc, succulence, FT_Adj)



##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
traitdat $S_initdiam<-scale(traitdat $init.diam,center=TRUE, scale=TRUE)

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)

traitdat $sLAR<-scale(traitdat $avg_LAR,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)


head(traitdat)


#Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Addition","Removal"))


#*******************************************************************************
####  Probability of reproduction for vegetative traits only  #####
#*******************************************************************************

repro_model <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+ 
                             Water*Herbivore*sSUC+ I(sSUC^2)+
                             Water*Herbivore*sSLA+ I(sSLA^2) 
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model,type="III") 
summary(repro_model)

simulationOutput <- simulateResiduals(fittedModel= repro_model, plot = T, re.form = NULL)

repro_model_noplantID <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                                  Water*Herbivore*sLAR+
                                  Water*Herbivore*sSUC+I(sSUC^2)+
                                  Water*Herbivore*sSLA+I(sSLA^2) 
                                +(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
                                
repro_model_nocb <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                                  Water*Herbivore*sLAR+
                                  Water*Herbivore*sSUC+I(sSUC^2)+
                                  Water*Herbivore*sSLA+I(sSLA^2) 
                                +(1|PlantID)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
                                

repro_model_nogeno <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                                    Water*Herbivore*sLAR+
                                    Water*Herbivore*sSUC+I(sSUC^2)+
                                    Water*Herbivore*sSLA+I(sSLA^2) 
                                  +(1|PlantID)+(1|Cage_Block),data=traitdat,family=binomial(link="logit"))
                                  
anova(repro_model,repro_model_noplantID)
anova(repro_model,repro_model_nocb)
anova(repro_model,repro_model_nogeno)



#Conversation of standardized traits back to raw units
-1*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
0*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
1*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
2*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
3*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)


sla_repro = visreg(repro_model,"sSLA", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
sla_repro_plot<-ggplot(sla_repro $fit, aes(sSLA, visregFit)) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdat, aes(sSLA, Reproduced),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Specific leaf area (cm2/g)")+ scale_y_continuous("Probability of reproduction")
sla_repro_plot

suc_repro = visreg(repro_model,"sSUC", by="Herbivore",overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
suc_repro_plot<-ggplot(suc_repro$fit, aes(sSUC, visregFit,group= Herbivore, colour= Herbivore, fill=factor(Herbivore))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Herbivore)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdat, aes(sSUC, Reproduced, color= Herbivore, shape=Herbivore),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf succulence (mg/cm2)")+ scale_y_continuous("Probability of reproduction")+scale_fill_manual(values = cols2, name = "Grasshopper treatment", labels = c("Addition","Removal"))+scale_colour_manual(values = cols2, name = "Grasshopper treatment", labels = c("Addition","Removal"))
suc_repro_plot

#Conversation of standardized trait values back to raw units
-1*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
0*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
1*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
2*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
3*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)


##To extract coefficients, we need to create the quadratic effects outside of the model
traitdat$sSLA2<-traitdat$sSLA*traitdat$sSLA
traitdat$sLAR2<-traitdat$sLAR*traitdat$sLAR
traitdat$sSUC2<-traitdat$sSUC*traitdat$sSUC


repro_model_quad <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+ 
                             Water*Herbivore*sSUC+ sSUC2 +
                             Water*Herbivore*sSLA+ sSLA2 
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model_quad,type="III") 
##Slope of selection surfaces 
emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSLA")
emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSLA2")
emtrends(repro_model_quad, specs = c("sSLA"), var = "sSLA")

emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSUC")
emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSUC2")

emtrends(repro_model_quad, specs = c("Herbivore"), var = "sSUC")
emtrends(repro_model_quad, specs = c("Herbivore"), var = "sSUC2")


####  Univariate models  #####

repro_model_SLA <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sSLA+ I(sSLA^2) 
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model_SLA,type="III") 

sla_pred_uni <- ggpredict(repro_model_SLA, terms = c("sSLA[all]"), type = "re", interval="confidence")
SLA_reproduction_uni <- plot(sla_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=FALSE,  facet=FALSE,dot_alpha=0.75)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")
SLA_reproduction_uni

repro_model_SUC <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sSUC+I(sSUC^2)+ 
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model_SUC,type="III") 
sSUC_pred_uni <- ggpredict(repro_model_SUC, terms = c("sSUC[all]","Herbivore"), type = "re", interval="confidence")

Succulence_reproduction_herbivore_uni <-plot(sSUC_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols2,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")
Succulence_reproduction_herbivore_uni


repro_model_LAR <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+ I(sLAR ^2)+
                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))                            
Anova(repro_model_LAR,type="III") 

repro_model_LAR_optim <- update(repro_model_LAR,
                control=glmmTMBControl(optimizer=optim,
                                       optArgs=list(method="BFGS")))
Anova(repro_model_LAR_optim,type="III") 
emtrends(repro_model_LAR_optim, specs = c("Herbivore","Water"), var = "I(sLAR^2)")


lar_univariate <- ggpredict(repro_model_LAR, terms = c("sLAR[all]"), type = "re", interval="confidence")
LAR_uni <- plot(lar_univariate, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE,dot_alpha=0.75, jitter=0.025)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of reproduction")
LAR_uni


#*******************************************************************************
####  Fecundity for all traits #####
#*******************************************************************************

traitdatRepro <- filter(traitdat, Reproduced == 1 ) 
##Rescale after removing the non-reproductive plants
traitdatRepro $sduration<-scale(traitdatRepro $flowering_duration,center=TRUE,scale=TRUE)
traitdatRepro $s_max_height<-scale(traitdatRepro $Max_height_flowering,center=TRUE,scale=TRUE)
traitdatRepro $sLAR<-scale(traitdatRepro $avg_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sSLA<-scale(traitdatRepro $SLA,center=TRUE,scale=TRUE)
traitdatRepro $sSUC<-scale(traitdatRepro $succulence,center=TRUE,scale=TRUE)
traitdatRepro $sFT<-scale(traitdatRepro $FT_Adj,center=TRUE,scale=TRUE)
traitdatRepro $S_initdiam <-scale(traitdatRepro $init.diam,center=TRUE,scale=TRUE)

fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ + Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+Water*Herbivore*I(sLAR^2)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits,type="III") 
simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)


fecund_modeltraits_nocb  <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ + Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+Water*Herbivore*I(sLAR^2)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

fecund_modeltraits_nogeno <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year + Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ + Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+Water*Herbivore*I(sLAR^2)+(1|Cage_Block),data=traitdatRepro,family=Gamma(link="log"))


anova(fecund_modeltraits,fecund_modeltraits_nocb)
anova(fecund_modeltraits,fecund_modeltraits_nogeno)

##To extract coefficients, we need to create the quadratic effects outside of the model
traitdatRepro $sSLA2<-traitdatRepro $sSLA* traitdatRepro $sSLA
traitdatRepro $sLAR2<-traitdatRepro $sLAR* traitdatRepro $sLAR
traitdatRepro $sSUC2<-traitdatRepro $sSUC* traitdatRepro $sSUC
traitdatRepro $sduration2<-traitdatRepro $sduration* traitdatRepro $sduration

fecund_modeltraits_quad <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*sduration2 +Herbivore*sduration2 +
                               Water*Herbivore* sSLA+ + Water*Herbivore*sSLA2 +
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+Water*Herbivore* sLAR2 +(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_quad,type="III")

emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sFT")
emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sduration")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sduration")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sduration2")
emtrends(fecund_modeltraits_quad, specs = c("s_max_height"), var = "s_max_height")

emtrends(fecund_modeltraits_quad, specs = c("Herbivore"), var = "sSUC")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sSLA")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sSLA2")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sLAR")
emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sLAR2")



#SLA
#Context dependent selection
SLA<-visregList(visreg(fecund_modeltraits,"sSLA", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(fecund_modeltraits,"sSLA", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
SLA_fecund_plot<-ggplot(SLA $fit, aes(sSLA, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sSLA, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("SLA (cm2/g)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
SLA_fecund_plot

-2*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
0*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
2*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
4*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)



#succulence
#Context dependent selection
SUC<-visregList(visreg(fecund_modeltraits,"sSUC", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(fecund_modeltraits,"sSUC", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


SUC_fecund_plot<-ggplot(SUC $fit, aes(sSUC, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sSUC, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf succulence (mg/cm2)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
SUC_fecund_plot

-2*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
-1*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
0*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
1*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
2*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
3*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)



#Context dependent selection

LAR<-visregList(visreg(fecund_modeltraits,"sLAR", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(fecund_modeltraits,"sLAR", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


LAR_fecund_plot<-ggplot(LAR $fit, aes(sLAR, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sLAR, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf damage from herbivores (proportion)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
LAR_fecund_plot

-1*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
0*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
1*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
2*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
3*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)


#Flowering phenology
#Context dependent selection
FT<-visregList(visreg(fecund_modeltraits,"sFT", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(fecund_modeltraits,"sFT", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)

FT_fecund_plot<-ggplot(FT $fit, aes(sFT, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sFT, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Flowering time (ordinal day of year)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
FT_fecund_plot



#Conversation of standardized traits back to raw units
-2*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
0*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
2*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
4*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)



#Duration
duration<-visregList(visreg(fecund_modeltraits,"sduration", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
		visreg(fecund_modeltraits,"sduration", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)

duration_fecund_plot<-ggplot(duration $fit, aes(sduration, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sduration, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Flowering duration (days)")+ scale_y_continuous("Fecundity (total fruit length, mm)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
duration_fecund_plot

#Conversation of standardized traits back to raw units
-1*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
0*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
1*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
2*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)



#height
height_fecund = visreg(fecund_modeltraits,"s_max_height", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
height_fecund_plot<-ggplot(height_fecund $fit, aes(s_max_height, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdatRepro, aes(s_max_height, Mature_length_siliques),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Tallest stem at flowering (cm)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))
height_fecund_plot

0*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)
2*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)
4*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)



### Univariate models####


###Flowering time####


fecund_modeltraits_FT <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+ Water*Herbivore*I(sFT^2)+(1|PlantID)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_FT,type="III") 
visreg(fecund_modeltraits_FT,"sFT", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Timing of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
visreg(fecund_modeltraits_FT,"sFT", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Timing of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")



#### Height ####


fecund_modeltraits_height <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*s_max_height +I(s_max_height^2)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_height,type="III") 
visreg(fecund_modeltraits_height,"s_max_height",  scale = "response", xlab="Height at flowering (m)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")



####Flowering duration####


fecund_modeltraits_duration <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore* sduration +Water*I(sduration^2)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_duration,type="III") 
visreg(fecund_modeltraits_duration,"sduration", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Duration of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
visreg(fecund_modeltraits_duration,"sduration", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Duration of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")


#### SLA ####


fecund_modeltraits_SLA <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                                Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_SLA,type="III") 
visreg(fecund_modeltraits_SLA,"sSLA", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Specific leaf area", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional",ylim=c(0,3000))
visreg(fecund_modeltraits_SLA,"sSLA", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Specific leaf area", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional",ylim=c(0,3000))


#### Succulence ####


fecund_modeltraits_SUC <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year + Water*Herbivore* sSUC+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_SUC,type="III") 
visreg(fecund_modeltraits_SUC,"sSUC", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
visreg(fecund_modeltraits_SUC,"sSUC", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")

visreg(fecund_modeltraits_SUC,"sSUC", by="Herbivore",overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional", ylim=c(0,3000))


#### LAR ####


fecund_modeltraits_LAR <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                                Water*Herbivore* sLAR+I(sLAR^2)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits_LAR,type="III") 
visreg(fecund_modeltraits_LAR,"sLAR", overlay=TRUE, scale = "response", xlab="Leaf damage from arthropod herbivores", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional", ylim=c(0,3000))




#*******************************************************************************
#### Hypothesis Three: Effect of water availability on fitness components and local adaptation #####
#*******************************************************************************

grasshopperFT <- filter(grasshopper, year != "2021")

#### Local adaptation analysis using repeated measures ####

hurdle_Model_LA <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))


emtrends(hurdle_Model_LA, specs = c("Water"), var = "S_elev",component="cond")

#To check residuals
simulationOutput <- simulateResiduals(fittedModel= hurdle_Model_LA, plot = T, re.form = NULL,allow.new.levels =TRUE)

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model_LA,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_LA,type="III", component="cond")


RM_fit <- subset(grasshopperFT, Mature_length_siliques>0) # for plotting from count part, excludes plants that did not reproduce


##Extract the coefficients using emtrends, which requires defining the quadratic effect prior to running the model
RM_fit $Selev2<-RM_fit $S_elev* RM_fit $S_elev


hurdle_Model_LA_quad <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                                +Water*Herbivore*S_elev
                                +Water*Herbivore*Selev2
                                +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= RM_fit, 
                                zi=~S_initdiam+year
                                +Water*Herbivore*S_elev
                                +Water*Herbivore*Selev2
                                +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))


visreg(hurdle_Model_LA_quad,"S_elev", by="Water", overlay = TRUE, partial = TRUE,scale="response",ylim=c(0,2000))

emtrends(hurdle_Model_LA_quad, specs = c("Water"), var = "S_elev",type="response", component="cond")
emtrends(hurdle_Model_LA_quad, specs = c("Water"), var = "Selev2",type="response", component="cond")


##Optimal source elevations when we just consider the water treatment:  0.2640449, 3242.76 (95% CI: 3142.12)
-  0.188 /(2*-0.356 ) #0.2640449
0.0378/(2* -0.656) #-0.02881098
-0.4141/(2*-0.0647) #3.200155


#Converting to elevation in m
0.2640449*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  
-0.02881098*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  
3.200155*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  



pred_hurdle <- ggpredict(hurdle_Model_LA, terms = c("S_elev[all]", "Water"), type = "random", interval="confidence")

hurdle_plot <-plot(pred_hurdle, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity")+ylim(0,2000)
hurdle_plot


##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
(2600-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)
(3100-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)
(3600-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)



fecundity_cline_water<-visreg(hurdle_Model_LA,"S_elev", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

fec_clinal_variation_water<-ggplot(fecundity_cline_water $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= RM_fit, aes(S_elev, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.564114,-0.1030306,1.358053))+ scale_y_continuous("Fecundity (length of longest fruit)",limits=c(0,2500))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Water)
fec_clinal_variation_water



##Test the random effects

hurdle_Model_LA_plant_count <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_LA_plant_zero <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_LA_cage_count <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_LA_cage_zero <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_LA_geno_count <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_LA_geno_zero <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block),family=ziGamma(link="log"))

anova(hurdle_Model_LA, hurdle_Model_LA_plant_count)
anova(hurdle_Model_LA, hurdle_Model_LA_plant_zero)
anova(hurdle_Model_LA, hurdle_Model_LA_cage_count)
anova(hurdle_Model_LA, hurdle_Model_LA_cage_zero)
anova(hurdle_Model_LA, hurdle_Model_LA_geno_count)
anova(hurdle_Model_LA, hurdle_Model_LA_geno_zero)



#### Local adaptation analysis using lifetime fitness ####

##Create lifetime fitness data file by summarizing fitness across years
repro_total <- aggregate(Reproduced ~ PlantID:Water:Herbivore, grasshopper, FUN=sum) 
fecundity_averages <- aggregate(Mature_length_siliques ~ PlantID:Water:Herbivore, grasshopper, FUN=sum) 
lifefit<-merge(repro_total,fecundity_averages, by=c("PlantID", "Water","Herbivore"),all=TRUE)
grasshopperMain<-subset(grasshopper,year=="2021")
grasshopperMaster <-grasshopperMain[c(1:14,16)]
lifefitness<-merge(grasshopperMaster, lifefit, by=c("PlantID", "Water","Herbivore"),all=TRUE)
head(lifefitness)

lifefitness $S_initdiam<-scale(lifefitness $init.diam,center=TRUE, scale=TRUE)
lifefitness $S_elev <-scale(lifefitness $elevation,center=TRUE, scale=TRUE)
hist(lifefitness$Reproduced)
hist(lifefitness$Mature_length_siliques)
lifefitness$Reproduction<-lifefitness$Reproduced
lifefitness$Reproduction[lifefitness$Reproduction == '2'] <- 1
hist(lifefitness$Reproduction)
##write.csv(lifefitness,"lifefitness.csv")
lifefitness <- read.csv("lifefitness.csv",stringsAsFactors=T)
life_fit <- drop_na(lifefitness, Mature_length_siliques) 

hurdle_Model_local_adapt <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore*I(S_elev^2)+(1|Cage_Block)+(1|Genotype), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore*I(S_elev^2)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
emtrends(hurdle_Model_local_adapt, specs = c("Water"), var = "S_elev",component="cond")

#To check residuals
simulationOutput <- simulateResiduals(fittedModel= hurdle_Model_local_adapt, plot = T, re.form = NULL,allow.new.levels =TRUE)

##This is the ANOVA table for the logistic regression part (probability of reproduction).
Anova(hurdle_Model_local_adapt,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_local_adapt,type="III", component="cond")


pred_hurdle <- ggpredict(hurdle_Model_local_adapt, terms = c("S_elev[all]", "Water"), type = "random", interval="confidence")

hurdle_plot <-plot(pred_hurdle, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity")+ylim(0,2000)
hurdle_plot

##Test the random effects
hurdle_Model_local_adapt_cage_count <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore*I(S_elev^2)+(1|Genotype), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore*I(S_elev^2)+(1|Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
hurdle_Model_local_adapt_cage_zero <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore*I(S_elev^2)+(1|Cage_Block)+(1|Genotype), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore*I(S_elev^2)+(1|Genotype),family=ziGamma(link="log"))

hurdle_Model_local_adapt_geno_count <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore*I(S_elev^2)+(1|Cage_Block), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore*I(S_elev^2)+(1|Genotype)+ (1| Cage_Block),family=ziGamma(link="log"))
hurdle_Model_local_adapt_geno_zero <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore*I(S_elev^2)+(1|Genotype)+(1|Cage_Block), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore*I(S_elev^2)+ (1| Cage_Block),family=ziGamma(link="log"))
anova(hurdle_Model_local_adapt, hurdle_Model_local_adapt_cage_count)
anova(hurdle_Model_local_adapt, hurdle_Model_local_adapt_cage_zero)
anova(hurdle_Model_local_adapt, hurdle_Model_local_adapt_geno_count)
anova(hurdle_Model_local_adapt, hurdle_Model_local_adapt_geno_zero)



##Extract the coefficients using emtrends, which requires defining the quadratic effect prior to running the model
life_fit $Selev2<-life_fit $S_elev* life_fit $S_elev


hurdle_Model_local_adapt_quad <- glmmTMB(Mature_length_siliques ~ S_initdiam+Water*Herbivore* S_elev+Water*Herbivore* Selev2
                                         +(1|Cage_Block)+(1|Genotype), data= life_fit, zi=~S_initdiam+Water*Herbivore*S_elev+Water*Herbivore* Selev2
                                         + (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
visreg(hurdle_Model_local_adapt_quad,"S_elev", by="Water", overlay = TRUE, partial = TRUE,scale="response",ylim=c(0,2000))

emtrends(hurdle_Model_local_adapt_quad, specs = c("Water"), var = "S_elev",type="response", component="cond")
emtrends(hurdle_Model_local_adapt_quad, specs = c("Water"), var = "Selev2",type="response", component="cond")
##Optimal source elevations when we just consider the water treatment:  0.4021352, 3267.04 (95% CI: 3136.13)
-  0.226 /(2*-0.281 ) 
-0.0195/(2* -0.520)
-0.432901/(2*-0.042)
#Converting to elevation in m
0.4021352*sd( lifefitness $elevation,na.rm = TRUE)+mean(lifefitness $elevation, na.rm = TRUE)  
0.01875*sd( lifefitness $elevation,na.rm = TRUE)+mean(lifefitness $elevation, na.rm = TRUE)  
5.153583*sd( lifefitness $elevation,na.rm = TRUE)+mean(lifefitness $elevation, na.rm = TRUE)  

pred_hurdle <- ggpredict(hurdle_Model_local_adapt_quad, terms = c("S_elev[all]", "Water"), type = "random", interval="confidence")

hurdle_plot <-plot(pred_hurdle, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Fecundity")+ylim(0,2500)
hurdle_plot



##To plot, we have to combine panels for the quadratic term for the supplemental watering treatment with the linear term from the restricted watering treatment. I did that by making these two plots and then using the quadratic supplemental and linear restricted in Illustrator along with recoding Source elevation fromo standardize units to raw units in meters

##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
(2600-mean(lifefitness $elevation, na.rm = TRUE))/sd( lifefitness $elevation,na.rm = TRUE)
(3100-mean(lifefitness $elevation, na.rm = TRUE))/sd( lifefitness $elevation,na.rm = TRUE)
(3600-mean(lifefitness $elevation, na.rm = TRUE))/sd( lifefitness $elevation,na.rm = TRUE)



fecundity_cline_water<-visreg(hurdle_Model_local_adapt,"S_elev", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

fec_clinal_variation_water<-ggplot(fecundity_cline_water $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= life_fit, aes(S_elev, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.551343,-0.08705503,1.377233))+ scale_y_continuous("Fecundity (length of longest fruit)",limits=c(0,2500))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Water)
fec_clinal_variation_water



fecundity_cline_water_restricted<-visreg(hurdle_Model_local_adapt_quad,"S_elev", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

fec_clinal_variation_water_restricted<-ggplot(fecundity_cline_water_restricted $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= life_fit, aes(S_elev, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.551343,-0.08705503,1.377233))+ scale_y_continuous("Fecundity (length of longest fruit)",limits=c(0,2500))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Water)
fec_clinal_variation_water_restricted



fecundity_cline_herb<-visreg(hurdle_Model_local_adapt,"S_elev", by="Herbivore", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

fec_clinal_variation_herb<-ggplot(fecundity_cline_herb $fit, aes(S_elev, visregFit,group= Herbivore, colour= Herbivore, fill=factor(Herbivore))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Herbivore)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= life_fit, aes(S_elev, Mature_length_siliques, color= Herbivore, shape=Herbivore), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.551343,-0.08705503,1.377233))+ scale_y_continuous("Fecundity (length of longest fruit)",limits=c(0,2500))+scale_fill_manual(values = cols2, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
fec_clinal_variation_herb



fecundity_cline_water_restricted<-visreg(hurdle_Model_local_adapt_quad,"S_elev", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

fec_clinal_variation_water_restricted<-ggplot(fecundity_cline_water_restricted $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= life_fit, aes(S_elev, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.551343,-0.08705503,1.377233))+ scale_y_continuous("Fecundity (length of longest fruit)",limits=c(0,2500))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Water)
fec_clinal_variation_water_restricted






####plotting figures####

#Figure3

damage_cline 

library(ggpubr)
Fig2<-ggarrange(damage_cline, damage_treatment,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
Fig2

Fig2<-ggarrange(damage_cline, damage_treatment,LAR_fecund_plot,
                    labels = c("A", "B","C"),
                    ncol = 3, nrow = 1)
Fig2


Figure2 <- ggarrange(damage_clinal_variation, 
  LAR_box,LAR_fecund_plot,
   labels = c("A", "B","C"),
  ncol = 3, nrow = 1, common.legend = TRUE, legend="right")
Figure2

Fig2 <- ggarrange(
  damage_clinal_variation, 
  LAR_box,LAR_fecund_plot,
   labels = c("A", "B"),
  ncol = 1, nrow = 3, common.legend = TRUE, legend="right")
Fig2



library(ggpubr)
figure3 <- ggarrange(
  SLA_clinal_variation, SLA_treatment, SLA_fecund_plot, 
  #labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 3, nrow = 1)
figure3




library(ggpubr)
figure3 <- ggarrange(
 SLA_clinal_variation, SLA_treatment, sla_repro_plot, SLA_fecund_plot,succ_clinal_variation, Suc_treatment, suc_repro_plot, SUC_fecund_plot , 
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 4, nrow = 2)
figure3


#Figure4

library(ggpubr)
figure4 <- ggarrange(
  phen_cline,
  FT_treatment,
  FT_fecund_plot,
  Duration_cline,
  Duration_treatment,
  duration_fecund_plot,
  height_cline,
  height_treatment,
  height_fecund_plot,
  labels = c("A", "B", "C","D", "E", "F","G", "H", "I"),
  ncol = 3, nrow = 3)
figure4

library(ggpubr)
figure4 <- ggarrange(
  phen_cline,
  FT_treatment,
  FT_fecund_plot,
  Fduration_clinal_variation,
  Duration_treatment,
  duration_fecund_plot,
  height_clinal_variation,
  height_treatment,
  height_fecund_plot,
  labels = c("A", "B", "C","D", "E", "F","G", "H", "I"),
  ncol = 3, nrow = 3)
figure4


figure4 <- ggarrange(
  phen_cline,
  FT_treatment,
  FT_fecund_plot,
  #Fduration_clinal_variation,
  #Duration_treatment,
  #duration_fecund_plot,
  #height_clinal_variation,
  #height_treatment,
  #height_fecund_plot,
  #labels = c("A", "B", "C","D", "E", "F","G", "H", "I"),
  ncol = 3, nrow = 2)
figure4





#Figure5

figure5 <- ggarrange(
  #repro_treatment,
  Local_adaptation,
  fecundity_treatment,
  
  
  labels = c("A", "B"),
  ncol = 1, nrow = 2)

figure5

###


