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

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")
#setwd("~/Documents/personnel/Jameel/grasshopper")

#this is where you specify the folder where you have the data on your computer


#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_updated_21March24.csv",stringsAsFactors=T)

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

##Let's correlate rosette and bolt leaf data to see if we can come up with composite figures
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

##Create composite leaf SLA, succuclence and lwc variables based on the regressions above. This gives us foliar trait data for the 67 plants for which we have bolt but not rosette collections
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


VWC_plot # Fig SX



### Hypothesis one #####

#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************

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
LAR_data <- drop_na(LAR_data,LAR_prop) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

LAR_data $y_beta<- (LAR_data $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

min(LAR_data $y_beta)

max(LAR_data $y_beta)




#Box_plot
LAR_box <-ggplot(LAR_data, aes(x = Herbivore, y = LAR_prop, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

LAR_box + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~year)


## model
LAR_Model<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)



#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>%  
  mutate(fit.m = predict(LAR_Model, se.fit=FALSE),              
         resid = residuals(LAR_Model))

##Convert coefficients to probabilities
newdf2 $predicted<-1/(1+exp(-(newdf2 $fit.m))) 
newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 

##Box plot
LAR_box <-ggplot(LAR_data, aes(x = Water, y = LAR_prop, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ 
  scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

damage_treatment<-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~year+Herbivore)



#cline, plotting predicted data
damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted, group= Water, 
                                 colour= Water))+geom_point(pch = 20) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + facet_grid(Herbivore~ year) + theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Water availability", labels = c("Restricted","Supplemental"))

damage_cline #fig 2


##Conversation of standardized elevation back to raw units
-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)




#env concatenates water and herbivore, needed for slope calculations

LAR_data $env<-interaction(LAR_data $Water, LAR_data $Herbivore)

LAR_data $env_year<-interaction(LAR_data $env, LAR_data $year)

LAR_ModelE<- gamlss (formula= y_beta ~S_elev* env_year-1+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(LAR_ModelE)

plot(LAR_Model)



#env concatenates water and herbivore. This model allows us to extract means and slopes for each combination of year, treatment and garden
LAR_data $env<-interaction(LAR_data $Water, LAR_data $Herbivore)
LAR_data $env_year<-interaction(LAR_data $env, LAR_data $year)
LAR_ModelE<- gamlss (formula= y_beta ~S_elev* env_year-1+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(LAR_ModelE)


#Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+S_elev*Water*year+S_elev* Herbivore*year+ Water* Herbivore*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)

#Subsets of models for drop 1
LAR_Model_two<- gamlss (formula= y_beta ~S_elev*Water+S_elev*Herbivore+ Water* Herbivore+year*Water+year*Herbivore+S_elev*year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)

LAR_Model_one<- gamlss (formula= y_beta ~S_elev+Water+Herbivore+ year+ random(census) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)


#*******************************************************************************
#### Specific Leaf Area #####
#*******************************************************************************

foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)


sla_model <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


SLA <-emmeans(sla_model, ~ Water*Herbivore, type="response", adjust = "sidak")
cld(SLA, details=TRUE)
Anova(sla_model, type = "III") 

coefficients_SLA <- emtrends(sla_model, specs = c("year"), var = "S_elev",type="response")
SLA_table<- as.data.frame(summary(coefficients_SLA))[c('year','S_elev.trend', 'SE')]
SLA_table <- SLA_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

sla_model_no_plantID <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))

sla_model_no_genotype <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))

sla_model_no_cageblock <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype), data = foliar, family= lognormal(link="log"))


anova(sla_model,sla_model_no_plantID)
anova(sla_model,sla_model_no_genotype)
anova(sla_model,sla_model_no_cageblock)

pred_sla <- ggpredict(sla_model, terms = c("S_elev[all]", "year"), type = "re", interval="confidence")


sla_cline <-plot(pred_sla, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = c("black","black"), facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Specific leaf area")


sla_cline


coefficients_SLA <- emtrends(sla_model, specs = c("year","Herbivore","Water"), var = "S_elev",type="response")
SLA_table<- as.data.frame(summary(coefficients_SLA))[c('year',"Herbivore","Water",'S_elev.trend', 'SE')]
SLA_table <- SLA_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

pred_sla <- ggpredict(sla_model, terms = c("S_elev[all]", "Water","Herbivore","year"), type = "re", interval="confidence")

sla_cline <-plot(pred_sla, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = c(cols), facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Specific leaf area")


sla_cline

##Box plot
#SLA_box <-ggplot(foliar, aes(x = Herbivore, y = SLA, fill = Water)) +
#  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Specific leaf area") +
#  scale_y_continuous(element_blank()) +
#  geom_point(pch = 21, position = position_jitterdodge(0.3))

#SLA_treatment<-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "Top")+ scale_x_discrete(labels=c("Grasshopper addition", "Grasshopper removal"))+ scale_fill_manual(values = c(cols), name = "Herbivore manipulation", labels = c("Water-restricted","Supplemental watering"))


SLA_box <-ggplot(foliar, aes(x = Water, y = SLA, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Specific leaf area") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

SLA_treatment<-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "Top")+ scale_x_discrete()+ scale_fill_manual(values = c(cols), name = "Water availability", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)


#23 cline, not significant
pred_sla23 <- ggpredict(sla_model, terms = c("S_elev[all]", "year[2023]"))

sla_cline23 <-plot(pred_sla23, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = c("black","black"), facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous(element_blank())+ scale_y_continuous(element_blank())


sla_cline23



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

succulence <-emmeans(succ_model, ~ Herbivore:year, type="response", adjust = "sidak")
cld(succulence, details=TRUE)


pred_succ <- ggpredict(succ_model, terms = c("S_elev[all]","Water","Herbivore"), type = "re", interval="confidence")


Succulence_cline <-plot(pred_succ, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Leaf succulence (mg)")

Succulence_cline 

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





####


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


coefficients_FT <- emtrends(FT_model, specs = c("Water","Herbivore"), var = "S_elev",type="response")
FT_table<- as.data.frame(summary(coefficients_FT))[c('Water',"Herbivore",'S_elev.trend', 'SE')]
FT_table <- FT_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

coefficients_FT <- emtrends(FT_model, specs = c("Water","year"), var = "S_elev",type="response")
FT_table<- as.data.frame(summary(coefficients_FT))[c('Water',"year",'S_elev.trend', 'SE')]
FT_table <- FT_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))




FT_model_no_geno<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
FT_model_no_cageblock<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1| PlantID),data= flowering,family=lognormal(link="log"))
FT_model_noplantID<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

anova(FT_model,FT_model_no_geno)
anova(FT_model,FT_model_no_cageblock)
anova(FT_model,FT_model_noplantID)



pred_FT <- ggpredict(FT_model, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_cline <-plot(pred_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")

phen_cline #main figure, fig1


FT_box <-ggplot(flowering, aes(x = Water, y = FT_Adj, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

FT_treatment <-FT_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)


pred_FT2 <- ggpredict(FT_model, terms = c("S_elev[all]", "Water","year"), type = "re", interval="confidence")
phen_cline_year <-plot(pred_FT2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")

####


phen_cline_year #supplemental figure, S5



#*******************************************************************************
#### Flowering duration #####
#*******************************************************************************


flowering_duration<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)

Anova(flowering_duration, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)


flowering_duration_nogeno<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Cage_Block)+(1| PlantID),data= flowering)

flowering_duration_nocageblock<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+(1|Genotype)+(1| PlantID),data= flowering)

flowering_duration_noplantid<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block),data= flowering)


anova(flowering_duration,flowering_duration_nogeno)
anova(flowering_duration,flowering_duration_nocageblock)
anova(flowering_duration,flowering_duration_noplantid)


duration <-emmeans(flowering_duration, ~ Herbivore:Water, type="response", adjust = "sidak")
cld(duration, details=TRUE)

pred_duration <- ggpredict(flowering_duration, terms = c("S_elev[all]","Water","Herbivore"), type = "re", interval="confidence")

Duration_cline <-plot(pred_duration, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(#text = element_text(size=10),
  axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering Duration (days)")


Duration_cline <-plot(pred_duration, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Flowering duration (days)") 

Duration_cline #no significant interaction

Duration_box <-ggplot(flowering, aes(x = Herbivore, y = flowering_duration, fill = Herbivore)) +
  geom_boxplot(outlier.shape = NA) +xlab(" Herbivore")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

Duration_treatment <-Duration_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_wrap(~Herbivore)



Duration_box <-ggplot(flowering, aes(x = Water, y = flowering_duration, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availibility")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

Duration_treatment <-Duration_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)

####

figureS6 <- ggarrange(
  Duration_cline,
  Duration_treatment,
  
  
  labels = c("A", "B"),
  ncol = 2, nrow = 1)
figureS6



#*******************************************************************************
#### Tallest bolt at flowering #####
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

coefficients_height <- emtrends(max_height, specs = c("Water","Herbivore"), var = "S_elev",type="response")
Height_table<- as.data.frame(summary(coefficients_height))[c('Water',"Herbivore",'S_elev.trend', 'SE')]
Height_table <- Height_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

pred_height <- ggpredict(max_height, terms = c("S_elev[all]", "Water", "Herbivore"), type = "re", interval="confidence")


height_cline <-plot(pred_height, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Height of tallest bolt at flowering")



height_cline #supplemental figure, S6

height_box <-ggplot(flowering, aes(x = Water, y = Max_height_flowering, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

height_treatment <-height_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)


height_box <-ggplot(flowering, aes(x = year, y = Max_height_flowering, fill = year)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

height_treatment <-height_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_wrap(~Herbivore)



### Hypothesis two ######

#*******************************************************************************
####  Probability of reproduction #####
#*******************************************************************************
grasshopperFT <- filter(grasshopper, year != "2021")

prob_repro <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))

Anova(prob_repro, type="III")
simulationOutput <- simulateResiduals(fittedModel= prob_repro, plot = T, re.form = NULL,allow.new.levels =TRUE)

testOutliers(simulationOutput, type = c(
  "bootstrap"), nBoot = 100, plot = T)

prob_repro_nogeno <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year+(1|PlantID) + (1|Cage_Block), data= grasshopperFT,family=binomial(link="logit"))

prob_repro_nopid <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year + (1|Cage_Block)+(1|Genotype), data= grasshopperFT,family=binomial(link="logit"))

prob_repro_nocageblock <- glmmTMB(Reproduced ~S_initdiam+ Water*Herbivore * S_elev *year+(1|PlantID) + (1|Genotype), data= grasshopperFT,family=binomial(link="logit"))


anova(prob_repro,prob_repro_nogeno)
anova(prob_repro,prob_repro_nopid)
anova(prob_repro,prob_repro_nocageblock)


repro<-emmeans(prob_repro, ~ Water:Herbivore:year, type="response", adjust = "sidak")
cld(repro, details=TRUE)



##Box plot
Prob_violin <-ggplot(grasshopperFT, aes(x = Water, y = Reproduced, fill = Water)) +
  geom_violin() +xlab("Water availibility")+ 
  stat_summary(fun.data = "mean_cl_normal")+
  scale_y_continuous("Probability of reproduction") +
  geom_point(pch = 21, position = position_jitterdodge())

repro_treatment<-test + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~year+Herbivore)




repro_herb<- ggplot(grasshopperFT, aes(x = Water, y = Reproduced, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.5, width = 0.33) +
  facet_grid(~ year+Herbivore) +
  theme_classic() +
  scale_fill_manual(values = cols) +
  #geom_point(size=0.5,fill = "white", position = position_jitterdodge()) +
  geom_jitter(shape = 21, fill = "NA", size = 2, alpha = 0.5, width = 0.3, height = 0)+
  #labs(y = "Probability of reproduction") + 
  labs(y = element_blank()) +
  labs(x = "Herbivore treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none")

#prepro_box <-ggplot(grasshopperFT, aes(x = Herbivore, y = Reproduced, fill = Water)) +
#  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
#  scale_y_continuous(element_blank()) +
#  geom_point(pch = 21, position = position_jitterdodge(0.3))

#prepro_treatment <-prepro_box + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Grasshopper addition", "Grasshopper removal")) +  scale_fill_manual(values = c(cols2), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))#+facet_wrap(~Herbivore)


repro_herb #supplemental figure




#*******************************************************************************
####  Fecundity #####
#*******************************************************************************
repro<-subset(grasshopperFT, Reproduced=="1")

#repro $Water1 <-factor(repro $Water, levels = c("Supplemental", "Restricted"))
repro $Water <-factor(repro $Water, levels = c("Restricted","Supplemental"))


fecundity <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*year* S_elev+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
Anova(fecundity,type="III")

summary(fecundity)



fecundity_nopid <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*year* S_elev + (1|Cage_Block)+(1|Genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

fecundity_nocb <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*year* S_elev+(1|PlantID) +(1|Genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

fecundity_nogeno <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*year* S_elev+(1|PlantID) + (1|Cage_Block), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(fecundity,fecundity_nopid)
anova(fecundity,fecundity_nocb)
anova(fecundity,fecundity_nogeno)


coefficients_fec <- emtrends(fecundity, specs = c("Water"), var = "S_elev",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('Water','S_elev.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table

Fecmeans<-emmeans(fecundity, ~ Water:Herbivore, type="response", adjust = "sidak")
cld(Fecmeans, details=TRUE)

plot(Fecmeans)


#pull the fitted values out and plot them in ggplot

newdf2 <- repro %>%  
  mutate(fit.m = predict(fecundity, se.fit=FALSE),              
         resid = residuals(fecundity))

##Convert coefficients to probabilities
newdf2 $predicted<-exp((newdf2 $fit.m)) 
newdf2 $resid_trans<-exp((newdf2 $resid)) 


pred_fec <- ggpredict(fecundity, terms = c("S_elev[all]","Water"), type = "re", interval="confidence")


Local_adaptation <-plot(pred_fec, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,2000)

Local_adaptation

fec_box <-ggplot(newdf2, aes(x = Water, y = predicted, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ 
  scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_x_continuous("Fecundity (mature fruit length)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

fecundity_treatment <-fec_box + theme_classic() + ylim(0,1000) +theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)


fec_box <-ggplot(repro, aes(x = Water, y = Mature_length_siliques, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ 
  scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_x_continuous("Fecundity (mature fruit length)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

fecundity_treatment <-fec_box + theme_classic() + ylim(0,2000) +theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete() +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)


####

fec_box_herb <-ggplot(repro, aes(x = Herbivore, y = Mature_length_siliques, fill = Herbivore)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  #scale_y_continuous("Flowering phenology (ordinal day of year)") +
  scale_y_continuous(element_blank()) +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

fecundity_treatment_herb <-fec_box_herb + theme_classic() + theme(axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Grasshopper addition", "Grasshopper removal")) +  scale_fill_manual(values = c(cols2), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~year)

pred_fec_year <- ggpredict(fecundity, terms = c("S_elev[all]","year"), type = "re", interval="confidence")

fec_cline_year <-plot(pred_fec_year, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = c("black","black"), facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation (m)")+ scale_y_continuous("Fecundity (length of mature siliques)")+ ylim(0,2000) #supplemental figure

coefficients_fec2 <- emtrends(fecundity, specs = c("year"), var = "S_elev",type="response")
fec_table2<- as.data.frame(summary(coefficients_fec2))[c('year','S_elev.trend', 'SE')]
fec_table2 <- fec_table2%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table2


#fecundity_treatment <- ggplot(repro, aes(x = Water, y = Mature_length_siliques, group = Water)) +
#  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
#  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
#  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
#  facet_grid(~ Herbivore) +
#  theme_light() +
#  scale_fill_manual(values = cols) +
#  labs(y = "Fecundity (length of mature siliques)") + 
#  labs(x = "Watering treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
#  theme(legend.position = "none")

#fecundity_treatment




###


figureS5 <- ggarrange(
  fec_cline_year,
  fecundity_treatment_herb,
  
  
  labels = c("A", "B"),
  ncol = 2, nrow = 1)
figureS5





### Hypothesis three ######

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced_updated
grasshopper_no2021<-subset(grasshopper, year!="2021")
# retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper,Flowered,FT_Adj,Max_height_flowering, Max_height_peak ,Season_survival, avg_leafnumber, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block, elevation,elev_km,  year,Stem_number_flowering, rosette_succulence,Ordinal_Date_flowering,Sum_height_flowering,Date_peak_flowering,Stem_number_peak,Mature_length_siliques,Reproduced,flowering_duration,days_until_mortality, SLA, lwc, succulence, rosette_SLA, rosette_lwc,FT_Adj,Snowmelt_FT_Adj,max_LAR)

##only plants that Reproduced_updated (no longer in dataset)
#Reproduced_updated<-subset(traitdat,Reproduced_updated=="1")
#ggplot(Reproduced_updated, aes(x= Mature_length_siliques_updated))+ geom_histogram()+ facet_grid(treat ~ .)
#ggplot(Reproduced_updated, aes(x= Mature_length_siliques_updated))+ geom_histogram(color="black", fill="white")


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

traitdat $sLAR_max<-scale(traitdat $max_LAR,center=TRUE,scale=TRUE)
traitdat $sLAR<-scale(traitdat $avg_LAR,center=TRUE,scale=TRUE)
traitdat $sSucc<-scale(traitdat $rosette_succulence,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sLWC<-scale(traitdat $lwc,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)
traitdat $sFT<-scale(traitdat $FT_Adj,center=TRUE,scale=TRUE)


head(traitdat)


#Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Addition","Removal"))


#*******************************************************************************
####  Probability of reproduction for vegetative traits only  #####
#*******************************************************************************


repro_model_full <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+
                             Water*Herbivore*sSUC+I(sSUC^2)+
                             Water*Herbivore*sSLA+I(sSLA^2) 
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
Anova(repro_model_full,type="III") 

simulationOutput <- simulateResiduals(fittedModel= repro_model_full, plot = T, re.form = NULL)
testOutliers(simulationOutput, type="bootstrap")

repro_model_noplantID <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+
                             Water*Herbivore*sSUC+I(sSUC^2)+
                             Water*Herbivore*sSLA+I(sSLA^2) 
                           +(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))


repro_model_full_nocb <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+
                             Water*Herbivore*sSUC+I(sSUC^2)+
                             Water*Herbivore*sSLA+I(sSLA^2) 
                           +(1|PlantID)+(1|Genotype),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

repro_model_full_nogeno <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                             Water*Herbivore*sLAR+
                             Water*Herbivore*sSUC+I(sSUC^2)+
                             Water*Herbivore*sSLA+I(sSLA^2) 
                           +(1|PlantID)+(1|Cage_Block),data=traitdat,family=binomial(link="logit"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

anova(repro_model_full,repro_model_noplantID)
anova(repro_model_full,repro_model_full_nocb)
anova(repro_model_full,repro_model_full_nogeno)




#SLA
emtrends(repro_model_full, specs = c("Water"), var = "sSLA")

repro <- emtrends(repro_model_full, specs = c("Water"), var = "sSLA",type="response")
repro_table<- as.data.frame(summary(repro))[c('Water','sSLA.trend', 'SE')]
repro_table <- repro_table%>% mutate(
  oddsratio = exp(sSLA.trend),
  Lower95 = oddsratio * exp(-1.96*SE),
  Upper95 = oddsratio * exp(1.96*SE),
  Quadslopes = (2*oddsratio),
  QuadSE = (2*exp(SE)),
  QLower95 = Quadslopes * (-1.96*QuadSE),
  QUpper95 = Quadslopes * (1.96*QuadSE))
repro_table


#succulence
emtrends(repro_model_full, specs = c("Herbivore"), var = "sSUC")

quad_repro <- emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sduration",type="response")
fec_table<- as.data.frame(summary(quad_fec))[c('Water',"Herbivore",'sduration.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  Quadslopes = (2*sduration.trend),
  QuadSE = (2*SE),
  Lower95 = Quadslopes * (-1.96*QuadSE),
  Upper95 = Quadslopes * (1.96*QuadSE))




sSUC_pred <- ggpredict(repro_model_full, terms = c("sSUC[all]","Herbivore"), type = "re", interval="confidence")

Succulence_reproduction_herbivore <-plot(sSUC_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols2,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")


sla_pred <- ggpredict(repro_model_full, terms = c("sSLA[all]", "Water"), type = "re", interval="confidence")

SLA_reproduction_water <- plot(sla_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")


#*******************************************************************************
####  Fecundity  #####
#*******************************************************************************



traitdatRepro <- filter(traitdat, Reproduced == 1 ) 
##Rescale once we have removed the non-reproductive plants
traitdatRepro $sduration<-scale(traitdatRepro $flowering_duration,center=TRUE,scale=TRUE)
traitdatRepro $sheight<-scale(traitdatRepro $Sum_height_flowering,center=TRUE,scale=TRUE)
traitdatRepro $s_max_height<-scale(traitdatRepro $Max_height_flowering,center=TRUE,scale=TRUE)
traitdatRepro $sLAR_max<-scale(traitdatRepro $max_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sLAR<-scale(traitdatRepro $avg_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sSLA<-scale(traitdatRepro $SLA,center=TRUE,scale=TRUE)
traitdatRepro $sSUC<-scale(traitdatRepro $succulence,center=TRUE,scale=TRUE)
traitdatRepro $sFT<-scale(traitdatRepro $FT_Adj,center=TRUE,scale=TRUE)

fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore* s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2) +
                               Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+
                               Water*Herbivore* sLAR+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits,type="III") 

simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)



fecund_modeltraits_nocb <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore* s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2) +
                               Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+
                               Water*Herbivore* sLAR+ (1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

fecund_modeltraits_nogeno <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
                               Water*Herbivore*sFT+
                               Water*Herbivore* s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2) +
                               Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+
                               Water*Herbivore* sLAR+ (1|Cage_Block),data=traitdatRepro,family=Gamma(link="log"))


anova(fecund_modeltraits,fecund_modeltraits_nocb)
anova(fecund_modeltraits,fecund_modeltraits_nogeno)




#SLA
emtrends(fecund_modeltraits, specs = c("Herbivore"), var = "I(sSLA^2)")

#Context dependent selection
standard_sSLA <- ggpredict(fecund_modeltraits, terms = c("sSLA[all]","Water","Herbivore"), type = "re", interval="confidence")

SLA_fecundity <-plot(standard_sSLA, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)


#succulence
emtrends(fecund_modeltraits, specs = c("Herbivore"), var = "sSUC")

#Context dependent selection

standard_sSUC <- ggpredict(fecund_modeltraits, terms = c("sSUC[all]","Water", "Herbivore"), type = "re", interval="confidence")

succulence_fecundity <-plot(standard_sSUC, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Leaf succulence")+ #scale_y_continuous("Fecundity (Summed length of siliques)")
scale_y_continuous(element_blank())+ ylim(0,2000)


#Flowering phenology
emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sFT")

#Context dependent selection
standard_FT <- ggpredict(fecund_modeltraits, terms = c("sFT[all]", "Water","Herbivore"), type = "re", interval="confidence")

FT_fecundity <-plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = c(cols), facet=TRUE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Flowering phenology")+ #scale_y_continuous("Fecundity (Summed length of siliques)")
scale_y_continuous(element_blank()) + ylim(0,2000)


#Duration
emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sduration")

quad_fec <- emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sduration",type="response")
fec_table<- as.data.frame(summary(quad_fec))[c('Water',"Herbivore",'sduration.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  Quadslopes = (2*sduration.trend),
  QuadSE = (2*SE),
  Lower95 = Quadslopes * (-1.96*QuadSE),
  Upper95 = Quadslopes * (1.96*QuadSE))

#Context dependent selection
standard_duration<- ggpredict(fecund_modeltraits, terms = c("sduration[all]", "Water","Herbivore"), type = "re", interval="confidence")

duration_fecundity <-plot(standard_duration, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = c(cols), facet=TRUE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Flowering duration")+ylim(0,2000)+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)



#height

standard_height<- ggpredict(fecund_modeltraits, terms = c("s_max_height[all]","Water","Herbivore"), type = "re", interval="confidence")


height_fecundity <-plot(standard_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="none")+scale_x_continuous("Height of tallest bolt at flowering")+#scale_y_continuous("Fecundity (Summed length of siliques)")
  scale_y_continuous(element_blank()) + ylim(0,2000)







##No selection
standard_sLAR <- ggpredict(fecund_modeltraits, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")

LAR_fecundity <-plot(standard_sLAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="none")+ ylim(0,2000)+scale_y_continuous("Fecundity (Summed length of siliques)")+scale_x_continuous("Leaf area removed by herbivores (Proportion)")




####### plotting figures #####





#Figure2

library(ggpubr)
figure2 <- ggarrange(
  damage_cline,
  damage_treatment,
  
  #LAR_fecundity,
  
  labels = c("A", "B"),
  ncol = 1, nrow = 2)
figure2


#Figure3

library(ggpubr)
figure3 <- ggarrange(
  sla_cline,
  SLA_treatment,
  SLA_reproduction_water,
  SLA_fecundity,
  
  Succulence_cline,
  Suc_treatment,
  Succulence_reproduction_herbivore,
  succulence_fecundity,
  
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 2, nrow = 4)
figure3


#Figure4

library(ggpubr)
figure4 <- ggarrange(
  phen_cline,
  FT_treatment,
  FT_fecundity,
  Duration_cline,
  Duration_treatment,
  duration_fecundity,
  height_cline,
  height_treatment,
  height_fecundity,
  labels = c("A", "B", "C","D", "E", "F","G", "H", "I"),
  ncol = 3, nrow = 3)
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



#####local trait means####

##### addition restricted
local_grasshopper_RA <- filter(grasshopper_no2021, Herbivore == "Addition", Water == "Restricted")
local_grasshopper_RA <- filter(local_grasshopper_RA, elev_dist >= -10, elev_dist <= 13)


#LAR 3.05586, 0.8239485
mean(local_grasshopper_RA$avg_LAR,na.rm = TRUE)
sd(local_grasshopper_RA$avg_LAR,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$avg_LAR[!is.na(local_grasshopper_RA$avg_LAR)]))


#SLA 196.9138, 17.03393
mean(local_grasshopper_RA$SLA,na.rm = TRUE)
sd(local_grasshopper_RA$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$SLA[!is.na(local_grasshopper_RA$SLA)]))

#Succulence 4.784538, 1.100558
mean(local_grasshopper_RA$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_RA$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$succulence[!is.na(local_grasshopper_RA$succulence)]))*1000


#day of flowering 173.6667, 2.603417
mean(local_grasshopper_RA$FT_Adj,na.rm = TRUE)
sd(local_grasshopper_RA$FT_Adj,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$FT_Adj[!is.na(local_grasshopper_RA$FT_Adj)]))


#duration 19.83333, 1.661659
mean(local_grasshopper_RA$flowering_duration,na.rm = TRUE)
sd(local_grasshopper_RA$flowering_duration,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$flowering_duration[!is.na(local_grasshopper_RA$flowering_duration)]))


#height 33.63333, 6.186904
mean(local_grasshopper_RA$Max_height_flowering,na.rm = TRUE)
sd(local_grasshopper_RA$Max_height_flowering,na.rm = TRUE) / sqrt(length(local_grasshopper_RA$Max_height_flowering[!is.na(local_grasshopper_RA$Max_height_flowering)]))




##### addition supplemental
local_grasshopper_AS <- filter(grasshopper_no2021, Herbivore == "Addition", Water == "Supplemental")
local_grasshopper_AS <- filter(local_grasshopper_AS, elev_dist >= -10, elev_dist <= 13)


#LAR 2.814267, 0.6896965
mean(local_grasshopper_AS$avg_LAR,na.rm = TRUE)
sd(local_grasshopper_AS$avg_LAR,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$avg_LAR[!is.na(local_grasshopper_AS$avg_LAR)]))


#SLA 307.3608, 28.70361
mean(local_grasshopper_AS$SLA,na.rm = TRUE)
sd(local_grasshopper_AS$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$SLA[!is.na(local_grasshopper_AS$SLA)]))

#Succulence 6.637142, 1.163733
mean(local_grasshopper_AS$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_AS$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$succulence[!is.na(local_grasshopper_AS$succulence)]))*1000


#day of flowering 177.5, 0.5
mean(local_grasshopper_AS$FT_Adj,na.rm = TRUE)
sd(local_grasshopper_AS$FT_Adj,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$FT_Adj[!is.na(local_grasshopper_AS$FT_Adj)]))


#duration 23, 1
mean(local_grasshopper_AS$flowering_duration,na.rm = TRUE)
sd(local_grasshopper_AS$flowering_duration,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$flowering_duration[!is.na(local_grasshopper_AS$flowering_duration)]))


#height 28.5, 5
mean(local_grasshopper_AS$Max_height_flowering,na.rm = TRUE)
sd(local_grasshopper_AS$Max_height_flowering,na.rm = TRUE) / sqrt(length(local_grasshopper_AS$Max_height_flowering[!is.na(local_grasshopper_AS$Max_height_flowering)]))




##### removal restricted
local_grasshopper_RR <- filter(grasshopper_no2021, Herbivore == "Removal", Water == "Restricted")
local_grasshopper_RR <- filter(local_grasshopper_RR, elev_dist >= -10, elev_dist <= 13)


#LAR 4.745702, 1.84663
mean(local_grasshopper_RR$avg_LAR,na.rm = TRUE)
sd(local_grasshopper_RR$avg_LAR,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$avg_LAR[!is.na(local_grasshopper_RR$avg_LAR)]))


#SLA 206.9137, 25.64905
mean(local_grasshopper_RR$SLA,na.rm = TRUE)
sd(local_grasshopper_RR$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$SLA[!is.na(local_grasshopper_RR$SLA)]))

#Succulence 5.540915, 0.797654
mean(local_grasshopper_RR$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_RR$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$succulence[!is.na(local_grasshopper_RR$succulence)]))*1000


#day of flowering 165.6667, 2.333333
mean(local_grasshopper_RR$FT_Adj,na.rm = TRUE)
sd(local_grasshopper_RR$FT_Adj,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$FT_Adj[!is.na(local_grasshopper_RR$FT_Adj)]))


#duration 17.66667, 3.480102
mean(local_grasshopper_RR$flowering_duration,na.rm = TRUE)
sd(local_grasshopper_RR$flowering_duration,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$flowering_duration[!is.na(local_grasshopper_RR$flowering_duration)]))


#height 18, 5.781868
mean(local_grasshopper_RR$Max_height_flowering,na.rm = TRUE)
sd(local_grasshopper_RR$Max_height_flowering,na.rm = TRUE) / sqrt(length(local_grasshopper_RR$Max_height_flowering[!is.na(local_grasshopper_RR$Max_height_flowering)]))






##### removal supplemental
local_grasshopper_RS <- filter(grasshopper_no2021, Herbivore == "Removal", Water == "Supplemental")
local_grasshopper_RS <- filter(local_grasshopper_RS, elev_dist >= -10, elev_dist <= 13)


#LAR 1.621326, 0.6246629
mean(local_grasshopper_RS$avg_LAR,na.rm = TRUE)
sd(local_grasshopper_RS$avg_LAR,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$avg_LAR[!is.na(local_grasshopper_RS$avg_LAR)]))


#SLA 232.9859, 19.25452
mean(local_grasshopper_RS$SLA,na.rm = TRUE)
sd(local_grasshopper_RS$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$SLA[!is.na(local_grasshopper_RS$SLA)]))

#Succulence 6.283142, 1.093187
mean(local_grasshopper_RS$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_RS$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$succulence[!is.na(local_grasshopper_RS$succulence)]))*1000


#day of flowering 172.6667, 6.765928
mean(local_grasshopper_RS$FT_Adj,na.rm = TRUE)
sd(local_grasshopper_RS$FT_Adj,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$FT_Adj[!is.na(local_grasshopper_RS$FT_Adj)]))


#duration 24.66667, 4.910307
mean(local_grasshopper_RS$flowering_duration,na.rm = TRUE)
sd(local_grasshopper_RS$flowering_duration,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$flowering_duration[!is.na(local_grasshopper_RS$flowering_duration)]))


#height 24.53333, 4.771559
mean(local_grasshopper_RS$Max_height_flowering,na.rm = TRUE)
sd(local_grasshopper_RS$Max_height_flowering,na.rm = TRUE) / sqrt(length(local_grasshopper_RS$Max_height_flowering[!is.na(local_grasshopper_RS$Max_height_flowering)]))




##### SLA restricted
local_grasshopper_R <- filter(grasshopper_no2021, Water == "Restricted")
local_grasshopper_R <- filter(local_grasshopper_R, elev_dist >= -10, elev_dist <= 13)

#SLA 201.7525, 14.97542
mean(local_grasshopper_R$SLA,na.rm = TRUE)
sd(local_grasshopper_R$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_R$SLA[!is.na(local_grasshopper_R$SLA)]))

##### SLA supplemental
local_grasshopper_S <- filter(grasshopper_no2021, Water == "Supplemental")
local_grasshopper_S <- filter(local_grasshopper_S, elev_dist >= -10, elev_dist <= 13)

#SLA 271.5506, 18.68434
mean(local_grasshopper_S$SLA,na.rm = TRUE)
sd(local_grasshopper_S$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper_S$SLA[!is.na(local_grasshopper_S$SLA)]))



##### succulence addition
local_grasshopper_A <- filter(grasshopper_no2021, Herbivore == "Addition")
local_grasshopper_A <- filter(local_grasshopper_A, elev_dist >= -10, elev_dist <= 13)

#Succulence 5.649086, 0.8042984
mean(local_grasshopper_A$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_A$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_A$succulence[!is.na(local_grasshopper_A$succulence)]))*1000

##### succulence removal
local_grasshopper_R <- filter(grasshopper_no2021, Herbivore == "Removal")
local_grasshopper_R <- filter(local_grasshopper_R, elev_dist >= -10, elev_dist <= 13)

#Succulence 5.885521, 0.654531
mean(local_grasshopper_R$succulence,na.rm = TRUE)*1000
sd(local_grasshopper_R$succulence,na.rm = TRUE) / sqrt(length(local_grasshopper_R$succulence[!is.na(local_grasshopper_R$succulence)]))*1000




m <- mean(traitdatRepro $avg_LAR, na.rm = TRUE)
sd <- sd(traitdatRepro $avg_LAR, na.rm = TRUE)

unscaled_x <- (-2 * sd) + m 

SLA_figure <- ggarrange(
  sla_cline,
  SLA_treatment,
  SLA_reproduction_water,
  SLA_fecundity,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2)

FDF_figure <- ggarrange(
  phen_cline,
  FT_treatment,
  FT_fecundity,
  labels = c("A", "B", "C"),
  ncol = 3, nrow = 1)


Duration_figure <- ggarrange(
  Duration_cline,
  Duration_treatment,
  duration_fecundity,
  labels = c("A", "B", "C"),
  ncol = 3, nrow = 1)


Succulence_figure<- ggarrange(
Succulence_cline,
Succulence_treatment_year,
Succulence_treatment_water,
succulence_fecundity,
Succulence_reproduction_herbivore,
labels = c("A", "B", "C", "D","E"),
ncol = 2, nrow = 3)

height_figure <- ggarrange(
  height_cline,
  Height_treatment,
  height_fecundity,
  labels = c("A", "B", "C"),
  ncol = 3, nrow = 1)


#testing

foliar_LAR<-drop_na(foliar,avg_LAR)
plot(foliar_LAR$SLA~foliar_LAR$avg_LAR)



