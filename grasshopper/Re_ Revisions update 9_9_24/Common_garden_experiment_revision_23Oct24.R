######## PROJECT: Common garden experiment: Examining clines and plasticity in response to water availability and grasshopper abundance
#### PURPOSE:Examine clines, trait values and fitness in response to water availability and grasshopper abundance in the field.
#### DATE LAST MODIFIED: 31 Oct 2024

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
library(lme4) # for running models



##this is where you specify the folder where you have the data on your computer


#setwd("~/Downloads/")

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/Grasshopper_manuscript_files/Grasshopper_manuscript_Submission/scripts_data/")


#setwd("~/Documents/personnel/Jameel/grasshopper")


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

##Create a dataframe for reproductive phenology traits which filters out the two plants that flowered without vernalization in the first season.
grasshopperFT <- filter(grasshopper, year != "2021")
##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(grasshopperFT, Ordinal_Date_flowering!="NA")

##set colors
cols=c("#CC79A7","lightblue") #for water treatment
cols2=c("#882255","#77DDAA") #for the grasshopper treatment

##This examines whether initial size varies as a function of source elevation, Table S1c
#initial_size <- lmer(init.diam ~ Water*Herbivore* elev_km +(1|Genotype), data = grasshopper)
#Anova(initial_size, type = "III")

## check residuals
#simulationOutput <- simulateResiduals(fittedModel= initial_size, plot = T, re.form = NULL,allow.new.levels =T)

## checking random effects
#initial_size_geno <- lm(init.diam ~ Water*Herbivore* elev_km, data = grasshopper)
#anova(initial_size, initial_size_geno)


#*******************************************************************************
#### Treatment Water effect ####
#*******************************************************************************

##read in data, in a separate file from the plant data
VWC <- read.csv("Common_garden_experiment_VWC.csv",stringsAsFactors=T) 
VWC $Year <-as.factor(VWC $Year)


## reformat datafile to analyze as repeated measures
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

## check to see that census and year are read as factors
VWC_data $census <-as.factor(VWC_data $census)
VWC_data $Year <-as.factor(VWC_data $Year)


## running model to analyze volumetic water content across water availability, herbivore abundance, and year
VWC_mod <- glmmTMB(VWC ~Water*Herbivore*Year+(1|census) + (1|Cage_Block), data= VWC_data,family=Gamma(link="log"))
Anova(VWC_mod,type="III")


##check residuals
#simulationOutput <- simulateResiduals(fittedModel= VWC_mod, plot = T, re.form = NULL,allow.new.levels =TRUE)

## Pairwise comparisons across waterng treatment and year
#VWC_means<-emmeans(VWC_mod, ~ Water*Year, type="response", adjust = "sidak")
#cld(VWC_means, details=TRUE)

## Pairwise comparisons across waterng treatment
#VWC_means<-emmeans(VWC_mod, ~ Water, type="response", adjust = "sidak")
#cld(VWC_means, details=TRUE)


##testing the random effects
#VWC_mod_no_census <- glmmTMB(VWC ~Water*Year*Herbivore+(1|census), data= VWC_data,family=Gamma(link="log"))

#VWC_mod_no_cageblock <- glmmTMB(VWC ~Water*Year*Herbivore + (1|Cage_Block), data= VWC_data,family=Gamma(link="log"))

#anova(VWC_mod,VWC_mod_no_cageblock)
#anova(VWC_mod,VWC_mod_no_census)

##Box plot
VWC_box <-ggplot(VWC_data, aes(x = Water, y = VWC, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous("Volumetric Water Content (percent)") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

VWC_plot <-VWC_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                             axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                             panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Year)

VWC_plot # Fig S3

#*******************************************************************************
####   Hypothesis one: Context dependency of clines and plasticity    #####
#*******************************************************************************

#*******************************************************************************
#### Herbivore resistance, analyzed as leaf area removed by herbivores #######
#*******************************************************************************

##* Censuses 1 and 2 occurred for all years. Census 3 was for 2021 and 2022 only

##reformat datafile for repeated measures
LAR_data_long_form<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3"),
                                                  names_to='census',
                                                  values_to='LAR')

LAR_data_long_form <- dplyr::select(LAR_data_long_form, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,census, year)

LAR_data_long_form$census[LAR_data_long_form$census == "LAR_1"] <- "1"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_2"] <- "2"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_3"] <-"3"

## check to see that census and year are read as factors
LAR_data_long_form $census <-as.factor(LAR_data_long_form $census)
LAR_data_long_form $year <-as.factor(LAR_data_long_form $year)

##Let's concatenate census and year to optimize terms in the model
LAR_data_long_form $census_year<-interaction(LAR_data_long_form$census, LAR_data_long_form$year,sep = "_")

## Convert leaf area removed by herbivores into proportion
LAR_data_long_form$LAR_prop<-LAR_data_long_form $LAR/100
hist(LAR_data_long_form$LAR_prop) #check to see if converstion worked


## model analyzing herbivore resistance across source elevation (cline), treatment levels (plasticity), and year.Need to run to generate Figure 2A, damage_clinal_variation
## Takes a moment to run

LAR_Model <- glmmTMB(LAR_prop ~Water*Herbivore*year+S_elev*Herbivore*year+S_elev*Water*year+S_elev*Water*Herbivore + (1| census_year)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), zi=~1, data=LAR_data_long_form, family=beta_family(link="logit"))

Anova(LAR_Model)

## Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= LAR_Model, plot = T, re.form = NULL,allow.new.levels =T)


## pairwise comparisons of the water availability treatment
#LAR <-emmeans(LAR_Model, ~ Water, type="response", adjust = "sidak")
#cld(LAR, details=TRUE)

## Slope of the relationship between source elevation, year, and herbivore abundance. The beta values need to be exponentiated
#emtrends(LAR_Model, specs = c("Herbivore","year"), var = "S_elev")

##Set tick marks on the X axis for these elevations: 2600, 3100, 3600. these breaks are in the code to generate damage_cline_variation
#(2600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3100-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)

##Conversation of standardized elevation back to raw units
#-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)

## creating plot for the cline using visreg, takes a moment to run. it will generate panels for all year by treatment levels
damage_cline<-visregList(
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2021"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2021"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2022"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2022"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2023"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2023"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

damfit<-data.frame(damage_cline$fit)
damfit$visregFITexp<-exp(damfit$visregFit)

damage_clinal_variation<-ggplot(damfit, aes(S_elev, visregFITexp,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= LAR_data_long_form, aes(S_elev, LAR_prop, color= Water, shape=Water), alpha=0.50, color="black")+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.557947,-0.1085093, 1.340929))+ scale_y_continuous("Leaf area removed by herbivores (proportion)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+ geom_line(aes(group=Water),linewidth=1)+scale_linetype_manual(values=c("solid", "dotted"))+facet_grid(Herbivore~year)

damage_clinal_variation #Fig 2A


##Box_plot,  significant effect of water
LAR_box <-ggplot(LAR_data_long_form, aes(x = Water, y = LAR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

LAR_treatment <-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)
LAR_treatment #Fig 2B


#*******************************************************************************
#### Specific Leaf Area #####
#*******************************************************************************

## model analyzing specific leaf area across source elevation (cline), treatment levels (plasticity), and year.Need to run to generate Figure 3A, SLA_clinal_variation
sla_model <- glmmTMB(SLA ~ Water*Herbivore*year+S_elev*Herbivore*year+S_elev*Water*year+S_elev*Water*Herbivore+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))

Anova(sla_model, type = "III")


## Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)


## pairwise comparisons of the water availability treatment
#SLA <-emmeans(sla_model, ~ Water, type="response", adjust = "sidak")
#cld(SLA, details=TRUE)

## Slope of the relationship between source elevation and SLA. The beta values need to be exponentiated
#emtrends(sla_model, specs = c("year"), var = "S_elev")

##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
#(2600-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)
#(3100-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)
#(3600-mean(foliar $elevation, na.rm = TRUE))/sd( foliar $elevation,na.rm = TRUE)

## cline for specific leaf area, only the cline in 2022 is significant
SLA_clinal<-visreg(sla_model,"S_elev", by="year", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
SLA_clinal_variation<-ggplot(SLA_clinal $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= foliar, aes(S_elev, SLA), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("SLA (cm2/g)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))+facet_wrap(~year)
SLA_clinal_variation

## Box_plot, y axis starts at 0

SLA_box <-ggplot(foliar, aes(x = Water, y = SLA, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous(limits = c(0,500),"SLA cm2/g") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

SLA_treatment <-SLA_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                  panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))
SLA_treatment


#*******************************************************************************
#### Succulence #####
#*******************************************************************************

succulence_data <- dplyr::select(foliar, succulence, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, Cage_Block, elev_km, S_elev,  year, elev_dist_km) #create separate data frame with relevant data columns for beta transformation

succulence_data <- drop_na(succulence_data,succulence) #remove missing data

#beta transformation
n<-nrow(succulence_data)
succulence_data $y_beta<- (succulence_data $succulence*(n-1) + 0.5)/n

## model analyzing leaf succulence across source elevation (cline), treatment levels (plasticity), and year.Need to run to generate Figure 3E, succulence_clinal_variation
succ_model <- glmmTMB(y_beta ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=beta_family())
Anova(succ_model,type="III")

succ_model <- glmmTMB(y_beta ~Water*Herbivore*year+S_elev*Herbivore*year+S_elev*Water*year+S_elev*Water*Herbivore+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=beta_family(link="logit"))

Anova(succ_model, type = "III")

##Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= succ_model, plot = T, re.form = NULL,allow.new.levels =T)

## pariwise comparisons across herbivore abundance and year
#succulence <-emmeans(succ_model, ~ Herbivore:year, type="response", adjust = "sidak")
#cld(succulence, details=TRUE)

## cline for succulence (not significant)
succ_clinal<-visreg(succ_model,"S_elev", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

Succ_clinal_variation<-ggplot(succ_clinal $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= succulence_data, aes(S_elev, y_beta), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("Leaf succulence (g/cm2)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))

Succ_clinal_variation #Fig 3E


##Box plot for herbivore by year interaction

Suc_box_herb <-ggplot(succulence_data, aes(x = Herbivore, y = y_beta, fill = Herbivore,shape=Herbivore)) +geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf succulence (g/cm2)") +geom_point(aes(shape=factor(Herbivore)), size = 2,position = position_jitterdodge(0.3))

Suc_treatment_herb <-Suc_box_herb + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Addition", "Removal")) +  scale_fill_manual(values = c(cols2), name = "Water treatment", labels = c("Herbivore Addition","Herbivore removal"))+scale_shape_manual(values=c(22,25), name = "Herbivore treatment", labels = c("Herbivore Addition","Herbivore removal"))+facet_wrap(~year)



Suc_treatment_herb #Fig 3F


#*******************************************************************************
#### Flowering phenology #####
#*******************************************************************************

## model analyzing flowering phenology as a function of source elevation (cline) and the interaction between treatment levels (plasticity) and year.

FT_model <-glmmTMB(FT_Adj ~ Water*Herbivore*year+S_elev*Herbivore*year+S_elev*Water*year+S_elev*Water*Herbivore+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(FT_model, type = "III") # 

##check residuals
simulationOutput <- simulateResiduals(fittedModel= FT_model, plot = T, re.form = NULL,allow.new.levels =TRUE)


##get slopes
#FT_table<- as.data.frame(summary(coefficients_FT))[c("Water","Herbivore",'S_elev.trend', 'SE')]
#FT_table <- FT_table%>% mutate(
#  slopes = exp(S_elev.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## pariwise comparisons across herbivore treatment levels
#FT <-emmeans(FT_model, ~ Herbivore*year, type="response", adjust = "sidak")
#cld(FT, details=TRUE)


##test random effects
#FT_model_no_geno<-glmmTMB(FT_Adj ~ S_elev+Water*Herbivore+year+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
#FT_model_no_cageblock<-glmmTMB(FT_Adj ~ S_elev+Water*Herbivore+year+(1|Genotype)+(1| PlantID),data= flowering,family=lognormal(link="log"))
#FT_model_noplantID<-glmmTMB(FT_Adj ~ S_elev+Water*Herbivore+year+(1|Genotype)+(1|Cage_Block),data= flowering,family=lognormal(link="log"))

#anova(FT_model,FT_model_no_geno)
#anova(FT_model,FT_model_no_cageblock)
#anova(FT_model,FT_model_noplantID)

##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
#(2600-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)
#(3100-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)
#(3600-mean(flowering $elevation, na.rm = TRUE))/sd( flowering $elevation,na.rm = TRUE)


##flowering phenology cline

flowering_time_cline<-visregList(visreg(FT_model,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                                 visreg(FT_model,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
phen_cline <-ggplot(flowering_time_cline $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= flowering, aes(S_elev, FT_Adj, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.630394,-0.1998366,  1.23072))+ scale_y_continuous("Flowering phenology (ordinal day of year)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)


phen_cline #Fig 4A

##Box_plot,  significant effect of herbivore abundance treatment by year


FT_box_herb <-ggplot(flowering, aes(x = Herbivore, y = FT_Adj, fill = Herbivore,shape=Herbivore)) +geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Flowering phenology (g)") +geom_point(aes(shape=factor(Herbivore)), size = 2,position = position_jitterdodge(0.3))

FT_treatment_herb <-FT_box_herb + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Addition", "Removal")) +  scale_fill_manual(values = c(cols2), name = "Water treatment", labels = c("Herbivore Addition","Herbivore removal"))+scale_shape_manual(values=c(22,25), name = "Herbivore treatment", labels = c("Herbivore Addition","Herbivore removal"))+facet_wrap(~year)



FT_treatment_herb #Fig 4B

#*******************************************************************************
#### Flowering duration #####
#*******************************************************************************

## model analyzing flowering duration as a funciton of source elevation (cline) and the interaction between treatment levels (plasticity) and year. Need to run to generate Figure 4D, Duration_cline
flowering_duration <-glmmTMB(flowering_duration ~ Water*Herbivore+S_elev*Herbivore+S_elev*Water+S_elev*year+Water*year+Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

Anova(flowering_duration, type = "III")

##Use the DHARMa package to examine the residuals, which are reasonable
#simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)

## pariwise comparisons across treatment levels
#Duration <-emmeans(flowering_duration, ~ Herbivore*Water, type="response", adjust = "sidak")
#cld(Duration, details=TRUE)

## pariwise comparisons across water treatment levels and year
#Duration_2 <-emmeans(flowering_duration, ~ year*Water, type="response", adjust = "sidak")
#cld(Duration_2, details=TRUE)


## cline plot, not significant
floweringduration_cline<-visreg(flowering_duration,"S_elev", partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

Duration_cline<-ggplot(floweringduration_cline $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= foliar, aes(S_elev, flowering_duration), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("Flowering duration (days)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))#+facet_wrap(~year)


Duration_cline #Fig 4D


##Box_plot for treatment effect
Duration_box <-ggplot(flowering, aes(x = Water, y = flowering_duration, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous("Flowering duration") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

Duration_treatment <-Duration_box_year + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)

Duration_treatment #Fig 4E



#*******************************************************************************
#### Height at flowering #####
#*******************************************************************************

## model analyzing height at flowering as a funciton of source elevation (cline) and the interaction between treatment levels (plasticity) and year. Need to run to generate Figure 4G, height_cline
max_height <-glmmTMB(Max_height_flowering ~ Water*Herbivore+S_elev*Herbivore+S_elev*Water+S_elev*year+Water*year+Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

Anova(max_height, type = "III") # # 

##Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= max_height, plot = T, re.form = NULL)


##get slopes that need to be exponentiated

#H_table<- as.data.frame(emtrends(max_height, specs = c("S_elev"), var = "S_elev"))
#H_table <- H_table%>% mutate(
#  slopes = exp(S_elev.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

##cline plot, significant cline
height_cline<-visreg(max_height,"S_elev", partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

height_clinal_variation<-ggplot(height_cline $fit, aes(S_elev, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= foliar, aes(S_elev, Max_height_flowering), alpha=0.5, color="black")+scale_linetype_manual(values=c("dashed","solid"))+scale_y_continuous("Tallest stem at flowering (cm)")+scale_x_continuous("Source elevation (m)",breaks=c(-1.563474, -0.1282131,1.307048))#+facet_wrap(~year)

height_clinal_variation #Figure 4G

##Box_plot no significant plasticity
height_box <-ggplot(flowering, aes(x = Water, y = Max_height_flowering, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ scale_y_continuous("Tallest stem at flowering (cm)") +geom_point(aes(shape=factor(Water)), size = 2,position = position_jitterdodge(0.3))

height_treatment <-height_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                                        axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                        panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Restricted", "Supplemental")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~Herbivore)
height_treatment #Figure 4H




#*******************************************************************************
####   Hypothesis Two: Selection via two fitness components on traits   #####
#*******************************************************************************

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")
##retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper_no2021,FT_Adj,Max_height_flowering, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block, elevation, year, Mature_length_siliques,Reproduced,flowering_duration, SLA, succulence, FT_Adj)



##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
traitdat $S_initdiam<-scale(traitdat $init.diam,center=TRUE, scale=TRUE)

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing trait values 

traitdat $sLAR<-scale(traitdat $avg_LAR,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)

head(traitdat)


##Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Addition","Removal"))


#*******************************************************************************
#### Selection via probability of reproduction for vegetative traits only  #####
#*******************************************************************************

repro_model <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
                        Water*Herbivore*sLAR+ 
                        Water*Herbivore*sSUC+ I(sSUC^2)+
                        Water*Herbivore*sSLA+ I(sSLA^2) +
                        +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model,type="III") 

#summary(repro_model)

##To check residuals
#simulationOutput <- simulateResiduals(fittedModel= repro_model, plot = T, re.form = NULL)


##Conversation of standardized traits back to raw units
#-1*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
#0*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
#1*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
#2*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)
#3*sd(traitdat $SLA,na.rm=TRUE)+mean(traitdat $SLA,na.rm=TRUE)

##plots for selection anaylses

##selection on Specific leaf area 
sla_repro = visreg(repro_model,"sSLA", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
sla_repro_plot<-ggplot(sla_repro $fit, aes(sSLA, visregFit)) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdat, aes(sSLA, Reproduced),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Specific leaf area (cm2/g)")+ scale_y_continuous("Probability of reproduction")

sla_repro_plot #figure 3C

##selection on leaf succulence
suc_repro = visreg(repro_model,"sSUC", by="Herbivore",overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")
suc_repro_plot<-ggplot(suc_repro$fit, aes(sSUC, visregFit,group= Herbivore, colour= Herbivore, fill=factor(Herbivore))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Herbivore)) +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdat, aes(sSUC, Reproduced, color= Herbivore, shape=Herbivore),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_shape_manual(values=c(22,25))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf succulence (g/cm2)")+ scale_y_continuous("Probability of reproduction")+scale_fill_manual(values = cols2, name = "Grasshopper treatment", labels = c("Addition","Removal"))+scale_colour_manual(values = cols2, name = "Grasshopper treatment", labels = c("Addition","Removal"))

suc_repro_plot #figure 3G

##Conversation of standardized trait values back to raw units
#-1*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
#0*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
#1*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
#2*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)
#3*sd(traitdat $succulence,na.rm=TRUE)+mean(traitdat $succulence,na.rm=TRUE)


## To extract coefficients, we need to create the quadratic effects outside of the model
#traitdat$sSLA2<-traitdat$sSLA*traitdat$sSLA
#traitdat$sLAR2<-traitdat$sLAR*traitdat$sLAR
#traitdat$sSUC2<-traitdat$sSUC*traitdat$sSUC


#repro_model_quad <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
#                             Water*Herbivore*sLAR+ 
#                             Water*Herbivore*sSUC+ sSUC2 +
#                             Water*Herbivore*sSLA+ sSLA2 
#                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
#Anova(repro_model_quad,type="III") 

##Slope of selection surfaces 
#emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSLA")
#emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSLA2")
#emtrends(repro_model_quad, specs = c("sSLA2"), var = "sSLA2")
#emtrends(repro_model_quad, specs = c("sSLA"), var = "sSLA")

#emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSUC")
#emtrends(repro_model_quad, specs = c("Water","Herbivore"), var = "sSUC2")

#emtrends(repro_model_quad, specs = c("Herbivore"), var = "sSUC")
#emtrends(repro_model_quad, specs = c("Herbivore"), var = "sSUC2")


####  Univariate models  #####

#repro_model_SLA <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
#                             Water*Herbivore*sSLA+ I(sSLA^2) 
#                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
#Anova(repro_model_SLA,type="III") 

#sla_pred_uni <- ggpredict(repro_model_SLA, terms = c("sSLA[all]"), type = "re", interval="confidence")
#SLA_reproduction_uni <- plot(sla_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=FALSE,  facet=FALSE,dot_alpha=0.75)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Probability of reproduction")
#SLA_reproduction_uni

#repro_model_SUC <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
#                             Water*Herbivore*sSUC+I(sSUC^2)+ 
#                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
#Anova(repro_model_SUC,type="III") 
#sSUC_pred_uni <- ggpredict(repro_model_SUC, terms = c("sSUC[all]","Herbivore"), type = "re", interval="confidence")

#Succulence_reproduction_herbivore_uni <-plot(sSUC_pred_uni, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols2,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")
#Succulence_reproduction_herbivore_uni


#repro_model_LAR <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore*year+
#                             Water*Herbivore*sLAR+ I(sLAR ^2)+
#                            +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))                            
#Anova(repro_model_LAR,type="III") 

#repro_model_LAR_optim <- update(repro_model_LAR,
#                control=glmmTMBControl(optimizer=optim,
#                                       optArgs=list(method="BFGS")))
#Anova(repro_model_LAR_optim,type="III") 
#emtrends(repro_model_LAR_optim, specs = c("Herbivore","Water"), var = "I(sLAR^2)")


#lar_univariate <- ggpredict(repro_model_LAR, terms = c("sLAR[all]"), type = "re", interval="confidence")
#LAR_uni <- plot(lar_univariate, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE,dot_alpha=0.75, jitter=0.025)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf area removed")+ scale_y_continuous("Probability of reproduction")
#LAR_uni

#*******************************************************************************
####  Selection via Seed set using all traits #####
#*******************************************************************************

##Create datafile with only reproductive plants
traitdatRepro <- filter(traitdat, Reproduced == 1 )

##Rescale after removing the non-reproductive plants
traitdatRepro $sduration<-scale(traitdatRepro $flowering_duration,center=TRUE,scale=TRUE)
traitdatRepro $s_max_height<-scale(traitdatRepro $Max_height_flowering,center=TRUE,scale=TRUE)
traitdatRepro $sLAR<-scale(traitdatRepro $avg_LAR,center=TRUE,scale=TRUE)
traitdatRepro $sSLA<-scale(traitdatRepro $SLA,center=TRUE,scale=TRUE)
traitdatRepro $sSUC<-scale(traitdatRepro $succulence,center=TRUE,scale=TRUE)
traitdatRepro $sFT<-scale(traitdatRepro $FT_Adj,center=TRUE,scale=TRUE)
traitdatRepro $S_initdiam <-scale(traitdatRepro $init.diam,center=TRUE,scale=TRUE)


fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore+year +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sLAR+Water*Herbivore*I(sLAR^2)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits,type="III") 



##To check residuals
#simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)

##To extract coefficients, we need to create the quadratic effects outside of the model
#traitdatRepro $sSLA2<-traitdatRepro $sSLA* traitdatRepro $sSLA
#traitdatRepro $sLAR2<-traitdatRepro $sLAR* traitdatRepro $sLAR
#traitdatRepro $sSUC2<-traitdatRepro $sSUC* traitdatRepro $sSUC
#traitdatRepro $sduration2<-traitdatRepro $sduration* traitdatRepro $sduration

#fecund_modeltraits_quad <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore+year +
#                                    Water*Herbivore*sFT+
#                                    Water*Herbivore*s_max_height +
#                                    Water*Herbivore* sduration +Water*sduration2 +Herbivore*sduration2 +
#                                    Water*Herbivore* sSLA+ + Water*Herbivore*sSLA2 +
#                                    Water*Herbivore* sSUC+ 
#                                    Water*Herbivore* sLAR+Water*Herbivore* sLAR2 +(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
#Anova(fecund_modeltraits_quad,type="III")

## extract selection gradients for flowering time
#emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sFT")
#FT
#coefficients_FT <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sFT",type="response")
#FT_table<- as.data.frame(summary(coefficients_FT))[c('sFT.trend', 'SE')]
#FT_table <- FT_table%>% mutate(
#  slopes = exp(sFT.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))


## extract selection gradients for flowering duration
#emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sduration")
#coefficients_dur <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sduration",type="response")
#dur_table<- as.data.frame(summary(coefficients_dur))[c('sduration.trend', 'SE')]
#dur_table <- dur_table%>% mutate(
#  slopes = exp(sduration.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract quadratic terms for flowering duration
#emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"), var = "sduration2")
#coefficients_dur <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sduration2",type="response")
#dur_table<- as.data.frame(summary(coefficients_dur))[c('sduration2.trend', 'SE')]
#dur_table <- dur_table%>% mutate(
#  slopes = exp(sduration2.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract selection gradients for height at flowering 
#coefficients_height <- emtrends(fecund_modeltraits_quad, specs = c("s_max_height"),var = "s_max_height",type="response")
#height_table<- as.data.frame(summary(coefficients_height))[c('s_max_height.trend', 'SE')]
#height_table <- height_table%>% mutate(
#  slopes = exp(s_max_height.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract selection gradients for succulence
#coefficients_succ <- emtrends(fecund_modeltraits_quad, specs = c("Herbivore"),var = "sSUC",type="response")
#succ_table<- as.data.frame(summary(coefficients_succ))[c("Herbivore",'sSUC.trend', 'SE')]
#succ_table <- succ_table%>% mutate(
#  slopes = exp(sSUC.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract selection gradients for specific leaf area
#coefficients_sla <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sSLA",type="response")
#sla_table<- as.data.frame(summary(coefficients_sla))[c("Water","Herbivore",'sSLA.trend', 'SE')]
#sla_table <- sla_table%>% mutate(
#  slopes = exp(sSLA.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract quadratic terms for specific leaf area
#coefficients_sla <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sSLA2",type="response")
#sla_table<- as.data.frame(summary(coefficients_sla))[c("Water","Herbivore",'sSLA2.trend', 'SE')]
#sla_table <- sla_table%>% mutate(
#  slopes = exp(sSLA2.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract selection gradients for herbivore resistance
#coefficients_LAR <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = "sLAR",type="response")
#LAR_table<- as.data.frame(summary(coefficients_LAR))[c("Water","Herbivore",'sLAR.trend', 'SE')]
#LAR_table <- LAR_table%>% mutate(
#  slopes = exp(sLAR.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))

## extract quadratic terms for herbivore resistance
#coefficients_LAR <- emtrends(fecund_modeltraits_quad, specs = c("Water","Herbivore"),var = #"sLAR2",type="response")
#LAR_table<- as.data.frame(summary(coefficients_LAR))[c("Water","Herbivore",'sLAR2.trend', 'SE')]
#LAR_table <- LAR_table%>% mutate(
#  slopes = exp(sLAR2.trend),
#  Lower95 = slopes * exp(-1.96*SE),
#  Upper95 = slopes * exp(1.96*SE))


##plots for selection anaylses

##SLA
##Context dependent selection
SLA<-visregList(visreg(fecund_modeltraits,"sSLA", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                visreg(fecund_modeltraits,"sSLA", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)
SLA_fecund_plot<-ggplot(SLA $fit, aes(sSLA, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sSLA, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("SLA (cm2/g)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
SLA_fecund_plot #figure 3D

##Conversation of standardized traits back to raw units

#-2*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
#0*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
#2*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)
#4*sd(traitdatRepro $SLA,na.rm=TRUE)+mean(traitdatRepro $SLA,na.rm=TRUE)



##Succulence
##Context dependent selection
SUC<-visregList(visreg(fecund_modeltraits,"sSUC", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                visreg(fecund_modeltraits,"sSUC", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


SUC_fecund_plot<-ggplot(SUC $fit, aes(sSUC, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sSUC, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf succulence (g/cm2)")+ scale_y_continuous("Fecundity (total fruit length)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
SUC_fecund_plot #figure 3H

##Conversation of standardized traits back to raw units

#-2*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
#-1*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
#0*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
#1*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
#2*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)
#3*sd(traitdatRepro $succulence,na.rm=TRUE)+mean(traitdatRepro $succulence,na.rm=TRUE)


## Leaf area removed
##Context dependent selection

LAR<-visregList(visreg(fecund_modeltraits,"sLAR", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                visreg(fecund_modeltraits,"sLAR", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


LAR_fecund_plot<-ggplot(LAR $fit, aes(sLAR, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+  geom_point(data= traitdatRepro, aes(sLAR, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Leaf damage from herbivores (proportion)")+ scale_y_continuous("Fecundity (total fruit length, mm)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
LAR_fecund_plot #figure 2C

##Conversation of standardized traits back to raw units

#-1*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
#0*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
#1*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
#2*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)
#3*sd(traitdatRepro $avg_LAR,na.rm=TRUE)+mean(traitdatRepro $avg_LAR,na.rm=TRUE)


##Flowering phenology
##no longer Context dependent selection

#FT
FT = visreg(fecund_modeltraits,"sFT", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

FT_fecund_plot<-ggplot(FT $fit, aes(sFT, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdatRepro, aes(sFT, Mature_length_siliques),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Flowering phenology (ordinal day of year")+ scale_y_continuous("Fecundity (total fruit length, mm)",limits=c(0,2000))
FT_fecund_plot #figure 4C


##Conversation of standardized traits back to raw units
#-2*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
#0*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
#2*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)
#4*sd(traitdatRepro $FT_Adj,na.rm=TRUE)+mean(traitdatRepro $FT_Adj,na.rm=TRUE)



#Duration
duration<-visregList(visreg(fecund_modeltraits,"sduration", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                     visreg(fecund_modeltraits,"sduration", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)

duration_fecund_plot<-ggplot(duration $fit, aes(sduration, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= traitdatRepro, aes(sduration, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Flowering duration (days)")+ scale_y_continuous("Fecundity (total fruit length, mm)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
duration_fecund_plot #figure 4F

##Conversation of standardized traits back to raw units

#-1*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
#0*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
#1*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)
#2*sd(traitdatRepro $flowering_duration,na.rm=TRUE)+mean(traitdatRepro $flowering_duration,na.rm=TRUE)



#height
height_fecund = visreg(fecund_modeltraits,"s_max_height", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

height_fecund_plot<-ggplot(height_fecund $fit, aes(s_max_height, visregFit)) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) + geom_line() +theme_bw()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_jitter(data= traitdatRepro, aes(s_max_height, Mature_length_siliques),width=0,height=0.025,size=2, alpha=0.75, color="black")+scale_x_continuous("Tallest stem at flowering (cm)")+ scale_y_continuous("Fecundity (total fruit length, mm)",limits=c(0,2000))
height_fecund_plot #figure 4I

##Conversation of standardized traits back to raw units

#0*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)
#2*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)
#4*sd(traitdatRepro $Max_height_flowering,na.rm=TRUE)+mean(traitdatRepro $Max_height_flowering,na.rm=TRUE)



### Univariate models####


###Flowering time###


#fecund_modeltraits_FT <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore+year +
#                             Water*Herbivore*sFT
#                     +(1|PlantID)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#Anova(fecund_modeltraits_FT,type="III") 

#visreg(fecund_modeltraits_FT,"sFT", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Timing of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
#visreg(fecund_modeltraits_FT,"sFT", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Timing of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")



#### Height ###


#fecund_modeltraits_height <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore+year +
#                               Water*Herbivore*s_max_height +I(s_max_height^2)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))

#Anova(fecund_modeltraits_height,type="III") 
#visreg(fecund_modeltraits_height,"s_max_height",  scale = "response", xlab="Height at flowering (m)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")



####Flowering duration###


#fecund_modeltraits_duration <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
#                               Water*Herbivore* sduration +Water*I(sduration^2)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
#Anova(fecund_modeltraits_duration,type="III") 
#visreg(fecund_modeltraits_duration,"sduration", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Duration of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
#visreg(fecund_modeltraits_duration,"sduration", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Duration of reproduction (ordinal day)", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")


#### SLA ###


#fecund_modeltraits_SLA <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
#                                Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
#                               (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
#Anova(fecund_modeltraits_SLA,type="III") 
#visreg(fecund_modeltraits_SLA,"sSLA", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Specific leaf area", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional",ylim=c(0,3000))
#visreg(fecund_modeltraits_SLA,"sSLA", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Specific leaf area", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional",ylim=c(0,3000))


#### Succulence ###


#fecund_modeltraits_SUC <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year + Water*Herbivore* sSUC+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
#Anova(fecund_modeltraits_SUC,type="III") 
#visreg(fecund_modeltraits_SUC,"sSUC", by="Water",cond=list("Herbivore"="Addition"), overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")
#visreg(fecund_modeltraits_SUC,"sSUC", by="Water",cond=list("Herbivore"="Removal"), overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=FALSE, type="conditional")

#visreg(fecund_modeltraits_SUC,"sSUC", by="Herbivore",overlay=TRUE, scale = "response", xlab="Leaf succulence", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional", ylim=c(0,3000))


#### LAR ###


#fecund_modeltraits_LAR <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore*year +
#                                Water*Herbivore* sLAR+I(sLAR^2)+ (1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
#Anova(fecund_modeltraits_LAR,type="III") 
#visreg(fecund_modeltraits_LAR,"sLAR", overlay=TRUE, scale = "response", xlab="Leaf damage from arthropod herbivores", ylab="Fecundity (total fruit length)", partial=TRUE, type="conditional", ylim=c(0,3000))




#*******************************************************************************
#### Hypothesis Three: Effect of water availability on fitness components and local adaptation #####
#*******************************************************************************

#### Local adaptation analysis using repeated measures ####

hurdle_Model_LA <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= grasshopperFT, 
                           zi=~S_initdiam+year
                           +Water*Herbivore*S_elev
                           +Water*Herbivore*I(S_elev^2)
                           +(1|PlantID)+ (1| Cage_Block)+(1|Genotype),family=ziGamma(link="log"))
RM_fit <- subset(grasshopperFT, Mature_length_siliques>0) # for plotting from count part, excludes plants that did not reproduce


##To check residuals
#simulationOutput <- simulateResiduals(fittedModel= hurdle_Model_LA, plot = T, re.form = NULL,allow.new.levels =TRUE)

##This is the ANOVA table for the logistic regression part (probability of reproduction).
#Anova(hurdle_Model_LA,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
#Anova(hurdle_Model_LA,type="III", component="cond")



##Extract the coefficients using emtrends, which requires defining the quadratic effect prior to running the model
#RM_fit $Selev2<-RM_fit $S_elev* RM_fit $S_elev

#hurdle_Model_LA_quad <- glmmTMB(Mature_length_siliques ~ S_initdiam+year
#                                +Water*Herbivore*S_elev
#                                +Water*Herbivore*Selev2
#                                +(1|PlantID)+(1|Cage_Block)+(1|Genotype), data= RM_fit,family=Gamma(link="log"))


#visreg(hurdle_Model_LA_quad,"S_elev", by="Water", overlay = TRUE, partial = TRUE,scale="response",ylim=c(0,2000))

#emtrends(hurdle_Model_LA_quad, specs = c("Water"), var = "S_elev",type="response", component="cond")
#emtrends(hurdle_Model_LA_quad, specs = c("Water"), var = "Selev2",type="response", component="cond")


##Optimal source elevations when we just consider the water treatment:  0.2640449, 3242.76 (95% CI: 3142.12)
#-  0.188 /(2*-0.356 ) #0.2640449
#0.0378/(2* -0.656) #-0.02881098
#-0.4141/(2*-0.0647) #3.200155


##Converting to elevation in m
#0.2640449*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  
#-0.02881098*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  
#3.200155*sd( RM_fit $elevation,na.rm = TRUE)+mean(RM_fit $elevation, na.rm = TRUE)  

##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
#(2600-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)
#(3100-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)
#(3600-mean(RM_fit $elevation, na.rm = TRUE))/sd( RM_fit $elevation,na.rm = TRUE)


## plotting fitness across water availability treatment

fecundity_cline_water<-visreg(hurdle_Model_LA,"S_elev", by="Water", overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response")

Local_adaptation<-ggplot(fecundity_cline_water $fit, aes(S_elev, visregFit,group= Water, colour= Water, fill=factor(Water))) + geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= RM_fit, aes(S_elev, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.564114,-0.1030306,1.358053))+ scale_y_continuous("Seedset (total fruit length, mm)",limits=c(0,2500))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Water)

##Figure 5
Local_adaptation



#### plotting main figures ####
## requires ggpubr to be loaded and code for the models and plots already run

## Figure 2

figure2 <- ggarrange(damage_clinal_variation, 
                     LAR_treatment,LAR_fecund_plot,
                     labels = c("A", "B","C"),
                     ncol = 3, nrow = 1, common.legend = TRUE, legend="none")
figure2

## Figure 3

figure3 <- ggarrange(
  SLA_clinal_variation, SLA_treatment, sla_repro_plot, SLA_fecund_plot,Succ_clinal_variation, Suc_treatment_herb, suc_repro_plot, SUC_fecund_plot , 
  labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  ncol = 4, nrow = 2)
figure3


## Figure 4

figure4 <- ggarrange(
  phen_cline,
  FT_treatment_herb,
  FT_fecund_plot,
  Duration_cline,
  Duration_treatment,
  duration_fecund_plot,
  height_clinal_variation,
  height_treatment,
  height_fecund_plot,
  labels = c("A", "B", "C","D", "E", "F","G", "H", "I"),
  ncol = 3, nrow = 3)
figure4


##Figure5

figure5 <- Local_adaptation

figure5

###


