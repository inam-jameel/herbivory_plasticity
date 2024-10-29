## checking LAR model random effects


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
library(MuMIn) #For model selection


##this is where you specify the folder where you have the data on your computer

##setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/Grasshopper_manuscript_files/Grasshopper_manuscript_Submission/scripts_data/")

setwd("~/Documents/personnel/Jameel/grasshopper")


##read in data 
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

##This examines whether initial size varies as a function of source elevation, Table S15
initial_size <- lmer(init.diam ~ Water*Herbivore* elev_km +(1|Genotype), data = grasshopper)
Anova(initial_size, type = "III") # 

simulationOutput <- simulateResiduals(fittedModel= initial_size, plot = T, re.form = NULL,allow.new.levels =T)

# checking random effects
#initial_size_geno <- lm(init.diam ~ Water*Herbivore* elev_km, data = grasshopper)
#anova(initial_size, initial_size_geno)


#*******************************************************************************
####   Hypothesis one: Context dependency of clines and plasticity    #####
#*******************************************************************************

#*******************************************************************************
#### Herbivore damage #######
#*******************************************************************************

##* Censuses 1 and 2 occurred for all years. Census 3 was for 2021 and 2022 only

##reformat datafile

LAR_data_long_form<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3"),
                                                  names_to='census',
                                                  values_to='LAR')

LAR_data_long_form <- dplyr::select(LAR_data_long_form, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,census, year)

LAR_data_long_form$census[LAR_data_long_form$census == "LAR_1"] <- "1"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_2"] <- "2"
LAR_data_long_form$census[LAR_data_long_form$census == "LAR_3"] <-"3"

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

##Conversation of standardized elevation back to raw units
#-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)


########################
########################
## glmmtmb
########################
########################

damage_model <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year + (1| census_year)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data_long_form, family=beta_family())

Anova(damage_model, type = "III")


simulationOutput <- simulateResiduals(fittedModel= damage_model, plot = T, re.form = NULL,allow.new.levels =T)


damage_model_no_plantID <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year 
                                   + (1| census_year)
                                   #+ (1| PlantID)
                                   + (1| Cage_Block)
                                   + (1| Genotype)
                                   , data=LAR_data_long_form, family=beta_family())


damage_model_no_genotype <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year 
                                    + (1| census_year)
                                    + (1| PlantID)
                                    + (1| Cage_Block)
                                    #+ (1| Genotype)
                                    , data=LAR_data_long_form, family=beta_family())

damage_model_no_cb <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year 
                              + (1| census_year)
                              + (1| PlantID)
                              # + (1| Cage_Block)
                              + (1| Genotype)
                              , data=LAR_data_long_form, family=beta_family())

damage_model_no_census <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year 
                                  #+ (1| census_year)
                                  + (1| PlantID)
                                  + (1| Cage_Block)
                                  + (1| Genotype)
                                  , data=LAR_data_long_form, family=beta_family())


anova(damage_model,damage_model_no_plantID)
anova(damage_model,damage_model_no_genotype)
anova(damage_model,damage_model_no_cb)
anova(damage_model,damage_model_no_census)

damage_model_three_a <- glmmTMB(y_beta ~S_elev*Water*Herbivore+year + (1| census_year)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data_long_form, family=beta_family())

Anova(damage_model_three_a, type = "III")

damage_model_three_b <- glmmTMB(y_beta ~S_elev+Water*Herbivore*year + (1| census_year)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data_long_form, family=beta_family())

Anova(damage_model_three_b, type = "III")


damage_model_two <- glmmTMB(y_beta ~S_elev+Water*Herbivore+year + (1| census_year)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data_long_form, family=beta_family())

Anova(damage_model_two, type = "III")


model.sel(damage_model,damage_model_three_a,damage_model_three_b,damage_model_two) # four way model is best, delta is 137.56

damage_cline<-visregList(
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2021"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2021"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2022"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2022"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2023"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2023"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

#damage_cline<-visregList(visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_glmmTMB_full,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

damfit<-data.frame(damage_cline$fit)
damfit$visregFITexp<-exp(damfit$visregFit)

damage_clinal_variation<-ggplot(damfit, aes(S_elev, visregFITexp,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= LAR_data_long_form, aes(S_elev, LAR_prop, color= Water, shape=Water), alpha=0.50, color="black")+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.557947,-0.1085093, 1.340929))+ scale_y_continuous("Leaf area removed by herbivores (proportion)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+ geom_line(aes(group=Water),linewidth=1)+scale_linetype_manual(values=c("solid", "dotted"))+facet_grid(Herbivore~year)

damage_clinal_variation #Fig 2A









########################
########################
## GAMLSS
########################
########################


#fourway_model

LAR_Model_four<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year
                         + random(census_year) 
                         + random(PlantID)
                         + random(Cage_Block)
                         +random(Genotype)
                         ,family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc =500))

plot(LAR_Model_four)
summary(LAR_Model_four)
drop1(LAR_Model_four)
save(LAR_Model_four, file='LAR_Model.rda')   


LAR_Model_four_no_census<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year
                         #+ random(census_year) 
                         + random(PlantID)
                         + random(Cage_Block)
                         +random(Genotype)
                         ,family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))


LAR_Model_four_no_plantID<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year
                         + random(census_year) 
                         #+ random(PlantID)
                         + random(Cage_Block)
                         +random(Genotype)
                         ,family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

LAR_Model_four_no_geno<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year
                         + random(census_year) 
                         + random(PlantID)
                         + random(Cage_Block)
                         #+random(Genotype)
                         ,family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

LAR_Model_four_no_cb<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year
                         + random(census_year) 
                         + random(PlantID)
                         #+ random(Cage_Block)
                         +random(Genotype)
                         ,family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))












##Threeway_model
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore+year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

plot(LAR_Model_three)
summary(LAR_Model_three)
drop1(LAR_Model_three)
#save(LAR_Model_three, file='LAR_Model_three.rda')  
load(LAR_Model_three)

##Alternative three=way_model
LAR_Model_threeB<- gamlss (formula= y_beta ~year*Water*Herbivore+ S_elev + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

plot(LAR_Model_threeB)
summary(LAR_Model_threeB)
drop1(LAR_Model_threeB)
#save(LAR_Model_threeB, file='LAR_Model_threeB.rda')  


LAR_Model_twob<- gamlss (formula= y_beta ~S_elev*Water
                         + Water* Herbivore
                         + year*Water
                         + year*Herbivore
                         + S_elev
                         + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_twob)





##Two-way model
LAR_Model_two<- gamlss (formula= y_beta ~year+Water*Herbivore+ S_elev + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))

plot(LAR_Model_two)
summary(LAR_Model_two)
drop1(LAR_Model_two)
#save(LAR_Model_two, file='LAR_Model_two.rda')  




##Compare the three- and four-way models. The four-way model has a lower AICc value and is therefore preferred.
AICc(LAR_Model_threeB) #-26355.79
AICc(LAR_Model_three) #-26355.79
AICc(LAR_Model_four) #-26515.16


##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
#(2600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3100-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)

## creating plot for the cline using visreg, takes a moment to run
damage_cline<-visregList(
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2021"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2021"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2022"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2022"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2023"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2023"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

#damage_cline<-visregList(visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

damfit<-data.frame(damage_cline$fit)
damfit$visregFITexp<-exp(damfit$visregFit)

damage_clinal_variation<-ggplot(damfit, aes(S_elev, visregFITexp,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= LAR_data_long_form, aes(S_elev, LAR_prop, color= Water, shape=Water), alpha=0.50, color="black")+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.557947,-0.1085093, 1.340929))+ scale_y_continuous("Leaf area removed by herbivores (proportion)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+ geom_line(aes(group=Water),linewidth=1)+scale_linetype_manual(values=c("solid", "dotted"))+facet_grid(Herbivore~year)

damage_clinal_variation #Fig 2A

##Box_plot
LAR_box <-ggplot(LAR_data_long_form, aes(x = Herbivore, y = LAR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = 1,position = position_jitterdodge(0.3))

LAR_box <-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~year)
LAR_box


## grid format
LAR_box <-ggplot(LAR_data_long_form, aes(x = Water, y = LAR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous(element_blank()) +geom_point(aes(shape=factor(Water)), size = 1.5,position = position_jitterdodge(0.3))

LAR_box <-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(Herbivore~year)
LAR_box




##env concatenates water and herbivore. This model allows us to extract means and slopes for each combination of treatment and elevation
LAR_data_long_form $env<-interaction(LAR_data_long_form $Water, LAR_data_long_form $Herbivore)
LAR_data_long_form $env_year<-interaction(LAR_data_long_form $env, LAR_data_long_form $year)


#LAR_Model_three<- gamlss (formula= y_beta ~S_elev* env+ year-1 +random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
#summary(LAR_Model_three)
#save(LAR_ModelE, file='LAR_ModelE.rda')   


LAR_ModelE_four<- gamlss (formula= y_beta ~S_elev* env_year-1 +random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
summary(LAR_ModelE_four)




##Subsets of models for drop 1
LAR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore
                          +S_elev*Water*year
                          +S_elev* Herbivore*year
                          + Water* Herbivore*year
                          + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_three)
#save(LAR_Model_three, file='LAR_Model_three.rda')   

##Subsets of models for drop 1
LAR_Model_two<- gamlss (formula= y_beta ~S_elev*Water
                        + S_elev*Herbivore
                        + Water* Herbivore
                        + year*Water
                        + year*Herbivore
                        + S_elev*year
                        + year
                        + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)
#save(LAR_Model_two, file='LAR_Model_two.rda')  


LAR_Model_one<- gamlss (formula= y_beta ~S_elev+Water+Herbivore+ year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)
#save(LAR_Model_one, file='LAR_Model_one.rda')







