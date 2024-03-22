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
cols2=c("#882255","lightgreen")


### Hypothesis one ######

#*******************************************************************************
#### Herbivore damage #####
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

ggplot(LAR_data, aes(x= LAR_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

LAR_data <- drop_na(LAR_data,LAR_prop) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

LAR_data $y_beta<- (LAR_data $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

min(LAR_data $y_beta)

max(LAR_data $y_beta)

## model
damage_model <- glmmTMB(y_beta ~S_elev*Water*Herbivore*year + (1| census)+ (1| PlantID)+ (1| Cage_Block)+ (1| Genotype), data=LAR_data, family=beta_family())

Anova(damage_model,type="III") 

simulationOutput <- simulateResiduals(fittedModel= damage_model, plot = T, re.form = NULL,allow.new.levels =TRUE)


## elevation plot
pred_dam <- ggpredict(damage_model, terms = c("S_elev[all]","Herbivore","year"), type = "re", interval="confidence")

damage_cline <-plot(pred_dam, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols2, facet=TRUE)+theme(text = element_text(size=10), axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Standardized source elevation (km)")+ scale_y_continuous("Leaf Area Removed") 

damage_cline

## treatment plot

damage_treatment <-ggplot(LAR_data, aes(x = Water, y = y_beta, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33,fill = "gray") +
  
  facet_grid(Herbivore~ year) +
  theme_bw() + ylim(0,.7)+
  scale_fill_manual(values = c("#CC79A7","lightblue")) +
  labs(y = "Damage by Herbivores (%)") + 
  labs(x = element_blank()) +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none")

damage_treatment

#*******************************************************************************
#### Specific Leaf Area #####
#*******************************************************************************

foliar<-subset(grasshopper,SLA>0)
foliar $S_elev<-scale(foliar $elevation,center=TRUE, scale=TRUE)

sla_model <- glmmTMB(SLA ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = foliar, family= lognormal(link="log"))
Anova(sla_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= sla_model, plot = T, re.form = NULL,allow.new.levels =T)

pred_sla <- ggpredict(sla_model, terms = c("S_elev[all]","year"), type = "re", interval="confidence")

sla_cline <-plot(pred_sla, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = c("black","black"), facet=TRUE)+theme(#text = element_text(size=20),
                                                                                                                              axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Specific leaf area")

sla_cline

SLA_treatment <- ggplot(SLA_df, aes(x = Water, y = SLA, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33,fill = "gray") +
  
  #facet_grid(Herbivore~ year) +
  theme_bw() + ylim(75,250)+
  scale_fill_manual(values = c("#CC79A7","lightblue")) +
  labs(y = "Specific leaf area") + 
  labs(x = element_blank()) +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

SLA_treatment

#*******************************************************************************
#### Succulence #####
#*******************************************************************************

foliar$succulence_mg<-foliar$succulence*1000
foliar$Suc_mg<- foliar$succulence_mg+0.01

succulence_data <- dplyr::select(foliar, succulence, elevation, Genotype, Cage, Water, Herbivore, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,  year, Suc_mg)

succulence_data <- drop_na(succulence_data,succulence) 

succ_model <- glmmTMB(Suc_mg ~ Water*Herbivore*S_elev*year+(1|PlantID)+(1|Genotype)+(1|Cage_Block), data = succulence_data, family=lognormal(link="log"))
Anova(succ_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= succ_model, plot = T, re.form = NULL,allow.new.levels =T)


pred_succ <- ggpredict(succ_model, terms = c("S_elev[all]","year"), type = "re", interval="confidence")

Succulence_cline <-plot(pred_succ, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = c("black","black"), facet=TRUE)+theme(
  #text = element_text(size=20),
  axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (standardized)")+ scale_y_continuous("Leaf succulence")

Succulence_cline


Succulence_treatment<- ggplot(succulence_data, aes(x = Water, y = Suc_mg, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = cols) +
  labs(y = "Leaf succulence") + 
  labs(x = "Watering Treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

Succulence_treatment #supplemental figure

Succulence_treatment_year<- ggplot(succulence_data, aes(x = Herbivore, y = Suc_mg, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = cols2) +
  labs(y = "Leaf succulence") + 
  labs(x = "Herbivore Treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

Succulence_treatment_year #supplemental figure

####

figureS4 <- ggarrange(
  Succulence_treatment,
  Succulence_treatment_year,
  
  
  labels = c("A", "B"),
  ncol = 2, nrow = 1)
figureS4


#*******************************************************************************
#### Flowering time #####
#*******************************************************************************

#filter out the two plants that flowered in the first season. This step is likely not necessary because those plants didn't have flowering times.
grasshopperFT <- filter(grasshopper, year != "2021")

##To enable simulate residuals to work, we have to exclude plants that did not flower
flowering<-subset(grasshopperFT, Ordinal_Date_flowering!="NA")

FT_model<-glmmTMB(FT_Adj ~ S_elev*Water*Herbivore*year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))
Anova(FT_model, type = "III") # 

pred_FT <- ggpredict(FT_model, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
phen_cline <-plot(pred_FT, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(#text = element_text(size=10),
  axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")

phen_cline #main figure, fig1


pred_FT2 <- ggpredict(FT_model, terms = c("S_elev[all]", "Water","year"), type = "re", interval="confidence")
phen_cline_year <-plot(pred_FT2, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering phenology (ordinal day of year)")+ geom_abline(intercept = 0, slope = 1, linetype = "dashed")


phen_cline_year #supplemental figure, S5


#FT_treatment<- ggplot(flowering, aes(x = Water, y = FT_Adj, group = Water)) +
#  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
#  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
#  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
#  facet_grid(~ Herbivore) +
#  theme_light() +
#  scale_fill_manual(values = cols) +
#  labs(y = "Flowering phenology (ordinal day of year)") + 
#  labs(x = element_blank()) +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
#  theme(legend.position = "none") 

#FT_treatment #not significant


#*******************************************************************************
#### Flowering duration #####
#*******************************************************************************


flowering_duration<-glmmTMB(flowering_duration ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering)

Anova(flowering_duration, type = "III") # elevation is significant, not year
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= flowering_duration, plot = T, re.form = NULL,allow.new.levels =TRUE)

pred_duration <- ggpredict(flowering_duration, terms = c("S_elev[all]", "Water","Herbivore"), type = "re", interval="confidence")
Duration_cline <-plot(pred_duration, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(#text = element_text(size=10),
  axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Flowering Duration (days)")


Duration_cline #no significant interaction


Duration_treatment<- ggplot(flowering, aes(x = Water, y = flowering_duration, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = cols) +
  labs(y = "Flowering duration") + 
  labs(x = element_blank()) +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 

Duration_treatment #not significant


#*******************************************************************************
#### Tallest bolt at flowering #####
#*******************************************************************************

max_height<-glmmTMB(Max_height_flowering ~ S_elev*Water*Herbivore+year+(1|Genotype)+(1|Cage_Block)+(1| PlantID),data= flowering,family=lognormal(link="log"))

Anova(max_height, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable
simulationOutput <- simulateResiduals(fittedModel= max_height, plot = T, re.form = NULL)

pred_height <- ggpredict(max_height, terms = c("S_elev[all]"), type = "re", interval="confidence")
height_cline <-plot(pred_height, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = "black", facet=TRUE)+theme(#text = element_text(size=10),
  axis.title.x=element_blank(),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="none")+scale_x_continuous("Source elevation (km)")+ scale_y_continuous("Height of tallest bolt at flowering")


height_cline #supplemental figure, S5


#Height_treatment<- ggplot(flowering, aes(x = Water, y = Max_height_flowering, group = Water)) +
#  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
#  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
#  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
#  facet_grid(~ Herbivore) +
#  theme_light() +
#  scale_fill_manual(values = cols) +
#  labs(y = "Height of tallest bolt at flowering") + 
#  labs(x = "Watering treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
#  theme(legend.position = "none") 

#Height_treatment #not significant


####

library(ggpubr)
figure1 <- ggarrange(damage_treatment,
                     damage_cline,
                     
                     SLA_treatment,
                     sla_cline,
                     
                     Succulence_treatment,
                     Succulence_cline,
                     
                     FT_treatment,
                     phen_cline,
                     
                     Duration_treatment,
                     Duration_cline,
                     
                     Height_treatment,
                     height_cline,
                     
                     labels = c("A", "B", "C", "D","E", "F", "G", "H",
                                "I", "J", "K", "L","M", "N", "O", "P"),
                     ncol = 2, nrow = 6)
figure1


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

repro_herb<- ggplot(grasshopperFT, aes(x = Herbivore, y = Reproduced, group = Herbivore)) +
  geom_violin(aes(fill = Herbivore), alpha = 0.95, trim = T) +
  #geom_jitter(shape = 21, fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = cols2) +
  labs(y = "Probability of reproduction") + 
  labs(x = "Herbivore treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none")

repro_herb


#*******************************************************************************
####  Fecundity #####
#*******************************************************************************
repro<-subset(grasshopperFT, Reproduced=="1")

repro $Water1 <-factor(repro $Water, levels = c("Supplemental", "Restricted"))
repro $Water <-factor(repro $Water, levels = c("Restricted","Supplemental"))


fecundity <- glmmTMB(Mature_length_siliques ~S_initdiam+Water*Herbivore*year* S_elev+(1|PlantID) + (1|Cage_Block)+(1|Genotype), data= repro,family=Gamma(link="log"),control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
Anova(fecundity,type="III")


summary(fecundity)

coefficients_fec <- emtrends(fecundity, specs = c("Water"), var = "S_elev",type="response")
fec_table<- as.data.frame(summary(coefficients_fec))[c('Water','S_elev.trend', 'SE')]
fec_table <- fec_table%>% mutate(
  slopes = exp(S_elev.trend),
  Lower95 = slopes * exp(-1.96*SE),
  Upper95 = slopes * exp(1.96*SE))

fec_table

Fecmeans<-emmeans(fecundity, ~ Water:Herbivore, type="response")
cld(Fecmeans, details=TRUE)

pred_fec <- ggpredict(fecundity, terms = c("S_elev[all]","Water"), type = "re", interval="confidence")
Local_adaptation <-plot(pred_fec, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=FALSE)+theme(#text = element_text(size=10),
  axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+
  #theme(legend.position="right")+
  theme(legend.position="none")+
  scale_x_continuous("Source elevation ")+ scale_y_continuous("Fecundity (length of mature siliques)")


Local_adaptation


fecundity_treatment <- ggplot(repro, aes(x = Water, y = Mature_length_siliques, group = Water)) +
  geom_violin(aes(fill = Water), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.1), fill = "white", size = .5, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ Herbivore) +
  theme_light() +
  scale_fill_manual(values = cols) +
  labs(y = "Fecundity (length of mature siliques)") + 
  labs(x = "Watering treatment") +
  #theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none")

fecundity_treatment


####

library(ggpubr)
figure2 <- ggarrange(
    repro_herb,
    fecundity_treatment,
    Local_adaptation,
    
                     
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1)
figure2



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
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Removal","Addition"))


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




sSUC_pred <- ggpredict(repro_model_full, terms = c("sSUC[all]", "Herbivore"), type = "re", interval="confidence")

Succulence_reproduction_herbivore <-plot(full_sSUC, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols2,facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Leaf succulence")+ scale_y_continuous("Probability of reproduction")


sla_pred <- ggpredict(repro_model_full, terms = c("sSLA[all]", "Water"), type = "re", interval="confidence")

SLA_reproduction_water <- plot(sla_pred, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous(element_blank())

####

library(ggpubr)
figure3 <- ggarrange(
  Succulence_reproduction_herbivore,
  SLA_reproduction_water,
  
  
  labels = c("A", "B"),
  ncol = 2, nrow = 1)
figure3

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

#SLA
emtrends(fecund_modeltraits, specs = c("Herbivore"), var = "I(sSLA^2)")

#Context dependent selection
standard_sSLA <- ggpredict(fecund_modeltraits, terms = c("sSLA[all]","Water","Herbivore"), type = "re", interval="confidence")

SLA_fecundity <- plot(standard_sSLA, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Specific leaf area")+ scale_y_continuous("Fecundity (Summed length of siliques)")+ ylim(0,2000)

#succulence
emtrends(fecund_modeltraits, specs = c("Herbivore"), var = "sSUC")

#Context dependent selection

standard_sSUC <- ggpredict(fecund_modeltraits, terms = c("sSUC[all]","Herbivore"), type = "re", interval="confidence")

succulence_fecundity <-plot(standard_sSUC, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols2, facet=FALSE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Succulence")+ scale_y_continuous("Fecundity")

#Flowering phenology
emtrends(fecund_modeltraits, specs = c("Water","Herbivore"), var = "sFT")

#Context dependent selection
standard_FT <- ggpredict(fecund_modeltraits, terms = c("sFT[all]", "Water","Herbivore"), type = "re", interval="confidence")
FT_fecundity <-plot(standard_FT, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering phenology")+ scale_y_continuous(element_blank())+ ylim(0,2000)


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

duration_fecundity <-plot(standard_duration, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("Flowering duration")+ scale_y_continuous(element_blank())+ ylim(0,2000)






##No selection
standard_sLAR <- ggpredict(fecund_modeltraits, terms = c("sLAR[all]", "Water","Herbivore"), type = "re", interval="confidence")
plot(standard_sLAR, show_data=TRUE, show_title =FALSE, show_legend=TRUE, colors = cols, facet=TRUE)+theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+theme(legend.position="right")+scale_x_continuous("LAR")+ scale_y_continuous("Fecundity") + ylim(0,2000)



####

library(ggpubr)
figure4 <- ggarrange(
  SLA_fecundity,
  FT_fecundity,
  duration_fecundity,
  
  
  labels = c("A", "B", "C"),
  ncol = 3, nrow = 1)
figure4

