######## PROJECT: grasshopper experiment: variation in herbivore damage due to water availability 
#### PURPOSE:Examine fitness, traits in response to water availability and herbivory .
#### AUTHORS: Jill Anderson
#### DATE LAST MODIFIED: 28 March 24

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

 #this is where you specify the folder where you have the data on your computer
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")
setwd("~/Documents/grasshopper")



#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_updated_21March24.csv",stringsAsFactors=T)

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper$year<-as.factor(grasshopper$year)

##Change the baseline for the water treatment
grasshopper $Water <-factor(grasshopper $Water, levels = c("Restricted", "Supplemental"))

##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)




#*******************************************************************************
#### 2.leaf damage across censuses #####
#*******************************************************************************
#####repeated measures with all damage data####

#reformat datafile

LAR_data<- grasshopper %>% pivot_longer(cols=c("LAR_1","LAR_2","LAR_3","LAR_4","LAR_5","LAR_6"),
                                        names_to='census',
                                        values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam,  Cage_Block,  S_elev,census, year)

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

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

LAR_data $y_beta<- (LAR_data $LAR_prop*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(LAR_data $y_beta)

max(LAR_data $y_beta)



##Box plot
LAR_box <-ggplot(LAR_data, aes(x = Herbivore, y = LAR_prop, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore Water")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~year)



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


LAR_fig= ggplot(newdf2, aes(x= S_elev,y= predicted, group= Water, colour= Water))+geom_point(size=2,aes(x= S_elev,y= LAR_prop)) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source Elevation")  
LAR_fig + theme_bw() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "top")  +scale_colour_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Restricted","Supplemental")) +geom_smooth(method = "glm", 
    method.args = list(family = "quasibinomial"),  se = FALSE) + facet_wrap(~ Herbivore:year)
    
    ##Conversation of standardized elevation back to raw units
    -1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
    0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
    1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)



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




##Box plot
LAR_box <-ggplot(LAR_data, aes(x = Herbivore, y = LAR_prop, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore manipulation")+ 
  scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

damage_treatment<-LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c(cols), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_wrap(~year)

#cline, plotting predicted data
damage_cline= ggplot(newdf2, aes(x= S_elev,y= predicted, group= Water, 
                                 colour= Water))+geom_point(pch = 20) + scale_y_continuous("Leaf area removed by herbivores (proportion)")+ scale_x_continuous("Source elevation")   + theme_classic() + facet_grid(Herbivore~ year) + theme(axis.title.y = element_text(size=10, colour = "black")) +theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),panel.grid.minor=element_blank(),legend.position = "Top")+geom_smooth(method = "glm", method.args = list(family = "quasibinomial"),  se = FALSE) +scale_colour_manual(values = cols, name = "Water availability", labels = c("Restricted","Supplemental"))

damage_cline #fig 2

