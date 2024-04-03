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

##Box plot
VWC_box<-ggplot(VWC_data, aes(x = Water, y = VWC, fill = Water)) +
  geom_boxplot(outlier.shape = NA) +xlab("Water availability")+ 
  scale_y_continuous("Volumetric water content") +
  geom_point(pch = 21, position = position_jitterdodge(0.3))

VWC_plot<-VWC_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none") +  scale_fill_manual(values = c(cols), name = "Water treatment")+facet_wrap(~Year)


VWC_plot # Fig SX

