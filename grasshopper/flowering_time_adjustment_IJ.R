######## PROJECT: grasshopper experiment: variation in herbivore damage due to water availability 
#### PURPOSE:Calculate adjusted flowering time values
#### AUTHORS: Inam Jameel & Jill Anderson
#### DATE LAST MODIFIED: 14 feb 24

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
library(MuMIn)

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper/files_from_jill/")
setwd("~/Documents/personnel/Jameel/grasshopper/summary")

 #this is where you specify the folder where you have the data on your computer
 
#read in data 
grasshopper <- read.csv("flowering_time_adjustment.csv",stringsAsFactors=T)
grasshopper$year<-as.factor(grasshopper$year)

grasshopper$Elapsed_days<-grasshopper$Census_2_OD -grasshopper$Census_1_OD

sapply(grasshopper,class)

##Restrict the dataset to cases where there were no more than 10 elapsed days

sub1<-subset(grasshopper, Elapsed_days<11)


final_model<-lm(Silique_length_Census2~Elapsed_days-1 , data= sub1)
summary(final_model)
Anova(final_model,type="III")
visreg(final_model, "Elapsed_days",  ylab="Length of siliques (mm)", xlab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(final_model)
plot(final_model)

mod1 <-lmer(Elapsed_days~ Silique_length_Census2-1 +(1|year), data= sub1)
summary(mod1)
Anova(mod1,type="III")
visreg(mod1, "Silique_length_Census2",  xlab="Length of siliques (mm)", ylab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod1)




mod2<-lmer(Silique_number_Census2 ~Elapsed_days-1 +(1|year), data= sub1)
summary(mod2)
Anova(mod2,type="III")
visreg(mod2, "Elapsed_days",  ylab="Number of siliques", xlab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod2)


##Could it be that both number and silique length are important predictors? No. (I also tried gamma and it didn't really change anything)
mod3<-lmer(Elapsed_days~Silique_number_Census2 +Silique_length_Census2-1 +(1|year), data= sub1, )
summary(mod3)
Anova(mod3,type="III")
visreg(mod3, "Silique_number_Census2",  xlab="Number of siliques", ylab="Time (elapsed days)", partial= TRUE , type="conditional")
visreg(mod3, "Silique_length_Census2",  xlab="Length of siliques (mm)", ylab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod3)



