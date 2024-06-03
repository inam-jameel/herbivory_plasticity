######## PROJECT: greenhouse experiment: local adaptation to herbivory
#### PURPOSE:Calculate adjusted flowering time values
#### AUTHORS:  Jill Anderson
#### DATE LAST MODIFIED: 8 April 24

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

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/greenhouse/data")

#setwd("~/Documents/personnel/Jameel/chapter3/Greenhouse_census_summary")

 #this is where you specify the folder where you have the data on your computer
 
#read in data 
greenhouse <- read.csv("flowering_time_adjustment_greenhouse.csv",stringsAsFactors=T)
greenhouse$Elapsed_days<-greenhouse$Census_2_OD -greenhouse$Census_1_OD
hist(greenhouse$Elapsed_days)
sapply(greenhouse,class)

##Restrict the dataset to cases where there were no more than 30 elapsed days



#Adjust flowering time based on final_model in flowering_time_adjustment.R
greenhouse1$FT_Adj<-round((greenhouse1$Date_flowering_exp - (greenhouse1$Silique_length_flowering /7)),0)

plot(greenhouse1$FT_Adj~greenhouse1$Date_flowering_exp)

grasshopper$Snowmelt_FT_Adj<-grasshopper$FT_Adj-grasshopper$Day_of_snowmelt
plot(grasshopper$Snowmelt_FT_Adj ~grasshopper$Snowmelt_Date_flowering)

hist(grasshopper$Ordinal_Date_flowering)
hist(grasshopper$FT_Adj)
#####

sub1<-subset(greenhouse, Elapsed_days<31)


mod1 <-lm(Elapsed_days~ Silique_length_Census2-1 , data= sub1)
summary(mod1)
Anova(mod1,type="III")
visreg(mod1, "Silique_length_Census2",  xlab="Length of siliques (mm)", ylab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod1)
plot(mod1)





mod2<-lm(Silique_number_Census2 ~Elapsed_days-1 , data= sub1)
summary(mod2)
Anova(mod2,type="III")
visreg(mod2, "Elapsed_days",  ylab="Number of siliques", xlab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod2)


##Could it be that both number and silique length are important predictors? 
mod3<-lm(Elapsed_days~Silique_number_Census2 +Silique_length_Census2-1 , data= sub1)
summary(mod3)
Anova(mod3,type="III")
visreg(mod3, "Silique_number_Census2",  xlab="Number of siliques", ylab="Time (elapsed days)", partial= TRUE , type="conditional")

visreg(mod3, "Silique_length_Census2",  xlab="Length of siliques (mm)", ylab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod3)

mod4<-lm(Silique_length_Census2~Elapsed_days-1 , data= sub1)
summary(mod4)
Anova(mod4,type="III")
visreg(mod4, "Elapsed_days",  ylab="Length of siliques (mm)", xlab="Time (elapsed days)", partial= TRUE , type="conditional")
r.squaredGLMM(mod4)
plot(mod4)






