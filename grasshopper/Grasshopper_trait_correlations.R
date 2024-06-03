#creating correlation tables for each trait x treatment

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

foliar$succulence_mg<-foliar$succulence*1000
foliar$Suc_mg<- foliar$succulence_mg+0.01


##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced_updated
grasshopper_no2021<-subset(grasshopper, year!="2021")
# retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper_no2021,FT_Adj,Max_height_flowering, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block,  year, Mature_length_siliques,Reproduced,flowering_duration, SLA, lwc, succulence, FT_Adj)





## Restricted, addition ####

Restricted_Addition <- filter(traitdat, Water == "Restricted", Herbivore == "Addition")

#avg_LAR

hist(Restricted_Addition$avg_LAR)

plot(Restricted_Addition$avg_LAR~Restricted_Addition$SLA)
cor.test(Restricted_Addition$avg_LAR,Restricted_Addition$SLA, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$avg_LAR~Restricted_Addition$succulence)
cor.test(Restricted_Addition$avg_LAR,Restricted_Addition$succulence, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$avg_LAR~Restricted_Addition$FT_Adj)
cor.test(Restricted_Addition$avg_LAR,Restricted_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$avg_LAR~Restricted_Addition$flowering_duration)
cor.test(Restricted_Addition$avg_LAR,Restricted_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$avg_LAR~Restricted_Addition$Max_height_flowering)
cor.test(Restricted_Addition$avg_LAR,Restricted_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")



#SLA

hist(Restricted_Addition$SLA)

#plot(Restricted_Addition$SLA~Restricted_Addition$avg_LAR)
plot(Restricted_Addition$SLA~Restricted_Addition$succulence)
cor.test(Restricted_Addition$SLA,Restricted_Addition$succulence, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$SLA~Restricted_Addition$FT_Adj)
cor.test(Restricted_Addition$SLA,Restricted_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$SLA~Restricted_Addition$flowering_duration)
cor.test(Restricted_Addition$SLA,Restricted_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$SLA~Restricted_Addition$Max_height_flowering)
cor.test(Restricted_Addition$SLA,Restricted_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#succulence

hist(Restricted_Addition$succulence)

plot(Restricted_Addition$succulence~Restricted_Addition$FT_Adj)
cor.test(Restricted_Addition$succulence,Restricted_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$succulence~Restricted_Addition$flowering_duration)
cor.test(Restricted_Addition$succulence,Restricted_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$succulence~Restricted_Addition$Max_height_flowering)
cor.test(Restricted_Addition$succulence,Restricted_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#FT_Adj

hist(Restricted_Addition$FT_Adj)

plot(Restricted_Addition$FT_Adj~Restricted_Addition$flowering_duration)
cor.test(Restricted_Addition$FT_Adj,Restricted_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Addition$FT_Adj~Restricted_Addition$Max_height_flowering)
cor.test(Restricted_Addition$FT_Adj,Restricted_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#flowering_duration

hist(Restricted_Addition$flowering_duration)

plot(Restricted_Addition$flowering_duration~Restricted_Addition$Max_height_flowering)
cor.test(Restricted_Addition$flowering_duration,Restricted_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#Max_height_flowering

hist(Restricted_Addition$Max_height_flowering)





## Supplemental addition ####

Supplemental_Addition <- filter(traitdat, Water == "Supplemental", Herbivore == "Addition")


#avg_LAR

hist(Supplemental_Addition$avg_LAR)

plot(Supplemental_Addition$avg_LAR~Supplemental_Addition$SLA)
cor.test(Supplemental_Addition$avg_LAR,Supplemental_Addition$SLA, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$avg_LAR~Supplemental_Addition$succulence)
cor.test(Supplemental_Addition$avg_LAR,Supplemental_Addition$succulence, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$avg_LAR~Supplemental_Addition$FT_Adj)
cor.test(Supplemental_Addition$avg_LAR,Supplemental_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$avg_LAR~Supplemental_Addition$flowering_duration)
cor.test(Supplemental_Addition$avg_LAR,Supplemental_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$avg_LAR~Supplemental_Addition$Max_height_flowering)
cor.test(Supplemental_Addition$avg_LAR,Supplemental_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")



#SLA

hist(Supplemental_Addition$SLA)

#plot(Supplemental_Addition$SLA~Supplemental_Addition$avg_LAR)
plot(Supplemental_Addition$SLA~Supplemental_Addition$succulence)
cor.test(Supplemental_Addition$SLA,Supplemental_Addition$succulence, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$SLA~Supplemental_Addition$FT_Adj)
cor.test(Supplemental_Addition$SLA,Supplemental_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$SLA~Supplemental_Addition$flowering_duration)
cor.test(Supplemental_Addition$SLA,Supplemental_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$SLA~Supplemental_Addition$Max_height_flowering)
cor.test(Supplemental_Addition$SLA,Supplemental_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#succulence

hist(Supplemental_Addition$succulence)

plot(Supplemental_Addition$succulence~Supplemental_Addition$FT_Adj)
cor.test(Supplemental_Addition$succulence,Supplemental_Addition$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$succulence~Supplemental_Addition$flowering_duration)
cor.test(Supplemental_Addition$succulence,Supplemental_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$succulence~Supplemental_Addition$Max_height_flowering)
cor.test(Supplemental_Addition$succulence,Supplemental_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#FT_Adj

hist(Supplemental_Addition$FT_Adj)

plot(Supplemental_Addition$FT_Adj~Supplemental_Addition$flowering_duration)
cor.test(Supplemental_Addition$FT_Adj,Supplemental_Addition$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Addition$FT_Adj~Supplemental_Addition$Max_height_flowering)
cor.test(Supplemental_Addition$FT_Adj,Supplemental_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#flowering_duration

hist(Supplemental_Addition$flowering_duration)

plot(Supplemental_Addition$flowering_duration~Supplemental_Addition$Max_height_flowering)
cor.test(Supplemental_Addition$flowering_duration,Supplemental_Addition$Max_height_flowering, method = 'pearson', use = "complete.obs")


#Max_height_flowering

hist(Supplemental_Addition$Max_height_flowering)


## Restricted removal ####

Restricted_Removal <- filter(traitdat, Water == "Restricted", Herbivore == "Removal")


#avg_LAR

hist(Restricted_Removal$avg_LAR)

plot(Restricted_Removal$avg_LAR~Restricted_Removal$SLA)
cor.test(Restricted_Removal$avg_LAR,Restricted_Removal$SLA, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$avg_LAR~Restricted_Removal$succulence)
cor.test(Restricted_Removal$avg_LAR,Restricted_Removal$succulence, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$avg_LAR~Restricted_Removal$FT_Adj)
cor.test(Restricted_Removal$avg_LAR,Restricted_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$avg_LAR~Restricted_Removal$flowering_duration)
cor.test(Restricted_Removal$avg_LAR,Restricted_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$avg_LAR~Restricted_Removal$Max_height_flowering)
cor.test(Restricted_Removal$avg_LAR,Restricted_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")



#SLA

hist(Restricted_Removal$SLA)

#plot(Restricted_Removal$SLA~Restricted_Removal$avg_LAR)
plot(Restricted_Removal$SLA~Restricted_Removal$succulence)
cor.test(Restricted_Removal$SLA,Restricted_Removal$succulence, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$SLA~Restricted_Removal$FT_Adj)
cor.test(Restricted_Removal$SLA,Restricted_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$SLA~Restricted_Removal$flowering_duration)
cor.test(Restricted_Removal$SLA,Restricted_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$SLA~Restricted_Removal$Max_height_flowering)
cor.test(Restricted_Removal$SLA,Restricted_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#succulence

hist(Restricted_Removal$succulence)

plot(Restricted_Removal$succulence~Restricted_Removal$FT_Adj)
cor.test(Restricted_Removal$succulence,Restricted_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$succulence~Restricted_Removal$flowering_duration)
cor.test(Restricted_Removal$succulence,Restricted_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$succulence~Restricted_Removal$Max_height_flowering)
cor.test(Restricted_Removal$succulence,Restricted_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#FT_Adj

hist(Restricted_Removal$FT_Adj)

plot(Restricted_Removal$FT_Adj~Restricted_Removal$flowering_duration)
cor.test(Restricted_Removal$FT_Adj,Restricted_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Restricted_Removal$FT_Adj~Restricted_Removal$Max_height_flowering)
cor.test(Restricted_Removal$FT_Adj,Restricted_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#flowering_duration

hist(Restricted_Removal$flowering_duration)

plot(Restricted_Removal$flowering_duration~Restricted_Removal$Max_height_flowering)
cor.test(Restricted_Removal$flowering_duration,Restricted_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#Max_height_flowering

hist(Restricted_Removal$Max_height_flowering)


## Supplemental removal ####

Supplemental_Removal <- filter(traitdat, Water == "Supplemental", Herbivore == "Removal")


#avg_LAR

hist(Supplemental_Removal$avg_LAR)

plot(Supplemental_Removal$avg_LAR~Supplemental_Removal$SLA)
cor.test(Supplemental_Removal$avg_LAR,Supplemental_Removal$SLA, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$avg_LAR~Supplemental_Removal$succulence)
cor.test(Supplemental_Removal$avg_LAR,Supplemental_Removal$succulence, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$avg_LAR~Supplemental_Removal$FT_Adj)
cor.test(Supplemental_Removal$avg_LAR,Supplemental_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$avg_LAR~Supplemental_Removal$flowering_duration)
cor.test(Supplemental_Removal$avg_LAR,Supplemental_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$avg_LAR~Supplemental_Removal$Max_height_flowering)
cor.test(Supplemental_Removal$avg_LAR,Supplemental_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")



#SLA

hist(Supplemental_Removal$SLA)

#plot(Supplemental_Removal$SLA~Supplemental_Removal$avg_LAR)
plot(Supplemental_Removal$SLA~Supplemental_Removal$succulence)
cor.test(Supplemental_Removal$SLA,Supplemental_Removal$succulence, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$SLA~Supplemental_Removal$FT_Adj)
cor.test(Supplemental_Removal$SLA,Supplemental_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$SLA~Supplemental_Removal$flowering_duration)
cor.test(Supplemental_Removal$SLA,Supplemental_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$SLA~Supplemental_Removal$Max_height_flowering)
cor.test(Supplemental_Removal$SLA,Supplemental_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#succulence

hist(Supplemental_Removal$succulence)

plot(Supplemental_Removal$succulence~Supplemental_Removal$FT_Adj)
cor.test(Supplemental_Removal$succulence,Supplemental_Removal$FT_Adj, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$succulence~Supplemental_Removal$flowering_duration)
cor.test(Supplemental_Removal$succulence,Supplemental_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$succulence~Supplemental_Removal$Max_height_flowering)
cor.test(Supplemental_Removal$succulence,Supplemental_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#FT_Adj

hist(Supplemental_Removal$FT_Adj)

plot(Supplemental_Removal$FT_Adj~Supplemental_Removal$flowering_duration)
cor.test(Supplemental_Removal$FT_Adj,Supplemental_Removal$flowering_duration, method = 'pearson', use = "complete.obs")

plot(Supplemental_Removal$FT_Adj~Supplemental_Removal$Max_height_flowering)
cor.test(Supplemental_Removal$FT_Adj,Supplemental_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#flowering_duration

hist(Supplemental_Removal$flowering_duration)

plot(Supplemental_Removal$flowering_duration~Supplemental_Removal$Max_height_flowering)
cor.test(Supplemental_Removal$flowering_duration,Supplemental_Removal$Max_height_flowering, method = 'pearson', use = "complete.obs")


#Max_height_flowering

hist(Supplemental_Removal$Max_height_flowering)



#avg_LAR

hist(traitdat$avg_LAR)

plot(traitdat$avg_LAR~traitdat$SLA)
cor.test(traitdat$avg_LAR,traitdat$SLA, method = 'pearson', use = "complete.obs")

plot(traitdat$avg_LAR~traitdat$succulence)
cor.test(traitdat$avg_LAR,traitdat$succulence, method = 'pearson', use = "complete.obs")

plot(traitdat$avg_LAR~traitdat$FT_Adj)
cor.test(traitdat$avg_LAR,traitdat$FT_Adj, method = 'pearson', use = "complete.obs")

plot(traitdat$avg_LAR~traitdat$flowering_duration)
cor.test(traitdat$avg_LAR,traitdat$flowering_duration, method = 'pearson', use = "complete.obs")

plot(traitdat$avg_LAR~traitdat$Max_height_flowering)
cor.test(traitdat$avg_LAR,traitdat$Max_height_flowering, method = 'pearson', use = "complete.obs")



#SLA

hist(traitdat$SLA)

#plot(traitdat$SLA~traitdat$avg_LAR)
plot(traitdat$SLA~traitdat$succulence)
cor.test(traitdat$SLA,traitdat$succulence, method = 'pearson', use = "complete.obs")

plot(traitdat$SLA~traitdat$FT_Adj)
cor.test(traitdat$SLA,traitdat$FT_Adj, method = 'pearson', use = "complete.obs")

plot(traitdat$SLA~traitdat$flowering_duration)
cor.test(traitdat$SLA,traitdat$flowering_duration, method = 'pearson', use = "complete.obs")

plot(traitdat$SLA~traitdat$Max_height_flowering)
cor.test(traitdat$SLA,traitdat$Max_height_flowering, method = 'pearson', use = "complete.obs")


#succulence

hist(traitdat$succulence)

#plot(traitdat$succulence~traitdat$SLA)
#plot(traitdat$succulence~traitdat$avg_LAR)
plot(traitdat$succulence~traitdat$FT_Adj)
cor.test(traitdat$succulence,traitdat$FT_Adj, method = 'pearson', use = "complete.obs")

plot(traitdat$succulence~traitdat$flowering_duration)
cor.test(traitdat$succulence,traitdat$flowering_duration, method = 'pearson', use = "complete.obs")

plot(traitdat$succulence~traitdat$Max_height_flowering)
cor.test(traitdat$succulence,traitdat$Max_height_flowering, method = 'pearson', use = "complete.obs")


#FT_Adj

hist(traitdat$FT_Adj)

#plot(traitdat$FT_Adj~traitdat$SLA)
#plot(traitdat$FT_Adj~traitdat$succulence)
#plot(traitdat$FT_Adj~traitdat$avg_LAR)
plot(traitdat$FT_Adj~traitdat$flowering_duration)
cor.test(traitdat$FT_Adj,traitdat$flowering_duration, method = 'pearson', use = "complete.obs")

plot(traitdat$FT_Adj~traitdat$Max_height_flowering)
cor.test(traitdat$FT_Adj,traitdat$Max_height_flowering, method = 'pearson', use = "complete.obs")


#flowering_duration

hist(traitdat$flowering_duration)

#plot(traitdat$flowering_duration~traitdat$SLA)
#plot(traitdat$flowering_duration~traitdat$succulence)
#plot(traitdat$flowering_duration~traitdat$FT_Adj)
#plot(traitdat$flowering_duration~traitdat$avg_LAR)
plot(traitdat$flowering_duration~traitdat$Max_height_flowering)
cor.test(traitdat$flowering_duration,traitdat$Max_height_flowering, method = 'pearson', use = "complete.obs")


#Max_height_flowering

hist(traitdat$Max_height_flowering)

#plot(traitdat$Max_height_flowering~traitdat$SLA)
#plot(traitdat$Max_height_flowering~traitdat$succulence)
#plot(traitdat$Max_height_flowering~traitdat$FT_Adj)
#plot(traitdat$Max_height_flowering~traitdat$flowering_duration)
#plot(traitdat$Max_height_flowering~traitdat$avg_LAR)





