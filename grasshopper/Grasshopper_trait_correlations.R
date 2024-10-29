#creating correlation tables for each trait x treatment

##this is where you specify the folder where you have the data on your computer

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


##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")

traitdat <- dplyr::select(grasshopper_no2021,FT_Adj,Max_height_flowering, avg_LAR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block,  year, Mature_length_siliques,Reproduced,flowering_duration, SLA, succulence, FT_Adj)





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


## all traits ####

local_grasshopper <- filter(grasshopper_no2021, elev_dist >= -10, elev_dist <= 13)

#SLA 234.24, 12.58
mean(local_grasshopper$SLA,na.rm = TRUE)
sd(local_grasshopper$SLA,na.rm = TRUE) / sqrt(length(local_grasshopper$SLA[!is.na(local_grasshopper$SLA)]))

#Flowering phenology 172.28, 1.97
mean(local_grasshopper$FT_Adj,na.rm = TRUE)
sd(local_grasshopper$FT_Adj,na.rm = TRUE) / sqrt(length(local_grasshopper$FT_Adj[!is.na(local_grasshopper$FT_Adj)]))

#Height at flowering 27.6, 3.36
mean(local_grasshopper$Max_height_flowering,na.rm = TRUE)
sd(local_grasshopper$Max_height_flowering,na.rm = TRUE) / sqrt(length(local_grasshopper$Max_height_flowering[!is.na(local_grasshopper$Max_height_flowering)]))


