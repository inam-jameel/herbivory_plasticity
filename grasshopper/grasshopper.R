######## PROJECT: grasshopper experiment: variation in herbivore damage due to water availability 
#### PURPOSE:Examine fitness, traits in response to water availability and herbivory .
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 4 jan 24

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

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/grasshopper")

 #this is where you specify the folder where you have the data on your computer


#read in data 
grasshopper <- read.csv("Grasshopper_fulldata_long_05Jan24.csv")

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Genotype<-as.factor(grasshopper$Genotype)
grasshopper$population<-as.factor(grasshopper$population)
grasshopper$Treatment<-as.factor(grasshopper$Treatment)
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Herbivore<-as.factor(grasshopper$Herbivore)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper $PlantID <-as.factor(grasshopper $PlantID)
grasshopper $Cage_Block <-as.factor(grasshopper $Cage_Block)

##Change the baseline for treatment
grasshopper $Treatment <-factor(grasshopper $Treatment, levels = c("Control", "watered"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)

#This rescales year
grasshopper $S_year<-scale(grasshopper $year,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

#Let's concatenate cage and block to reduce nesting necessary for random effects.
grasshopper $block<-interaction(grasshopper$Cage, grasshopper$Block,sep = "_")

#Let's concatenate herbivore and watering treatments, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Treatment,sep = "_")


grasshopper <- filter(grasshopper, Exclude == "Include")

#if trying only in herbivore treatment

#herbivore_only <- filter(grasshopper, Herbivore == "Grasshopper")

#*******************************************************************************
#### watering treatment effect #####
#*******************************************************************************

#just looking at cage/block. no plants

VWC_gg<-ggplot(grasshopper, aes(x = Treatment, y = avg_vwc, fill = Treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Watering treatment")+ scale_y_continuous("Volumetric Water Content") +
  geom_point(pch = 21, position = position_jitterdodge())
VWC_gg + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                               panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("Ample water","Water-restricted")) +  scale_fill_manual(values = c( "lightblue","#CC79A7"))


VWC_repeated <- lmer(avg_vwc~ Treatment*Herbivore*S_year + (1|Cage_Block), data=grasshopper )

visreg(VWC_repeated, 'Treatment', type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Volumetric Water Content", partial=TRUE) 


Anova(VWC_repeated, type = "III")


#pull the fitted values out and plot them in ggplot

vwc_df <- grasshopper %>% 
  
  mutate(fit.m = predict(VWC_repeated, re.form = NA),
         
         fit.c = predict(VWC_repeated, re.form = NULL), #all random effects
         
         resid = residuals(VWC_repeated))

vioplot(fit.c ~ Treatment, data= vwc_df, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,ylim=c(0,20),colMed = "black",colMed2 = c("lightblue","#CC79A7"), col=c("lightblue","#CC79A7"), ylab="Volumetric water content", xlab="Watering treatment") +stripchart(fit.c ~ Treatment, data= vwc_df, col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)

ggplot(vwc_df, aes(x = Treatment, y = fit.c, group = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.95, trim = T) +
  geom_jitter(shape = 21, position = position_jitter(0.15), fill = "gray", size = 2, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", size = 0.5, width = 0.33) +
  
  facet_grid(~ year) +
  theme_light() +
  scale_fill_manual(values = c("cadetblue2", "brown1")) +
  labs(y = "Volumetric Water Content") +
  labs(x = "") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  theme(legend.position = "none") 



#*******************************************************************************
#### 1.Fitness #####
#*******************************************************************************

hurdle_Model <- glmmTMB(Mature_length_siliques ~Treatment*S_elev*Herbivore+S_year + (1|PlantID)+ (1| population)+ (1| Cage_Block), data=grasshopper, zi=~Treatment*S_elev*Herbivore+S_year + (1|PlantID)+ (1| population)+ (1| Cage_Block),family=ziGamma(link="log"))


##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows significiant source elevation.
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model,type="III", component="cond")

library(ggeffects)

plot(ggpredict(hurdle_Model, terms = "Herbivore","Treatment"))

result <- ggpredict(hurdle_Model, c("Herbivore","Treatment"))
plot(result)



require(AICcmodavg)
require(performance)
require(DHARMa)
diagnose(hurdle_Model)

plotQQunif(hurdle_Model)
plotResiduals(hurdle_Model)

summary(hurdle_Model)









## Probability of reproduction ####


grasshopper_pf= ggplot(grasshopper, aes(x= elev_km,y= Reproduced, group= Herbivore, 
                       colour= Herbivore))+geom_point(size=5) + scale_y_continuous("Probability of flowering")+ scale_x_continuous("Source Elevation")  
grasshopper_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ Treatment, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))


####treatment####

mod_repro <-glmmTMB(Reproduced~S_elev*Treatment*Herbivore+year+(1|PlantID),data=grasshopper,family=binomial(link="logit"))

Anova(mod_repro, type="III") #sig interaction btw mat_treat and elev

visreg(mod_repro_REavgLARa,"elev", by="mat_treat", overlay=FALSE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg(mod_repro_REavgLAR, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Probability of flowering", partial=FALSE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255"))) 




#*******************************************************************************
#### 1.leaf damage across treatments #####
#*******************************************************************************

#repeated measures with average damage data across years 

grasshopper1 <- drop_na(grasshopper,avg_LAR) 


##convert LAR (measured as % of leaf area removed by herbivores) to proportion
grasshopper1 $LAR_prop<- grasshopper1 $avg_LAR/100
##Check that the conversion worked
hist(grasshopper1 $LAR_prop)
head(grasshopper1)  
sapply(grasshopper1,class)    


#Let's concatenate water and herbivore treatment so that we can visualize all treatment levels in the same visreg model.
grasshopper1 $treat_herb<-interaction(grasshopper1 $Treatment, grasshopper1 $Herbivore,sep = "_")

head(grasshopper1)

#First, we will use ggplot to look at the data:

LAR_RM = ggplot(grasshopper1, aes(x= elevation,y=LAR_prop, group= treat_herb))+geom_point(aes(colour= treat_herb ),size=5) +scale_x_continuous("elevation")+ scale_y_continuous("LAR_prop")
LAR_RM + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank())+geom_smooth(aes(colour=treat_herb), method="glm",size=1.6)



LAR_RM_box<-ggplot(grasshopper1, aes(x = Herbivore, y = LAR_prop, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))

#basic lmer model. This framework isn't great for proportional data, but it's the most straightforward approach. 

modA<- lmer (LAR_prop~ Treatment*Herbivore*elev_km+ (1|Cage_Block)+(1|Genotype),data= grasshopper1)
Anova(modA,type="III")
##You can see here that the residuals aren't great
plot(modA)

#Despite the heteroscedasticity in the residuals, hte plot here pretty much demonstrates the major patterns
visreg(modA, overlay = TRUE, "elev_km", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


##This model allows us to visualize all treatments on the same figure
modB <- lmer (LAR_prop~treat_herb*elev_km+ (1|Cage_Block)+(1|Genotype),data= grasshopper1)
Anova(modB,type="III")
visreg(modB, overlay = TRUE, "elev_km", by="treat_herb", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

## The most appropriate model for these proportional data is a zero-inflated beta regression. This is because we have values of 0, but no values of one (or we would need a zero-one inflated beta regression). The GAMLSS package in R can execute these models. This code fits data on a [0,1) interval, which is what we have.

##Gamlss models fail if there are any NA values in the entire dataset. So, I exclude NAs here

grasshopper2 <- dplyr::select(grasshopper1, LAR_prop, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, elev_km, S_elev,block, treat_herb,Cage_Block,avg_vwc)



##We can proceed to the zero inflated beta regression.

gamlss.model<- gamlss (formula=LAR_prop~Treatment* S_elev*Herbivore+ random(Cage_Block)+random(Genotype), 
                       sigma.formula=LAR_prop~Treatment* S_elev*Herbivore+ random(Cage_Block)+random(Genotype), nu.formula=LAR_prop~Treatment* S_elev*Herbivore+ random(Cage_Block)+random(Genotype),family=BEZI, data= grasshopper2)
summary(gamlss.model)

#has not converged



visreg(gamlss.model, overlay = TRUE, "S_elev", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Elevation (KM)", ylab="Leaf Area Herbivorized (Scaled)", partial=TRUE,
       fill=list(col=grey
                 (c(0.99), alpha=0)
       ), band = TRUE, gg = TRUE,
       
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

dropterm(gamlss.model)

mod2<-stepGAIC(gamlss.model) # suggests LAR_prop ~ Treatment + S_elev + Herbivore + random(Cage_Block) +    Treatment:S_elev 

mod2$anova

summary(mod2)


# beta transformation

ggplot(grasshopper1, aes(x= LAR_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(year ~  .)

n<-nrow(grasshopper1)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

grasshopper1 $y_beta<- (grasshopper1 $LAR_prop*(n-1) + 0.5)/n

hist(grasshopper1 $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(grasshopper1 $y_beta)

max(grasshopper1 $y_beta)


#Then, the analysis for Treatment/Herbivore

library(betareg)
beta_modela<- betareg ( y_beta ~S_elev*Treatment*Herbivore+year, data=grasshopper1)
Anova(beta_modela,type="III")
visreg(beta_modela, "elev_km", scale="response")

visreg(beta_modela, 'Treatment', by= "Herbivore", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=TRUE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=.2,col=c("#6699cc","#882255"))) 

summary(beta_modela)

plot(beta_modela)

grasshopper3 <- dplyr::select(grasshopper1, LAR_prop, y_beta, year, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, elev_km, S_elev,block, treat_herb,Cage_Block,avg_vwc)


gamlss_moda<- gamlss (formula= y_beta ~S_elev*Treatment*Herbivore+year + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=grasshopper3,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_moda)
plot(gamlss_moda)
Anova(gamlss_moda)


#pull the fitted values out and plot them in ggplot

newdf1 <- grasshopper3 %>% 
  
  mutate(fit.m = predict(gamlss_moda, re.form = NA),
         
         fit.c = predict(gamlss_moda, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_moda))

##Convert fit.m and fit.c back to the proportional scale

newdf1 $fit.m_trans<-1/(1+exp(-(newdf1 $fit.m))) 

newdf1 $fit.c_trans<-1/(1+exp(-(newdf1 $fit.m))) 

newdf1 $resid_trans<-1/(1+exp(-(newdf1 $resid))) 



LAR_fig =ggplot(newdf1,aes(x= elev_km,y= fit.c_trans,shape= treat_herb, linetype= treat_herb,color= treat_herb, group= treat_herb)) + 
  
  geom_point(aes(shape= treat_herb),size=4)+scale_shape_manual(values = c(0,2,15,17)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +scale_linetype_manual(values=c("dotdash", "dotdash","solid","solid"))+
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","lightblue","darkred"))

LAR_fig


LAR_RM_box<-ggplot(newdf1, aes(x = Herbivore, y = fit.c_trans, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_RM_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))






####repeated measures with all damage data####

#needs the file with herbivore data
#read in data 
grasshopperLAR <- read.csv("Grasshopper_damage.csv")

sapply(grasshopperLAR,class)

grasshopperLAR$Genotype<-as.factor(grasshopperLAR$Genotype)
grasshopperLAR$population<-as.factor(grasshopperLAR$population)
grasshopperLAR$Treatment<-as.factor(grasshopperLAR$Treatment)
grasshopperLAR$Block<-as.factor(grasshopperLAR$Block)
grasshopperLAR$Herbivore<-as.factor(grasshopperLAR$Herbivore)
grasshopperLAR$Cage<-as.factor(grasshopperLAR$Cage)
grasshopperLAR $PlantID <-as.factor(grasshopperLAR $PlantID)

grasshopperLAR $Treatment <-factor(grasshopperLAR $Treatment, levels = c("watered", "Control"))
grasshopperLAR $S_elev<-scale(grasshopperLAR $elevation,center=TRUE, scale=TRUE)
grasshopperLAR $S_initdiam<-scale(grasshopperLAR $init.diam,center=TRUE, scale=TRUE)
grasshopperLAR$elev_km<-grasshopperLAR $elevation/1000
grasshopperLAR $block<-interaction(grasshopperLAR$Cage, grasshopperLAR$Block,sep = "_")
grasshopperLAR $treat<-interaction(grasshopperLAR$Herbivore, grasshopperLAR$Treatment,sep = "_")

#reformat datafile
LAR_data<- grasshopperLAR %>% pivot_longer(cols=c("LAR_census1_2021","LAR_census2_2021","LAR_census3_2021","LAR_census3_2022","LAR_census4_2022","LAR_census6_2022","LAR_census8_2022","LAR_census11_2022","LAR_census14_2022","LAR_census5_2023","LAR_census7_2023"),
                                 names_to='exposure',
                                 values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, elev_km, S_elev, exposure,block, treat)

LAR_data$census<-LAR_data$exposure
LAR_data$census[LAR_data$census == "LAR_census1_2021"] <- "1"
LAR_data$census[LAR_data$census == "LAR_census2_2021"] <- "2"
LAR_data$census[LAR_data$census == "LAR_census3_2021"] <-"3"
LAR_data$census[LAR_data$census == "LAR_census3_2022"] <- "4"
LAR_data$census[LAR_data$census == "LAR_census4_2022"] <- "5"
LAR_data$census[LAR_data$census == "LAR_census6_2022"] <- "6"
LAR_data$census[LAR_data$census == "LAR_census8_2022"] <- "7"
LAR_data$census[LAR_data$census == "LAR_census11_2022"] <-"8"
LAR_data$census[LAR_data$census == "LAR_census14_2022"] <- "9"
LAR_data$census[LAR_data$census == "LAR_census5_2023"] <- "10"
LAR_data$census[LAR_data$census == "LAR_census7_2023"] <- "11"

LAR_data $census <-as.numeric(LAR_data $census)

LAR_data$year<-LAR_data$exposure
LAR_data$year[LAR_data$year == "LAR_census1_2021"] <- "2021"
LAR_data$year[LAR_data$year == "LAR_census2_2021"] <- "2021"
LAR_data$year[LAR_data$year == "LAR_census3_2021"] <-"2021"
LAR_data$year[LAR_data$year == "LAR_census3_2022"] <- "2022"
LAR_data$year[LAR_data$year == "LAR_census4_2022"] <- "2022"
LAR_data$year[LAR_data$year == "LAR_census6_2022"] <- "2022"
LAR_data$year[LAR_data$year == "LAR_census8_2022"] <- "2022"
LAR_data$year[LAR_data$year == "LAR_census11_2022"] <-"2022"
LAR_data$year[LAR_data$year == "LAR_census14_2022"] <- "2022"
LAR_data$year[LAR_data$year == "LAR_census5_2023"] <- "2023"
LAR_data$year[LAR_data$year == "LAR_census7_2023"] <- "2023"

LAR_data $year <-as.numeric(LAR_data $year)


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


#Then, the analysis for Treatment/Herbivore

library(betareg)
beta_modelb<- betareg ( y_beta ~S_elev*Treatment*Herbivore+census, data=LAR_data)
Anova(beta_modelb,type="III")
visreg(beta_modelb, "S_elev", scale="response")

visreg(beta_modelb, 'S_elev', by= "Herbivore", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=TRUE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=.5,col=c("#6699cc","#882255"))) 

visreg(beta_modelb, 'Herbivore', by= "Treatment", overlay = FALSE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=TRUE)
       
       # cex.lab = 1.5,cex.axis = 1.5,
       #fill=list(col=grey(c(1), alpha=0.1)),
       #line=list(col=c("#6699cc","#882255")),
       #points=list(cex=.5,col=c("#6699cc","#882255"))) 

summary(beta_modelb)

plot(beta_modelb)

gamlss_modb<- gamlss (formula= y_beta ~S_elev*Treatment*Herbivore+census + random(PlantID)+ random(Block)+random(population),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_modb)
plot(gamlss_modb)
Anova(gamlss_modb)

# from Jill's email on 1/4/2024
#using block is cage/block concatenated 

#This asks whether patterns differ across years and controls for censuses with a random effect
LAR_data $census1 <-as.factor(LAR_data $census)


Mod1<- gamlss (formula= y_beta ~S_elev*Treatment*Herbivore*year+ random(census1) + random(PlantID)+ random(block)+random(population),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(Mod1)


# from Jill's email on 1/4/2024
#This simply treats year as a potential nuisance factor – ideally, we would nest census within year, but I don’t think nesting is possible in gamlss

##### LAR Model to use says Jill ####

Mod2<- gamlss (formula= y_beta ~S_elev*Treatment*Herbivore+year+ random(census1) + random(PlantID)+ random(block)+random(population),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500))
summary(Mod2)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(Mod2, re.form = NA),
         
         fit.c = predict(Mod2, re.form = NULL), #all random effects
         
         resid = residuals(Mod2))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 


LAR_box <-ggplot(newdf2, aes(x = Herbivore, y = fit.c_trans, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())

LAR_box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))

LAR_fig =ggplot(newdf2,aes(x= elev_km,y= fit.c_trans,shape= Herbivore, linetype= Herbivore,color= Herbivore, group= Herbivore)) + 
  
  geom_point(aes(shape= Herbivore),size=4)+scale_shape_manual(values = c(0,17,20)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed (%)") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))

LAR_fig



####repeated measures in 2023####

LAR_2023<-filter(LAR_data, census > 9)

hist(LAR_2023$LAR_prop)

ggplot(LAR_2023, aes(x= LAR_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

#LAR_2023 <- drop_na(LAR_2023,LAR_prop) 

#n<-nrow(LAR_2023)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

#LAR_2023 $y_beta<- (LAR_2023 $LAR_prop*(n-1) + 0.5)/n

#hist(LAR_2023 $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

#min(LAR_2023 $y_beta)

#max(LAR_2023 $y_beta)


#Then, the analysis for Treatment/Herbivore

library(betareg)
beta_modelb<- betareg ( y_beta ~elev_km*avg_vwc_2023*Herbivore+census, data=LAR_2023)
Anova(beta_modelb,type="III")
visreg(beta_modelb, "elev_km", scale="response")

visreg(beta_modelb, 'elev_km', by= "Herbivore", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=TRUE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=.5,col=c("#6699cc","#882255"))) 

summary(beta_modelb)

plot(beta_modelb)
gamlss_modb<- gamlss (formula= y_beta ~elev_km*avg_vwc_2023*Herbivore+census + random(PlantID)+ random(Block)+random(population),family=BE(mu.link = "logit"), data=LAR_2023,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_modb)
plot(gamlss_modb)
Anova(gamlss_modb)


#pull the fitted values out and plot them in ggplot

newdf3 <- LAR_2023 %>% 
  
  mutate(fit.m = predict(gamlss_modb, re.form = NA),
         
         fit.c = predict(gamlss_modb, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_modb))

##Convert fit.m and fit.c back to the proportional scale

newdf3 $fit.m_trans<-1/(1+exp(-(newdf3 $fit.m))) 

newdf3 $fit.c_trans<-1/(1+exp(-(newdf3 $fit.m))) 

newdf3 $resid_trans<-1/(1+exp(-(newdf3 $resid))) 



LAR_fig =ggplot(newdf3,aes(x= elev_km,y= LAR,shape= Herbivore, linetype= Herbivore,color= Herbivore, group= Herbivore)) + 
  
  geom_point(aes(shape= Herbivore),size=4)+scale_shape_manual(values = c(0,17,20)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))

LAR_fig



##plot the partial residuals (instead of raw data) 


predictedProb_fig =ggplot(newdf3,aes(x= elev_km,y= fit.m_trans + resid_trans,shape= Herbivore, linetype= Herbivore,color= Herbivore, group= Herbivore)) + 
  
  #geom_line(aes(y= fit.m_trans, lty= Herbivore), size=0.8) +
  
  geom_smooth(method = "lm", se = TRUE,formula=y~x) + 
  
  theme_bw()  

predictedProb_fig



#*******************************************************************************
#### 2.leaf traits #####
#*******************************************************************************

#repeated measures, have to make seperate data frams for SLA, LAR, VWC and combine 

rosette_SLA_data<- grasshopper %>% pivot_longer(cols=c("rosette_SLA_2023","rosette_SLA_2022"),
                                        names_to='year',
                                        values_to='SLA')
r_SLA_data <- dplyr::select(rosette_SLA_data, SLA, PlantID, year,elev_km)


r_SLA_data$year[r_SLA_data$year == "rosette_SLA_2022"] <- "1"
r_SLA_data$year[r_SLA_data$year == "rosette_SLA_2023"] <- "2"

r_SLA_data$year<-as.numeric(r_SLA_data$year)
ggplot(r_SLA_data, aes(x= SLA))+ geom_histogram(color="black", fill="white")+ facet_grid(year ~  .)


#VWC for 2022 and 2023
VWC_data_SLA<- grasshopper %>% pivot_longer(cols=c("avg_vwc_2023","avg_vwc_2022"),
                                            names_to='year',
                                            values_to='VWC')

VWC_data_SLA <- dplyr::select(VWC_data_SLA, VWC, PlantID,year)


VWC_data_SLA$year[VWC_data_SLA$year == "avg_vwc_2022"] <- "1"
VWC_data_SLA$year[VWC_data_SLA$year == "avg_vwc_2023"] <- "2"

VWC_data_SLA$year<-as.numeric(VWC_data_SLA$year)


#LAR for 2022 and 2023
LAR_data_SLA<- grasshopper %>% pivot_longer(cols=c("avg_LAR_2023","avg_LAR_2022"),
                                            names_to='year',
                                            values_to='LAR')

LAR_data_SLA <- dplyr::select(LAR_data_SLA, LAR, PlantID,year)

LAR_data_SLA$year[LAR_data_SLA$year == "avg_LAR_2022"] <- "1"
LAR_data_SLA$year[LAR_data_SLA$year == "avg_LAR_2023"] <- "2"

LAR_data_SLA$year<-as.numeric(LAR_data_SLA$year)
LAR_data_SLA$LAR_prop<-LAR_data_SLA $LAR/100
hist(LAR_data_SLA$LAR_prop)



LAR_VWC <- merge(VWC_data_SLA,LAR_data_SLA, by=c("PlantID","year"))

SLA_LAR_VWC <- merge(LAR_VWC_test,r_SLA_data, by=c("PlantID","year"))

plant_data <- dplyr::select(grasshopper, elevation, Genotype, population, Cage, Treatment, Herbivore, Block, PlantID, init.diam, S_initdiam, block, S_elev,treat)


repeated_measures <- merge(SLA_LAR_VWC_test,plant_data, by=c("PlantID"))

# damage and VWC 
SLA_RM <- lmer(SLA ~ VWC*LAR_prop*elev_km+year+(1|PlantID)+(1|population)+(1|Cage) +(1|Block:Cage),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = repeated_measures)

plot(SLA_RM)

Anova(SLA_RM, type = "III") # elevation is significant, not year

visreg(SLA_RM,"elev_km", by="VWC", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))


# treatment and herbivore
SLA_RMa <- lmer(SLA ~ Treatment*Herbivore*elev_km+year+(1|PlantID) +(1|population)+(1|Cage) +(1|Block:Cage),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = repeated_measures)

plot(SLA_RMa)

Anova(SLA_RMa, type = "III") # elevation, year is significant

visreg(SLA_RMa,"elev_km", by="Treatment", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

#### SLA individual for 2022 ####

grasshopper$avg_LAR_2022_prop<-grasshopper $avg_LAR_2022/100
hist(grasshopper$avg_LAR_2022_prop)

SLA_2022 =ggplot(grasshopper, aes(x= elevation,y= rosette_SLA_2022, color = treat))+ geom_point(size=3.5)+ ggtitle("Rosette SLA 2022 by Source Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_2022

SLA_2022b =ggplot(grasshopper, aes(x= elevation,y= bolt_SLA_2022, color = treat))+ geom_point(size=3.5)+ ggtitle("Bolt SLA 2022 by Source Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_2023b


#with treatment
SLA_2022_mod <- lmer(rosette_SLA_2022 ~ Treatment*Herbivore*elev_km +(1|population)+(1|Cage) +(1|Block:Cage),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = grasshopper)

plot(SLA_2022_mod)

Anova(SLA_2022_mod, type = "III")

visreg(SLA_2022_mod,"elev_km", by="Treatment", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

#with avg damage/mat avg damage
SLA_2022_moda <- lmer(rosette_SLA_2022 ~ avg_vwc_2022*avg_LAR_2022_prop*elev_km +(1|Genotype)+(1|Block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = grasshopper)
plot(SLA_2022_moda)

Anova(SLA_2022_moda, type = "III")

summary(SLA_2022_moda)

visreg2d(SLA_2022_moda,"elev_km","avg_LAR_2022_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")


#### SLA individual for 2023 ####


hist(grasshopper$avg_LAR_2023_prop)

SLA_2023 =ggplot(grasshopper, aes(x= elevation,y= rosette_SLA_2023, color = treat))+ geom_point(size=3.5)+ ggtitle("Rosette SLA 2023 by Source Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_2023

SLA_2023b =ggplot(grasshopper, aes(x= elevation,y= bolt_SLA_2023, color = treat))+ geom_point(size=3.5)+ ggtitle("Bolt SLA 2023 by Source Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_2023b


#with treatment
SLA_2023_mod <- lmer(rosette_SLA_2023 ~ Treatment*Herbivore*elev_km +(1|population)+(1|Cage) +(1|Block:Cage),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = grasshopper)

plot(SLA_2023_mod)

Anova(SLA_2023_mod, type = "III")

visreg(SLA_2023_mod,"elev_km", by="Treatment", overlay=TRUE,   scale = "response", xlab="Source elevation (Km)", ylab="SLA CM3/g ", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#56B4E9","darkred")), points=list(col=c("#56B4E9","darkred")),fill=list(col=grey(c(1), alpha=0.4)))

#with avg damage/avg_vwc_2023
SLA_2023_moda <- lmer(rosette_SLA_2023 ~avg_vwc_2023*avg_LAR_2023_prop*elev_km +(1|population)+(1|Cage) +(1|Block:Cage),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = grasshopper)
plot(SLA_2023_moda)

Anova(SLA_2023_moda, type = "III")

summary(SLA_2023_moda)

visreg2d(SLA_2023_moda,"elev_km","avg_LAR_2023_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")

