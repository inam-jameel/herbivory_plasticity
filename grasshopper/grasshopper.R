######## PROJECT: grasshopper experiment: water availability on herbivory 
#### PURPOSE:Examine fitness, traits in response to water availibility and herbivory .
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 27 Nov 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for plotting
require(glmmTMB) # for running survival model, have to load twice
require(gamlss) # for running phenology model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(vioplot) #for violin plots

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field/2023_files/census_files")

setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field/grasshopper")
 #this is where you specify the folder where you have the data on your computer

#for macbook air
#setwd("/Users/inam/OneDrive - University of Georgia/Inam_experiments/field_2023_files/census_files")

#read in data 
grasshopper <- read.csv("Grasshopper_damage.csv")

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Genotype<-as.factor(grasshopper$Genotype)
grasshopper$population<-as.factor(grasshopper$population)
grasshopper$Treatment<-as.factor(grasshopper$Treatment)
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Herbivore<-as.factor(grasshopper$Herbivore)
grasshopper$Cage<-as.factor(grasshopper$Cage)
grasshopper $PlantID <-as.factor(grasshopper $PlantID)

#grasshopper$Include<-as.factor(grasshopper$Include)

##Change the baseline for treatment
grasshopper $Treatment <-factor(grasshopper $Treatment, levels = c("watered", "Control"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)


#This rescales source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

#Let's concatenate cage and block to reduce nesting necessary for random effects.
grasshopper $block<-interaction(grasshopper$Cage, grasshopper$Block,sep = "_")

#Let's concatenate herbivore and watering treatments, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Treatment,sep = "_")


#trying only in treatment

grasshopper <- filter(grasshopper, Herbivore == "Grasshopper")

#*******************************************************************************
#### 1.leaf damage  across treatments #####
#*******************************************************************************

#repeated measures 

LAR_data<- grasshopper %>% pivot_longer(cols=c("LAR_census1_2021","LAR_census2_2021","LAR_census3_2021","LAR_census3_2022","LAR_census4_2022","LAR_census6_2022","LAR_census8_2022","LAR_census11_2022","LAR_census14_2022","LAR_census5_2023","LAR_census7_2023"),
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



#Then, the analysis for treatment/mat_treat

library(betareg)
beta_modela<- betareg ( y_beta ~elev_km*Treatment+census, data=LAR_data)
Anova(beta_modela,type="III")
visreg(beta_modela, "elev_km", scale="response")

visreg(beta_modela, 'elev_km', by= "Treatment", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 

summary(beta_modela)

plot(beta_modela)
gamlss_moda<- gamlss (formula= y_beta ~elev_km*Treatment+census + random(PlantID)+ random(Block)+random(population),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_moda)
plot(gamlss_moda)
Anova(gamlss_moda)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(gamlss_moda, re.form = NA),
         
         fit.c = predict(gamlss_moda, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_moda))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 



LAR_fig =ggplot(newdf2,aes(x= elev_km,y= LAR,shape= Herbivore, linetype= Herbivore,color= Herbivore, group= Herbivore)) + 
  
  geom_point(aes(shape= Herbivore),size=4)+scale_shape_manual(values = c(0,17,20)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))

LAR_fig




##And if you wanted to plot the partial residuals (instead of raw data) – again, you could customize this to your heart’s delight:


predictedProb_fig =ggplot(newdf2,aes(x= elev,y= fit.m_trans + resid_trans,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  geom_smooth(method = "lm", se = TRUE,formula=y~x) + 
  
  theme_bw()  

predictedProb_fig



##In these analyses, we find that leaf damage increases with source elevation (as expected) more in the grasshopper treatment than the control.

##Analysis of leaf damage (LAR = leaf area removed) as a function of source elevation by watering treament by grasshopper treatment. We are going to use leaf damage data from the most recent census. Utimately, we will do a repeated measures analysis, but not right now.
grasshopper1 <- drop_na(census5,LAR) 
grasshopper2 <- drop_na(census7,LAR) 
  
         
##convert LAR (measured as % of leaf area removed by herbivores) to proportion
grasshopper1 $LAR<- grasshopper1 $LAR/100
grasshopper2 $LAR<- grasshopper2 $LAR/100
      ##Check that the conversion worked
       hist(grasshopper1 $LAR)
       hist(grasshopper2 $LAR)


#Let's concatenate water and herbivore treatment so that we can visualize all treatment levels in the same visreg model.
      grasshopper1 $treat_herb<-interaction(grasshopper1 $Treatment, grasshopper1 $Herbivore,sep = "_")
      grasshopper2 $treat_herb<-interaction(grasshopper2 $Treatment, grasshopper2 $Herbivore,sep = "_")
      
head(grasshopper1)

#First, we will use ggplot to look at the data:

LAR_prop1 = ggplot(grasshopper1, aes(x= elevation,y=LAR, group= treat_herb))+geom_point(aes(colour= treat_herb ),size=5) +scale_x_continuous("elevation")+ scale_y_continuous("LAR")
LAR_prop1 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank())+geom_smooth(aes(colour=treat_herb), method="glm",size=1.6)

  
  
  LAR_prop1<-ggplot(grasshopper1, aes(x = Herbivore, y = LAR, fill = Treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
  
LAR_prop1 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))

LAR_prop2<-ggplot(grasshopper2, aes(x = Herbivore, y = LAR, fill = Treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
LAR_prop2 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","#CC79A7"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))
                           
 
 
 
 #Let's concatenate cage and Block so that we don't have to worry about nesting 
      grasshopper1 $cage_block<-interaction(grasshopper1 $Cage, grasshopper1 $Block,sep = "_")
            grasshopper1 $trtherb<-interaction(grasshopper1 $Treatment, grasshopper1 $Herbivore,sep = "_")
            
            grasshopper2 $cage_block<-interaction(grasshopper2 $Cage, grasshopper2 $Block,sep = "_")
            grasshopper2 $trtherb<-interaction(grasshopper2 $Treatment, grasshopper2 $Herbivore,sep = "_")


##Now, we will do a basic lmer model. This frameowrk isn't great for proportional data, but it's the most straightforward approach. 

    modA<- lmer (LAR~ Treatment*Herbivore*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper1)
 Anova(modA,type="III")
 
 modB<- lmer (LAR~ Treatment+Herbivore*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper2)
 Anova(modB,type="III")
 
 modB<- lmer (LAR~elev_km*Treatment*Herbivore+(1|cage_block)+(1|Genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=grasshopper2)
 Anova(modB)
 
 visreg(modB,overlay=T,"Treatment" ,by="Herbivore")
 plot(predictorEffects(modB, ~ elev_km), type="response",partial.residuals=TRUE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Leaf area herbivorized",line=list(multiline=TRUE, lty=1:2,col=c("lightblue","darkred")))
 
 
 
 ##You can see here that the residuals aren't great
 plot(modB)
 
 #Despite the heteroscedasticity in the residuals, hte plot here pretty much demonstrates the major patterns
  visreg(modA, overlay = TRUE, "elev_km", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))
  
  visreg(modB, overlay = TRUE, "elev_km", by="Treatment", type="conditional", #scale = "response", 
         xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
         fill=list(col="light grey"
                   #(c(0.99), alpha=0)
         ), band = FALSE,
         #line=list(col=grey(c(0.2,0.6))),
         points=list(cex=0.65,  pch=(19)))
       

##This model allows us to visualize all treatments on the same figure
 modC <- lmer (LAR~trtherb*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper1)
 Anova(modC,type="III")
   visreg(modC, overlay = TRUE, "elev_km", by="trtherb", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))
   
   modD <- lmer (LAR~trtherb*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper2)
   Anova(modD,type="III")
   visreg(modD, overlay = TRUE, "elev_km", by="trtherb", type="conditional", scale = "response", 
          xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
          fill=list(col="light grey"
                    #(c(0.99), alpha=0)
          ), band = TRUE,
          line=list(col=c("#009E73","#E69F00","#0072B2","#CC79A7")),
          points=list(col=c("#009E73","#E69F00","#0072B2","#CC79A7"),cex=0.65,  pch=(19)))

## The most appropriate model for these proportional data is a zero-inflated beta regression. This is because we have values of 0, but no values of one (or we would need a zero-one inflated beta regression). The GAMLSS package in R can execute these models. This code fits data on a [0,1) interval, which is what we have.

##Gamlss models fail if there are any NA values in the entire dataset. So, I exclude NAs here
grasshopper1a <- grasshopper1[c(1,4,5,7:13,40,97:104)]

grasshopper2a <- grasshopper2[c(1,4,5,7:13,40,97:104)]


##We can proceed to the zero inflated beta regression. The issue is that these models are complicated to explain. Given that the result is nearly the same as the lmer model above, I suggest that Lisa presents the results from the lmer model      
      gamlss.model<- gamlss (formula=LAR~trtherb* elev_km+ random(cage_block)+random(Genotype), 
         sigma.formula=LAR~trtherb* elev_km*Herbivore+ random(cage_block)+random(Genotype), nu.formula=LAR~trtherb* elev_km*Herbivore+ random(cage_block)+random(Genotype),family=BEZI, data= grasshopper1a)
      summary(gamlss.model)
      


    
          
visreg(gamlss.model, overlay = TRUE, "elev_km", by="trtherb", type="conditional", #scale = "response", 
       xlab="Elevation (KM)", ylab="Leaf Area Herbivorized (Scaled)", partial=TRUE,
       fill=list(col=grey
                 (c(0.99), alpha=0)
       ), band = TRUE, gg = TRUE,

       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


gamlss.model2<- gamlss (formula=LAR~trtherb* elev_km+ random(cage_block)+random(Genotype), 
                       sigma.formula=LAR~trtherb* elev_km*Herbivore+ random(cage_block)+random(Genotype), nu.formula=LAR~trtherb* elev_km*Herbivore+ random(cage_block)+random(Genotype),family=BEZI, data= grasshopper2a)
summary(gamlss.model2)





visreg(gamlss.model2, overlay = TRUE, "elev_km", by="trtherb", type="conditional", #scale = "response", 
       xlab="Elevation (KM)", ylab="Leaf Area Herbivorized (Scaled)", partial=TRUE,
       fill=list(col=grey
                 (c(0.99), alpha=0)
       ), band = TRUE, gg = TRUE,
       
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

pred<-predict(gamlss.model2, type="response", what="mu")
plot_data <- data.frame( predicted_data = pred, actual_data= grasshopper2a$LAR )

grasshopper2b <- cbind(grasshopper2a,pred)


ggplot(grasshopper2b, aes(x = elev_km, y = pred, color = trtherb))+
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Source elevation", y = "Leaf Area Removed", color = "trtherb") + 
  theme_bw() + scale_color_manual(values=c("#009E73","#E69F00","#0072B2","#CC79A7"))


                
