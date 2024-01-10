##### This R script contains pilot analyses of data from the 2022 growing season for the grasshopper experiment. We do not have final fitness data, so fecundity analyses should be treated as prpreliminary.


#libraries
library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
library(tidyr)
library(broom)
library(car)
library(MASS)
library(emmeans)
library(visreg)
library(gamlss)
library(effects)
    library(mgcv)
library(glmmTMB)
library(bbmle)

setwd("~/Documents/Jill/postdoc/RMBL/2022/Gothic")


#read in data 
grasshopper <- read.csv("grasshopper_summary26July22.csv")

sapply(grasshopper,class)
##Some  variables are being read as characters not factors. Let's fix that
grasshopper$Genotype<-as.factor(grasshopper$Genotype)
grasshopper$population<-as.factor(grasshopper$population)
grasshopper$Treatment<-as.factor(grasshopper$Treatment)
grasshopper$Block<-as.factor(grasshopper$Block)
grasshopper$Herbivore<-as.factor(grasshopper$Herbivore)
grasshopper$Cage<-as.factor(grasshopper$Cage)
##Change the baseline for treatment
grasshopper $Treatment <-factor(grasshopper $Treatment, levels = c("watered", "Control"))


##This standardizes source elevation to a mean of 0 and standard deviation of 1, which is often necessary for model convergence
grasshopper $S_elev<-scale(grasshopper $elevation,center=TRUE, scale=TRUE)

##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
grasshopper $S_initdiam<-scale(grasshopper $init.diam,center=TRUE, scale=TRUE)


#This resacles source elevation from meters to km
grasshopper$elev_km<-grasshopper $elevation/1000

#Let's concatenate cage and block to reduce nesting necessary for random effects.
grasshopper $block<-interaction(grasshopper$Cage, grasshopper$Block,sep = "_")

#Let's concatenate herbivore and watering treatments, which is helpful for some models.
grasshopper $treat<-interaction(grasshopper$Herbivore, grasshopper$Treatment,sep = "_")


###########################################
#############Analysis 1: Fitness ##########
###########################################

##We will take a hurdle model approach to fitness, first analyzing the probability of reproduction amongst individuals that successfully overwintered, and then fecundity amongst individuals that reproduced.

### Analysis 1a: Logistic regression of the probability of reproduction. This model includes a covariate for size at planting
logistic_model<-glmer(reproduction ~ S_initdiam + elev_km* Herbivore* Treatment+ (1|block) +(1|population), data= grasshopper, family=binomial(link=logit), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
  Anova(logistic_model , type="III")

## There is some indication of a treatment effect and a treatment by source elevation effect, but the p-values are marginal, and the patterns are opposite of predictions. (Trend for increasing probability of reproduction with source elevation under dry conditions, and decreasing prob. of repro. with source elevation under water conditions)
 ##Plotting the relationship between the probability of reproduction and source elevation under different environments.
    logmod <-predictorEffect("elev_km",  partial.residuals=TRUE, logistic_model)
  plot(logmod, lwd=2,xlab="Source elevation (km)", ylab="Probability of reproduction", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),
       partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1.0))

###Analysis 1b: Fruit production amongst individuals that reproduced. This analysis needs to be treated cautiously because the plants haven't finished maturing fruits. Here we are using the number of fruits as a fitness proxy because we do not yet have mature fruit lengths. The length of mature fruits is much more highly correlated with seed set.

fecund<- glmer(silique_number ~  S_initdiam +elev_km* Herbivore* Treatment+ (1|block) +(1|population), data= grasshopper,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecund,type="III")
plot(fecund)

#Plotting the significant source elevation by herbivore treatment by watering treatment interaction.

##As expected, there is a decline in fecundity with source elevation under dry conditions with grasshopper enrichment (Treatment = Control; Herbivore = Grasshopper). This decline indicates that low elevation genotypes have the greatest fecundity in dry/herbivore enriched environments. This is likely the only significant slope in the 4 plots below.
   fecundmod <-predictorEffect("elev_km",  partial.residuals=TRUE, fecund)
  plot(fecundmod, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),
       partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,60))
       
       

## This is just another approach to visualizing these patterns:
fecundB<- glmer(silique_number ~  S_initdiam +elev_km* treat+ (1|block) +(1|population), data= grasshopper,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecundB,type="III")

  visreg(fecundB,"elev_km", by="treat", overlay=FALSE, scale = "response", xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))



##Analysis 1c: We can also do a hurdle model in one step in the glmmTMB package. First, we'll check a couple of diferent distributions. Running all of these models takes a while
hurdle_Model = glmmTMB (total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), zi=~ S_initdiam +elev_km* Herbivore* Treatment  + (1|block)+ (1| population),data = grasshopper ,family=truncated_nbinom2)

fit_hnbinom1 = glmmTMB (total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), zi=~ S_initdiam +elev_km* Herbivore* Treatment  + (1|block)+ (1| population),data = grasshopper ,family=truncated_nbinom1)

fit_trunc_poisson = glmmTMB (total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), zi=~ S_initdiam +elev_km* Herbivore* Treatment  + (1|block)+ (1| population),data = grasshopper ,family=truncated_poisson)

fit_ziGamma = glmmTMB (total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), zi=~ S_initdiam +elev_km* Herbivore* Treatment  + (1|block)+ (1| population),data = grasshopper ,family=ziGamma)


AICtab(hurdle_Model, fit_hnbinom1, fit_trunc_poisson, fit_ziGamma)

##This is the ANOVA table for the logistic regression part (probability of reproduction). It shows a signifciant source elevation by treatment interaction, similar to Analysis 1a.
Anova(hurdle_Model,type="III", component="zi")
##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). It shows a a marginally significant source elevation by herbivore by treatment interaction.
Anova(hurdle_Model,type="III", component="cond")

##Analysis 1d: We can also do a hurdle model in one step in the glmmTMB package. First, we'll check a couple of diferent distributions. Running all of these models takes a while
model_poisson <- glmmTMB(total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), data=grasshopper, ziformula = ~1, family = poisson)


model_nb2 <- glmmTMB(total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), data=grasshopper,  ziformula = ~1, family = nbinom2)

model_nb1 <- glmmTMB(total_fitness ~ S_initdiam +elev_km* Herbivore* Treatment + (1|block)+ (1| population), data=grasshopper,  ziformula = ~1, family = nbinom1)

AICtab(model_nb2, model_nb1, model_poisson)

Anova(model_nb2,type="III")

  zeroinflated <-predictorEffect("elev_km",  partial.residuals=TRUE, model_nb2)
  plot(zeroinflated, lwd=2,xlab="Source elevation (km)", ylab="Fecundity", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),
       partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,60))

###########################################
#############Analysis 2: Leaf damage ##########
###########################################
##In these analyses, we find that leaf damage increases with source elevation (as expected) more in the grasshopper treatment than the control.

##Analysis of leaf damage (LAR = leaf area removed) as a function of source elevation by watering treament by grasshopper treatment. We are going to use leaf damage data from the most recent census. Utimately, we will do a repeated measures analysis, but not right now.
grasshopper1 <- drop_na(grasshopper,LAR_census14) 

         
##convert LAR (measured as % of leaf area removed by herbivores) to proportion
grasshopper1 $LAR_prop14<- grasshopper1 $LAR_census14/100
      ##Check that the conversion worked
       hist(grasshopper1 $LAR_prop14)
head(grasshopper1)  
sapply(grasshopper1,class)    


#Let's concatenate water and herbivore treatment so that we can visualize all treatment levels in the same visreg model.
      grasshopper1 $treat_herb<-interaction(grasshopper1 $Treatment, grasshopper1 $Herbivore,sep = "_")
 
head(grasshopper1)

#First, we will use ggplot to look at the data:

LAR_prop14 = ggplot(grasshopper1, aes(x= elevation,y=LAR_prop14, group= treat_herb))+geom_point(aes(colour= treat_herb ),size=5) +scale_x_continuous("elevation")+ scale_y_continuous("LAR_prop14")
LAR_prop14 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank())+geom_smooth(aes(colour=treat_herb), method="glm",size=1.6)



LAR_prop143<-ggplot(grasshopper1, aes(x = Herbivore, y = LAR_prop14, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
LAR_prop143 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(),legend.position = "bottom")
                           
  
  
  LAR_prop143<-ggplot(grasshopper1, aes(x = Herbivore, y = LAR_prop14, fill = Treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
LAR_prop143 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("grasshopper exclosure", "grasshopper enclosure")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Watering treatment", labels = c("Ample water","Water-restricted"))
                           
 
 
 
 #Let's concatenate cage and Block so that we don't have to worry about nesting 
      grasshopper1 $cage_block<-interaction(grasshopper1 $Cage, grasshopper1 $Block,sep = "_")
            grasshopper1 $trtherb<-interaction(grasshopper1 $Treatment, grasshopper1 $Herbivore,sep = "_")


##Now, we will do a basic lmer model. This frameowrk isn't great for proportional data, but it's the most straightforward approach. 

    modA<- lmer (LAR_prop14~ Treatment*Herbivore*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper1)
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
 modB <- lmer (LAR_prop14~trtherb*elev_km+ (1|cage_block)+(1|Genotype),data= grasshopper1)
 Anova(modB,type="III")
   visreg(modB, overlay = TRUE, "elev_km", by="trtherb", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf damage from insect herbivores (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

## The most appropriate model for these proportional data is a zero-inflated beta regression. This is because we have values of 0, but no values of one (or we would need a zero-one inflated beta regression). The GAMLSS package in R can execute these models. This code fits data on a [0,1) interval, which is what we have.

##Gamlss models fail if there are any NA values in the entire dataset. So, I exclude NAs here
grasshopper1 <- grasshopper1[c(1,4,5,7:13,35:36,85,87:88)]
 grasshopper1 $LAR_prop14<- grasshopper1 $LAR_census14/100

##We can proceed to the zero inflated beta regression. The issue is that these models are complicated to explain. Given that the result is nearly the same as the lmer model above, I suggest that Lisa presents the results from the lmer model      
      gamlss.model<- gamlss (formula=LAR_prop14~Treatment* elev_km*Herbivore+ random(block)+random(Genotype), 
         sigma.formula=LAR_prop14~Treatment* elev_km*Herbivore+ random(block)+random(Genotype), nu.formula=LAR_prop14~Treatment* elev_km*Herbivore+ random(block)+random(Genotype),family=BEZI, data= grasshopper1)
      summary(gamlss.model)

    
          
visreg(gamlss.model, overlay = TRUE, "elev_km", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR_prop14", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

   

                
