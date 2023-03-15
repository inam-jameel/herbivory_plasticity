#libraries
library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
packageVersion("lme4")
library(tidyr)
library(dplyr)
library(broom)
library(car)
library(MASS)
library(emmeans)
library(visreg)
library(fitdistrplus)
library(gamlss)
library(rjags)
library(effects)
library(glmmTMB)


setwd("~/Desktop/Anderson_data/herbivory/data/field/2022_herbivory/")

#### read in data ####
data<- read.csv("2022_summary_incomplete.csv")



data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)


data$Genotype<-as.factor(data$Genotype)
data$treatment<-as.factor(data$Treatment)
data$block<-as.factor(data$block)
data$cohort<-as.factor(data$cohort)

data$elev_dist<-data$elevation-data$garden_elevation

#data$elev_dist<- abs(data$elev_dist)
data$abs_elev_dist<- abs(data$elev_dist)

data $S_elev_dist_abs<-scale(data$abs_elev_dist,center=TRUE, scale=TRUE) #needed to be scaled to avoid issues in model

data$elev_dist_km<- data$elev_dist/1000
data $S_elev_dist<-scale(data$elev_dist,center=TRUE, scale=TRUE)



# Fitness models: Fitness ~ garden x treatment x elev
#survival, reproduction, fecundity 

#survival

mod_surv<- glmer(survival~elev_km*treatment*garden+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #when i add cohort to the model, it doesnt work.. 
Anova(mod_surv)
plot(mod_surv)

#no residuals, multiline
survived <-predictorEffect("elev_km",  partial.residuals=FALSE, mod_surv)
plot(survived, lwd=2,xlab="Source elevation (Km)", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#residuals
survived <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_surv)
plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

#concatenate garden and treatment to plot from visreg
data $garden_treat<-interaction(data $garden, data $treatment,sep = "_")

mod_2<- glmer(survival~garden_treat*elev_km+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, but still nothing
Anova(mod_2) 

#no residuals, multiline

visreg(mod_2, overlay = FALSE, "elev_km", by="garden_treat", type="conditional", scale = "response", 
       xlab="Source elevation (Km)", ylab="Fitness (Prob Survival)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       #line=list(col=grey(c(0,0.8))),
       #points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = TRUE)





#prob reproduction

mod_reproduction<- glmer(Flowered~treatment*elev_km*garden+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_reproduction,type="III")  #significance of garden
plot(mod_reproduction)

visreg(mod_reproduction, overlay = TRUE, "elev_km", by="treatment", type="conditional", scale = "response", 
       xlab="Source Elevation (Km)", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = TRUE)
       

reproduction <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_reproduction)
plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (Prob Reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

#subset only surviving plants

surv <- subset(data, survival=="1")

mod_reproduction1<- glmer(Flowered~treatment*elev_km*garden+ (1|block)+(1|Genotype), data = surv, control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_reproduction1,type="III")  
plot(mod_reproduction1)


reproduction <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_reproduction1)
plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

mod_reproduction2<- glmer(Flowered~garden_treat*elev_km+ (1|block)+(1|Genotype), data = surv, control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_reproduction2,type="III")  
plot(mod_reproduction2)

visreg(mod_reproduction2, overlay = TRUE, "elev_km", by="garden_treat", type="conditional", scale = "response", 
       xlab="Source elevation", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       #points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = FALSE)


#num fruits
#not working... low numbers? only had a handful reproduce in the first place


fecundA<- glmer(Mature_length_siliques ~  init_size +elev_km* garden* treatment+ (1|block) +(1|population), data= filter(surv, Mature_length_siliques > 0),family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecundA,type="III")

fruit <-predictorEffect("elev_km",  partial.residuals=TRUE,fecundA)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1200))




library(glmmTMB)
library(glmmTMB) #load it twice 

hurdle_Model = glmmTMB (Mature_silique_number ~ treatment * elev_km *garden  + (1| block)+(1|Genotype), zi=~ treatment* elev_km *garden  + (1| block)+(1|Genotype),data = surv ,family=truncated_nbinom2)

diagnose(hurdle_Model)

check_zeroinflation(hurdle_Model)

check_overdispersion(hurdle_Model)

Anova(hurdle_Model,type="III", component="zi")

summary(hurdle_Model)

Anova(hurdle_Model,component="cond")

Anova(hurdle_Model,component="zi")


fruit <-predictorEffect("elev_km",  partial.residuals=TRUE,hurdle_Model)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,100))

##


# local adaptation models: fitness ~ T x elevational distance (source – transplant)
#survival, reproduction, fecundity 

hurdle_ModelLA = glmmTMB (Mature_length_siliques ~ treatment * elev_dist_km+I(elev_dist_km^2)*treatment  + (1| block)+(1|Genotype), zi=~ treatment* elev_dist_km+I(elev_dist_km^2)*treatment *garden  + (1| block)+(1|Genotype),data = subset(surv, Mature_length_siliques > 0) ,family=Gamma(link=log))

Anova(fecundA,type="III")

summary(hurdle_ModelLA)

fruit <-predictorEffect("elev_dist_km",  partial.residuals=FALSE,hurdle_ModelLA)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1000))

fecundA<- glmer(Mature_length_siliques ~  init_size +treatment * elev_dist_km+I(elev_dist_km^2)*treatment+ (1|block) +(1|population), data= filter(surv, Mature_length_siliques > 0),family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecundA,type="III")

fruitA <-predictorEffect("elev_dist_km",  partial.residuals=TRUE,fecundA)
plot(fruitA, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1500))

flowered <- subset(surv, Flowered=="1")


#mod_fruit<- lmer( Mature_silique_number ~ treatment* S_elev_dist_abs+I(S_elev_dist_abs^2)*treatment +cohort+(1|Genotype) ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = flowered)

Anova(mod_fruit)

fruit <-predictorEffect("S_elev_dist_abs",  partial.residuals=TRUE,mod_fruit)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,20))

#model for RMBL talk
visreg(mod_fruit, overlay = FALSE, "S_elev_dist_abs", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (survival)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

mod_fruit1<- lmer( Mature_silique_number ~ treatment* S_elev_dist+I(S_elev_dist^2)*treatment +cohort+(1|block)+(1|Genotype) +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = flowered)

Anova(mod_fruit1)

fruit1 <-predictorEffect("S_elev_dist",  partial.residuals=TRUE,mod_fruit1)
plot(fruit1, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,20))

visreg(mod_fruit1, overlay = FALSE, "S_elev_dist", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (survival)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))




##Then, we examine local adaptation using the probability of flowering as the fitness component. Here, values of 0 include individuals that failed to flower along with those that died

mod_repro<- glmer(Flowered~treatment* elev_dist_km+I(elev_dist_km^2)*treatment +cohort+(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) 

Anova(mod_repro)

pred_repro <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_repro)

plot(pred_repro, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))

#model for RMBL talk
visreg(mod_repro, overlay = FALSE, "elev_dist_km", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))



#survival

mod_surv<- glmer(survival~treatment* elev_dist_km+I(elev_dist_km^2)*treatment +cohort+(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) 

Anova(mod_surv)

mod_surv1 <-predictorEffect("elev_dist_km",  partial.residuals=TRUE, mod_surv)

plot(mod_surv1, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=TRUE, pch=19, col="black"))


#model for RMBL talk
visreg(mod_surv, overlay = FALSE, "elev_dist_km", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (survival)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

### couldnt figure this out
mod_surv<- glmer(survival~treatment+S_elev_dist +(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, gardenXelevation interaction only

Anova(mod_surv)

mod_surv1 <-predictorEffect("S_elev_dist",  partial.residuals=TRUE, mod_surv)

plot(mod_surv1, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=TRUE, pch=19, col="black"))



mod_surv<- glmer(survival~treatment+elev_dist +(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, gardenXelevation interaction only

Anova(mod_surv)
plot(mod_surv)

survived <-predictorEffect("elev_dist",  partial.residuals=TRUE, mod_surv)
plot(survived, lwd=2,xlab="Distance from Source Elevation", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

survived2 <-predictorEffect("elev_dist",  partial.residuals=FALSE, mod_surv)
plot(survived2, lwd=2,xlab="Distance from Source Elevation", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

#### LAR for BT effectiveness ####

LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=TRUE)
LSmeans_LAR <-na.omit(LSmeans_LAR)

LSmeans_LAR$emmean100 <- LSmeans_LAR$emmean/100
#had to delete the genotypes that did not have both c and p 

mod <- lmer (emmean~treatment+elev_km+garden+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)

Anova(mod)

visreg(mod, overlay = TRUE, "garden", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))))  


#estess 
Estess <- subset(LSmeans_LAR,garden=="Estess")

estess_mod <- lmer (emmean~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Estess)
plot(estess_mod)
Anova(estess_mod,type="III")

E1 <- visreg(estess_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 12),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  


#trying gamlss model
E_gam<-na.omit(Estess) #get rid of NAs
E_gam$genotype <- as.factor(E_gam$genotype) #random needs to be a factor
gamlss.model_estess<- gamlss (formula=emmean~elevation*treatment+random(genotype), family=BEZI, data= Estess) #response out of range


#gothic
Gothic <- subset(LSmeans_LAR,garden=="Gothic")

gothic_mod <- lmer (emmean~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Gothic)
plot(gothic_mod)
Anova(gothic_mod,type="III")

G1 <- visreg(gothic_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 12),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  

#schofield

Scho <- subset(LSmeans_LAR,garden=="Schofield")
scho_mod <- lmer (emmean~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Scho)
plot(scho_mod)
Anova(scho_mod,type="III")


S1 <- visreg(scho_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 12),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  



Estess_com <- subset(data_LAR,garden=="Estess")
E_LAR = ggplot(Estess_com, aes(x= elevation,y=LAR))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = treatment), size = 4)
E <- E_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)





####LAR,#untested/outdated with summary file #### 

#need to update with specific censuses to
data$LAR_prop<- data$LAR/100

##Check that the conversion worked

hist(data$LAR_prop)

##creating new dataset for 2021 cohort

cohort2021 <- subset(data,cohort=="2021")


##Check to make sure that you have the necessary columns.

head(cohort2021)      

#Let's concatenate block and garden so that we don't have to worry about nesting block within garden

cohort2021 $block_garden<-interaction(cohort2021 $block, cohort2021 $garden,sep = "_")

##some variables are being read as chr not factors

cohort2021$garden<-as.factor(cohort2021$garden)
cohort2021$genotype<-as.factor(cohort2021$genotype)
cohort2021$treatment<-as.factor(cohort2021$treatment)
cohort2021$block<-as.factor(cohort2021$block)

#remove NA values for gamlss

c2021<-na.omit(cohort2021)

##Okay, it looks like there were no NAs in that column, so we can proceed to the zero-one inflated negative binomial regression. What I've written out here is the full model, which has convergence issues.

#LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype)
gamlss.modelA<- gamlss (formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), tau.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), family=BEINF, data= c2021)


##Perhaps we don't want to interact everything. Please keep in mind that garden should be a fixed effect, not random (you only have 2 gardens), but maybe the patterns don't differ much across gardens:



#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.model<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), tau.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), family=BEINF, data= c2021)

#many many warnings

gamlss.model
plot(gamlss.model)
summary(gamlss.model) #nu only has garden as significant

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.model)
mod2<-stepGAIC(gamlss.model)
mod2$anova

summary(mod2) # there are significant effects of treatment, elevation of origin and their interaction for mu (the probability of being damaged). For nu (amount of damage on plants with >0 and <1 damage; i.e., the beta component), there is a significant effect of elevation of origin and garden, but not treatment or elevation by treatment). For tau (probability of 100% damage – which you likely don’t have), there are no significant effects.


#So, we might be able to remove the tau component from the model, and run a zero-inflated model:
  
  gamlss.modelB<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), family=BEINF, data= c2021)

plot(gamlss.modelB)

summary(gamlss.modelB)
  
  
visreg(gamlss.modelB, overlay = TRUE, "S_elev", by="treatment", type="conditional", #scale = "response", 
         xlab="Elev", ylab="LAR", partial=TRUE,
         fill=list(col="light grey"
                   #(c(0.99), alpha=0)
         ), band = FALSE,
         #line=list(col=grey(c(0.2,0.6))),
          points=list(cex=0.65,  pch=(19)))  

c2021 $garden_treat<-interaction(c2021 $treatment, c2021 $garden,sep = "_")
  
gamlss.modelE<- gamlss (formula=LAR_prop~S_elev* garden_treat + random(block_garden)+random(genotype), sigma.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(genotype), nu.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(genotype), family=BEINF, data= c2021)
  
  plot(gamlss.modelE)
  
  summary(gamlss.modelE)
  
  visreg(gamlss.modelE, overlay = FALSE, "S_elev", by="garden_treat", type="conditional", 
         #scale = "response",   
         xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = FALSE)    


  #You can create a datafile with the predicted data from the nu component using this code for plotting:
  pred<-predict(gamlss.modelB, newdata=c2021, type="response", what="nu")

  str(gamlss.modelB)
  
  pred <- as.data.frame(pred)  
  View(pred)
  
  pred_data<- bind_cols(c2021, pred)
  
  View(c2021)
  
  LAR_2021 = ggplot(subset(pred_data,garden=="Gothic"), aes(x= elevation,y=pred))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = factor(treatment)), size = 4)
  P <-LAR_2021 + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                           panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)
P



## linear model, this is the incorrect fit  but nothing else is working
##Now to get LSMEANS
modB<- lmer (LAR_prop~genotype*treatment*garden+(1|Block_Garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)


fam_avg_linear<-emmeans(modB, ~genotype:treatment:garden)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_LAR.csv")


LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=FALSE)
LSmeans_LAR <-na.omit(LSmeans_LAR)



  
  