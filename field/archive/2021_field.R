#looks like i had actual anaylses witht eh 2021 field data here. I use the 2021 summary file. should refer back to this for the code if needed.


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
library(gamlss)
library(effects)

setwd("~/Desktop/Anderson_data/herbivory/data/field/2021_herbivory/")

#### read in data ####
data<- read.csv("Field2021_summary.csv")



data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)

data$Genotype<-as.factor(data$Genotype)
data$treatment<-as.factor(data$Treatment)
data$block<-as.factor(data$block)
data$cohort<-as.factor(data$cohort)

data$elev_dist<-data$elevation-data$garden_elevation

#data$elev_dist<- abs(data$elev_dist)
#data$abs_elev_dist<- abs(data$elev_dist)

#data $S_elev_dist_abs<-scale(data$abs_elev_dist,center=TRUE, scale=TRUE) #needed to be scaled to avoid issues in model

data$elev_dist_km<- data$elev_dist/1000
data $S_elev_dist<-scale(data$elev_dist,center=TRUE, scale=TRUE)


data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $S_init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)

#filtering for commitee meeting
data $block_garden<-interaction(data $block, data $garden,sep = "_")
data $block_treatment<-interaction(data $block, data $treatment,sep = "_")
data $cohort_garden<-interaction(data $cohort, data $garden,sep = "_")


data <- filter(data, garden != "Estess")
data <- filter(data, cohort_garden != "2020_Gothic")


#### Fitness models: Fitness ~ garden x treatment x elev ####
#survival, reproduction, fecundity 
#removed cohort from these models

#surivival

mod_surv<- glmer(survival~garden*elev_km+treatment+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, gardenXelevation interaction only
Anova(mod_surv)
plot(mod_surv)

#concatenate garden and treatment to plot from visreg
data $garden_treat<-interaction(data $garden, data $treatment,sep = "_")


#no residuals, multiline
survived <-predictorEffect("elev_km",  partial.residuals=FALSE, mod_surv)
plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#residuals
survived <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_surv)
plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=3:1, col="black"),
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#prob reproduction

#mod_reproduction<- glmer(flowered~garden+treatment+elev_km+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)


mod_reproduction<- glmer(flowered~garden+treatment+elev_km+ (1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)


mod_reproduction<- glm(flowered~garden+treatment+elev_km, data = data, family=binomial)


Anova(mod_reproduction,type="III")  #significance of garden
plot(mod_reproduction)

visreg(mod_reproduction, overlay = TRUE, "garden", by="treatment", type="conditional", scale = "response", 
       xlab="Garden", ylab="Fitness (Prob Flowering)", partial=F, cex.lab = 1.5,cex.axis = 1.5,
       #fill=list(col="blue"),
       #line=list(col=grey(c(0,0.8))),
       #points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = TRUE)
       

reproduction <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_reproduction)
plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (Prob Reproduction)", type="response",lines=list(multiline=FALSE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, col="black"), ylim=c(0,1))

#subset only surviving plants

surv <- subset(data, survival=="1")

mod_reproduction1<- glmer(flowered~garden+treatment+elev_km+ (1|block)+(1|Genotype), data = surv, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_reproduction1,type="III")  #only trending signif of treatment and elev_km
plot(mod_reproduction1)

visreg(mod_reproduction, overlay = TRUE, "elev_km", by="garden", type="conditional", scale = "response", 
       xlab="Source elevation", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       #fill=list(col="blue"),
       line=list(col=grey(c(0,0.8,0.5))),
       points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = FALSE)

reproduction <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_reproduction1)
plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#num fruits
#not working... low numbers? only had a handful reproduce in the first place


library(glmmTMB)
library(glmmTMB) #load it twice 

hurdle_Model = glmmTMB (Mature_silique_number ~ treatment*elevation  + (1| block)+(1|Genotype), zi=~ treatment  + (1| block)+(1|Genotype),data = surv ,family=truncated_nbinom2)

diagnose(hurdle_Model)

check_zeroinflation(hurdle_Model)

check_overdispersion(hurdle_Model)

Anova(hurdle_Model,type="III", component="zi")

summary(hurdle_Model)

Anova(hurdle_Model,component="cond")

Anova(hurdle_Model,component="zi")


fruit <-predictorEffect("elevation",  partial.residuals=TRUE,hurdle_Model)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

##


#### local adaptation models: fitness ~ T x elevational distance (source – transplant)####
#survival, reproduction, fecundity 


##Then, we examine local adaptation using the probability of flowering as the fitness component. Here, values of 0 include individuals that failed to flower along with those that died

mod_repro<- glmer(flowered~treatment* elev_dist_km+I(elev_dist_km^2)*treatment +cohort+(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) 

Anova(mod_repro)

pred_repro <-predictorEffect("elev_dist_km",  partial.residuals=FALSE, mod_repro)

plot(pred_repro, lwd=2,xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", pch=19, type="response",lines=list(multiline=FALSE, lty=1, col="black"),ylim=c(0,0.25), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))

#model for RMBL talk
visreg(mod_repro, overlay = FALSE, "elev_dist_km", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (probability of flowering)", partial=FALSE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
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

#fecundity

flowered <- subset(surv, flowered=="1")


mod_fruit<- lmer( Mature_silique_number ~ treatment* elev_dist_km+I(elev_dist_km^2)*treatment +cohort+(1|block)+(1|Genotype) +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = flowered)

Anova(mod_fruit)

fruit <-predictorEffect("elev_dist_km",  partial.residuals=TRUE,mod_fruit)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,20))

#model for RMBL talk
visreg(mod_fruit, overlay = FALSE, "elev_dist_km", by="treatment", type="contrast", 
       scale = "response",
       xlab="Elevational transfer distance", ylab="Fitness (survival)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

 



### couldnt figure this out
mod_surv<- glmer(survival~treatment+S_elev_dist +(1|block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, gardenXelevation interaction only

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

mod <- lmer (emmean~treatment+elev_km*garden+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)

Anova(mod)

visreg(mod, overlay = TRUE, "elev_km", by="treatment", type="conditional", scale = "response", 
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

#need to update with specific censuses. going with LAR_7 right now
data$LAR_prop<- data$LAR_7/100

##Check that the conversion worked

hist(data$LAR_prop)

#Let's concatenate block and garden so that we don't have to worry about nesting block within garden

data $block_garden<-interaction(data $block, data $garden,sep = "_")


data_LAR<-data %>% dplyr::select(treatment, S_elev, elev_km,garden,block_garden,Genotype,LAR_prop)

#remove NA values for gamlss

data_LAR<-na.omit(data_LAR)


LSmeans_LAR <- read.csv("apr1623LSmeans_LAR.csv", stringsAsFactors=FALSE)


##Okay, it looks like there were no NAs in that column, so we can proceed to the zero-one inflated negative binomial regression. What I've written out here is the full model, which has convergence issues.

#LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype)
gamlss.modelA<- gamlss (formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(Genotype), sigma.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(Genotype), tau.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)


##Perhaps we don't want to interact everything. Please keep in mind that garden should be a fixed effect, not random (you only have 2 gardens), but maybe the patterns don't differ much across gardens:



#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.model<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), tau.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)

#many many warnings

gamlss.model
plot(gamlss.model)
summary(gamlss.model) #nu only has garden as significant

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.model) #what does this do?
mod2<-stepGAIC(gamlss.model) #50 or more warnings
mod2$anova

summary(mod2) #Mu, probability of being damaged: there are significant effects of treatment, elevation of origin and their interaction. For nu (amount of damage on plants with >0 and <1 damage; i.e., the beta component), there is a significant effect of garden, but not treatment or elevation by treatment. For tau (probability of 100% damage – which you likely don’t have), there are no significant effects.what is sigma?


#So, we might be able to remove the tau component from the model, and run a zero-inflated model:
  
  gamlss.modelB<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)

plot(gamlss.modelB)

summary(gamlss.modelB)
  
  
visreg(gamlss.modelB, "S_elev", by="treatment", type="conditional", #scale = "response", 
         xlab="Elev", ylab="LAR", partial=TRUE,
         fill=list(col="light grey"
                   #(c(0.99), alpha=0)
         ), band = FALSE,
         #line=list(col=grey(c(0.2,0.6))),
          points=list(cex=0.65,  pch=(19)))  

data_LAR $garden_treat<-interaction(data_LAR $treatment, data_LAR $garden,sep = "_")
  
gamlss.modelE<- gamlss (formula=LAR_prop~S_elev* garden_treat + random(block_garden)+random(Genotype), sigma.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)
  
  plot(gamlss.modelE)
  
  summary(gamlss.modelE)
  
  visreg(gamlss.modelE, overlay = FALSE, "S_elev", by="garden_treat", type="conditional", 
         #scale = "response",   
         xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = FALSE)    

##this code doenst work!
  #You can create a datafile with the predicted data from the nu component using this code for plotting:
  pred<-predict(gamlss.modelE, newdata=data_LAR, type="response", what="nu")

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
modB<- lmer (LAR_prop~Genotype*treatment*garden+(1|block_garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)


fam_avg_linear<-emmeans(modB, ~Genotype:treatment:garden)
fam_means_linear<-as.data.frame(fam_avg_linear)

elev <- data_LAR[c("Genotype","elev_km")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_LAR <- merge(fam_means_linear,elev,by="Genotype") #merge the dataframes


write.csv(fam_means_linear,file="apr1623LSmeans_LAR.csv")


####commitee meeting 2023 at this point i abandonded recalculating LSmeans and went to the means i previously calculated last year.####
LSmeans_LAR <- read.csv("LSmeans_LAR3.csv", stringsAsFactors=FALSE)

#remove estess
LSmeans_LAR <- filter(LSmeans_LAR, garden != "Estess")

LSmeans_LAR$Genotype<-as.factor(LSmeans_LAR$genotype)

#had to go in and manually remove the genos that didnt have values for both treatments
# 194_2, schofield
# 267_13A, schofield
# 270_3, schofield

modC<- lmer (emmean~elev*treatment+elev*garden+(1|Genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)

modC<- lm (emmean~elev_km*garden*treatment, data=LSmeans_LAR)

Anova(modC)

visreg(modC,overlay=F,"elev", by ="garden")


#boxplot
LARgga<-ggplot(LSmeans_LAR, aes(x = garden, y = emmean, fill = treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Leaf Area Herbivorized",limits =c(0, 0.15)) +# this sets the scale
  geom_point(pch = 21, position = position_jitterdodge())
LARgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic (Low)", "Schofield (High)")) +  scale_fill_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide"))
                                                                                                                                                                           
# scatter plot
p_lar= ggplot(LSmeans_LAR, aes(x= elev,y= emmean, #group= treatment, 
                               colour= treatment))+geom_point(size=5) + scale_y_continuous("Leaf Area Herbivorized",limits =c(0, 0.15))+ scale_x_continuous("Source Elevation")  
p_lar + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide"))



  
  