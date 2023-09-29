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
library(effects)
library(glmmTMB)



setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/field/2023_files/census_files")
#### read in data ####
data<- read.csv("inam_schofield_092623.csv")

#data<- filter(data,include == "1")


data $S_elev<-scale(data$Elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$Elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)


data$Genotype<-as.factor(data$Genotype)
data$treatment<-as.factor(data$treatment)
data$Block<-as.factor(data$Block)
data$cohort<-as.factor(data$cohort)

data$elev_dist<-data$Elevation-3133

#data$elev_dist<- abs(data$elev_dist)
data$abs_elev_dist<- abs(data$elev_dist)

data $S_elev_dist_abs<-scale(data$abs_elev_dist,center=TRUE, scale=TRUE) #needed to be scaled to avoid issues in model

data$elev_dist_km<- data$elev_dist/1000
data $S_elev_dist<-scale(data$elev_dist,center=TRUE, scale=TRUE)

#concatenate garden and treatment to plot from visreg
#data $garden_treat<-interaction(data $garden, data $treatment,sep = "_")

#Let's concatenate block and garden so that we don't have to worry about nesting block within garden

#data $block_garden<-interaction(data $block, data $garden,sep = "_")

#census 2 and 4 have the damage census

census2 <- filter(data, Census == "2")
census4 <- filter(data, Census == "4") #incomplete census

#census 7 is final census
census7 <- filter(data, Census == "7")


#### Traits ####

##SLA 

hist(data$rosette_SLA)


SLAgg= ggplot(data, aes(x= elevation,y= rosette_SLA, #group= treatment, 
                               colour= treatment))+geom_point(size=5) + scale_y_continuous("SLA CM3/g")+ scale_x_continuous("Source Elevation")  
SLAgg + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide"))



#specifying the concatenated treatment 
SLA<- lmer(rosette_SLA ~ treatment*elev_km*garden #+I(elev_km^2)*treatment
           +init_size +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = data)
Anova(SLA)
plot (SLA)

visreg(SLA, overlay = FALSE, "elev_km", by="garden", #type="conditional", #scale = "response", 
       xlab="Elevation (m)", ylab="Specific leaf area", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(0.7,0.9))),
       #line=list(col=grey(c(0,0.5))),
       # points=list(cex=1.5,col=grey(c(0.2,0.8)))
)  

SLAmod <-predictorEffect("elev_km",  partial.residuals=TRUE, SLA)
plot(SLAmod, lwd=2,xlab="Source elevation (km)", ylab="SLA", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=20, col="black"), ylim=c(0,250))

#### SLA genotypic means ####
##Now to get LSMEANS
SLAmod<- lmer (rosette_SLA~Genotype*treatment*garden+(1|block_garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data)


fam_avg_linear<-emmeans(SLAmod, ~Genotype:treatment:garden)
fam_means_linear<-as.data.frame(fam_avg_linear)

elev <- data[c("Genotype","elev_km")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_SLA <- merge(fam_means_linear,elev,by="Genotype") #merge the dataframes

#LSmeans_LAR$elev<-LSmeans_LAR$elevation/1000
write.csv(LSmeans_SLA,file="LSmeans_SLA.csv")

#had to go in and manually remove the genos that didnt have values for both treatments
# 194_2, schofield
# 267_13A, schofield
# 270_3, schofield
# 273_18, gothic

LSmeans_SLA <- read.csv("LSmeans_SLA.csv", stringsAsFactors=FALSE)

LSmeans_SLA$Genotype<-as.factor(LSmeans_SLA$Genotype)

m_SLA<- lmer (emmean~elev_km+treatment+garden+(1|Genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_SLA)

modC<- lm (emmean~elev_km*garden*treatment, data=LSmeans_LAR)

Anova(m_SLA)

visreg(m_SLA,overlay=F,"elev_km", by ="garden")

#boxplot
SLAgga<-ggplot(LSmeans_SLA, aes(x = garden, y = emmean, fill = treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("SLA CM3/g") +
  geom_point(pch = 21, position = position_jitterdodge())
SLAgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic (Low)", "Schofield (High)")) +  scale_fill_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide")
                            )


SLAgga= ggplot(LSmeans_SLA, aes(x= elev_km,y= emmean, #group= treatment, 
                         colour= treatment))+geom_point(size=5) + scale_y_continuous("SLA CM3/g")+ scale_x_continuous("Source Elevation")  
SLAgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide"))


##LMDC 

hist(data$rosette_LDMC)

ldmcgg = ggplot(data, aes(x= elevation,y=rosette_LDMC, group= treatment))+geom_point(aes(colour= treatment ),size=5) +scale_x_continuous("elevation")+ scale_y_continuous("SLA")
ldmcgg + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank())+geom_smooth(aes(colour=treatment), method="glm",size=1.6)



ldmcgg<-ggplot(data, aes(x = treatment, y = rosette_LDMC, fill = treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("SLA CM3/g") +
  geom_point(pch = 21, position = position_jitterdodge())
SLAgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Control", "Pesticide")) +  scale_fill_manual(values = c( "lightblue","darkred")) #name = "Watering treatment", labels = c("Ample water","Water-restricted"))


#specifying the concatenated treatment 
LDMC<- lmer(rosette_LDMC ~ treatment+elev_km*garden+
           +init_size +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = data)
Anova(LDMC)
plot (LDMC)

visreg(LDMC, overlay = FALSE, "elev_km", by="garden", #type="conditional", #scale = "response", 
       xlab="Elevation (m)", ylab="LDMC", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(0.7,0.9))),
       #line=list(col=grey(c(0,0.5))),
       # points=list(cex=1.5,col=grey(c(0.2,0.8)))
)  

LDMCmod <-predictorEffect("elev_km",  partial.residuals=TRUE, LDMC)
plot(LDMCmod, lwd=2,xlab="Source elevation (km)", ylab="LDMC", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=20, col="black"))#, ylim=c(0,002))



#### Fitness models: Fitness ~ garden x treatment x elev ####
#survival, reproduction, fecundity 

#survival

mod_surv<- glmer(Survival~elev_km*treatment + cohort+init_size+ (1|Block)+(1|Genotype), data = census7, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #when i add cohort to the model, it doesnt work.. 
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

mod_reproduction<- glmer(Reproduction~treatment*elev_km+ cohort+init_size+(1|Block)+(1|Genotype), data = census7, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
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

surv <- subset(census7, Survival=="1")
surv2022 <- subset(census7, cohort=="2022")


mod_reproduction1<- glmer(Reproduction~treatment*elev_km+init_size (1|Block)+(1|Genotype), data = surv2022, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)
Anova(mod_reproduction1,type="III")  
plot(mod_reproduction1)

visreg(mod_reproduction1, overlay = TRUE, "cohort", type="conditional", scale = "response", 
       xlab="Source Elevation (Km)", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = TRUE)



reproduction <-predictorEffect("elev_km",  partial.residuals=TRUE, mod_reproduction1)
plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))

mod_reproduction2<- glmer(Flowered~garden_treat*elev_km+init_size+ (1|block)+(1|Genotype), data = surv, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)
Anova(mod_reproduction2,type="III")  
plot(mod_reproduction2)

visreg(mod_reproduction2, overlay = TRUE, "elev_km", by="garden_treat", type="conditional", scale = "response", 
       xlab="Source elevation", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       #fill=list(col="blue"),
       #line=list(col=grey(c(0,0.8))),
       #points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = FALSE)


#num fruits
#not working... low numbers? only had a handful reproduce in the first place


fecundA<- glmer(Mature_length_siliques ~  init_size +elev_km* treatment+elev_km*garden+ (1|block) +(1|population), data= filter(surv, Mature_length_siliques > 0),family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecundA,type="III")

fruit <-predictorEffect("elev_km",  partial.residuals=TRUE,fecundA)
plot(fruit, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduction)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1200))


visreg(fecundA, overlay = TRUE, "treatment", by="treatment", type="conditional", scale = "response", 
       xlab="Source elevation", ylab="Fitness (Prob Flowering)", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       #fill=list(col="blue"),
       #line=list(col=grey(c(0,0.8))),
       #points=list(cex=1.5,col=grey(c(0,0.8))),
       jitter = FALSE)



library(glmmTMB)
library(glmmTMB) #load it twice 

hurdle_Model = glmmTMB (Mature_silique_number ~ treatment * elev_km + garden  + (1| block)+(1|Genotype), zi=~ treatment* elev_km *garden  + (1| block)+(1|Genotype),data = surv ,family=truncated_nbinom2)

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


#### local adaptation models: fitness ~ T x elevational distance (source – transplant) ####
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


mod_fruit<- lmer( Mature_silique_number ~ treatment* S_elev_dist_abs+I(S_elev_dist_abs^2)*treatment +cohort+(1|Genotype) ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = flowered)

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

mod_repro<- glmer(Reproduction~treatment* elev_dist_km+I(elev_dist_km^2)*treatment +(1|Block)+(1|Genotype), data = census7, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) 

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

mod_surv<- glmer(Survival~treatment* elev_dist_km+I(elev_dist_km^2)*treatment +cohort+(1|Block)+(1|Genotype), data = census7, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) 

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
mod_surv<- glmer(Survival~treatment+S_elev_dist +(1|Block)+(1|Genotype), data = data, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial) #best model, gardenXelevation interaction only

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

#need to update with specific censuses
census2$LAR_prop<- census2$LAR/100

##Check that the conversion worked

hist(census2$LAR_prop)


#remove NA values for gamlss


data_LAR<-census2 %>% dplyr::select(treatment, S_elev, elev_km,Genotype,LAR_prop,Block)

data_LAR<-na.omit(data_LAR)


##Okay, it looks like there were no NAs in that column, so we can proceed to the zero-one inflated negative binomial regression. What I've written out here is the full model, which has convergence issues.

#LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype)
gamlss.modelA<- gamlss (formula=LAR_prop~treatment*S_elev+ random(Block)+random(Genotype), sigma.formula=LAR_prop~treatment*S_elev+ random(Block)+random(Genotype), nu.formula=LAR_prop~treatment*S_elev+ random(Block)+random(Genotype), tau.formula=LAR_prop~treatment*S_elev+ random(Block)+random(Genotype), family=BEINF, data= data_LAR)


##Perhaps we don't want to interact everything. Please keep in mind that garden should be a fixed effect, not random (you only have 2 gardens), but maybe the patterns don't differ much across gardens:



#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.model<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), tau.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)

#many many warnings

gamlss.modelA
plot(gamlss.modelA)
summary(gamlss.modelA) #nu only has garden as significant

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.modelA)
mod2<-stepGAIC(gamlss.modelA)
mod2$anova

summary(mod2) # there are significant effects of treatment, elevation of origin and their interaction for mu (the probability of being damaged). For nu (amount of damage on plants with >0 and <1 damage; i.e., the beta component), there is a significant effect of elevation of origin and garden, but not treatment or elevation by treatment). For tau (probability of 100% damage – which you likely don’t have), there are no significant effects.


#So, we might be able to remove the tau component from the model, and run a zero-inflated model:
  
  gamlss.modelB<- gamlss (formula=LAR_prop~treatment+S_elev+ random(Block)+random(Genotype), 
                          sigma.formula=LAR_prop~treatment+S_elev+ random(Block)+random(Genotype), 
                          nu.formula=LAR_prop~treatment+S_elev+ random(Block)+random(Genotype), family=BEINF, data= data_LAR)

plot(gamlss.modelB)

summary(gamlss.modelB)
  
  
visreg(gamlss.modelB, "treatment", type="conditional", #scale = "response", 
         xlab="Elev", ylab="LAR", partial=TRUE,
         fill=list(col="light grey"
                   #(c(0.99), alpha=0)
         ), band = FALSE,
         #line=list(col=grey(c(0.2,0.6))),
          points=list(cex=0.65,  pch=(19)))  

c2021 $garden_treat<-interaction(c2021 $treatment, c2021 $garden,sep = "_")
  
gamlss.modelE<- gamlss (formula=LAR_prop~S_elev+garden_treat + random(block_garden)+random(Genotype), sigma.formula=LAR_prop~S_elev+ garden_treat+ random(block_garden)+random(Genotype), nu.formula=LAR_prop~S_elev+ garden_treat+ random(block_garden)+random(Genotype), family=BEINF, data= data_LAR)
  
  plot(gamlss.modelE)
  
  summary(gamlss.modelE)
  
  visreg(gamlss.modelE, overlay = FALSE, "S_elev", by="garden_treat", type="conditional", 
         #scale = "response",   
         xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = FALSE)    


  #You can create a datafile with the predicted data from the nu component using this code for plotting:
  pred<-predict(gamlss.modelB, newdata=data_LAR, type="response", what="nu"
                )
  
  pred<-predict(gamlss.modelB, type="response", what="mu")
  data_LAR <- cbind(data_LAR,pred)
  

  str(gamlss.modelB)
  
  pred <- as.data.frame(pred)  
  View(pred)
  
  pred_data<- bind_cols(data_LAR, pred)
  
  data_LAR$LAR <- (data_LAR$pred*100)
  
  LAR_schofield = ggplot(data_LAR, aes(x= elev_km,y=LAR))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Percent Leaf area herbivorized") + geom_point(aes(colour = factor(treatment)), size = 4)
  P <-LAR_schofield + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                           panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="lm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="lm",size=1.6, se=T)
P

ggplot(data_LAR, aes(x = elev_km, y = LAR, color = treatment))+
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Source elevation", y = "Leaf Area Removed", color = "treatment") + 
  theme_bw() + scale_color_manual(values=c("#009E73","#CC79A7"))


#gamlss models are not converging!


#### linear model, this is the incorrect fit  but nothing else is working ####
##Now to get LSMEANS
modB<- lmer (LAR_prop~Genotype*treatment+cohort+init_size+(1|Block),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=census2)

plot(modB)
Anova(modB)

fam_avg_linear<-emmeans(modB, ~Genotype:treatment)
fam_means_linear<-as.data.frame(fam_avg_linear)

elev <- data[c("Genotype","Elevation")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_LAR <- merge(fam_means_linear,elev,by="Genotype") #merge the dataframes

LSmeans_LAR$elev<-LSmeans_LAR$Elevation/1000
write.csv(LSmeans_LAR,file="LSmeans_LAR.csv")


write.csv(fam_means_linear,file="LSmeans_LAR.csv")

#had to go in and manually remove the genos that didnt have values for both treatments
#185_2, gothic
# 194_2, schofield
# 267_13A, schofield
# 270_3, schofield

LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=FALSE)
LSmeans_LAR <-na.omit(LSmeans_LAR)

LSmeans_LAR $S_elev<-scale(LSmeans_LAR$elevation,center=TRUE, scale=TRUE)


LSmeans_LAR$Genotype<-as.factor(LSmeans_LAR$Genotype)

modC<- lmer (emmean~elev*treatment+(1|Genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)


Anova(modC)

visreg(modC,overlay=T,"elev", by ="treatment")


#boxplot
LARgga<-ggplot(LSmeans_LAR, aes(x = garden, y = emmean, fill = treatment,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Leaf Area Herbivorized",limits =c(0, 0.15)) +
  geom_point(pch = 21, position = position_jitterdodge())
LARgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic (Low)", "Schofield (High)")) +  scale_fill_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide")
                                                                                                                                                                           )
# scatter plot
p_lar= ggplot(LSmeans_LAR, aes(x= elev_km,y= emmean, #group= treatment, 
                               colour= treatment))+geom_point(size=5) + scale_y_continuous("Leaf Area Herbivorized",limits =c(0, 0.15))+ scale_x_continuous("Source Elevation")  
p_lar + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ garden, scales="free_x") +scale_colour_manual(values = c( "lightgreen","pink"), name = "Herbivory treatment", labels = c("Control","Pesticide"))
  

