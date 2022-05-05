library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)
packageVersion("lme4")
library(tidyr)
library(broom)
library(car)
library(MASS)
library(emmeans)
library(visreg)
library(fitdistrplus)
library(gamlss)
library(effects)

setwd("~/Desktop/Anderson_data/herbivory/data/feeding_assay/Feeding_trial_spring2021")

#import data
Feeding_trial <- read.csv("Feeding_trialNov2021a.csv")

##Convert elevation to km (it helps with model convergence)
Feeding_trial$elev<-Feeding_trial$elevation/1000
Feeding_trial$S_weight<-scale(Feeding_trial$ini_insect_weight.mcg.,center=TRUE, scale=TRUE)
Feeding_trial$time<-scale(Feeding_trial$elapsed_time,center=TRUE, scale=TRUE) #need to include elapsed time in the models

#Change treatment, genotype, batch to factor
Feeding_trial$mat_treat<-as.factor(Feeding_trial$mat_treat)
Feeding_trial$genotype <-as.factor(Feeding_trial$genotype)
Feeding_trial$batch<-as.factor(Feeding_trial$batch)
Feeding_trial$mat_exp_ID <-as.factor(Feeding_trial$mat_exp_ID) #need to include this as random effect since multiple reps per mat line

#checking the distribution 
#shapiro.test(Feeding_trial$RGR)
#qqnorm(Feeding_trial$RGR)
#hist(Feeding_trial$RGR_min)

#hist(Feeding_trial$LAR)
#shapiro.test(Feeding_trial$LAR) #significant, distribution is not normal


#### LAR ####
#looking at LAR by elevation
plot(Feeding_trial$elev, Feeding_trial$LAR, xlab="elev", ylab="LAR")


#gamlss model, need to remove NAs
data_LAR <- drop_na(Feeding_trial,LAR) #this removes the rows without LAR values
data_LAR <-na.omit(data_LAR)


gamlss_mod1<- gamlss (formula=LAR~S_weight+elev*mat_treat+I(elev^2)*mat_treat+ random(genotype)+ random(batch)+random(mat_exp_ID), family=BEINF, data= data_LAR)

visreg(gamlss_mod1, 'elev', by= "mat_treat", overlay = FALSE, type="conditional", 
       #scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8)))) 

gamlss.modelB<- gamlss (formula= LAR ~ S_weight + time+elev*mat_treat+I(elev^2)*mat_treat+ random(genotype)+ random(batch)+random(mat_exp_ID), nu.formula=LAR ~ S_weight + time+elev*mat_treat+I(elev^2)*mat_treat+random(genotype)+random(batch)+random(mat_exp_ID), tau.formula=LAR ~ S_weight + time+elev*mat_treat+I(elev^2)*mat_treat+ random(batch)+random(genotype)+random(mat_exp_ID), family=BEINF, data= data_LAR) #BEINF family 

Anova(gamlss_mod1)

visreg(gamlss.modelB, overlay = TRUE, "elev", by="mat_treat", type="conditional", scale = "linear", 
       xlab="Source Elevation (KM)", ylab="Probability of reproduction", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(lty=1, col=grey(c(0.2,0.6))),
       points=list(cex=0.65, col=grey(c(0.3,0.6)),  pch=(19)))

modelB <-predictorEffect("elev",  partial.residuals= TRUE, gamlss.modelB)
plot(modelB, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))

plot(gamlss.modelB)
summary(gamlss.modelB)

data_LAR$elev2<-(data_LAR$elev)*(data_LAR$elev)

gamlss.modelC<- gamlss (formula= LAR ~ S_weight + time+elev*mat_treat+ elev2*mat_treat+ random(genotype)+random(mat_exp_ID),nu.formula=LAR ~ S_weight + time+elev*mat_treat+ elev2*mat_treat+random(genotype)+random(mat_exp_ID), tau.formula=LAR ~ S_weight + time+elev*mat_treat+ elev2*mat_treat+ random(genotype)+random(mat_exp_ID), family=BEINF, data= data_LAR)

plot(gamlss.modelC)
summary(gamlss.modelC)



#The scatter is very high if you just plot the relationship with ggplot:
  
p_lar= ggplot(data_LAR, aes(x= elev,y= LAR, group= mat_treat, colour= mat_treat))+geom_point(size=5) +scale_x_continuous("elevation")+ scale_y_continuous("Leaf damage") 
p_lar + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ mat_treat, scales="free_x")

#LSMEANS: dont have to run again, already done
  
#modB<- lmer(LAR ~ genotype*mat_treat+S_weight+time +(1|mat_exp_ID)+(1|batch),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = Feeding_trial)

#fam_avg_linear<-emmeans(modB, ~genotype:mat_treat)
#fam_means_linear<-as.data.frame(fam_avg_linear)

#write.csv(fam_means_linear,file="LSmeans_LAR.csv")
#LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=TRUE)

#now we have lsmeans, but need to add corresponding elevation
#elev <- data_LAR[c("genotype","elevation")] #make dataframe of genotypes and corresponding elev
#elev <- unique(elev) #calls unique rows 
#LSmeans_LAR <- merge(LSmeans_LAR,elev,by="genotype") #merge the dataframes

#LSmeans_LAR$elev<-LSmeans_LAR$elevation/1000
write.csv(LSmeans_LAR,file="LSmeans_LAR.csv")
LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=TRUE)

LSmeans_LAR1 <-na.omit(LSmeans_LAR)

gamlss_mod1<- gamlss (formula=emmean~mat_treat*elev+mat_treat*I(elev^2)+ random(genotype), family=BE, data= LSmeans_LAR1)

visreg(gamlss_mod1, 'elev', by= "mat_treat", overlay = FALSE, type="conditional", 
       #scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0))),
       points=list(cex=1.5,col=grey(c(0)))
       ) 

modC<- lmer (emmean~mat_treat*elev+mat_treat*I(elev^2)
             +(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)

Anova(modC,type="III")
plot(modC)

model1 <-predictorEffect("elev",  partial.residuals= TRUE, gamlss_mod1)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))


##For a moment, let's treat LAR as a normal variable and run a linear mixed model

modA<- lmer (LAR~elev*mat_treat+I(elev^2)*mat_treat+S_weight+(1|batch)+(1|genotype/mat_exp_ID),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)

Anova(modA,type="III")
plot(modA)

model1 <-predictorEffect("elev",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))


# Generalized additive model
library(mgcv)
mod_1a <- gam(LAR ~mat_treat+ini_insect_weight.mcg.+s(elev, by=mat_treat)+s(genotype,bs="re")+s(batch,bs="re"),data= data_LAR, method="REML")

plot(mod_1a,pages=1,scheme=1,residuals=TRUE,unconditional=TRUE) 
plot(mod_1a, shade = TRUE, pages = 1, scale = 0)
summary(mod_1a) #batch is significant here.. 

visreg(mod_1a, overlay = FALSE, "elev", by="mat_treat", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

mod_1b <- gam(emmean ~mat_treat+s(elev_km, by=mat_treat)+s(genotype,bs="re"),data= LSmeans_LAR , method="REML")
#error term, model has more coeff than data


#### RGR ####

#trying to transform RGR, did not work out that well
Feeding_trial$t_RGR=sqrt(Feeding_trial$RGR+4.701560096)
hist(Feeding_trial$t_RGR)

##Insect mortality, this analysis was skipped for the fall 2021 assay since all insects survived
#mort <- glmer(insect_alive~elev*mat_treat + (1|genotype), data= Feeding_trial, family=binomial(link=logit)) 
##Anova(mort,type="III")
#  visreg(mort,"elev", by="mat_treat", overlay=TRUE, scale = "response", xlab=" Source elevation", ylab="Probability of insect survival")
  
##Relative growth rate
# initial insect weight was significant
mod_1<- lmer(RGR ~ elev+mat_treat + ini_insect_weight.mcg. + (1|genotype) +(1|batch), control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = Feeding_trial)
Anova(mod_1,type="III")
plot(mod_1)


##Generalized additive model
library(mgcv)
mod_1a <- gam(RGR ~mat_treat+ini_insect_weight.mcg.+s(elev, by=mat_treat)+s(genotype,bs="re")+s(batch,bs="re"),data= Feeding_trial, method="REML")
## use alternative plotting scheme, and way intervals include
## smoothing parameter uncertainty...
plot(mod_1a,pages=1,scheme=1,unconditional=TRUE) 
plot(mod_1a,pages=1,scheme=1,residuals=TRUE,unconditional=TRUE) 
summary(mod_1a)
plot(mod_1a, shade = TRUE, pages = 1, scale = 0)


##lsmeans for RGR
modF<- lmer(RGR ~ genotype*mat_treat+ini_insect_weight.mcg. +(1|batch),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = Feeding_trial)

fam_avg_linear<-emmeans(modF, ~genotype:mat_treat)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_RGR.csv")
#manually added the elevation to the file for the genos

LSmeans_RGR <- read.csv("LSmeans_RGR.csv", stringsAsFactors=TRUE)
LSmeans_RGR$elev_km<-LSmeans_RGR$elevation/1000
LSmeans_RGR <- subset(LSmeans_RGR, mat_treat!="neg_cont" & mat_treat!="post_cont")#gets rid of the postive and neg control rows

modG<- lm(emmean ~ elevation,data=LSmeans_RGR)

visreg(modG,"elev_km", type="contrast", scale = "response", xlab=" Source Elevation (km)", ylab="Insect Relative Growth Rate", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)

modA<- lmer (emmean~elev_km*mat_treat+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_RGR)
Anova(modA,type="III")
summary(modA)
plot(modA)

model1 <-predictorEffect("elev_km",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))


#Ploting RGR
library(ggpubr)


RGR2 = ggplot(LSmeans_RGR, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Insect Relative Growth Rate") + geom_point(aes(colour = mat_treat), size = 4)
R <- RGR2 + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = mat_treat),method="glm",size=2.0, se=F)+geom_smooth(aes(group = mat_treat, colour = mat_treat),method="lm",size=1.6, se=F)
print(R)



### final insect weight ####
mod_2<- lmer(final_insect_weight ~ elev*mat_treat + ini_insect_weight.mcg. + (1|genotype) +(1|batch), control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = Feeding_trial)
Anova(mod_1,type="III")
plot(mod_1)


visreg(mod_1, "elev", by="mat_treat")
transgen<-predictorEffects(mod_1, partial.residuals=T)
plot(transgen,type="response")
transgen <-predictorEffect("elev",  partial.residuals=TRUE, mod_1)
plot(transgen, lwd=2,xlab="Source elevation (KM)", ylab="Insect relative growth rate", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))
##lsmeans for final weight
modH<- lmer(final_insect_weight ~ genotype*mat_treat+ini_insect_weight.mcg. +(1|batch),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = Feeding_trial)

fam_avg_linear<-emmeans(modG, ~genotype:mat_treat)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_FW.csv")
#manually added the elevation to the file for the genos

LSmeans_FW <- read.csv("LSmeans_FW.csv", stringsAsFactors=TRUE)
LSmeans_FW$elev_km<-LSmeans_FW$elevation/1000
LSmeans_FW <- subset(LSmeans_FW, mat_treat!="neg_cont" & mat_treat!="post_cont")#gets rid of the postive and neg control rows

modI<- lm(emmean ~ mat_treat*elev_km,data=LSmeans_FW)

visreg(modI,"elev_km", type="contrast", scale = "response", xlab=" Source Elevation (km)", ylab="Final insect weight (mcg)", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)


modA<- lmer (emmean~elev_km*mat_treat+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_FW)
Anova(modA,type="III")
summary(modA)
plot(modA)

model1 <-predictorEffect("elev_km",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))


#Ploting FW
library(ggpubr)


FW2 = ggplot(LSmeans_FW, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Final insect weight (mcg)") + geom_point(aes(colour = mat_treat), size = 4)
F <- FW2 + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                         panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = mat_treat),method="glm",size=2.0, se=F)+geom_smooth(aes(group = mat_treat, colour = mat_treat),method="lm",size=1.6, se=F)
print(F)








