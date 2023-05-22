#script for the 2021 greenhouse experiment for summary file

#libraries
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
library(rjags)
library(gamlss)
library(effects)
 
setwd("~/Desktop/Anderson_data/herbivory/data/greenhouse/")


#read in data 
gh2 <- read.csv("GHSummary2021.csv")
#gh2 <- drop_na(gh2,LAR)
gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)

##convert LAR to proportion
  gh2$LAR_prop<- gh2$LAR13/100
  hist(gh2$LAR13)
  
#Check that the conversion worked
  hist(gh2$LAR_prop)
  #head(gh2)  

##Some of your variables are being read as characters not factors. Let's fix that
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line

#plasticity

#SLA





modB<- lmer(SLA ~ genotype:mat_treat:treatment +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = gh2)
fam_avg_linear<-emmeans(modB, ~genotype:mat_treat:treatment)
fam_means_linear<-as.data.frame(fam_avg_linear)
write.csv(fam_means_linear,file="LSmeans_SLA.csv")



LSmeans_SLA <- read.csv("LSmeans_SLA.csv", stringsAsFactors=TRUE)
#now we have lsmeans, but need to add corresponding elevation
#elev <- gh2[c("genotype","elevation")] #make dataframe of genotypes and corresponding elev
#elev <- unique(elev) #calls unique rows 
#LSmeans_SLA <- merge(LSmeans_SLA,elev,by="genotype") #merge the dataframes

#LSmeans_SLA$elev<-LSmeans_SLA$elevation/1000
#LSmeans_SLA$emmean<-LSmeans_SLA$emmean*100

#write.csv(LSmeans_SLA,file="LSmeans_SLA.csv")




hist(LSmeans_SLA$emmean)

SLAgga<-ggplot(LSmeans_SLA, aes(x = mat_treat, y = emmean, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Maternal treatment")+ scale_y_continuous("SLA CM3/g") +
  geom_point(pch = 21, position = position_jitterdodge())
SLAgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("No Herbivores", "Herbivorized", "Parent Control")) +  scale_fill_manual(values = c( "gold","purple"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))



gh_sla= ggplot(LSmeans_SLA, aes(x= elev,y= emmean, group= treatment, 
                                colour= treatment))+geom_point(size=5) + scale_y_continuous("SLA CM3/g")+ scale_x_continuous("Source Elevation")  
gh_sla + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "gold","purple","grey"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))




modB<- lm(emmean ~ elevation*treatment+mat_treat ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = LSmeans_SLA)
Anova(modB)


modB<- lm(emmean ~ elevation*treatment*mat_treat, data = LSmeans_SLA)
Anova(modB)

##use this
visreg(modB,"elevation", by="treatment", overlay=TRUE, scale = "response", xlab="Source elevation", ylab="Specific leaf area", line=list(col="black"),partial=FALSE,points=list(cex=1.2, col="black"))

visreg(modB, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="leaf area", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(0.7,0.9))),
       line=list(col=grey(c(0,0.5))),
       points=list(cex=1.5,col=grey(c(0.2,0.8))))  

SLA_1<- lmer(SLA~elev*treatment*mat_treat+ (1|block)+(1|genotype), data = gh2)
Anova(SLA_1,type="III")

visreg(SLA_1, overlay = TRUE, "elev", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="leaf area", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(0.7,0.9))),
       line=list(col=grey(c(0,0.5))),
       points=list(cex=1.5,col=grey(c(0.2,0.8))))  

SLA_p1 <-predictorEffect("elev",  partial.residuals=FALSE, SLA_1)

plot(SLA_p1, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))



#plasticity figure
ggplot(
  #subset(
    LSmeans_SLA,
   # treatment == "Herbivorized"),
  aes(x=treatment, y=emmean, color = elevation)) + # Change fill to color
  theme_classic(base_size = 22) + 
  #geom_point() + 
  stat_summary(fun=mean, position = "dodge") + 
  #stat_summary(
  #  geom="errorbar", 
  #  fun.data= mean_cl_boot,
  #  width = 0.1, size = 0.2, col = "grey57"
  #  ) + 
  # Lines by species using grouping
  #scale_color_viridis(discrete = FALSE) +
  #scale_fill_manual(values = cbp2)+
  stat_summary(aes(group = genotype), geom = "line", size=1.5, fun = mean) +
  ylab("SLA") +
  xlab("maternal trt")

herb <- subset(gh2,mat_treat=="herb")
cont <- subset(gh2,mat_treat=="cont")
parent <- subset(gh2,mat_treat=="Parent")

ggplot(parent,aes(x=treatment, y=SLA, color = elevation)) + # Change fill to color
  theme_classic(base_size = 22) + 
  #geom_point() + 
  stat_summary(fun=mean, position = "dodge") + 
  #stat_summary(
  #  geom="errorbar", 
  #  fun.data= mean_cl_boot,
  #  width = 0.1, size = 0.2, col = "grey57"
  #  ) + 
  # Lines by species using grouping
  #scale_color_viridis(discrete = FALSE) +
  #scale_fill_manual(values = cbp2)+
  stat_summary(aes(group = elevation), geom = "line", size=1.5, fun = mean) +
  ylab("Leaf Area Removed") +
  xlab("maternal trt")


##use this
visreg(modB,"elevation", by="treatment", overlay=FALSE, scale = "response", xlab="Source elevation", ylab="Probability of reproduction", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))

plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     +           partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))



#subset parents
p_SLA <- subset(LSmeans_SLA, mat_treat=="Parent")

modC<- lm(emmean ~ elev+treatment ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = p_SLA)
Anova(modC)

#leaf area
modLA<- lmer(leaf_area ~ genotype:mat_treat:treatment +(1|block),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = gh2)
fam_avg_linear<-emmeans(modLA, ~genotype:mat_treat:treatment)
fam_means_linear<-as.data.frame(fam_avg_linear)
write.csv(fam_means_linear,file="LSmeans_leafarea.csv")
LSmeans_leafarea<- read.csv("LSmeans_leafarea.csv", stringsAsFactors=TRUE)
#now we have lsmeans, but need to add corresponding elevation
elev <- gh2[c("genotype","elevation")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_leafarea <- merge(LSmeans_leafarea,elev,by="genotype") #merge the dataframes
LSmeans_leafarea$elev<-LSmeans_leafarea$elevation/1000

write.csv(LSmeans_leafarea,file="LSmeans_leafarea.csv")

modF<- lm(emmean ~ elevation*treatment*mat_treat ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data =LSmeans_leafarea)

Anova(modF)

visreg(modF, overlay = FALSE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="leaf area", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="grey"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0.2,0.8,5))))  

#subset 
h_LA <- subset(LSmeans_leafarea, mat_treat=="herb")

modG<- lm(emmean ~ elevation+treatment ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data =h_LA)
Anova(modG)

visreg(modG, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8)))) 

c_LA <- subset(LSmeans_leafarea, mat_treat=="cont")

modH<- lm(emmean ~ elevation+treatment ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data =c_LA)
Anova(modH)

visreg(modH, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8)))) 

p_LA <- subset(LSmeans_leafarea, mat_treat=="Parent")

modI<- lm(emmean ~ elevation+treatment ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data =p_LA)
Anova(modI)

visreg(modH, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8)))) 


p_leafarea= ggplot(gh2, aes(x= treatment,y= leaf_area, group= mat_treat, colour= mat_treat))+geom_point(size=5)+ #+scale_x_continuous("elevation")
 scale_y_continuous("leaf_area") 
p_leafarea + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ mat_treat, scales="free_x")







#lsmeans for prob flowering, not working

modB<- glmer(flowered~ genotype*mat_treat*treatment+ (1|block), data = gh2, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(modB,type="III")

fam_avg_linear<-emmeans(modB, ~genotype:mat_treat:treatment)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_flowered.csv")
LSmeans_flowered <- read.csv("LSmeans_flowered.csv", stringsAsFactors=TRUE)

#now we have lsmeans, but need to add corresponding elevation
elev <- gh2[c("genotype","elevation")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_flowered <- merge(LSmeans_flowered,elev,by="genotype") #merge the dataframes

LSmeans_flowered$elev<-LSmeans_flowered$elevation/1000

LSmeans_flowered <- data.frame(LSmeans_flowered)

mod_pr<- glmer(emmeans~ elev*mat_treat+treatment+ (1|genotype), data = LSmeans_flowered, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_pr,type="III")

#gh2 <- drop_na(gh2,survival) #removes the plants that died before the experiment #manually removed them. had master note

#survival

mod_surv<- glmer(survival~treatment*S_elev+I(S_elev^2)*treatment +(1|block)+(1|genotype), data = gh2, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_surv,type="III") #scaled initial size is significant
plot(mod_surv) # wacky residuals

survived <-predictorEffect("S_elev",  partial.residuals=FALSE, mod_surv)
plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
          partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))





mod_surv1<- glmer(survival~mat_treat+S_init_size+ (1|block)+(1|genotype), data = gh2, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_surv1,type="III") 

#no parent
mod_surv2<- glmer(survival~S_elev*mat_treat+I(S_elev^2)+S_init_size+ (1|block)+(1|genotype), data = subset(gh2, mat_treat!="Parent"), control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial) 

Anova(mod_surv2,type="III") #S_elev squared is significant, but does not look real

survived <-predictorEffect("S_elev",  partial.residuals=TRUE, mod_surv2)
plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


#Prob_reproduction
mod_pr<- glmer(flowered~elev*mat_treat*treatment+ (1|block)+(1|genotype), data = gh2, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_pr,type="III")

reproduction <-predictorEffect("elev",  partial.residuals=FALSE, mod_pr)

plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

visreg(mod_pr,"elev", by="mat_treat", overlay=FALSE, scale = "response", xlab="Source elevation", ylab="Probability of reproduction", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))


prob_repro <- glmer(flowered ~ treatment*mat_treat*elev  + (1| block)+(1|genotype),data = gh2, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

Anova(prob_repro,type="III")

visreg(prob_repro,"elev", by="mat_treat", overlay=FALSE, scale = "response", xlab="Source elevation", ylab="Probability of reproduction", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))

#redoing rosemary figure with most recent census

#treatment is added rather than having the interaction in model since the interactions with it were not significant
#no parents
mod_pr1<- glmer(flowered~S_init_size +elev*mat_treat+treatment+(1|block)+(1|genotype), data = subset(gh2, mat_treat!="Parent"), control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_pr1,type="III")

reproduction6 <-predictorEffect("elev",  partial.residuals=TRUE, mod_pr1)

plot(reproduction6, lwd=2,xlab="Source Elevation (KM)", ylab="Probability of reproduction", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

drop1(mod_pr1, ~elev*mat_treat, test="Chisq")

mod_pr2<- glm(flowered~ elevation*mat_treat, data = subset(gh2, mat_treat!="Parent"))
Anova(mod_pr2,type="III")


reproduction7 <-predictorEffect("elevation",  partial.residuals=FALSE, mod_pr2)

plot(reproduction7, lwd=2,xlab="Source Elevation (KM)", ylab="Probability of reproduction", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

prob_repro <- glmer(repro ~ treatment*mat_treat*elev+S_init_size  + (1| block)+(1|genotype),data = gh2, family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

Anova(prob_repro,type="III")



visreg(prob_repro,"elev", by="mat_treat", overlay=FALSE, scale = "response", xlab="Source elevation", ylab="Probability of reproduction", line=list(col="black"),partial=TRUE,points=list(cex=1.2, col="black"))



#siliques
mod_slq<- lmer (Mature_silique_number~treatment*mat_treat*elev+S_init_size+(1|block)+(1|genotype/mat_exp_ID),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=gh2)
plot(mod_slq) #indicated poor fit, zero inflated

hist(gh2$Mature_silique_number)
gh_success<-subset(gh2,Mature_silique_number!="0")
hist(gh_success$Mature_silique_number)
Anova(mod_slq,type="III")

##jill email about hurdle model


##Hurdle model with nbin

library(glmmTMB)
hurdle_Model = glmmTMB (Mature_silique_number ~ treatment*mat_treat*elev+S_init_size  + (1| block)+(1|genotype/mat_exp_ID), zi=~ treatment*mat_treat*elev+S_init_size  + (1| block)+(1|genotype/mat_exp_ID),data = gh2 ,family=truncated_nbinom2)

diagnose(hurdle_Model)

check_zeroinflation(hurdle_Model)

check_overdispersion(hurdle_Model)

Anova(hurdle_Model,type="III", component="zi")

summary(hurdle_Model)

Anova(hurdle_Model,component="cond")

Anova(hurdle_Model,component="zi")

##

visreg(mod_slq, "treatment",xlab="Treatment", ylab="LAR", scale = "linear", type="conditional")

slq1 <-predictorEffect("treatment",  partial.residuals=FALSE, mod_slq)
plot(mod_slq, lwd=2,xlab="Source Elevation (KM)", ylab="Siliqe Number", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

visreg(mod_slq, overlay = TRUE, "elev", by="mat_treat", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Silique Number", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))))  
library(sjPlot)
plot_model(mod_slq, type = "pred", terms = c("elev", "mat_treat", "treatment"))


#LAR


#linear
mod_LAR<- lmer (LAR_prop~treatment+mat_treat+elev+(1|block)+(1|genotype/mat_exp_ID),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)
plot(mod_LAR)

Anova(mod_LAR,type="III")
visreg(mod_LAR, "treatment",xlab="Treatment", ylab="LAR", scale = "linear", type="conditional")



#gamlss model, need to remove NAs
data_LAR <- drop_na(c_13,LAR) #this removes the rows without LAR values
data_LAR <-na.omit(data_LAR)
sum(c_13$status == "dead")

gamlss.modelLAR<- gamlss (formula= LAR_prop ~ init_size + elev+mat_treat+treatment + random(genotype)+ random(block)+random(mat_exp_ID), nu.formula=LAR ~ init_size + elev+mat_treat+treatment + random(genotype)+ random(block)+random(mat_exp_ID), tau.formula=LAR ~ init_size + elev+mat_treat+treatment + random(genotype)+ random(block)+random(mat_exp_ID), family=BEINF, data= data_LAR) #BEINF family 

plot(gamlss.modelLAR)
summary(gamlss.modelLAR)

visreg(gamlss.modelLAR, "elev", by="treatment",xlab="Sourch Elevation", ylab="LAR", scale = "response", type="conditional")


visreg(gamlss.modelB, overlay = FALSE, "elev", by="mat_treat", type="conditional", scale = "linear", 
       xlab="Source Elevation (KM)", ylab="Probability of reproduction", partial=FALSE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(lty=1, col=grey(c(0.2,0.6))),
       points=list(cex=0.65, col=grey(c(0.3,0.6)),  pch=(19)))

data_LAR$elev2<-(data_LAR$elev)*(data_LAR$elev)

gamlss.modelC<- gamlss (formula= LAR ~ init_size + elev*mat_treat+ elev2*mat_treat+ treatment+ random(genotype)+random(mat_exp_ID),nu.formula=LAR ~ init_size +elev*mat_treat+ elev2*mat_treat+ treatment+ random(genotype)+random(mat_exp_ID), tau.formula=LAR ~ init_size + elev*mat_treat+ elev2*mat_treat+ treatment+ random(genotype)+random(mat_exp_ID), family=BEINF, data= data_LAR)

plot(gamlss.modelC)
summary(gamlss.modelC)


##### back to original 


c_7<-subset(gh2, census=="7")
c_7h<-subset(c_7, treatment=="Herbivorized")
c_7c<-subset(c_7, treatment=="Control")


mod_7c2<- glmer(reproduction~init_size +S_elev*mat_treat*treatment+I(S_elev^2)+ (1|block)+(1|genotype), data = c_7, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7c2,type="III")

reproduction2 <-predictorEffect("S_elev",  partial.residuals=TRUE, mod_7c2)

plot(reproduction2, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


mod_7c3<- glmer(reproduction~init_size +S_elev*mat_treat+treatment+(1|block)+(1|genotype), data = c_7, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7c3,type="III")

reproduction3 <-predictorEffect("S_elev",  partial.residuals=FALSE, mod_7c3)

plot(reproduction3, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


mod_7c4<- glmer(reproduction~mat_treat*treatment+(1|block)+(1|genotype), data = c_7, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7c4,type="III")

reproduction4 <-predictorEffect("treatment",  partial.residuals=FALSE, mod_7c4)

plot(reproduction4, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", type="response", ylim=c(0,1))



##No parentals
c_7b<-subset(c_7, mat_treat!="Parent")
mod_7c5<- glmer(reproduction~init_size +S_elev*mat_treat*treatment+ (1|block)+(1|genotype), data = c_7b, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7c5,type="III")

reproduction5 <-predictorEffect("S_elev",  partial.residuals=TRUE, mod_7c5)

plot(reproduction5, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


###THIS IS THE FIGURE IN THE ROSEMARY


mod_7c6<- glmer(reproduction~init_size +elev*mat_treat+treatment+(1|block)+(1|genotype), data = c_7b, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7c6,type="III")

reproduction6 <-predictorEffect("elev",  partial.residuals=FALSE, mod_7c6)

plot(reproduction6, lwd=2,xlab="Source Elevation (KM)", ylab="Probability of reproduction", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))

drop1(mod_7c6, ~S_elev*mat_treat, test="Chisq")

visreg(mod_7c6, overlay = FALSE, "elev", by="mat_treat", type="contrast", scale = "response", 
       xlab="Source Elevation (KM)", ylab="Probability of reproduction", partial=FALSE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(lty=1, col=grey(c(0.2,0.6))),
       points=list(cex=0.65, col=grey(c(0.3,0.6)),  pch=(19)))

mod_7c7<- glm(reproduction~ elev*mat_treat, data = c_7b)
Anova(mod_7c7,type="III")

visreg(mod_7c7, overlay = FALSE, "elev", by="mat_treat", type="contrast", scale = "response", 
       xlab="Source Elevation (KM)", ylab="Probability of reproduction", partial=FALSE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(lty=1, col=grey(c(0.2,0.6))),
       points=list(cex=0.65, col=grey(c(0.3,0.6)),  pch=(19)))

reproduction7 <-predictorEffect("elev",  partial.residuals=FALSE, mod_7c7)

plot(reproduction7, lwd=2,xlab="Source Elevation (KM)", ylab="Probability of reproduction", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


#just parents
c_7p<-subset(c_7, mat_treat=="Parent")
c_7p<-subset(c_7, mat_treat=="Parent")


mod_7p<- glmer(reproduction~init_size +S_elev*treatment+ (1|block)+(1|genotype), data = (c_7p<-subset(c_7, treatment=="Control")), control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7p,type="III")

mod_7p1c<- glmer(reproduction~init_size +S_elev + I(S_elev^2)+ (1|block)+(1|genotype), data = (c_7p<-subset(c_7, treatment=="Control")), control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_7p1c,type="III")

reproduction <-predictorEffect("S_elev",  partial.residuals=TRUE, mod_7p1c)

plot(reproduction, lwd=2,xlab="Elevation of origin", ylab="Fitness (reproduced)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))



#survival for just herb treatment
mod_7b<- glmer(survival~S_elev*mat_treat*treatment+ (1|block)+(1|genotype), data = c_7, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)

survived <-predictorEffect("S_elev",  partial.residuals=TRUE, mod_7b)

plot(survived, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1))


#example code for optimizers
#flower<-glmer(cbind(avg_flowered,avg_surv-avg_flowered)~ Ts_avg_height +Ts_avg_leaves + Ts_avg_shootmass+T_RDPI_HEIGHT+ T_RDPI_LEAVES+T_RDPI_SHOOTMASS+(1|pop), data= temp, family=binomial,   control=glmerControl(optimizer = "nloptwrap",optCtrl=list(maxfun=2e7)))
#Anova(flower,type="III")
#flower_min<-glmer(cbind(avg_flowered,avg_surv-avg_flowered)~T_RDPI_HEIGHT+ Ts_avg_height+(1|pop), data= temp, family=binomial, glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb"))) 
#Anova(flower_min,type="III")

#firths correction

summary(mod_7b)
Anova(mod_7b,type = "III")

visreg(mod_7b, overlay = TRUE, "S_elev", by="mat_treat", type="conditional", scale = "response", 
       xlab="elev_km", ylab="Prob_surivival", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

mod_7ah<- glmer(survival~S_elev*mat_treat+ (1|block)+(1|genotype), data = c_7h, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)

summary(mod_7ah)
Anova(mod_7ah,type = "III")

visreg(mod_7ah, overlay = TRUE, "S_elev", by="mat_treat", type="conditional", scale = "response", 
       xlab="elev_km", ylab="Prob_surivival", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

mod_7ac<- glmer(survival~S_elev*mat_treat+ (1|block)+(1|genotype), data = c_7c, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)

summary(mod_7ac)
Anova(mod_7ac,type = "III")

visreg(mod_7ac, overlay = TRUE, "S_elev", by="mat_treat", type="conditional", scale = "response", 
       xlab="elev_km", ylab="Prob_surivival", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


#LAR
c_1<-subset(c_1, treatment=="Herbivorized")
c_1<- drop_na(c_1,LAR_prop)
c_1 <- drop_na(c_1,LAR_prop) 

plot(c_1$elev_km, c_1$LAR, xlab="elev", ylab="LAR")

LAR = ggplot(c_1, aes(x= elevation,y=LAR, group=mat_treat))+geom_point(aes(colour=mat_treat)) +scale_x_continuous("Elev")+ scale_y_continuous("LAR")
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                         panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6)
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                         panel.grid.minor=element_blank(), legend.position = "bottom")+ geom_smooth(method="glm",size=1.6
                                                                                                    , formula=y~poly(x,2)) 

gamlss.model_1a<- gamlss (formula=LAR_prop~elev_km*mat_treat+ random(block)+random(genotype), family=BEZI, data= c_1)
summary(gamlss.model_1a)
plot(gamlss.model_1a)


gamlss.model_1b<- gamlss (formula=LAR_prop~ mat_treat + random(block)+random(genotype), 
                        nu.formula=LAR_prop~ mat_treat + random(block)+random(genotype),family=BEZI, data= c_1)
summary(gamlss.model_1b)
plot(gamlss.model_1b)


mod_1a<- lmer(LAR_prop~elev_km+mat_treat+ (1|block)+(1|genotype), data = c_1)
summary(mod_1a)
Anova(mod_1a,type = "III")


    
 visreg(gamlss.model_1a, overlay = TRUE, "elev_km", by="mat_treat", type="conditional", #scale = "response", 
       xlab="Treatment", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))

####make data set for third census####
 c_3<-subset(gh2, census=="3")
 c_3<-subset(c_3, treatment=="Herbivorized")
 c_3<- drop_na(c_3,LAR_prop)
 c_3 <- drop_na(c_3,LAR_prop) 
 
 plot(c_3$elev_km, c_3$LAR, xlab="elev", ylab="LAR")
 
 LAR = ggplot(c_3, aes(x= elevation,y=LAR, group=mat_treat))+geom_point(aes(colour=mat_treat)) +scale_x_continuous("Elev")+ scale_y_continuous("LAR")
 LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6)
 LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                          panel.grid.minor=element_blank(), legend.position = "bottom")+ geom_smooth(method="glm",size=1.6
                                                                                                     , formula=y~poly(x,2)) 
 
 gamlss.model_3a<- gamlss (formula=LAR_prop~elev_km+mat_treat+ random(block)+random(genotype), family=BEZI, data= c_3)
 summary(gamlss.model_3a)
 plot(gamlss.model_3a)
 
 
 gamlss.model_3b<- gamlss (formula=LAR_prop~ mat_treat + random(block)+random(genotype), 
                         nu.formula=LAR_prop~ mat_treat + random(block)+random(genotype),family=BEZI, data= c_3)
 summary(gamlss.model_3b)
 plot(gamlss.model_3b)
 
 
 mod3a<- lmer(LAR_prop~elev_km+mat_treat+ (1|block)+(1|genotype), data = c_3)
 summary(mod3a)
 Anova(modA,type = "III")
 
 
 
 visreg(gamlss.model_3a, overlay = TRUE, "elev_km", by="mat_treat", type="conditional", #scale = "response", 
        xlab="Treatment", ylab="LAR", partial=TRUE,
        fill=list(col="light grey"
                  #(c(0.99), alpha=0)
        ), band = FALSE,
        #line=list(col=grey(c(0.2,0.6))),
        points=list(cex=0.65,  pch=(19)))


 
##################################
##### Beta regression with beta transformation. No random effects are possible.
       
       library(betareg)

n<-nrow(grasshopper_census3)
#this is the beta transformation, which transforms all values of 0 to a small value.
grasshopper_census3$y_LAR_beta<- (grasshopper_census3$LAR_prop*(n-1) + 0.5)/n

hist(grasshopper_census3$y_LAR_beta)

mod1 <- betareg(y_LAR_beta ~ Treatment*Herbivore* elev_km, data = grasshopper_census3)
summary(mod1)


mod2 <- betareg(y_LAR_beta ~ Treatment*Herbivore+ elev_km, data = grasshopper_census3)
summary(mod2)


visreg(mod2, overlay = TRUE, "Treatment", by="Herbivore", type="conditional", #scale = "response", 
       xlab="Treatment", ylab="Leaf area removed (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))
       
  ##For visualizing     
mod3 <- betareg(y_LAR_beta ~ treat_herb+elev_km, data = grasshopper_census3)
summary(mod3)
visreg(mod3, overlay = TRUE, "elev_km", by="treat_herb", type="conditional", #scale = "response", 
       xlab="Source elevation (km)", ylab="Leaf area removed (proportion)", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))


       
