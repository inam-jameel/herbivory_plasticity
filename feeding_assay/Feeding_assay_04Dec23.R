######## PROJECT: feeding assay experiment: Fitness and phenotypes in response to herbivory
#### PURPOSE:Examine LAR in response to maternal herbivory treatment.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 6 Jul 24

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require("ggplot2") #for plotting 
require("visreg") # for plotting
require("car") # to run ANOVA on model output
require("plyr") # for data wrangling
require("dplyr") # for data wrangling
require("tidyr") # for data wrangling
require("effects") # for plotting
require("emmeans") #for plotting
require("glmmTMB") # for running survival model, have to load twice
require("gamlss") # for running phenology model
require("multcomp") #for pairwise comparisons
require("vioplot") #for violin plots


# set working directory

setwd("~/OneDrive - University of Georgia/Inam_experiments/Herbivory_data/feeding_assay/Feeding_trial_spring2021")  #this is where you specify the folder where you have the data on your computer

#import data
Feeding_trial <- read.csv("Feeding_trialNov2021a.csv")

##Convert elevation to km (it helps with model convergence)
Feeding_trial$elev<-Feeding_trial$elevation/1000
Feeding_trial$mat_avgLAR<-Feeding_trial$mat_avgLAR/100
Feeding_trial$S_weight<-scale(Feeding_trial$ini_insect_weight.mcg.,center=TRUE, scale=TRUE)
Feeding_trial$time<-scale(Feeding_trial$elapsed_time,center=TRUE, scale=TRUE) #need to include elapsed time in the models

#Change treatment, genotype, batch to factor
Feeding_trial$mat_treat<-as.factor(Feeding_trial$mat_treat)
Feeding_trial$Exp_ID <-as.factor(Feeding_trial$Exp_ID)
Feeding_trial$genotype <-as.factor(Feeding_trial$genotype)
Feeding_trial$batch<-as.factor(Feeding_trial$batch)
Feeding_trial$mat_exp_ID <-as.factor(Feeding_trial$mat_exp_ID) #need to include this as random effect since multiple reps per mat line

#set colors
cols=c("#882255","#56B4E9")

#LAR plots
plot(Feeding_trial$elev, Feeding_trial$LAR, xlab="elev", ylab="LAR")

plot(Feeding_trial$ini_insect_weight.mcg., Feeding_trial$LAR, xlab=" inital weight", ylab="LAR")

plot(Feeding_trial$final_insect_weight, Feeding_trial$LAR, xlab="final weight", ylab="LAR")

plot(Feeding_trial$final_insect_weight, Feeding_trial$ini_insect_weight.mcg., xlab="final weight", ylab="inital weight")

plot(Feeding_trial$time_minutes, Feeding_trial$LAR, xlab="Time", ylab="LAR") 

plot(Feeding_trial$time_started, Feeding_trial$LAR, xlab="Time", ylab="LAR") #lower value indicates earlier in the day

plot(Feeding_trial$time_end, Feeding_trial$LAR, xlab="Time", ylab="LAR") #lower value indicates earlier in the day

#RGR plots

plot(Feeding_trial$RGR, Feeding_trial$LAR, xlab="RGR", ylab="LAR")

print(cor(Feeding_trial$RGR, Feeding_trial$LAR, use="complete.obs"))
test <- cor.test(Feeding_trial$RGR, Feeding_trial$LAR)
test

plot(Feeding_trial$elev, Feeding_trial$RGR, xlab="elev", ylab="RGR")

plot(Feeding_trial$ini_insect_weight.mcg., Feeding_trial$LAR, xlab=" inital weight", ylab="RGR")

plot(Feeding_trial$time_minutes, Feeding_trial$LAR, xlab="Time", ylab="RGR") 

plot(Feeding_trial$time_started, Feeding_trial$LAR, xlab="Time", ylab="RGR") #lower value indicates earlier in the day

plot(Feeding_trial$time_end, Feeding_trial$LAR, xlab="Time", ylab="RGR") #lower value indicates earlier in the day



#### LAR ####
#looking at LAR by elevation

data_LAR <- dplyr::select(Feeding_trial, LAR, elev, genotype, mat_treat, Exp_ID, time, batch,mat_exp_ID,S_weight,ini_insect_weight.mcg.,mat_avgLAR)
data_LAR <- drop_na(data_LAR,LAR) #this removes the rows without LAR values


ggplot(data_LAR, aes(x= LAR))+ geom_histogram(color="black", fill="white")


#The scatter is very high if you just plot the relationship with ggplot:

assayLAR =ggplot(data_LAR, aes(x= elev,y= LAR,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")

assayLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))+
  geom_smooth(method="lm",se=FALSE,  size=1,formula=y~x)+facet_wrap(~ mat_treat, scales="free_x")

# quadratic
assayLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255"))+
  geom_smooth(method="lm",se=FALSE,  size=1,formula=y~poly(x,2))+facet_wrap(~ mat_treat, scales="free_x")


#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

n<-nrow(data_LAR)

data_LAR $y_beta<- (data_LAR $LAR*(n-1) + 0.5)/n

hist(data_LAR $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(data_LAR $y_beta)

max(data_LAR $y_beta)



#Then, the analysis is:
library(betareg)

beta_model<- betareg ( y_beta ~elev*mat_avgLAR+batch+ini_insect_weight.mcg.+time, data=data_LAR)
Anova(beta_model,type="III")
plot(beta_model)


visreg(beta_model, "elev", by="mat_avgLAR", overlay=TRUE, scale="response")

visreg(beta_model, 'elev', by= "mat_avgLAR", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 


beta_modelquad<- betareg ( y_beta ~elev*mat_avgLAR+I(elev^2)*mat_avgLAR+batch+ini_insect_weight.mcg.+time, data=data_LAR)
Anova(beta_model,type="III")
plot(beta_model)


visreg(beta_modelquad, "elev", by="mat_avgLAR", overlay=TRUE, scale="response")

visreg(beta_modelquad, 'elev', by= "mat_avgLAR", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 


#####

#Box_plot
LAR_box <-ggplot(data_LAR, aes(x = mat_treat, y = y_beta, fill = mat_treat)) +
  geom_boxplot(outlier.shape = NA) +xlab("Maternal herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

LAR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Herbivory", "Naïve")) +  scale_fill_manual(values = cols, name = "Maternal herbivory treatment", labels = c("Herbivory","Naïve"))#+facet_grid(~Season+mat_treat)




## model
LAR_Model<- gamlss (formula= y_beta ~elev*mat_treat*ini_insect_weight.mcg.+ random(batch)+random(mat_exp_ID) +random(genotype),family=BE(mu.link = "logit"), data=data_LAR,control = gamlss.control(n.cyc = 500))

plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)




## model
LAR_Model<- gamlss (formula= y_beta ~elev*mat_treat+S_weight+ random(batch)+random(mat_exp_ID) +random(genotype),family=BE(mu.link = "logit"), data=data_LAR,control = gamlss.control(n.cyc = 500))

plot(LAR_Model)
summary(LAR_Model)
drop1(LAR_Model)



####

gamlss_mod1<- gamlss (formula=y_beta~S_weight+time+elev+mat_treat+random(batch)+random(genotype)+ random(mat_exp_ID),family=BE, data=data_LAR,control = gamlss.control(n.cyc = 500))



plot(gamlss_mod1)

summary(gamlss_mod1)


visreg(gamlss_mod1, 'S_weight', overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Source Elevation (Km)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255"))) 


gamlss.modelB<- gamlss (formula= LAR ~ time+S_weight +elev*mat_treat+ random(genotype)+ random(batch) +random(mat_exp_ID), 
                        nu.formula=LAR ~ time+S_weight +elev*mat_treat+random(genotype)+random(batch) + random(mat_exp_ID), 
                        tau.formula=LAR ~ time+S_weight +elev*mat_treat+ random(batch)+random(genotype)+random(mat_exp_ID), family=BEINF, data= data_LAR) 
plot(gamlss.modelB)
summary(gamlss.modelB)

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.modelB)

mod2<-stepGAIC(gamlss.modelB)

mod2$anova

summary(mod2)

#mod2:  gamlss(formula = LAR ~ S_weight + random(batch) + random(mat_exp_ID),  nu.formula = LAR ~ time + S_weight + elev * mat_treat +          random(genotype) + random(batch) + random(mat_exp_ID),  tau.formula = LAR ~ time + S_weight + elev * mat_treat +          random(batch) + random(genotype) + random(mat_exp_ID),  family = BEINF, data = data_LAR) 

visreg(gamlss.modelB, overlay = TRUE, "elev", by="mat_treat", type="conditional", scale = "linear", 
       xlab="Source Elevation (KM)", ylab="Probability of reproduction", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(lty=1, col=grey(c(0.2,0.6))),
       points=list(cex=0.65, col=grey(c(0.3,0.6)),  pch=(19)))


gamlss.modelC<- gamlss (formula= LAR ~ S_weight +elev*mat_treat+I(elev^2)*mat_treat+ random(genotype)+ random(batch)+random(mat_exp_ID), nu.formula=LAR ~ S_weight +elev*mat_treat+I(elev^2)*mat_treat+random(genotype)+random(batch)+random(mat_exp_ID), tau.formula=LAR ~ S_weight +elev*mat_treat+I(elev^2)*mat_treat+ random(batch)+random(genotype)+random(mat_exp_ID), family=BEINF(mu.link = "identity", sigma.link = "log"), data=data_LAR,control = gamlss.control(n.cyc = 500))

plot(gamlss.modelC)
summary(gamlss.modelC)

dropterm(gamlss.modelC)

mod3<-stepGAIC(gamlss.modelC)

mod3$anova

summary(mod3)

#mod3:  gamlss(formula = LAR ~ S_weight + elev + mat_treat +      I(elev^2) + random(genotype) + random(batch) +  elev:mat_treat + mat_treat:I(elev^2), nu.formula = LAR ~      S_weight + elev * mat_treat + I(elev^2) * mat_treat +   random(genotype) + random(batch) + random(mat_exp_ID),      tau.formula = LAR ~ S_weight + elev * mat_treat +   I(elev^2) * mat_treat + random(batch) + random(genotype) +          random(mat_exp_ID), family = BEINF, data = data_LAR) 


visreg(gamlss.modelC, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Source Elevation (Km)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255"))) 







#lmer model

modA<- lmer (LAR~elev*mat_treat+I(elev^2)*mat_treat+S_weight+(1|batch)+(1|genotype/mat_exp_ID),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)

Anova(modA,type="III")
plot(modA)

model1 <-predictorEffect("elev",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))


# Generalized additive model
library(mgcv)
mod_1a <- gam(LAR ~mat_treat+ini_insect_weight.mcg.+s(elev, by=mat_avgLAR )+s(genotype,bs="re")+s(batch,bs="re"),data= data_LAR, method="REML")

plot(mod_1a,pages=1,scheme=1,residuals=TRUE,unconditional=TRUE) 
plot(mod_1a, shade = TRUE, pages = 1, scale = 0)
summary(mod_1a) #batch is significant here.. 

visreg(mod_1a, overlay = FALSE, "elev", by="mat_avgLAR", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))



#### RGR ####

#looking at RGR by elevation

data_RGR <- dplyr::select(Feeding_trial, LAR, RGR, elev, genotype, mat_treat, Exp_ID, time, batch,mat_exp_ID,S_weight,ini_insect_weight.mcg., final_insect_weight,mat_avgLAR)

data_RGR <- drop_na(data_RGR,RGR) #this removes the rows without RGR values
data_RGR <- drop_na(data_RGR,LAR) #this removes the rows without LAR values


ggplot(data_RGR, aes(x= RGR))+ geom_histogram(color="black", fill="white")

#gamlss model
gamlss_rgr<- gamlss (formula=RGR~elev+mat_treat+random(batch)+random(genotype)+ random(mat_exp_ID), data=data_RGR,control = gamlss.control(n.cyc = 500))

plot(gamlss_rgr)
summary(gamlss_rgr)
drop1(gamlss_rgr)

hist(data_RGR$RGR)

#lmer model
RGR <-lmer (RGR~ elev+mat_avgLAR +time
            +(1|genotype)
            +(1|batch)
            +(1|mat_exp_ID)
            , data = data_RGR,
            control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

#lmer model
RGR <-lmer (RGR~ elev*mat_treat
            +(1|genotype)
            +(1|batch)
            +(1|mat_exp_ID)
            , data = data_RGR,
            control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

plot(RGR)

Anova(RGR,type="III")


##Generalized additive model
library(mgcv)
mod_1a <- gam(RGR ~mat_treat+time+S_weight+s(elev, by=mat_treat)+s(genotype,bs="re")+s(batch,bs="re"),data= data_RGR, method="REML")
## use alternative plotting scheme, and way intervals include
## smoothing parameter uncertainty...
plot(mod_1a,pages=1,scheme=1,unconditional=TRUE) 
plot(mod_1a,pages=1,scheme=1,residuals=TRUE,unconditional=TRUE) 
summary(mod_1a)
plot(mod_1a, shade = TRUE, pages = 1, scale = 0)



#Ploting RGR
library(ggpubr)


RGR2 = ggplot(LSmeans_RGR, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Insect Relative Growth Rate") + geom_point(aes(colour = mat_treat), size = 4)
R <- RGR2 + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = mat_treat),method="glm",size=2.0, se=F)+geom_smooth(aes(group = mat_treat, colour = mat_treat),method="lm",size=1.6, se=F)
print(R)


RGR_box <-ggplot(data_RGR, aes(x = mat_treat, y = RGR, fill = mat_treat)) +
  geom_boxplot(outlier.shape = NA) +xlab("Maternal herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

RGR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                       axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                       panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Herbivory", "Naïve")) +  scale_fill_manual(values = cols, name = "Maternal herbivory treatment", labels = c("Herbivory","Naïve"))#+facet_grid(~Season+mat_treat)



### final insect weight ####
mod_2<- lmer(final_insect_weight ~ elev*mat_treat + time + ini_insect_weight.mcg.+ (1|genotype) +(1|batch), control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = data_RGR)
Anova(mod_2,type="III")
plot(mod_2)


final_weight_model <- glmmTMB(final_insect_weight ~ mat_treat*elev+(1|Exp_ID)+(1|genotype)+(1|batch), data = data_RGR, family= lognormal(link="log"))
Anova(final_weight_model, type = "III") # 
#Use the DHARMa package to examine the residuals, which are reasonable

simulationOutput <- simulateResiduals(fittedModel= final_weight_model, plot = T, re.form = NULL,allow.new.levels =T)
library(DHARMa)


#Box_plot
final_weight <-ggplot(data_RGR, aes(x = mat_treat, y = final_insect_weight, fill = mat_treat)) +
  geom_boxplot(outlier.shape = NA) +xlab("Maternal herbivore treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +
  geom_point(pch = 21, size = .5,position = position_jitterdodge(0.3))

final_weight + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                  axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                  panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Herbivory", "Naïve")) +  scale_fill_manual(values = cols, name = "Maternal herbivory treatment", labels = c("Herbivory","Naïve"))#+facet_grid(~Season+mat_treat)



visreg(mod_2, "elev", by="mat_treat")
transgen<-predictorEffects(mod_2, partial.residuals=T)
plot(transgen,type="response")
transgen <-predictorEffect("elev",  partial.residuals=TRUE, mod_2)
plot(transgen, lwd=2,xlab="Source elevation (KM)", ylab="Insect relative growth rate", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=FALSE, pch=19, col="black"))









