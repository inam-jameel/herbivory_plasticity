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

#status of this file is unknown. there is code at the very bottem (250+) that i used to make figures for BT effectiveness. Additionally, the code to get the lsmeans (Lsmeans_LAR) is found at 187

#### merging data ####
#schofield <- read.csv("Field2021_schofield.csv")
#s_data <- read.csv("schofield_plant_data.csv")

#schofield_a<- merge(schofield,s_data,by="plantID")

#write.csv(schofield_a,"~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/Field2021_schofield_merge.csv")

#g_2020 <- read.csv("Field2021_gothic2020.csv")
#g_data <- read.csv("gothic_2020_plant_data.csv")

#gothic_a<- merge(g_2020,g_data,by="plantID")

#write.csv(gothic_a,"~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/Field2021_gothic2020_merge.csv")


#### read in data ####
data<- read.csv("Field2021_merged.csv")


data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)

data$genotype<-as.factor(data$genotype)
data$treatment<-as.factor(data$Treatment)
data$block<-as.factor(data$block)
data$cohort<-as.factor(data$cohort)

data$elev_dist<-data$garden_elevation-data$elevation


# creating a df for LAR
data_LAR <- drop_na(data,LAR) #this removes the rows without LAR values

##convert LAR to proportion
data_LAR $LAR_prop<- data_LAR $LAR/100
##Check that the conversion worked
hist(data_LAR $LAR_prop)
hist(data_LAR $LAR)

data_LAR <-na.omit(data_LAR)

#Let's concatenate Block and Garden so that we don't have to worry about nesting Block within Garden
      data_LAR $Block_Garden<-interaction(data_LAR $garden, data_LAR $block,sep = "_")

 
gamlss.model1<- gamlss (formula=LAR_prop~ cohort+S_elev*treatment*garden+random(Block_Garden)+random(genotype), 
nu.formula=LAR_prop~ cohort+S_elev*treatment*garden+random(Block_Garden)+random(genotype),
family=BEZI, data = data_LAR)
summary(gamlss.model1)
plot(gamlss.model1)






#looking at LAR v elevation ggplot
LAR <- plot(data_LAR$elevation, data_LAR$LAR, xlab="elev", ylab="LAR")

LARfig = ggplot(LAR, aes(x= elevation,y=LAR))+geom_point(size=5) +scale_x_continuous("LAR")+ scale_y_continuous("elevation")
LARfig + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6)
LARfig + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                           panel.grid.minor=element_blank(), legend.position = "bottom")+ geom_smooth(method="glm",size=1.6
                                                                                                      , formula=y~poly(x,2)) 
                                                                                                      
                                                                                                      
                                                                                                      
  ##Concatenate treatment and garden so that we can plot from visreg
      data_LAR $garden_treat<-interaction(data_LAR $garden, data_LAR $treatment,sep = "_")

gamlss.modelA <- gamlss (formula=LAR_prop~ cohort+S_elev* garden_treat +random(Block_Garden)+random(genotype), 
nu.formula=LAR_prop~ cohort+S_elev* garden_treat +random(Block_Garden)+random(genotype),
family=BEZI, data = data_LAR)
summary(gamlss.modelA)
plot(gamlss.modelA)

Anova(gamlss.modelA,type="III") #error regarding aliased coefficients

  
  visreg(gamlss.modelA, overlay = FALSE, "S_elev", by="garden_treat", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))                                                                                                    
                                                                                                      


gamlss.modelA <- gamlss (formula=LAR_prop~ cohort+S_elev* garden_treat +random(Block_Garden)+random(genotype), 
nu.formula=LAR_prop~ cohort+S_elev* garden_treat +random(Block_Garden)+random(genotype),
family=BEZI, data = data_LAR)
summary(gamlss.modelA)
plot(gamlss.modelA)



##We could try a beta regression with a beta transformation (so that there are no 0 values) 
sample_size<-nrow(data_LAR)
data_LAR$lar_beta<-(data_LAR$LAR_prop * (sample_size-1) + 0.5)/sample_size
min(data_LAR$lar_beta)
#Now that there are no more 0 values, we can simplify our gamlss model
gamlss.modelC <- gamlss (formula= lar_beta ~genotype*treatment*garden+random(Block_Garden), 
family=BE, data = data_LAR)


##The error message about smoothing seems to go away if we remove the random effect for block
gamlss.modelC <- gamlss (formula= lar_beta ~genotype*treatment*garden, 
family=BE, data = data_LAR)

summary(gamlss.modelC)
plot(gamlss.modelC)
 fam_avg<-emmeans(gamlss.modelC, ~genotype:treatment:garden, type="response")
fam_LSMEANS<-as.data.frame(fam_avg)

fam_avg<-emmeans(gamlss.modelC, lar_beta~genotype:treatment:garden, what = c("mu"), type="response")

modC<-gamlss(formula = lar_beta~ genotype*treatment*garden, sigma.formula=lar_beta~genotype*treatment*garden, nu.formula = lar_beta~ genotype*treatment*garden, family=BEINF, data = data_LAR)
summary(modC)

fam_avg<-emmeans(modC, ~genotype:treatment:garden, type="response")

write.csv(LSmeans_estess,file="~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/LSmeans_LAR.csv")
write.csv(LSmeans_estessp,file="~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/LSmeans_LARp.csv")


LSmeans_LAR <- read.csv("~/Desktop/Anderson_data/Herb-plasticity/herb-plasticity/data/2021_herbivory/LSmeans_LAR.csv", stringsAsFactors=FALSE)

estess <- ggplot(filter(LSmeans_LAR), aes(elevation, response)) + geom_point() 
+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "Elevation") + ggtitle("Comparsion of two treatments after 6 hr exposure")
print(estess)

LAR = ggplot(LSmeans_LAR, aes(x= elevation,y=response))+geom_point(size=5) +scale_x_continuous("LAR")+ scale_y_continuous("elevation")
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                         panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6)
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                         panel.grid.minor=element_blank(), legend.position = "bottom")+ geom_smooth(method="glm",size=1.6)
                                                                                                    

##beta reg model

library(betareg)
beta_model<-betareg(lar_beta ~genotype*treatment*garden, data = data_LAR)
 fam_avgB<-emmeans(beta_model, ~genotype:treatment:garden, type="response")
fam_LSMEANSB<-as.data.frame(fam_avgB)


##For a moment, let's treat LAR as a normal variable and run a mixed model, even though this is clearly not appropriate

modA<- lmer (LAR~S_elev*treatment*garden+(1|Block_Garden)+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)
Anova(modA,type="III")
plot(modA)

model1 <-predictorEffect("S_elev",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
partial.residuals=list(smooth= TRUE, pch=19, col="black"))


##Now to get LSMEANS
modB<- lmer (LAR_prop~genotype*treatment*garden+(1|Block_Garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)
  

fam_avg_linear<-emmeans(modB, ~genotype:treatment:garden)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_LAR.csv")


LSmeans_LAR <- read.csv("LSmeans_LAR.csv", stringsAsFactors=FALSE)
LSmeans_LAR <-na.omit(LSmeans_LAR)



CtreatmentByelev <- ggplot(filter(LSmeans_LAR,treatment == "C"), aes(elevation, emmean)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "treatment") + ggtitle("control treatment across all gardens")
print(CtreatmentByelev)

PtreatmentByelev <- ggplot(filter(LSmeans_LAR,treatment == "P"), aes(elevation, emmean)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "treatment") + ggtitle("control treatment across all gardens")
print(PtreatmentByelev)

modG<- lm(emmean ~ elevation,data=filter(LSmeans_LAR,treatment == "herbivory"))
visreg(modG,"elev_km", type="contrast", scale = "response", xlab=" Source Elevation (km)", ylab="Relative Leaf Area Removed (%)", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)


estess <- ggplot(filter(LSmeans_LAR), aes(elevation, response)) + geom_point() 
+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "Elevation") + ggtitle("Comparsion of two treatments after 6 hr exposure")
print(estess)

LAR = ggplot(LSmeans_LAR, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous("LAR")+ scale_y_continuous("elevation") + geom_point(aes(colour = factor(treatment)), size = 4)
LAR + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                         panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=1.6)



modA<- lmer (emmean~elev_km*treatment*garden+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)
Anova(modA,type="III")
plot(modA)

model1 <-predictorEffect("elev_km",  partial.residuals= TRUE, modA)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=FALSE,  col="black"), 
     partial.residuals=list(smooth= TRUE, pch=19, col="black"))




#Ploting LAR 
library(ggpubr)

ggarrange(E, G, S, 
          #labels = c("A", "B", "C"), 
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1)



Estess <- subset(LSmeans_LAR,garden=="Estess")
E_LAR = ggplot(Estess, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = treatment), size = 4)
E <- E_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)

Gothic <- subset(LSmeans_LAR,garden=="Gothic")
G_LAR = ggplot(Gothic, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = factor(treatment)), size = 4)
G <- G_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)

Scho <- subset(LSmeans_LAR,garden=="Schofield")
S_LAR = ggplot(Scho, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = factor(treatment)), size = 4)
S <-S_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                         axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                         panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)


LAR = ggplot(LSmeans_LAR, aes(x= elevation,y=emmean))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Leaf area removed %") + geom_point(aes(colour = factor(treatment)), size = 4)
LAR + theme_bw() + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                                axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                                panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)

#####RMBL update figures here 

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

#Leaf area removed ~Elevation*maternal env+ Elevation2 *Maternal env + Scaled weight of insect +(1|batch)+(1|genotype/mat_exp_ID)


#estess 
Estess <- subset(LSmeans_LAR,garden=="Estess")

estess_mod <- lmer (emmean100~treatment*elevation*garde+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Estess)
plot(estess_mod)
Anova(estess_mod,type="III")

E1 <- visreg(mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))))  

model1 <-predictorEffect("elevation",  partial.residuals= FALSE, mod)
plot(model1, lwd=2,xlab="Elevation", ylab="LAR", pch=19, type="response",lines=list(multiline=TRUE),confint=list(style="auto"), 
     partial.residuals=list(smooth= TRUE))

model1 <-predictorEffect("elevation",  partial.residuals=TRUE, mod)
plot(model1, lwd=2,xlab="Elevation of origin", ylab="Fitness (survival)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1))


#trying gamlss model
E_gam<-na.omit(Estess) #get rid of NAs
E_gam$genotype <- as.factor(E_gam$genotype) #random needs to be a factor
gamlss.model_estess<- gamlss (formula=emmean~elevation*treatment+random(genotype), family=BEZI, data= Estess) #response out of range
   

#gothic
Gothic <- subset(LSmeans_LAR,garden=="Gothic")

gothic_mod <- lmer (emmean100~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Gothic)
plot(gothic_mod)
Anova(gothic_mod,type="III")

G1 <- visreg(gothic_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  


#schofield

Scho <- subset(LSmeans_LAR,garden=="Schofield")
scho_mod <- lmer (emmean100~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Scho)
plot(scho_mod)
Anova(scho_mod,type="III")

visreg(scho_mod, overlay = TRUE, "elev_km", by="treatment", type="conditional", scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = TRUE,
       line=list(col=grey(c(0.2,0.6))),
       points=list(col=grey(c(0.2,0.6))))  

S1 <- visreg(scho_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  



Estess_com <- subset(data_LAR,garden=="Estess")
E_LAR = ggplot(Estess_com, aes(x= elevation,y=LAR))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = treatment), size = 4)
E <- E_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)






