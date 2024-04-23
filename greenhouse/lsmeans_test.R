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

####
library(lme4)
library(emmeans)

##This code runs the model

modelA<-lme(SLA~ Garden*Treatment +(1|Garden_Block)+(1|genotype), data=dat)

summary(modelA)

Anova(modelA)



##This code extracts the LSMEANs (the average values with SE for each garden,treatment combination)

lsmeans_modelA <-emmeans (modelA, Garden~Treatment, type=”response”)

##this code should generate a plot for you. You will likely need to change the order of some of the gardens and treatments to get a plot similar to what I wrote (with the gardens and treatments on the X axis in order of decreasing aridity

emmip(modelA _full, Garden ~  Treatment, CIs=TRUE)+theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+ 
  
  ylab("SLA")

####


emmeans(modD, list(pairwise ~ elevation), adjust = "tukey")


# using lsmeans to plot things
modF<- lmer(LAR ~ treatment:elev_km:num_leaves_ +(1|block),data=herbivory)

LSmeans_LAR <- emmeans(modF,~elev_km:treatment:num_leaves_)
write.csv(LSmeans_LAR,file="Desktop/Anderson/Transgen_plasticity/LSmeans_LAR.csv")
LSmeans_LAR <- read.csv("Desktop/Anderson/Transgen_plasticity/LSmeans_LAR.csv", stringsAsFactors=FALSE)

treatmentByelev <- ggplot(filter(LSmeans_LAR,treatment == "herbivory"), aes(elev_km, emmean)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "treatment") + ggtitle("Comparsion of two treatments after 6 hr exposure")
print(treatmentByelev)

modG<- lm(emmean ~ elev_km,data=filter(LSmeans_LAR,treatment == "herbivory"))
visreg(modG,"elev_km", type="contrast", scale = "response", xlab=" Source Elevation (km)", ylab="Relative Leaf Area Removed (%)", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)



        