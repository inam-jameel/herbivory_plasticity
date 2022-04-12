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

herbivory<- read.csv("~/Desktop/Anderson/Transgen_plasticity/old_versions/Herbivory_genotypes_022020_ija.csv", stringsAsFactors=FALSE)

#makeing sure that columns have proper designations
herbivory$treatment <- as.factor(herbivory$treatment)
herbivory$population <- as.factor(herbivory$population)
herbivory$block <- as.factor(herbivory$block)
herbivory$elev_km<-herbivory$elevation/1000
herbivory$elev_km <- as.factor(herbivory$elev_km)

#making subsets
herb <- filter(herbivory, treatment == "herbivory")
cont <- filter(herbivory, treatment == "control")

#basic plots
treatmentVcontrol<- ggplot(herbivory, aes(treatment, LAR)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "treatment") + theme_classic(base_size = 25) # big, big text
print(treatmentVcontrol)

treatmentByelev <- ggplot(herb, aes(elev_km, LAR)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+ labs(y="leaf area removed (%)", x = "treatment") + ggtitle("Comparsion of two treatments after 6 hr exposure")
print(treatmentByelev)

modA<- lmer(LAR ~ treatment + elev_km + (1|block),data=herbivory)
modA<- lmer(LAR ~ elev_km + (1|block),data=herb)


modB<- lmer(LAR ~ elev_km + (1|block),data=herb)


modC<- lmer(LAR ~ treatment (1|block),data=herbivory)
modD<- lmer(LAR_2 ~ elev_km + (1|geno)+(1|block),data=herb)
modE<- lmer(num_leaves_ ~ elev_km + (1|geno)+(1|block),data=herbivory)

Anova(modA)
Anova(modB)

summary(modA)
t.test(herb$LAR, cont$LAR, var.equal = TRUE, paired = FALSE, alternative = "greater")
visreg(modA,"elev_km", by="treatment", type="contrast", scale = "response", xlab="num_leaves", ylab="Leaf Area Removed", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)
visreg(modA,"elev_km", type="contrast", scale = "response", xlab="elevation", ylab="Leaf Area Removed", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)

visreg(modE,"elev_km", type="contrast", scale = "response", xlab="elevation", ylab="Number of leaves", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)
visreg(modB,"inital_size_mm", by ="elevation", type="contrast", scale = "response", xlab="R1tot", ylab="Viability", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)
visreg(modB,"elev_km", type="contrast", scale = "response", xlab="Treatment", ylab="Leaf Area Removed", line=list(col="black"), points=list(cex=1, col="black", pch=16), partial=TRUE)



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







####keep this section for later####
## see email from jill from 2/18/2020 "Quick question about lsmeans"
modC<-gamlss(formula = LAR_2 ~ treatment*elevation, sigma.formula=LAR_2 ~ treatment*elevation, nu.formula = LAR_2 ~ treatment*elevation, family=BEINF, data = na.omit(herbivory))
# get this error Warning message: In RS() : Algorithm RS has not yet converged
sapply(herbivory,class)
#trying to add random effects:
modD<-gamlss(formula = LAR_2 ~ treatment*elevation, sigma.formula=LAR_2 ~ treatment*elevation, nu.formula = LAR_2 ~ treatment*elevation, family=BEINF, data = na.omit(herbivory), random(population,))

        