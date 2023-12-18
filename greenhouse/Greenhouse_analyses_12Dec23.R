######## PROJECT: Greenhouse experiment: Fitness and phenotypes in response to herbivory
#### PURPOSE:Examine fitness and traits in response to herbivory treatment in the greenhouse.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 12 Dec 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(tidyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for plotting
require(glmmTMB) # for running survival model, have to load twice
require(gamlss) # for running phenology model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(vioplot) #for violin plots
library(RColorBrewer) #for 2d visreg plots

# set working directory

setwd("~/OneDrive - University of Georgia/Inam_experiments/Herbivory_data/greenhouse/")  #this is where you specify the folder where you have the data on your computer

#read in data 
gh2 <- read.csv("gen2_GHSummary2022.csv")

gh2 <- filter(gh2, Include == "yes" ) #removed plants with unknown genotypes

gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)



# convert LAR  to proportion

#maternal leaf damage

gh2$mat_avgLAR_pzero<- gh2$mat_avgLAR_pzero/100
hist(gh2$mat_avgLAR_pzero) 

#summary
gh2$LAR_avg_prop<- gh2$LAR_avg/100
hist(gh2$LAR_avg_prop) 

gh2$LAR_max_prop<- gh2$LAR_max/100 
hist(gh2$LAR_max_prop) 

#season 2 #for repeated measures
gh2$S2_LAR1_prop<- gh2$S2_LAR1/100
hist(gh2$S2_LAR1_prop) 

gh2$S2_LAR2_prop<- gh2$S2_LAR2/100 #most recent
hist(gh2$S2_LAR2_prop) 


#season 1

gh2$S1_LAR1_prop<- gh2$S1_LAR1/100 
hist(gh2$S1_LAR1_prop)

gh2$S1_LAR3_prop<- gh2$S1_LAR3/100 #most damage
hist(gh2$S1_LAR3_prop) 

gh2$S1_LAR13_prop<- gh2$S1_LAR13/100 #last in S1
hist(gh2$S1_LAR13_prop)

##variables from characters to factors
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$exp_ID<-as.factor(gh2$exp_ID)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line
gh2 $treat_mat<-interaction(gh2 $treatment, gh2 $mat_treat,sep = "_")
gh2 $S_LAR_avg_prop<-scale(gh2$LAR_avg_prop,center=TRUE, scale=TRUE)


#*******************************************************************************
#### 1.SLA for greenhouse experiment #####
#*******************************************************************************

### plotting data

SLA_S1 =ggplot(gh2, aes(x= elevation,y= S1_SLA, color = treat_mat))+ geom_point(size=3.5)+ ggtitle("SLA season 1 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S1

SLA_S2 =ggplot(gh2, aes(x= elevation,y= S2_SLA, color = mat_avgLAR_pzero))+ geom_point(size=3.5)+ ggtitle("SLA season 2 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S2

#converting wide to long for repeated measures, this filters out the SLA that came from too small leaves

SLA_data<- gh2 %>% pivot_longer(cols=c('S2_SLA_include', 'S1_SLA_include'),
                                names_to='year',
                                values_to='SLA')

SLA_data <- dplyr::select(SLA_data, SLA, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, year,LAR_avg_prop, mat_avgLAR_pzero)

ggplot(SLA_data, aes(x= SLA))+ geom_histogram(color="black", fill="white")+ facet_grid(treatment ~  .)

#mat_avgLAR_pzero is the maternal damage, with the grandparental treatment listed as a zero

SLA_RM <- lmer(SLA ~ LAR_avg_prop*mat_avgLAR_pzero*elev+year + (1|block)+(1|genotype)+(1|exp_ID),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)

plot(SLA_RM)

Anova(SLA_RM)

visreg(SLA_RM,"elev", by="LAR_avg_prop", overlay=TRUE,   scale = "response", xlab="Leaf area removed by herbivores", ylab="SLA CM3/g ", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("darkred","orange","#56B4E9")), points=list(col=c("darkred","orange","#56B4E9")),fill=list(col=grey(c(1), alpha=0.4)))

visreg2d(SLA_RM,"elev","LAR_avg_prop",xlab="Source elevation (Km)", ylab="Leaf area herbivorized")
visreg2d(SLA_RM,"elev","LAR_avg_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(10),zlim=c(200,300),ylim=c(0,.25))


# want to do selection, but unsure if I need to do LSmeans? 

LSmeans_SLA<- lmer(SLA ~ treatment*mat_treat*genotype+year + (1|block)+(1|genotype)+(1|exp_ID),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)

LSmeans_SLA <- emmeans(LSmeans_SLA,~treatment*mat_treat*genotype)


elev <- SLA_data[c("genotype","elev")] #make dataframe of genotypes and corresponding elev
elev <- unique(elev) #calls unique rows 
LSmeans_SLA <- merge(LSmeans_SLA,elev,by="genotype") #merge the dataframes

write.csv(LSmeans_SLA,file="LSmeans_SLA.csv")
LSmeans_SLA <- read.csv("LSmeans_SLA.csv", stringsAsFactors=FALSE)

#*******************************************************************************
#### 2.LAR for greenhouse experiment #####
#*******************************************************************************


herb <- subset(gh2, treatment == "Herbivorized")

#histograms of avg damage and max damage
ggplot(herb, aes(x= LAR_avg_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)
ggplot(herb, aes(x= LAR_max_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)


#avg
elev_avgLAR =ggplot(herb, aes(x= elevation,y= LAR_avg_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")
elev_avgLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_avgLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

#repeated measures 

LAR_data<- herb %>% pivot_longer(cols=c('S1_LAR1_prop', 'S1_LAR3_prop', 'S1_LAR13_prop','S2_LAR1_prop', 'S2_LAR2_prop'),
                                 names_to='exposure',
                                 values_to='LAR')

LAR_data <- dplyr::select(LAR_data, LAR, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, exposure,mat_avgLAR_pzero)

LAR_data$census<-LAR_data$exposure
LAR_data$census[LAR_data$census == "S1_LAR1_prop"] <- "1"
LAR_data$census[LAR_data$census == "S1_LAR3_prop"] <- "2"
LAR_data$census[LAR_data$census == "S1_LAR13_prop"] <-"3"
LAR_data$census[LAR_data$census == "S2_LAR1_prop"] <- "4"
LAR_data$census[LAR_data$census == "S2_LAR2_prop"] <- "5"

ggplot(LAR_data, aes(x= LAR))+ geom_histogram(color="black", fill="white")+ facet_grid(census ~  .)

LAR_data <- drop_na(LAR_data,LAR) 

n<-nrow(LAR_data)

#this is the beta transformation, which transforms all values of 0 to a small value.

#Reference: Smithson M, Verkuilen J (2006). "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71.

LAR_data $y_beta<- (LAR_data $LAR*(n-1) + 0.5)/n

hist(LAR_data $y_beta)

#You can check that you no longer have zero values (or values of 1, but I doubt that would be true with your data.

min(LAR_data $y_beta)

max(LAR_data $y_beta)



#Then, the analysis is:
library(betareg)
beta_model<- betareg ( y_beta ~elev*mat_avgLAR_pzero+census, data=LAR_data)
Anova(beta_model,type="III")
visreg(beta_model, "elev", by="mat_avgLAR_pzero", overlay=TRUE, scale="response")

visreg(beta_model, 'elev', by= "mat_avgLAR_pzero", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Damage by herbivores (beta transformation)", partial=FALSE,# cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255",""))) 



summary(beta_model)

plot(beta_model)
gamlss_mod<- gamlss (formula= y_beta ~elev*mat_treat+census + random(exp_ID)+ random(block)+random(genotype),family=BE(mu.link = "logit"), data=LAR_data,control = gamlss.control(n.cyc = 500)) #have to use logit link, dont need to specify, but i did here

summary(gamlss_mod)
plot(gamlss_mod)
Anova(gamlss_mod)


#pull the fitted values out and plot them in ggplot

newdf2 <- LAR_data %>% 
  
  mutate(fit.m = predict(gamlss_mod, re.form = NA),
         
         fit.c = predict(gamlss_mod, re.form = NULL), #all random effects
         
         resid = residuals(gamlss_mod))

##Convert fit.m and fit.c back to the proportional scale

newdf2 $fit.m_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $fit.c_trans<-1/(1+exp(-(newdf2 $fit.m))) 

newdf2 $resid_trans<-1/(1+exp(-(newdf2 $resid))) 



LAR_fig =ggplot(newdf2,aes(x= elev,y= LAR,shape= mat_treat, linetype= mat_treat,color= mat_treat, group= mat_treat)) + 
  
  geom_point(aes(shape= mat_treat),size=4)+scale_shape_manual(values = c(0,17,20)) +scale_x_continuous("Elevation (km)")+ scale_y_continuous("Leaf area removed by herbivores") +  geom_smooth(method = "lm", se = TRUE,formula=y~x) +
  
  #geom_line(aes(y= fit.m_trans, lty= mat_treat), size=0.8) +
  
  theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255","grey"))

LAR_fig

#*******************************************************************************
#### 3.Fitness models: Probability of Reproduction ~ mat treatment x treatment x elev #####
#*******************************************************************************

# Probability of reproduction 

gh_pf= ggplot(gh2, aes(x= elev,y= Overall_flowered, group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Probability of flowering")+ scale_x_continuous("Source Elevation")  
gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))


#repeated measures
repro_data<- gh2 %>% pivot_longer(cols=c('S1_reproduced', 'S2_reproduced'),
                                  names_to='season',
                                  values_to='reproduced')

repro_data <- dplyr::select(repro_data, reproduced, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, LAR_max_prop,mat_avgLAR_pzero)
repro_data$season[repro_data$season == "S1_reproduced"] <- "1"
repro_data$season[repro_data$season == "S2_reproduced"] <- "2"
repro_data$season<-as.numeric(repro_data$season)

mod_repro_REavgLAR<-glmmTMB(reproduced~elev*mat_avgLAR_pzero*LAR_avg_prop+season+(1|block)+(1|genotype)+(1|exp_ID),data=repro_data,family=binomial(link="logit"))

Anova(mod_repro_REavgLAR) #sig interaction btw mat_treat and elev

visreg(mod_repro_REavgLAR,"elev", by="mat_avgLAR_pzero", overlay=FALSE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg(mod_repro_REavgLAR, 'elev', by= "mat_avgLAR_pzero", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Source Elevation (Km)", ylab="Probability of flowering", partial=FALSE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col=grey(c(1), alpha=0.1)),
       line=list(col=c("#6699cc","#882255")),
       points=list(cex=1.5,col=c("#6699cc","#882255"))) 

visreg2d(mod_repro_REavgLAR,"elev","LAR_avg_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(10),zlim=c(0,1),ylim=c(0,.25))

#*******************************************************************************
#### 5.Fitness models: Fecundity ~ garden x treatment x elev #####
#*******************************************************************************

repro <- filter(gh2, Overall_Mature_length_siliques > 0 )

gh_pf= ggplot(repro, aes(x= elev,y= round(Overall_Mature_length_siliques), group= treatment, 
                         colour= treatment))+geom_point(size=5) + scale_y_continuous("Fecundity(summed silique length)")+ scale_x_continuous("Source Elevation")  

gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))

vioplot(round(Overall_Mature_length_siliques) ~ treatment, data= repro, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255","grey"), col=c("#6699cc","#882255","grey"), ylab="Fecundity(summed silique length)", xlab="Herbivore treatment")+stripchart(round(Overall_Mature_length_siliques) ~ mat_treat, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)



#repeated measures
fruit_data<- gh2 %>% pivot_longer(cols=c('S1_Mature_length_siliques', 'S2_Mature_length_siliques'),
                                  names_to='season',
                                  values_to='fruit_length')

fruit_data <- dplyr::select(fruit_data, fruit_length, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, season,LAR_avg_prop, LAR_max_prop,mat_avgLAR_pzero)
fruit_data$season[fruit_data$season == "S1_Mature_length_siliques"] <- "1"
fruit_data$season[fruit_data$season == "S2_Mature_length_siliques"] <- "2"
fruit_data$season<-as.numeric(fruit_data$season)

#hurdle model in one step in the glmmTMB package. 

hist(gh2$Overall_Mature_length_siliques)

#repeated measures
hurdle_Model_repeated <- glmmTMB(fruit_length ~LAR_avg_prop*elev*mat_avgLAR_pzero+season + (1|block)+ (1| genotype)+ (1| exp_ID), data=fruit_data, zi=~LAR_avg_prop*elev*mat_avgLAR_pzero+season + (1|block)+ (1| genotype)+ (1| exp_ID),family=ziGamma(link="log"))


##This is the ANOVA table for the logistic regression part (probability of reproduction)
Anova(hurdle_Model_repeated,type="III", component="zi")

##This is the ANOVA table for the count part (fecundity amongst individuals that reproduced). 
Anova(hurdle_Model_repeated,type="III", component="cond")

require(performance)
require(DHARMa)

testDispersion(hurdle_Model_repeated)
plotQQunif(hurdle_Model_repeated)
plotResiduals(hurdle_Model_repeated)

summary(hurdle_Model_repeated)


visreg2d(hurdle_Model_repeated,"elev","LAR_avg_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,400, by=50))


zeroinflated <-predictorEffect("elev",  partial.residuals=FALSE, hurdle_Model_repeated)
plot(zeroinflated, lwd=2,xlab="Source elevation (km)", ylab="Fecundity", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,1500))



visreg(hurdle_Model_repeated,"elev", by="LAR_avg_prop", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Fecundity", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255", "grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))


# individual fecundity anaylsis
only_repro <- filter(fruit_data, fruit_length > 0 )

hist(only_repro$fruit_length)

mod_fecundity_avgLAR <- glmmTMB (round(fruit_length) ~ LAR_avg_prop*elev*mat_avgLAR_pzero+season  + (1| block)+(1|genotype)+(1|exp_ID), family=Gamma(link="log"), data = only_repro)

Anova(mod_fecundity_avgLAR)

