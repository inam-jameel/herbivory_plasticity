######## PROJECT: Greenhouse experiment: Fitness and phenotypes in response to herbivory
#### PURPOSE:Examine fitness and traits in response to herbivory treatment.
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 9 Nov 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for plotting
require(glmmTMB) # for running survival model, have to load twice
require(gamlss) # for running phenology model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(vioplot) #for violin plots


# set working directory

setwd("~/OneDrive - University of Georgia/Inam_experiments/Herbivory_data/greenhouse/")  #this is where you specify the folder where you have the data on your computer

#read in data 
gh2 <- read.csv("gen2_GHSummary2022.csv")

gh2 <- filter(gh2, Include == "yes" )

#gh2 <- drop_na(gh2,LAR)
gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)



# convert LAR from census 3 (most damage) and census 13 (most recent) to proportion
gh2$LAR_prop2<- gh2$S2_LAR2/100
hist(gh2$LAR_prop2)  

##Some of your variables are being read as characters not factors. Let's fix that
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line
#gh2$trichomes <- gh2$avg_trichomes +0.001 #added to do gamma glmer


# Relative contributions

#### SLA individual for full dataset ####
#SLA for second season is S2_SLA
#also have to filter for the indiv that we have appropriate values for SLA 
gh2=filter(gh2,S2_include_SLA == 1)


hist(gh2$S2_SLA)

SLAgga<-ggplot(gh2, aes(x = mat_treat, y = S2_SLA, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) +xlab("Maternal treatment")+ scale_y_continuous("SLA CM3/g") +
  geom_point(pch = 21, position = position_jitterdodge())
SLAgga + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                            panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("No Herbivores", "Herbivorized", "Parent Control")) +  scale_fill_manual(values = c( "gold","purple"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))



gh_sla= ggplot(gh2, aes(x= elev,y= S2_SLA, group= treatment, 
                        colour= treatment))+geom_point(size=5) + scale_y_continuous("SLA CM3/g")+ scale_x_continuous("Source Elevation")  
gh_sla + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                            panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="lm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "gold","purple","grey"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))





RC_SLA<- lmer(S2_SLA ~ treatment*elev*mat_treat+S_init_size +(1|block)+(1|genotype),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = gh2)
Anova(RC_SLA)
plot (RC_SLA)

#Testing the random effect of genotype

RC_SLAnogeno <- lmer(S2_SLA ~ treatment*elev*mat_treat+S_init_size 
                     +(1|block)
                     #+(1|genotype)
                     ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = gh2)


anova(RC_SLA, RC_SLAnogeno)

##Testing the random effect of block (tray)
RC_SLAnoblock <- lmer(S2_SLA ~ treatment*elev*mat_treat+S_init_size 
                      #+(1|block)
                      +(1|genotype)
                      ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = gh2)

anova(RC_SLA, RC_SLAnoblock)



visreg(RC_SLA,"elev", by="treatment", overlay=TRUE,   scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Final size (leaf number)", partial=TRUE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))


vioplot(S2_SLA ~ treatment, data= gh2, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255"), col=c("#6699cc","#882255"), ylab="SLA CM3/g", xlab="Herbivore treatment")+stripchart(S2_SLA ~ treatment, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 
                                                                                                                                                                                                                                                      ,vertical = TRUE, add = TRUE)

summary(glht(RC_SLA, linfct = mcp(treatment = "Tukey"))) #not significant? 



#### SLA individual for full dataset ####




