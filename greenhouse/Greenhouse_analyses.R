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
require(tidyr) # for data wrangling
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

gh2 $S_elev<-scale(gh2$elevation,center=TRUE, scale=TRUE)
gh2$elev<-gh2$elevation/1000
gh2 $S_init_size<-scale(gh2$ini_size,center=TRUE, scale=TRUE)



# convert LAR  to proportion

#summary
gh2$LAR_avg_prop<- gh2$LAR_avg/100
hist(gh2$LAR_avg_prop) 

gh2$LAR_max_prop<- gh2$LAR_max/100 
hist(gh2$LAR_max_prop) 

gh2$LAR_avg_S1_prop<- gh2$LAR_avg_S1/100
hist(gh2$LAR_avg_S1_prop) 

gh2$LAR_max_S1_prop<- gh2$LAR_max_S1/100
hist(gh2$LAR_max_S1_prop)

gh2$LAR_avg_S2_prop<- gh2$LAR_avg_S2/100
hist(gh2$LAR_avg_S2_prop) 

gh2$LAR_max_S2_prop<- gh2$LAR_max_S2/100 
hist(gh2$LAR_max_S2_prop) 



#season 2
gh2$S2_LAR1_prop<- gh2$S2_LAR1/100
hist(gh2$S2_LAR1_prop) 

gh2$S2_LAR2_prop<- gh2$S2_LAR2/100 #most recent
hist(gh2$S2_LAR2_prop) 


#season 1

gh2$S1_LAR1_prop<- gh2$S1_LAR1/100 
hist(gh2$S1_LAR1_prop)

gh2$S1_LAR3_prop<- gh2$S1_LAR3/100 #most damage
hist(gh2$S1_LAR3_prop) 

gh2$S1_LAR5_prop<- gh2$S1_LAR5/100 
hist(gh2$S1_LAR5_prop) 

gh2$S1_LAR13_prop<- gh2$S1_LAR13/100 #last in S1
hist(gh2$S1_LAR13_prop) 



##Some of your variables are being read as characters not factors. Let's fix that
gh2$genotype<-as.factor(gh2$genotype)
gh2$treatment<-as.factor(gh2$treatment)
gh2$block<-as.factor(gh2$block)
gh2$mat_treat<-as.factor(gh2$mat_treat)
gh2$mat_exp_ID <-as.factor(gh2$mat_exp_ID) #need to include this as random effect since multiple reps per mat line
gh2 $treat_mat<-interaction(gh2 $treatment, gh2 $mat_treat,sep = "_")
#gh2$trichomes <- gh2$avg_trichomes +0.001 #added to do gamma glmer


#*******************************************************************************
#### 1.SLA for greenhouse experiment #####
#*******************************************************************************

#### SLA individual for full dataset ####

SLA_S1 =ggplot(gh2, aes(x= elevation,y= S1_SLA, color = treat_mat))+ geom_point(size=3.5)+ ggtitle("SLA season 1 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S1

SLA_S2 =ggplot(gh2, aes(x= elevation,y= S2_SLA, color = treat_mat))+ geom_point(size=3.5)+ ggtitle("SLA season 2 by Souce Elevation")+scale_x_continuous("Source Elevation")+ scale_y_continuous("SLA") +theme_bw()+theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="lm", formula=y~x,  se=TRUE, size=1.6)
SLA_S2

#SLA for first season is S1_SLA
#also have to filter for the indiv that we have appropriate values for SLA 
SLA_data_S1=filter(gh2,S1_include_SLA == 1)

hist(SLA_data_S1$S1_SLA)

simple_SLA<- lmer(S1_SLA ~ treatment*mat_treat+elev+(1|block)+(1|genotype),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data_S1)
Anova(simple_SLA)
plot (simple_SLA)

visreg(simple_SLA,"mat_treat", by="treatment", overlay=TRUE,   scale = "response", xlab="Maternal environment", ylab="SLA", partial=TRUE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))

summary(glht(simple_SLA, linfct = mcp(mat_treat = "Tukey")))  

SLAmeans <- gh2 %>% 
  group_by(treatment) %>%
  get_summary_stats(S2_SLA, type="mean_se")
SLAmeans

vioplot(S1_SLA ~ treatment, data= SLA_data_S1, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255"), col=c("#6699cc","#882255"), ylab="SLA CM3/g", xlab="Herbivore treatment")+stripchart(S2_SLA ~ treatment, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)

vioplot(S1_SLA ~ mat_treat, data= SLA_data_S1, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255","grey"), col=c("#6699cc","#882255","grey"), ylab="SLA CM3/g", xlab="Maternal herbivore treatment")+stripchart(S1_SLA ~ mat_treat, data= SLA_data_S1,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)


emmip(simple_SLA, ~ treatment, type="response", CIs=TRUE)+theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank())+geom_point(size=3)+  ylab("Probability of survival")+ xlab("Minimum temperature in treatment")+scale_colour_grey()+ ylim(150,250)+ theme(legend.position = c(0.90, 0.85))



#SLA for second season is S2_SLA
#also have to filter for the indiv that we have appropriate values for SLA 
SLA_data_S2=filter(gh2,S2_include_SLA == 1)


hist(SLA_data_S2$S2_SLA)

simple_SLA<- lmer(S2_SLA ~ treatment*mat_treat+elev+(1|block)+(1|genotype),control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data_S2)
Anova(simple_SLA)
plot (simple_SLA)

visreg(simple_SLA,"elev", by="treatment", overlay=TRUE,   scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Final size (leaf number)", partial=TRUE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))

summary(glht(simple_SLA, linfct = mcp(treatment = "Tukey"))) 

SLAmeans <- gh2 %>% 
  group_by(treatment) %>%
  get_summary_stats(S2_SLA, type="mean_se")
SLAmeans

vioplot(S2_SLA ~ treatment, data= SLA_data_S2, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255","grey"), col=c("#6699cc","#882255","grey"), ylab="SLA CM3/g", xlab="herbivore treatment")+stripchart(S2_SLA ~ treatment, data= SLA_data_S1,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)





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


#converting wide to long for repeated measures

SLA_data<- gh2 %>% pivot_longer(cols=c('S2_SLA', 'S1_SLA'),
                 names_to='year',
                 values_to='SLA')

SLA_data <- dplyr::select(SLA_data, SLA, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block, year)

ggplot(SLA_data, aes(x= SLA))+ geom_histogram(color="black", fill="white")+ facet_grid(treatment ~  .)

SLA_RM <- lmer(SLA ~ treatment*elev*mat_treat*year 
                     +(1|block)
                     +(1|genotype)
                    +(1|exp_ID)
                     ,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)), data = SLA_data)

plot(SLA_RM)

Anova(SLA_RM)

visreg(SLA_RM,"elev", by="treatment", overlay=TRUE,   scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Final size (leaf number)", partial=FALSE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))



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

elev_avgLARbox<-ggplot(herb, aes(x = treatment, y = LAR_avg_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_avgLARbox + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                               panel.grid.minor=element_blank(),legend.position = "top")+ #scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +  
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))



#max
elev_maxLAR =ggplot(herb, aes(x= elevation,y= LAR_max_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")
elev_maxLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_maxLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

elev_maxLARbox<-ggplot(herb, aes(x = treatment, y = LAR_max_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_maxLARbox + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                               axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                               panel.grid.minor=element_blank(),legend.position = "top")+# scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))


#season 1

#histograms of avg damage and max damage
ggplot(herb, aes(x= LAR_avg_S1_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)
ggplot(herb, aes(x= LAR_max_S1_prop))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)


#avg
elev_avgS1LAR =ggplot(herb, aes(x= elevation,y= LAR_avg_S1_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")
elev_avgS1LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_avgS1LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

elev_avgS1LARbox<-ggplot(herb, aes(x = treatment, y = LAR_avg_S1_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_avgLARbox + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                    axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                    panel.grid.minor=element_blank(),legend.position = "top")+ #scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +  
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))



#max
elev_maxS1LAR =ggplot(herb, aes(x= elevation,y= LAR_max_S1_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")

elev_maxLAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_maxS1LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

elev_maxLARS1box<-ggplot(herb, aes(x = treatment, y = LAR_max_S1_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_maxLARbox + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                    axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                    panel.grid.minor=element_blank(),legend.position = "top")+# scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))

#season 2

#histograms of avg damage and max damage
ggplot(herb, aes(x= LAR_avg))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)
ggplot(herb, aes(x= LAR_max))+ geom_histogram(color="black", fill="white")+ facet_grid(mat_treat ~ .)


#avg
elev_avgS2LAR =ggplot(herb, aes(x= elevation,y= LAR_avg_S2_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")
elev_avgS2LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_avgS2LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

elev_avgS2LARbox<-ggplot(herb, aes(x = treatment, y = LAR_avg_S2_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_avgS2LARbox + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                    axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                    panel.grid.minor=element_blank(),legend.position = "top")+ #scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +  
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))



#max
elev_maxS2LAR =ggplot(herb, aes(x= elevation,y= LAR_max_S2_prop,shape= mat_treat, color= mat_treat,linetype= mat_treat))+geom_point(aes(shape= mat_treat),size=4)+scale_x_continuous("Source Elevation")+scale_y_continuous("Leaf area herbivorized")

elev_maxS2LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",se=FALSE,  size=1)

elev_maxS2LAR +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#6699cc","#882255", "grey"))+geom_smooth(method="lm",  size=1,formula=y~poly(x,2))

elev_maxLARS2box<-ggplot(herb, aes(x = treatment, y = LAR_max_S2_prop, fill = mat_treat,)) +
  geom_boxplot(outlier.shape = NA) +xlab("Herbivore treatment")+ scale_y_continuous("Leaf area removed (%)") +
  geom_point(pch = 21, position = position_jitterdodge())
elev_maxLARS2box + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                    axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                    panel.grid.minor=element_blank(),legend.position = "top")+# scale_x_discrete(labels=c("no herbivore", "herbivore addition")) +
  scale_fill_manual(values = c( "#6699cc","#882255", "grey"), name = "Maternal environment", labels = c("no herbivore", "herbivore addition","control"))


#analysis with LAR averaged across 2 seasons

LAR_data <- drop_na(herb,LAR_avg_prop) #this removes the rows without LAR values

LAR_data <- dplyr::select(LAR_data, LAR_avg_prop, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block)

#using BE since there are no zeros / ones for avg LAR
gamlss_mod<- gamlss (formula=LAR_avg_prop~elev+mat_treat+I(elev^2)+mat_treat +random(block)+random(genotype),family=BE(mu.link = "identity", sigma.link = "log"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(gamlss_mod)
plot(gamlss_mod)

visreg(gamlss_mod, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Source Elevation (Km)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255","grey"))) 

visreg(gamlss_mod, 'elev', type="conditional",
       #scale = "response", 
       xlab="Source elevation (Km)", ylab="Percentage leaf area herbivorized", partial=FALSE, cex.lab = 1.5,cex.axis = 1.5)


#lmer model
modA<- lmer(log(LAR_avg_prop)~ mat_treat+elev+I(elev^2)+ (1|block)+(1|genotype),data= LAR_data)
Anova(modA,type="III")
##You can see here that the residuals aren't great
plot(modA)

visreg(modA, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Source Elevation (Km)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255","grey"))) 

#proceeding with Max LAR

LAR_data <- drop_na(herb,LAR_max_prop) #this removes the rows without LAR values

LAR_data <- dplyr::select(LAR_data, LAR_max_prop, elev, genotype, treatment, mat_treat, exp_ID, S_init_size, block)

#using BE since there are no zeros / ones for MAX LAR
gamlss_mod<- gamlss (formula=LAR_max_prop~elev+mat_treat+I(elev^2)+mat_treat +random(block)+random(genotype),family=BE(mu.link = "identity", sigma.link = "log"), data=LAR_data,control = gamlss.control(n.cyc = 500))

summary(gamlss_mod)
plot(gamlss_mod)

visreg(gamlss_mod, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255","grey")))


visreg(gamlss_mod, 'elev', type="conditional",
       #scale = "response", 
       xlab="Source elevation (Km)", ylab="Percentage leaf area herbivorized", partial=FALSE, cex.lab = 1.5,cex.axis = 1.5)

#lmer model
modB<- lmer(LAR_max_prop~ mat_treat+elev+I(elev^2)+ (1|block)+(1|genotype),data= LAR_data)
Anova(modB,type="III")
##You can see here that the residuals aren't great
plot(modB)

visreg(modB, 'elev', by= "mat_treat", overlay = TRUE, type="conditional", 
       #scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255","grey")))




#*******************************************************************************
#### 3.Fitness models: Fitness ~ garden x treatment x elev #####
#*******************************************************************************

# Probability of reproduction 

gh_pf= ggplot(gh2, aes(x= elev,y= Overall_flowered, group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Probability of flowering")+ scale_x_continuous("Source Elevation")  
gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))




# Basic data cleaning and visualization packages
require(tidyverse)
require(ggplot2)
# Package for multinomial logistic regression
require(nnet)
# To transpose output results
require(broom)
# To create html tables
require(knitr)
require(kableExtra)
# To exponentiate results to get RRR
require(gtsummary)
# To plot predicted values
require(ggeffects)
# To get average marginal effect
require(marginaleffects)

ggplot(gh2 %>% filter(Flowered_seasons == 0 | 
                        Flowered_seasons == 1),
       aes(x = as.numeric(elevation), y = as.numeric(Flowered_seasons))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic()

ggplot(gh2 %>% filter(Flowered_seasons == 1 | 
                        Flowered_seasons == 2),
       aes(x = as.numeric(elevation), y = as.numeric(Flowered_seasons))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic()

fit_basic <- multinom(Flowered_seasons ~ elev*treatment*mat_treat, data = gh2)
summary(fit_basic)
Anova(fit_basic)

ggeffect(fit_basic, terms = "treatment") %>%
  plot()

multi <- lmer(Flowered_seasons~ elev+treatment*mat_treat + (1|block)+(1|genotype), data = gh2)
Anova(multi)
plot(multi)

visreg(multi,"treatment",  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255","grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))


#looking at overall flowered with treatment
mod_flower<- glmer(Overall_flowered~ elev*treatment*mat_treat+ S_init_size + (1|block)+(1|genotype), data = gh2, control=glmerControl(optimizer="optimx", tolPwrss=1e-3, optCtrl=list(method="nlminb")), family=binomial)
Anova(mod_flower) # all seperate factors are sig, sig interaction btw mat_treat and elev

visreg(mod_flower,"elev", by="treatment", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:2,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))



flower_model<-glmmTMB(Overall_flowered~elev*treatment*mat_treat + S_init_size+(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(flower_model,type="III")

visreg(flower_model,"elev", by="treatment", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))


#looking at overall flowered with average LAR
mod_flower_avgLAR<-glmmTMB(Overall_flowered~elev*LAR_avg_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_avgLAR) #sig interaction btw mat_treat and elev

visreg(mod_flower_avgLAR,"elev", by="LAR_avg_prop", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=TRUE,type="conditional",line=list(lty=1:3,col="black"), points=list(col="black"),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(mod_flower_avgLAR,"elev","LAR_avg_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(10),zlim=c(0,1),ylim=c(0,.25))

# average LAR x elev 2
mod_flower_avgLARquad<-glmmTMB(Overall_flowered~elev*LAR_max_prop*mat_treat+ I(elev^2)*LAR_max_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_avgLARquad) # slight significance

#looking at overall flowered with max LAR 

mod_flower_maxLAR<-glmmTMB(Overall_flowered~elev*LAR_max_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_maxLAR,type="III") #significant elev x lar

visreg(mod_flower_maxLAR,"elev", by="LAR_max_prop", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255")), points=list(col=c("#6699cc","#882255")),fill=list(col=grey(c(0.8), alpha=0.4)))

library(RColorBrewer)
visreg2d(mod_flower_maxLAR,"elev","LAR_max_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(10),zlim=c(0,1),ylim=c(0,1))

#looking at Season 1 flowered with season 1 avg  LAR
mod_flower_S1avgLAR<-glmmTMB(S1_flowered~elev*LAR_avg_S1_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S1avgLAR) #sig interaction btw mat_treat and elev

visreg(mod_flower_S1avgLAR,"elev", by="mat_treat", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255","grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(mod_flower_S1avgLAR,"elev","mat_treat", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,1, by=.05))

# average LAR x elev 2
mod_flower_S1avgLARquad<-glmmTMB(S1_flowered~elev*LAR_avg_S1*mat_treat+ I(elev^2)*LAR_avg_S1*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_avgLARquad) # individual factors are significant

#looking at Season 1 flowered with season 1 max LAR
mod_flower_S1maxLAR<-glmmTMB(S1_flowered~elev*LAR_max_S1_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S1maxLAR) #sig interaction btw mat_treat and elev

visreg(mod_flower_S1maxLAR,"elev", by="mat_treat", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255","grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(mod_flower_S1maxLAR,"elev","mat_treat", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,1, by=.05))

# max LAR x elev 2
mod_flower_S1maxLARquad<-glmmTMB(S1_flowered~elev*LAR_max_S1*mat_treat+ I(elev^2)*LAR_max_S1*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S1maxLARquad) # individual factors are significant

#looking at Season 2 flowered with season 2 avg  LAR
mod_flower_S2avgLAR<-glmmTMB(S2_flowered~elev*LAR_avg_S2_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S2avgLAR) #slight sig interaction btw mat_treat and elev

visreg(mod_flower_S2avgLAR,"elev", by="mat_treat", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Probability of flowering", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255","grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(mod_flower_S1avgLAR,"elev","mat_treat", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,1, by=.05))

# average LAR x elev 2
mod_flower_S2avgLARquad<-glmmTMB(S2_flowered~elev*LAR_avg_S2_prop*mat_treat+ I(elev^2)*LAR_avg_S2_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_avgLARquad) # individual factors are significant/ slightly significant

#looking at Season 2 flowered with season 2 max LAR
mod_flower_S2maxLAR<-glmmTMB(S2_flowered~elev*LAR_max_S2_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S2maxLAR) #nothing really

# max LAR x elev 2
mod_flower_S2maxLARquad<-glmmTMB(S2_flowered~elev*LAR_max_S2_prop*mat_treat+ I(elev^2)*LAR_max_S2_prop*mat_treat +(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(mod_flower_S2maxLARquad) # nothing really



# Probability of survival 
#nothing interesting here.. 

gh_pf= ggplot(gh2, aes(x= elev,y= S2_survival, group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Probability of survivial")+ scale_x_continuous("Source Elevation")

gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))

vioplot(S2_survival ~ treatment, data= gh2, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255"), col=c("#6699cc","#882255"), ylab="Survival", xlab="Herbivore treatment")+stripchart(S2_survival ~ treatment, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)
                                                                                                                                                      


mod_surv<- glmer(S2_survival~elev*treatment*mat_treat + S_init_size+ (1|block)+(1|genotype), data = gh2, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), family=binomial)
Anova(mod_surv) #only treatment is sig
plot(mod_surv)

survival_model<-glmmTMB(S2_survival~elev*treatment*mat_treat + S_init_size+(1|block)+(1|genotype),data=gh2,family=binomial(link="logit"))

Anova(survival_model,type="III")



# Fecundity

repro <- filter(gh2, Overall_Mature_length_siliques > 0 )

gh_pf= ggplot(repro, aes(x= elev,y= round(Overall_Mature_length_siliques), group= treatment, 
                       colour= treatment))+geom_point(size=5) + scale_y_continuous("Fecundity(summed silique length)")+ scale_x_continuous("Source Elevation")  

gh_pf + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                           axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(),
                           panel.grid.minor=element_blank(),legend.position = "top")+geom_smooth(method="glm",size=1.6, formula=y~x)+facet_wrap(~ mat_treat, scales="free_x") +scale_colour_manual(values = c( "#6699cc","#882255"), name = "Herbivore treatment", labels = c("No Herbivores","Herbivorized"))

vioplot(round(Overall_Mature_length_siliques) ~ treatment, data= repro, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#6699cc","#882255","grey"), col=c("#6699cc","#882255","grey"), ylab="Fecundity(summed silique length)", xlab="Herbivore treatment")+stripchart(round(Overall_Mature_length_siliques) ~ mat_treat, data= gh2,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)


mod_fecundity_avgLAR <- glmmTMB (Overall_Mature_length_siliques ~ LAR_avg_prop*elev*mat_treat  + (1| block)+(1|genotype), family=Gamma(link="log"),data = repro)

diagnose(mod_fecundity_avgLAR)

Anova(hurdle,type="III")

summary(hurdle)

## This is just another approach to visualizing these patterns:
fecundB<- glmer(Overall_Mature_length_siliques ~   LAR_avg_prop*elev*mat_treat+ (1|block) +(1|genotype), data= repro,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(fecundB,type="III")

fecundmod <-predictorEffect("elev",  partial.residuals=TRUE, fecundB)
plot(fecundmod, lwd=2,xlab="Source elevation (km)", ylab="Fecundity (number of fruits)", pch=19, type="response",lines=list(multiline=TRUE, lty=3:1, col="black"),
     partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,600))

visreg2d(fecundB,"elev","LAR_avg_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,750, by=50))


fruit <-predictorEffect("elev",  partial.residuals=TRUE, hurdle)

plot(fruit, lwd=2,xlab="Source Elevation (KM)", ylab="Fruit", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1,
                                                                                                 col="black"), partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1500))

visreg(hurdle,"mat_treat", by="LAR_max_prop", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Fecundity", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255", "grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg(hurdle,"elev", by="LAR_max_prop", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Fecundity", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255", "grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))

visreg2d(hurdle,"elev","LAR_max_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,750, by=50))



#library(glmmTMB)

#looking at overall

hurdle <- glmmTMB (Overall_Mature_length_siliques ~ LAR_max_prop*elev*mat_treat  + (1| block)+(1|genotype), zi=~ LAR_max_prop*elev*mat_treat  + (1| block)+(1|genotype),data = gh2)

hurdle <- glmmTMB (Overall_Mature_length_siliques ~ LAR_max_prop*elev*mat_treat  + LAR_max_prop*I(elev^2)*mat_treat +(1| block)+(1|genotype), zi=~ LAR_max_prop*elev*mat_treat  + LAR_max_prop*I(elev^2)*mat_treat (1| block)+(1|genotype),data = gh2)


diagnose(hurdle)

Anova(hurdle,type="III", component="zi")

summary(hurdle) #zero inflated model is significant, but not sure if i am plotting it?


fruit <-predictorEffect("elev",  partial.residuals=TRUE, hurdle,component="zi")

plot(fruit, lwd=2,xlab="Source Elevation (KM)", ylab="Fruit", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1,
                                                                                                 col="black"), partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,1500))

visreg(hurdle,"elev", by="mat_treat", overlay=TRUE,  scale = "response", xlab="Source elevation (Km)", ylab="Fecundity", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#6699cc","#882255", "grey")), points=list(col=c("#6699cc","#882255","grey")),fill=list(col=grey(c(0.8), alpha=0.4)))


library(RColorBrewer)
visreg2d(hurdle,"elev","LAR_max_prop", scale = "response", xlab="Source elevation (Km)", ylab="Leaf area herbivorized",col = colorRampPalette(brewer.pal(9,"Blues"))(20),zlim=c(0,500, by=50))




library(glmmTMB)

#looking at overll silique length


hurdle2 <- glmmTMB (Overall_Mature_length_siliques ~ treatment*mat_treat*elev+S_init_size  + (1| block)+(1|genotype), zi=~ treatment*mat_treat*elev+S_init_size  + (1| block)+(1|genotype),data = gh2 ,family=truncated_nbinom2)

diagnose(hurdle2)

Anova(hurdle2,type="III", component="zi")

summary(hurdle2)


fruit2 <-predictorEffect("elev",  partial.residuals=TRUE, hurdle2,component="zi")

plot(fruit2, lwd=2,xlab="Source Elevation (KM)", ylab="Fruit", pch=19, type="response",lines=list(multiline=FALSE, #lty=2:1,
                                                                                                  col="black"), partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,500))

plot(fruit2, lwd=2,xlab="Source Elevation (KM)", ylab="Fruit", pch=19, type="response",lines=list(multiline=TRUE, lty=2:1, col="black"), partial.residuals=list(smooth=FALSE, pch=19, col="black"), ylim=c(0,500))





modC<- lmer(Overall_Mature_length_siliques~ LAR_max_prop+mat_treat+elev+ (1|block)+(1|genotype),data= herb)
Anova(modC,type="III")
##You can see here that the residuals aren't great
plot(modA)

hist(herb$Overall_Mature_length_siliques)


visreg(modC, 'LAR_max_prop', by= "mat_treat", overlay = TRUE, type="conditional", 
       scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=c("#6699cc","#882255","grey")),
       points=list(cex=1.5,col=c("#6699cc","#882255","grey"))) 




reproLAR =filter(LAR_data, Overall_flowered== 1)


#converting wide to long for repeated measures

test %>% pivot_longer(cols=c('year1', 'year2'),
                    names_to='year',
                    values_to='points')






