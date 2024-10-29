#*******************************************************************************
#### Herbivore resistance #######
#*******************************************************************************

##* Censuses 1 and 2 occurred for all years. Census 3 was for 2021 and 2022 only

##reformat datafile

HR_data_long_form<- grasshopper %>% pivot_longer(cols=c("HR_1","HR_2","HR_3"),
                                                  names_to='census',
                                                  values_to='HR')

HR_data_long_form <- dplyr::select(HR_data_long_form, HR, elevation, Genotype, population, Cage, Water, Herbivore, Block, PlantID, init.diam, S_initdiam, Cage_Block, elev_km, S_elev,census, year)

HR_data_long_form$census[HR_data_long_form$census == "HR_1"] <- "1"
HR_data_long_form$census[HR_data_long_form$census == "HR_2"] <- "2"
HR_data_long_form$census[HR_data_long_form$census == "HR_3"] <-"3"

HR_data_long_form $census <-as.factor(HR_data_long_form $census)

HR_data_long_form $year <-as.factor(HR_data_long_form $year)
##Let's concatenate census and year
HR_data_long_form $census_year<-interaction(HR_data_long_form$census, HR_data_long_form$year,sep = "_")


HR_data_long_form$HR_prop<-HR_data_long_form $HR/100
hist(HR_data_long_form$HR_prop)
HR_data_long_form <- drop_na(HR_data_long_form,HR_prop) 

n<-nrow(HR_data_long_form)

##this is the beta transformation, which transforms all values of 0 to a small value.

HR_data_long_form $y_beta<- (HR_data_long_form $HR_prop*(n-1) + 0.5)/n

hist(HR_data_long_form $y_beta)

min(HR_data_long_form $y_beta)

max(HR_data_long_form $y_beta)

##Conversation of standardized elevation back to raw units
#-1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#0*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)
#1*sd(grasshopper $elevation,na.rm=TRUE)+mean(grasshopper $elevation,na.rm=TRUE)


#fourway_model
HR_Model_four<- gamlss (formula= y_beta ~S_elev*Water*Herbivore*year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=HR_data_long_form,control = gamlss.control(n.cyc = 500))
plot(HR_Model_four)
summary(HR_Model_four)
drop1(HR_Model_four)

save(HR_Model_four, file='LAR_Model.rda')   


##Set tick marks on the X axis for these elevations: 2600, 3100, 3600
#(2600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3100-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)
#(3600-mean(LAR_data_long_form $elevation, na.rm = TRUE))/sd( LAR_data_long_form $elevation,na.rm = TRUE)

## creating plot for the cline using visreg, takes a moment to run
damage_cline_H<-visregList(
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2021"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2021"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2022"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2022"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition",year="2023"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), 
  visreg(HR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal",year="2023"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

#damage_cline<-visregList(visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE), visreg(LAR_Model_four,"S_elev", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE),collapse=TRUE)

damfit_H<-data.frame(damage_cline_H$fit)
damfit_H$visregFITexp<-1-(exp((damfit_H$visregFit)*-1)) # just exponentiating the visreg values made the the values >100, so multiplying them by -1 first, and then subtracting by 1

damage_clinal_variation_H <-ggplot(damfit_H, aes(S_elev, visregFITexp,group= Water, colour= Water, fill=factor(Water))) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+ geom_point(data= HR_data_long_form, aes(S_elev, HR_prop, color= Water, shape=Water), alpha=0.50, color="black")+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_x_continuous("Source elevation (m)",breaks=c(-1.557947,-0.1085093, 1.340929))+ scale_y_continuous("Herbivore resistance (proportion)")+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+ geom_line(aes(group=Water),linewidth=1)+scale_linetype_manual(values=c("solid", "dotted"))+facet_grid(Herbivore~year)

damage_clinal_variation_H #Fig 2A

##Box_plot
HR_box <-ggplot(HR_data_long_form, aes(x = Herbivore, y = HR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous("Leaf area removed by herbivores (proportion)") +geom_point(aes(shape=factor(Water)), size = 1,position = position_jitterdodge(0.3))

HR_box <-HR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(~year)
HR_box


## grid format
HR_box <-ggplot(HR_data_long_form, aes(x = Water, y = HR_prop, fill = Water,shape=Water)) +geom_boxplot(outlier.shape = NA) +xlab("Grasshopper treatment")+ scale_y_continuous("Herbivore resistance (proportion)")+geom_point(aes(shape=factor(Water)), size = 1.5,position = position_jitterdodge(0.3))

HR_box <-HR_box + theme_classic() + theme(text = element_text(size=10),axis.line.x = element_line(colour = "black"), 
                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                            panel.grid.minor=element_blank(),legend.position = "none")+ scale_x_discrete(labels=c("grasshopper addition", "grasshopper removal")) +  scale_fill_manual(values = c("#CC79A7","lightblue"), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+scale_shape_manual(values=c(21,24), name = "Water treatment", labels = c("Water-restricted","Supplemental watering"))+facet_grid(Herbivore~year)#
HR_box




##env concatenates water and herbivore. This model allows us to extract means and slopes for each combination of treatment and elevation
HR_data_long_form $env<-interaction(HR_data_long_form $Water, HR_data_long_form $Herbivore)
HR_data_long_form $env_year<-interaction(HR_data_long_form $env, HR_data_long_form $year)


#LAR_Model_three<- gamlss (formula= y_beta ~S_elev* env+ year-1 +random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
#summary(LAR_Model_three)
#save(LAR_ModelE, file='LAR_ModelE.rda')   


LAR_ModelE_four<- gamlss (formula= y_beta ~S_elev* env_year-1 +random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=LAR_data_long_form,control = gamlss.control(n.cyc = 500))
summary(LAR_ModelE_four)




##Subsets of models for drop 1
HR_Model_three<- gamlss (formula= y_beta ~S_elev*Water*Herbivore
                          +S_elev*Water*year
                          +S_elev* Herbivore*year
                          + Water* Herbivore*year
                          + year
                          + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=HR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(HR_Model_three)
#save(LAR_Model_three, file='LAR_Model_three.rda')   

##Subsets of models for drop 1
HR_Model_two<- gamlss (formula= y_beta ~S_elev*Water
                        + S_elev*Herbivore
                        + Water* Herbivore
                        + year*Water
                        + year*Herbivore
                        + S_elev*year
                        + year
                        + random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=HR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_two)
#save(LAR_Model_two, file='LAR_Model_two.rda')  


HR_Model_one<- gamlss (formula= y_beta ~S_elev+Water+Herbivore+ year+ random(census_year) + random(PlantID)+ random(Cage_Block)+random(Genotype),family=BE(mu.link = "logit"), data=HR_data_long_form,control = gamlss.control(n.cyc = 500))
drop1(LAR_Model_one)
#save(LAR_Model_one, file='LAR_Model_one.rda')


#*******************************************************************************
####   Hypothesis Two: Selection via two fitness components on traits   #####
#*******************************************************************************

##Exclude 2021 because we only have LAR data from that year and only 2 plants Reproduced
grasshopper_no2021<-subset(grasshopper, year!="2021")
##retain only those traits to be included in the models;
colnames(grasshopper);

traitdat <- dplyr::select(grasshopper_no2021,FT_Adj,Max_height_flowering, avg_HR, Genotype, Water, Herbivore, PlantID, init.diam, Cage_Block,  year, Mature_length_siliques,Reproduced,flowering_duration, SLA, succulence, FT_Adj)



##This standardizes initial plant size (measured as diameter in mm)  to a mean of 0 and standard deviation of 1 for use as a covariate in fitness models
traitdat $S_initdiam<-scale(traitdat $init.diam,center=TRUE, scale=TRUE)

## Many quantitative genetic models have convergence issues (or run very slowly) using raw data because traits and fitness components are measured on different scales. For example, phenology could be measured in days, whereas egg or seed mass is measured in mg. It is generally useful to standardize traits to a mean of 0 and standard deviation of 1. Below is code for standardizing flowering phenology (the resulting standardized variable is sFP, for standardized flowering phenology) and other phenotypic traits. For leaf damage, the standardized variable is called sLAR (which uses our field abbreviation of LAR for leaf area removed by herbivores)

traitdat $sHR<-scale(traitdat $avg_HR,center=TRUE,scale=TRUE)
traitdat $sSLA<-scale(traitdat $SLA,center=TRUE,scale=TRUE)
traitdat $sSUC<-scale(traitdat $succulence,center=TRUE,scale=TRUE)


head(traitdat)


##Change baseline for plotting purposes
traitdat $Water<-factor(traitdat $Water, levels = c("Restricted","Supplemental"))
traitdat $Herbivore<-factor(traitdat $Herbivore, levels = c("Addition","Removal"))


#*******************************************************************************
####  Selection using Probability of reproduction for vegetative traits only  #####
#*******************************************************************************

repro_model <-glmmTMB(Reproduced~S_initdiam+Water*Herbivore+year+
                        Water*Herbivore*sHR+ 
                        Water*Herbivore*sSUC+ I(sSUC^2)+
                        Water*Herbivore*sSLA+ I(sSLA^2) 
                      +(1|PlantID)+(1|Cage_Block)+(1|Genotype),data=traitdat,family=binomial(link="logit"))
Anova(repro_model,type="III") 
#summary(repro_model)

#*******************************************************************************
####  Selection using Seed production for all traits #####
#*******************************************************************************

##Create datafile with only reproductive plants
traitdatRepro <- filter(traitdat, Reproduced == 1 )

##Rescale after removing the non-reproductive plants
traitdatRepro $sduration<-scale(traitdatRepro $flowering_duration,center=TRUE,scale=TRUE)
traitdatRepro $s_max_height<-scale(traitdatRepro $Max_height_flowering,center=TRUE,scale=TRUE)
traitdatRepro $sHR<-scale(traitdatRepro $avg_HR,center=TRUE,scale=TRUE)
traitdatRepro $sSLA<-scale(traitdatRepro $SLA,center=TRUE,scale=TRUE)
traitdatRepro $sSUC<-scale(traitdatRepro $succulence,center=TRUE,scale=TRUE)
traitdatRepro $sFT<-scale(traitdatRepro $FT_Adj,center=TRUE,scale=TRUE)
traitdatRepro $S_initdiam <-scale(traitdatRepro $init.diam,center=TRUE,scale=TRUE)

fecund_modeltraits <-glmmTMB(Mature_length_siliques~ S_initdiam+Water*Herbivore+year +
                               Water*Herbivore*sFT+
                               Water*Herbivore*s_max_height +
                               Water*Herbivore* sduration +Water*I(sduration^2) +Herbivore*I(sduration^2)+
                               Water*Herbivore* sSLA+ Water*Herbivore* I(sSLA^2)+
                               Water*Herbivore* sSUC+ 
                               Water*Herbivore* sHR+Water*Herbivore*I(sHR^2)+(1|Cage_Block)+(1|Genotype),data=traitdatRepro,family=Gamma(link="log"))
Anova(fecund_modeltraits,type="III") 

##To check residuals
simulationOutput <- simulateResiduals(fittedModel= fecund_modeltraits, plot = T, re.form = NULL,allow.new.levels =TRUE)


##To check residuals
simulationOutput <- simulateResiduals(fittedModel= repro_model, plot = T, re.form = NULL)


## Leaf area removed
##Context dependent selection

HR<-visregList(visreg(fecund_modeltraits,"sHR", by="Water",cond=list("Herbivore"="Addition"), overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"),
                visreg(fecund_modeltraits,"sHR", by="Water",cond=list("Herbivore"="Removal"),overlay = FALSE, partial = FALSE, rug = FALSE,plot=FALSE,scale="response"), collapse=TRUE)


HR_fecund_plot<-ggplot(HR $fit, aes(sHR, visregFit,group= Water, colour= Water, fill=factor(Water))) +
  geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), alpha=0.15, linetype=0) +
  geom_line(aes(group=Water)) +theme_classic()+theme(text = element_text(size=10), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "none")+  geom_point(data= traitdatRepro, aes(sHR, Mature_length_siliques, color= Water, shape=Water), alpha=0.75, color="black")+scale_shape_manual(values=c(21,24))+scale_linetype_manual(values=c("dashed","solid"))+scale_x_continuous("Herbivore resistance (proportion)")+ scale_y_continuous("Fecundity (total fruit length, mm)",limits=c(0,2000))+scale_fill_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+scale_colour_manual(values = cols, name = "Water treatment", labels = c("Restricted","Supplemental"))+facet_wrap(~Herbivore)
HR_fecund_plot

figure2 <- ggarrange(damage_clinal_variation_H, 
                     HR_box,HR_fecund_plot,
                     labels = c("A", "B","C"),
                     ncol = 3, nrow = 1, common.legend = TRUE, legend="none")
figure2
