####LAR,#untested/outdated with summary file #### 
#plan is to do run the model for each of the cohort/gardens

data<- read.csv("2021_LAR_060222.csv")


#need to update with specific censuses to
data$LAR_prop<- data$LAR/100

##Check that the conversion worked

hist(data$LAR_prop)

##Check to make sure that you have the necessary columns.

head(data)      

#Let's concatenate block and garden so that we don't have to worry about nesting block within garden

data $block_garden<-interaction(data $block, data $garden,sep = "_")

##some variables are being read as chr not factors

data$garden<-as.factor(data$garden)
data$genotype<-as.factor(data$Genotype)
data$treatment<-as.factor(data$Treatment)
data$block<-as.factor(data$block)
data $S_elev<-scale(data$elevation,center=TRUE, scale=TRUE)
data $elev_km<-data$elevation/1000
data $init_size<-scale(data$initial_size,center=TRUE, scale=TRUE)

#remove NA values for gamlss

data<-na.omit(data)

##Okay, it looks like there were no NAs in that column, so we can proceed to the zero-one inflated negative binomial regression. What I've written out here is the full model, which has convergence issues.

#LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype)
gamlss.modelA<- gamlss (formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), tau.formula=LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype), family=BEINF, data= data)


##Perhaps we don't want to interact everything. Please keep in mind that garden should be a fixed effect, not random (you only have 2 gardens), but maybe the patterns don't differ much across gardens:



#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.model<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), tau.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), family=BEINF, data= data)

#many many warnings

gamlss.model
plot(gamlss.model)
summary(gamlss.modelA) #nu only has garden as significant

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.modelA)
mod2<-stepGAIC(gamlss.modelA)
mod2$anova

summary(mod2) # there are significant effects of treatment, elevation of origin and their interaction for mu (the probability of being damaged). For nu (amount of damage on plants with >0 and <1 damage; i.e., the beta component), there is a significant effect of elevation of origin and garden, but not treatment or elevation by treatment). For tau (probability of 100% damage – which you likely don’t have), there are no significant effects.


#So, we might be able to remove the tau component from the model, and run a zero-inflated model:

gamlss.modelB<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype), family=BEINF, data= data)

plot(gamlss.modelB)

summary(gamlss.modelB)


visreg(gamlss.modelB, overlay = TRUE, "S_elev", by="treatment", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))  

data$garden_treat<-interaction(data $treatment, data $garden,sep = "_")

gamlss.modelE<- gamlss (formula=LAR_prop~S_elev* garden_treat + random(block_garden)+random(genotype), sigma.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(genotype), nu.formula=LAR_prop~S_elev* garden_treat+ random(block_garden)+random(genotype), family=BEINF, data= data)

plot(gamlss.modelE)

summary(gamlss.modelE)

visreg(gamlss.modelE, overlay = FALSE, "S_elev", by="garden_treat", type="conditional", 
       #scale = "response",   
       xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = FALSE)    


#You can create a datafile with the predicted data from the nu component using this code for plotting:
pred<-predict(gamlss.modelB, newdata=data, type="response", what="nu")

str(gamlss.modelB)

pred <- as.data.frame(pred)  
View(pred)

pred_data<- bind_cols(c2021, pred)

View(c2021)

LAR_2021 = ggplot(subset(pred_data,garden=="Gothic"), aes(x= elevation,y=pred))+geom_point(size=5) +scale_x_continuous("Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = factor(treatment)), size = 4)
P <-LAR_2021 + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                            axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                            panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)
P



## linear model, this is the incorrect fit  but nothing else is working
##Now to get LSMEANS
modB<- lmer (LAR_prop~genotype*treatment*garden+(1|block_garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data)


fam_avg_linear<-emmeans(modB, ~genotype:treatment:garden)
fam_means_linear<-as.data.frame(fam_avg_linear)

write.csv(fam_means_linear,file="LSmeans_LAR3.csv")


LSmeans_LAR <- read.csv("LSmeans_LAR3.csv", stringsAsFactors=FALSE)

#now we have lsmeans, but need to add corresponding elevation
#elev <- data[c("genotype","elevation")] #make dataframe of genotypes and corresponding elev
#elev <- unique(elev) #calls unique rows 
#LSmeans_LAR <- merge(LSmeans_LAR,elev,by="genotype") #merge the dataframes
#LSmeans_LAR$elev<-LSmeans_LAR$elevation/1000

#write.csv(LSmeans_leafarea,file="LSmeans_LAR3.csv")

LSmeans_LAR$emmean100 <- LSmeans_LAR$emmean*100

mod <- lmer (emmean~treatment*elev*garden+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=LSmeans_LAR)

Anova(mod)

visreg(mod, overlay = TRUE, "elev", by="treatment", type="conditional", scale = "response", 
       xlab="Elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
       fill=list(col="blue"),
       line=list(col=grey(c(0,0.8))),
       points=list(cex=1.5,col=grey(c(0,0.8))))  


#estess 
Estess <- subset(LSmeans_LAR,garden=="Estess")

estess_mod <- lmer (emmean100~treatment*elev+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Estess)
plot(estess_mod)
Anova(estess_mod,type="III")

E1 <- visreg(estess_mod, overlay = TRUE, "elev", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 30),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  


#trying gamlss model
E_gam<-na.omit(Estess) #get rid of NAs
Estess$genotype <- as.factor(Estess$genotype) #random needs to be a factor
gamlss.model_estess<- gamlss (formula=emmean~elev*treatment+random(genotype), family=BEZI, data= Estess) #response out of range

summary(gamlss.model_estess)


#gothic
Gothic <- subset(LSmeans_LAR,garden=="Gothic")

gothic_mod <- lmer (emmean100~treatment+elev+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Gothic)
plot(gothic_mod)
Anova(gothic_mod,type="III")

G1 <- visreg(gothic_mod, overlay = TRUE, "elev", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 12),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  

#schofield

Scho <- subset(LSmeans_LAR,garden=="Schofield")
scho_mod <- lmer (emmean100~treatment+elevation+(1|genotype),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=Scho)
plot(scho_mod)
Anova(scho_mod,type="III")


S1 <- visreg(scho_mod, overlay = TRUE, "elevation", by="treatment", type="conditional", scale = "response", 
             xlab="Source elevation (m)", ylab="Percentage leaf area herbivorized", partial=TRUE, cex.lab = 1.5,cex.axis = 1.5,
             fill=list(col="blue"),
             ylim=c(0, 12),
             line=list(col=grey(c(0,0.8))),
             points=list(cex=1.5,col=grey(c(0,0.8))))  



Estess_com <- subset(data_LAR,garden=="Estess")
E_LAR = ggplot(Estess_com, aes(x= elevation,y=LAR))+geom_point(size=5) +scale_x_continuous(" Source elevation")+ scale_y_continuous("Percentage leaf area herbivorized") + geom_point(aes(colour = treatment), size = 4)
E <- E_LAR + scale_color_grey() + theme_classic() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                                          axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                                          panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(aes(group = treatment),method="glm",size=2.0, se=F)+geom_smooth(aes(group = treatment, colour = treatment),method="glm",size=1.6, se=F)


