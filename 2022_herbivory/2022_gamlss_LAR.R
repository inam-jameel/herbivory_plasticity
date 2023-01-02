####LAR,#untested/outdated with summary file #### 

library(visreg)
library(ggplot2)

#plan is to do run the model for each of the cohort/gardens

data_LAR<- read.csv("2022_incomplete_LARsummary.csv")

hist(data_LAR$LAR_3)

#need to update with specific censuses to
data_LAR$LAR_prop<- data_LAR$LAR_3/100

##Check that the conversion worked

hist(data_LAR$LAR_prop)


p_lar= ggplot(data_LAR, aes(x= elevation,y= LAR_prop, group= Treatment, colour= Treatment))+geom_point(size=5) +scale_x_continuous("elevation")+ scale_y_continuous("Leaf damage") 
p_lar + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+geom_smooth(method="glm",size=1.6, formula=y~poly(x,2))+facet_wrap(~ Treatment, scales="free_x")


LAR_prop1<-ggplot(data_LAR, aes(x = garden, y = LAR_prop, fill = Treatment,)) +
        geom_boxplot(outlier.shape = NA) +xlab("Garden")+ scale_y_continuous("Leaf area removed (%)") +
        geom_point(pch = 21, position = position_jitterdodge())
LAR_prop1 + theme_bw() + theme(text = element_text(size=20),axis.line.x = element_line(colour = "black"), 
                                 axis.line.y = element_line(colour = "black"), panel.border = element_blank(), panel.grid.major =element_blank(), 
                                 panel.grid.minor=element_blank(),legend.position = "top")+ scale_x_discrete(labels=c("Gothic", "Schofield")) +  scale_fill_manual(values = c( "lightblue","darkred"), name = "Herbivore treatment", labels = c("Pesticide","Control"))

##Check to make sure that you have the necessary columns.

head(data_LAR)      

#Let's concatenate block and garden so that we don't have to worry about nesting block within garden

data_LAR $block_garden<-interaction(data_LAR $block, data_LAR $garden,sep = "_")

##some variables are being read as chr not factors

data_LAR$garden<-as.factor(data_LAR$garden)
data_LAR$genotype<-as.factor(data_LAR$Genotype)
data_LAR$treatment<-as.factor(data_LAR$Treatment)
data_LAR$block<-as.factor(data_LAR$block)
data_LAR $S_elev<-scale(data_LAR$elevation,center=TRUE, scale=TRUE)
data_LAR $elev_km<-data_LAR$elevation/1000
data_LAR $init_size<-scale(data_LAR$initial_size,center=TRUE, scale=TRUE)

#remove NA values for gamlss

data<-na.omit(data)

##Okay, it looks like there were no NAs in that column, so we can proceed to the zero-one inflated negative binomial regression. What I've written out here is the full model, which has convergence issues.

#LAR_prop~treatment*S_elev*garden+ random(block_garden)+random(genotype)
gamlss.modelA<- gamlss (formula=LAR_prop~treatment*S_elev*garden+ random(block)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev*garden+ random(block)+random(genotype), nu.formula=LAR_prop~treatment*S_elev*garden+ random(block)+random(genotype), tau.formula=LAR_prop~treatment*S_elev*garden+ random(block)+random(genotype), family=BEINF, data= data_LAR)


##Perhaps we don't want to interact everything. Please keep in mind that garden should be a fixed effect, not random (you only have 2 gardens), but maybe the patterns don't differ much across gardens:



#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.model<- gamlss (formula=LAR_prop~treatment*S_elev+garden+ random(block)+random(genotype), sigma.formula=LAR_prop~treatment*S_elev+garden+ random(block)+random(genotype), nu.formula=LAR_prop~treatment*S_elev+garden+ random(block)+random(genotype), tau.formula=LAR_prop~treatment*S_elev+garden+ random(block)+random(genotype), family=BEINF, data= data_LAR)

#many many warnings

gamlss.model
plot(gamlss.model)
plot(gamlss.modelA)
summary(gamlss.modelA)

#Stepwise analyis to find the best model (https://rdrr.io/cran/gamlss/man/stepGAIC.html). This process takes a couple of minutes.

dropterm(gamlss.modelA)
mod2<-stepGAIC(gamlss.modelA)
mod2$anova


#LAR_prop~treatment*S_elev+garden+ random(block_garden)+random(genotype)
gamlss.modelB<- gamlss (formula=LAR_prop~treatment+S_elev+garden+ random(block)+random(genotype), sigma.formula=LAR_prop~treatment+S_elev+garden+ random(block)+random(genotype), nu.formula=LAR_prop~treatment+S_elev+garden+ random(block)+random(genotype), tau.formula=LAR_prop~treatment+S_elev+garden+ random(block)+random(genotype), family=BEINF, data= data_LAR)

summary(gamlss.modelB) # there are significant effects of treatment, elevation of origin and their interaction for mu (the probability of being damaged). For nu (amount of damage on plants with >0 and <1 damage; i.e., the beta component), there is a significant effect of elevation of origin and garden, but not treatment or elevation by treatment). For tau (probability of 100% damage – which you likely don’t have), there are no significant effects.


#So, we might be able to remove the tau component from the model, and run a zero-inflated model:

gamlss.modelC<- gamlss (formula=LAR_prop~elev_km+treatment+garden+ random(block)+random(genotype), sigma.formula=LAR_prop~elev_km+treatment+garden+ random(block)+random(genotype), nu.formula=LAR_prop~elev_km+treatment+garden+ random(block)+random(genotype), family=BEINF, data= data_LAR)

plot(gamlss.modelC)

summary(gamlss.modelC)


visreg(gamlss.modelC,overlay = FALSE, "elev_km", by="treatment", type="conditional", #scale = "response", 
       xlab="Elev", ylab="LAR", partial=TRUE,
       fill=list(col="light grey"
                 #(c(0.99), alpha=0)
       ), band = FALSE,
       #line=list(col=grey(c(0.2,0.6))),
       points=list(cex=0.65,  pch=(19)))  

#concatenate for visreg
data_LAR$garden_treat<-interaction(data_LAR $treatment, data_LAR $garden,sep = "_")

gamlss.modelE<- gamlss (formula=LAR_prop~elev_km* garden_treat + random(block_garden)+random(genotype), sigma.formula=LAR_prop~elev_km* garden_treat+ random(block_garden)+random(genotype), nu.formula=LAR_prop~elev_km* garden_treat+ random(block_garden)+random(genotype), family=BEINF, data= data_LAR)

plot(gamlss.modelE)

summary(gamlss.modelE)

Anova(gamlss.modelE)

visreg(gamlss.modelE, overlay = FALSE, "elev_km", by="garden_treat", type="conditional", 
       #scale = "response",   
       xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = TRUE)    


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


## linear model, this is the incorrect fit

modB<- lmer (LAR_prop~elev_km*garden*treatment+(1|block_garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)


Anova(modA)

LAR <-predictorEffect("elev_km",  partial.residuals=TRUE, modA)
plot(LAR, lwd=2,xlab="Source Elevation (Km)", ylab="Leaf Area Herbivorized", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,0.3))

LAR <-predictorEffect("elev_km",  partial.residuals=TRUE, modB)
plot(LAR, lwd=2,xlab="Source Elevation (Km)", ylab="Leaf Area Herbivorized", pch=19, type="response",lines=list(multiline=FALSE, lty=2:1, col="black"), 
     partial.residuals=list(smooth=TRUE, pch=19, col="black"), ylim=c(0,0.3))


visreg(modB, overlay = FALSE, "elev_km", by="garden_treat", type="conditional", 
       #scale = "response",   
       xlab="Source elevation", ylab="Leaf area removed", partial=TRUE,  band = TRUE) 


##Now to get LSMEANS
modB<- lmer (LAR_prop~genotype*treatment*garden+(1|block_garden),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)), data=data_LAR)


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


