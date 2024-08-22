
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/greenhouse/data")


#read in data 
matseeds <- read.csv("maternal_effects_seeds.csv",stringsAsFactors=T)

sapply(matseeds,class)
##Some  variables are being read as characters not factors. Let's fix that
matseeds$exp_ID<-as.factor(matseeds$exp_ID)

#This rescales source elevation from meters to km
matseeds$elev_km<-matseeds $elevation/1000

##Change the base line for offspring
matseeds $treatment <-factor(matseeds $treatment, levels = c("Herbivory", "Naïve"))

#set colors
cols=c("#882255","#56B4E9")


mod1 <- glmmTMB(weight ~ treatment*elev_km+(1|exp_ID)+(1|block)+(1|genotype), data = matseeds)
simulationOutput <- simulateResiduals(fittedModel= mod1, plot = T, re.form = NULL,allow.new.levels =T)



#mod1 <- lmer(weight~ treatment*elev_km + (1|exp_ID), data=matseeds )
Anova(mod1, type = "III")
plot(mod1)


mod2 <- glmmTMB(weight ~ avg_dam_s2*elev_km+(1|exp_ID)+(1|block)+(1|genotype), data = matseeds, family=lognormal(link="log"))



pred_seed <- ggpredict(mod1, terms = c("elev_km[all]"), type = "re", interval="confidence")

seed_cline <-plot(pred_seed, show_data=TRUE, show_title =FALSE, show_legend=FALSE, colors = cols, facet=TRUE)+theme_classic()+theme(legend.position="right")+scale_x_continuous("Source elevation")+ scale_y_continuous("Seed Weight")


simulationOutput <- simulateResiduals(fittedModel= mod1, plot = T, re.form = NULL,allow.new.levels =TRUE)


Anova(mod2, type = "III")


weight_means<-emmeans(mod1, ~ treatment*elev_km, type="response", adjust = "sidak")

plot(weight_means, comparisons = TRUE)





