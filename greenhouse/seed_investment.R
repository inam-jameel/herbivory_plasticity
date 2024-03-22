
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/Inam_experiments/Herbivory_data/greenhouse")


#read in data 
matseeds <- read.csv("maternal_effects_seeds.csv",stringsAsFactors=T)

sapply(matseeds,class)
##Some  variables are being read as characters not factors. Let's fix that
matseeds$exp_ID<-as.factor(matseeds$exp_ID)

#This rescales source elevation from meters to km
matseeds$elev_km<-matseeds $elevation/1000


mod1 <- lmer(weight~ LAR_max_S2*elev_km + (1|exp_ID), data=matseeds )
Anova(mod1, type = "III")
plot(mod1)


mod2 <- glmmTMB(weight ~LAR_max_S2*elev_km + (1|exp_ID), data= matseeds)


simulationOutput <- simulateResiduals(fittedModel= mod2, plot = T, re.form = NULL,allow.new.levels =TRUE)


Anova(mod2, type = "III")


weight_means<-emmeans(mod2, ~ treatment*elev_km, type="response", adjust = "sidak")

plot(weight_means, comparisons = TRUE)





