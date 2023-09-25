## ---------------------------
##
## Script name: 02_SEM_GlobNut_1.0
##
## Purpose of script: Calculating piecewise structural equation models
##
## Author: Leonardo H. Teixeira
##
## Date Created: 2023-09-22
##
## Email: leonardo.htp@gmail.com
##
## ---------------------------
##
## Notes: Structural equation models are calculated using the R-package piecewiseSEM
##  
## References: Lefcheck, J.S. (2016). PiecewiseSEM: Piecewise structural equation modeling in R for ecology, evolution, and systematics. Methods in Ecology and Evolution 7: 573–579. DOI:10.1111/2041-210X.12512
##
## ---------------------------


#calling Packages

library(tidyr)
library(plyr) # plyr has to be loaded before dplyr to avoid malfunctioning of dplyr functions
library(dplyr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggplot2)
library(piecewiseSEM)
library(nlme)
#library(sesem)
#library(semPlot)
#library(lavaan)
#library(lavaanPlot)

# 1. opening & filtering data----
# data nutrients
data.np <- read.csv("Data/Raw/globnut_SEM_bio.csv", sep=",", na.strings = "NA")

# replacing spaces by underline (in case needed)
#data.np$species <- gsub(" ","_",data.np$species)

# 2. opening species richness data----
data.richness <- read.csv("Data/Raw/richness.csv", sep=",")
#rownames(traits) <- traits$species

# 4. opening climate data----
data.clim <- read.csv("Data/Raw/climate.csv",sep=",", na.strings = "NA")

# 4. opening nitrogen deposition data----
data.n_dep <- read.csv("Data/Raw/N_deposition.csv",sep=",", na.strings = "NA")

# 5. Subset data nutrients for plots present in data richness----
data.sem <- data.np[which(data.np$plot_ID %in% data.richness$plot_ID),]

# 6. checking for duplicates in data.sem----
data.sem.test <- unique(data.sem)
duplicated(data.sem$plot_ID)
#duplicated(data.sem[1500:4259,2])

# 7. Subset data climate for plots present in data richness----
data.clim1 <- data.clim[which(data.clim$plot_ID %in% data.richness$plot_ID),]

# 8. checking for duplicates in data.clim1----
data.clim.test <- unique(data.clim1)
duplicated(data.clim1$plot_ID)

# 7. Subset data nitrogen deposition for plots present in data richness----
data.dep1 <- data.n_dep[which(data.n_dep$plot_ID %in% data.richness$plot_ID),]

# 8. checking for duplicates in data.dep1----
data.dep.test <- unique(data.dep1)
duplicated(data.dep1$plot_ID)

# 11. adding species richness, evenness, mean temperature, precipitation & nitrogen deposition to data.sem
data.sem1 <- merge(data.sem,data.richness, all=T)
data.sem2 <- merge(data.sem1,data.clim1, all=T)
data.sem3 <- merge(data.sem2,data.dep1, all=T)

# 12. Structural equation modeling----
#To avoid issues with SEM, you might want make explanatory vars from factors to integers*

# preparing variables
str(data.sem3)

## create species richness (with correct name)
data.sem3$species_richness <- data.sem3$species_richnnes

#data.sem$Country <- as.factor(data.sem$Country)
### I think it would be nice to have the column for countries back in the data, so we can use it as a random factor

### transforming species richness to numeric
data.sem3$data.sem3$species_richness <- as.numeric(data.sem3$data.sem3$species_richness)

#data.sem$Harvest_Area_.cm2. <- as.numeric(data.sem$Harvest_Area_.cm2.)
### I think it would be nice to have the column for the harvested area back in the data, so we can use it as a random factor

### you might want to check the data distribution with histograms
# Source functions ----
source("R/zzz_functions.R") # functions to automate some plots

#Check the distribution of multiple continuous variables------
normality_check(data.sem3, variables=c("biomass","N","C","P","K","NP","NK","KP","species_richness","shannon_diversity",
                                       "shannon_evenness","worldclim_bio_1","worldclim_bio_12","mean_5yr","mean_10yr"))
####
##

### creating dataset to run the SEM models
data.sem.model <- select(data.sem3, plot_ID, biomass, N, C, P, K, NP, NK, KP, species_richness, shannon_diversity,
                         shannon_evenness, worldclim_bio_1, worldclim_bio_12, mean_5yr, mean_10yr)
as.data.frame(data.sem.model)


#Fit piecewise model with random effects (using package nlme)

#To avoid issues in piecewiseSEM, create column with transformed response variable or transform variable directly in dataframe
#example:
#data.no.contr$biomass.exoT <- (-1/(data.no.contr$biomass.exo+1))
###
##

#
###
# 6.1 - piecewiseSEM having temperature, precipitation, N-deposition, N, P, K & N:P, K:P and N:K ratios as predictors----

###removing NA values prior to fitting the model
data.sem.model <- na.omit(data.sem.model)

### log-transforming evenness
#data.sem.model$log_evenness <- log(data.sem.model$shannon_evenness+1) ### not working in the SEM

# Create component models and store in list
SEM.list.nutrients <- list(
  
  #Predicting N values 
  nitro.sem <- lme(N ~ worldclim_bio_1 + worldclim_bio_12 + mean_5yr + mean_10yr,
                 random = ~1|plot_ID,
                 data = data.sem.model),
  
  #Predicting K values 
  potas.sem <- lme(K ~ worldclim_bio_1 + worldclim_bio_12,
                   random = ~1|plot_ID,
                   data = data.sem.model),
  
  #Predicting P values 
  phosp.sem <- lme(P ~ worldclim_bio_1 + worldclim_bio_12,
                   random = ~1|plot_ID,
                   data = data.sem.model),
  
  #Predicting N:K ratio 
  NK.sem <- lme(NK ~ N + K,
                   random = ~1|plot_ID,
                   data = data.sem.model),
  
  #Predicting N:P ratio 
  NP.sem <- lme(NP ~ N + P,
                random = ~1|plot_ID,
                data = data.sem.model),
  
  #Predicting K:P ratio 
  KP.sem <- lme(KP ~ K + P,
                random = ~1|plot_ID,
                data = data.sem.model),
  
  #Predicting biomass 
  bio.sem <- lme(log(biomass) ~ worldclim_bio_1 + worldclim_bio_12 + N + P + K + NK + NP + species_richness,
                 random = ~1|plot_ID,
                 data = data.sem.model),
  
  #Predicting species richness 
  richness.sem <- lme(log(species_richness) ~ biomass + N + NP + shannon_evenness,
                      random = ~1|plot_ID,
                      data = data.sem.model))

  #Predicting evenness 
#  evenness.sem <- lme(shannon_evenness ~ NP,
#                    random = ~1|plot_ID,
#                    data = data.sem.model))

#Run psem with new syntax
nutrients.psem <- as.psem(SEM.list.nutrients)

summ.nutrients <- summary(nutrients.psem)

coef.nutrients <- summ.nutrients$coefficients

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect

write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_v1.csv")
####

#### plotting SEM
plot(nutrients.psem, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

####
####
#### end of code
####
####


#
###
# 6.2 - piecewiseSEM having N:P, K:P and N:K ratios as predictors----

# Create component models and store in list
SEM.list.nutrients.2 <- list(
  
  #Predicting biomass 
  bio.sem.2 <- lme(biomass_g.m2 ~ NP + KP + NK,
                   random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                   na.action = na.omit, data = data.sem.model),
  
  #Predicting species richness 
  richness.sem.2 <- lme(log(spp.richness) ~ biomass_g.m2 + NP + KP + NK,
                        random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                        na.action = na.omit, data = data.sem.model))

#Run psem with new syntax
nutrients.psem.2 <- as.psem(SEM.list.nutrients.2)

summ.nutrients.2 <- summary(nutrients.psem.2)

coef.nutrients.2 <- summ.nutrients.2$coefficients

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect

write.csv(coef.nutrients.2, "Outputs/Tables/results_piecewiseSEM_GlobNut-v77.csv")
####

#### plotting SEM
plot(nutrients.psem.2, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
###
# 6.3 - piecewiseSEM having the percentage of N, P and K as predictors----

# Create component models and store in list
SEM.list.nutrients.3 <- list(
  
  #Predicting biomass 
  bio.sem.3 <- lme(biomass_g.m2 ~ + N_percent + P_percent + K_percent,
                   random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                   na.action = na.omit, data = data.sem.model),
  
  #Predicting species richness 
  richness.sem.3 <- lme(log(spp.richness) ~ biomass_g.m2 + N_percent + P_percent + K_percent,
                        random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                        na.action = na.omit, data = data.sem.model))

#Run psem with new syntax
nutrients.psem.3 <- as.psem(SEM.list.nutrients.3)

summ.nutrients.3 <- summary(nutrients.psem.3)

coef.nutrients.3 <- summ.nutrients.3$coefficients

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect

write.csv(coef.nutrients.3, "Outputs/Tables/results_piecewiseSEM_GlobNut-v4.csv")
####

#### plotting SEM
plot(nutrients.psem.3, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
###
# 6.4 - piecewiseSEM having N, P, K (g.m2) as predictors of N:P, K:P and N:K ratios----

# Create component models and store in list
SEM.list.nutrients.4 <- list(
  
  #Predicting NP ratio 
  bio.np.4 <- lme(NP ~ N_g.m2 + P_g.m2, random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                  na.action = na.omit, data = data.sem.model),
  
  #Predicting KP ratio 
  bio.kp.4 <- lme(KP ~ K_g.m2 + P_g.m2, random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                  na.action = na.omit, data = data.sem.model),
  
  #Predicting KP ratio 
  bio.nk.4 <- lme(NK ~ N_g.m2 + K_g.m2, random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                  na.action = na.omit, data = data.sem.model),
  
  #Predicting biomass 
  bio.sem.4 <- lme(biomass_g.m2 ~ N_g.m2 + P_g.m2 + K_g.m2 + NP + KP + NK,
                   random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                   na.action = na.omit, data = data.sem.model),
  
  #Predicting species richness 
  richness.sem.4 <- lme(log(spp.richness) ~ biomass_g.m2 + N_g.m2 + P_g.m2 + K_g.m2 + NP + KP + NK,
                        random = list(~1|Country, ~1|Harvest_Area_.cm2.),
                        na.action = na.omit, data = data.sem.model))

#Run psem with new syntax
nutrients.psem.4 <- as.psem(SEM.list.nutrients.4)

summ.nutrients.4 <- summary(nutrients.psem.4)

coef.nutrients.4 <- summ.nutrients.4$coefficients

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect

write.csv(coef.nutrients.4, "Outputs/Tables/results_piecewiseSEM_GlobNut-v77.csv")
####

#### plotting SEM
plot(nutrients.psem.4, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##