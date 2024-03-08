## ---------------------------
##
## Script name: 04_SEM_GlobNut_1.0_v1
##
## Purpose of script: Calculating piecewise structural equation models
##
## Author: Leonardo H. Teixeira
##
## Date Created: 22-11-2023
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
library(dplyr) # handling datasets and creating tibbles (e.g., read_csv & left_join)
library(tidyverse)
library(purrr) # for easily merging multiple datasets simultaneously (using reduce)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(ggplot2)
library(piecewiseSEM)
library(nlme)
library(mgcv)
library(readr) # for the function read_csv
#library(sesem)
#library(semPlot)
#library(lavaan)
#library(lavaanPlot)

#
### 1. opening & filtering data----

# 1.1. opening data nutrients and biomass----
data.nut <- readRDS("Data/Raw/02_GlobNut.rds")
#data.nut <- read_csv("Data/Raw/globnut_SEM_bio.csv", col_names = T, na = "NA")

# replacing spaces by underline (in case needed)
#data.nut$species <- gsub(" ","_",data.nut$species)

# 1.2. opening species richness data----
data.richness <- read_csv("Data/Raw/richness.csv", col_names = T, na = "NA")
#rownames(traits) <- traits$species

### keeping only plot_id and shannon evenness in data.richness:
d.rich <- data.richness[,c(1,4)]

# 1.3. opening climate data----
#data.clim <- read_csv("Data/Raw/climate.csv", col_names = T, na = "NA")

# 1.4. opening nitrogen deposition data----
#data.n_dep <- read_csv("Data/Raw/N_deposition.csv", col_names = T, na = "NA")

# 1.5. opening soil age data----
#data.soil_age <- read_csv("Data/Raw/soil_age.csv", col_names = T, na = "NA")
####
##

#
### 2. merging datasets according to the plots present in data.nut----
data.sem <- purrr::reduce(list(data.nut, d.rich),
                          dplyr::left_join, by = 'plot_ID')

#
### 3. Checking data normality----
# checking data distribution with histograms
# Source functions to automate plots:
source("R/zzz_functions.R")

#Check the distribution of multiple continuous variables------
normality_check(data.sem, variables=c("biomass","N","C","P","K","NP","NK","KP","species_richness","shannon_diversity",
                                      "shan_evenness","pielou_evenness","MAT","MAP",
                                      "mean_5yr","mean_10yr", "soil_age"))

#
### 3.1. Checking collinearity among variables----

### response variables
cor_num <- cor(data.sem[,c(2:13,16,27,54,58,64)],use="pairwise.complete.obs") # Correlation matrix
cor_num

# Plot correlation matrix
# Option 1
library(corrplot)
corrplot(cor_num, 
         method="color", 
         type = "upper",
         addCoef.col = "black")

corrplot(cor_num, 
         method="color", 
         type = "lower",
         addCoef.col = "black")

# Option 2
pairs(data.sem[,c(2:13,16,27,54,58,64)], lower.panel = NULL, use="pairwise.complete.obs")

# Option 3
library(psych)
pairs.panels(data.sem[,c(2:13,16,27,54,58,64)],
             density = FALSE, 
             ellipses = FALSE)

#
### 4. Structural equation modeling----
#To avoid issues with SEM, you might want make explanatory vars from factors to integers*

# preparing variables
str(data.sem)

## create species richness (with correct name)
data.sem$species_richness <- data.sem$spec_ric

## create pielou evenness (with correct name)
data.sem$pielou_evenness <- data.sem$pilou_eve

## create N:K variable
data.sem$NK <- data.sem$N/data.sem$K

## create K:P variable
data.sem$KP <- data.sem$K/data.sem$P

data.sem$country <- as.factor(data.sem$country) # transforming country to factor
### I think it would be nice to have the column for countries back in the data, so we can use it as a random factor
data.sem$plot_ID <- as.integer(data.sem$plot_ID) # transforming plot ID to integer
data.sem$plot_size <- as.integer(data.sem$plot_size) # transforming plot size to integer
data.sem$lith_simp <- as.factor(data.sem$lith_simp) # transforming lithology to factor
data.sem$lim <- as.factor(data.sem$lim) # transforming limitation type to factor
data.sem$species_richness <- as.numeric(data.sem$species_richness) # transforming richness to numeric

###Obs.
# MAT = mean annual temperature
# MAP = mean annual precipitation

### 4.1. creating dataset to run the SEM models (soil age as integer)----
data.sem.model <- select(data.sem, plot_ID, country, cell, plot_size, lat, lon, sample_year, biomass, N, P, K, NP, NK, KP, species_richness, shan_evenness, pielou_evenness, MAT, MAP, PET, ndep, soil_age, lith_simp, lim)

str(data.sem.model)
data.sem.model <- as.data.frame(data.sem.model)
str(data.sem.model)

### 4.2. creating quadratic variables to use in the SEM----
# this needs to be done prior to running the SEM, otherwise variables won't be recognized
data.sem.model$MAT_quad <- data.sem.model$MAT^2
data.sem.model$MAP_quad <- data.sem.model$MAP^2
data.sem.model$species_richness_quad <- data.sem.model$species_richness^2

### 4.3. creating inverted ratios to use in the SEM----
# this is to check if partial correlations will remain the same with reversed nutrients ratio (i.e. from NP to PN)
data.sem.model$PN <- (data.sem.model$P)/(data.sem.model$N)
data.sem.model$KN <- (data.sem.model$K)/(data.sem.model$N)
data.sem.model$PK <- (data.sem.model$P)/(data.sem.model$K)

### 4.4. creating proxy for nutrient availability via multiplying by biomass-----
data.sem.model$Nava <- (data.sem.model$N)*(data.sem.model$biomass)
data.sem.model$Pava <- (data.sem.model$P)*(data.sem.model$biomass)
data.sem.model$Kava <- (data.sem.model$K)*(data.sem.model$biomass)

### 4.5. removing NA values prior to fitting the model----
data.sem.model <- na.omit(data.sem.model)

#
### 5.1. creating dataset to run the SEM models (acidic soils)----
data.sem.model.acid <- subset(data.sem.model, lith_simp=="acid")

#Fit piecewise model with random effects (using package nlme)

#To avoid issues in piecewiseSEM, create column with transformed response variable or transform variable directly in dataframe
#example:
#data.no.contr$biomass.exoT <- (-1/(data.no.contr$biomass.exo+1))
###
##

#
### 5.2. - piecewiseSEM having temperature, precipitation, N-deposition, soil age, N, P, K & N:P, K:P and N:K ratios as predictors----
#Fit piecewise model with random effects (using package nlme)

# Create component models and store in list
SEM.list.nutrients.acid <- list(
  
  #Predicting N values 
  nitro.sem.acid <- lme(Nava ~ MAT_quad + MAP_quad + PET + log(ndep) + soil_age,
                        random = ~1|cell,
                        data = data.sem.model.acid),
  
  #Predicting K values 
  potas.sem.acid <- lme(Kava ~ MAT + MAP + PET + soil_age,
                        random = ~1|cell,
                   data = data.sem.model.acid),
  
  #Predicting P values 
  phosp.sem.acid <- lme(Pava ~ MAT + MAP + PET + soil_age,
                    random = ~1|cell,
                   data = data.sem.model.acid),
  
  #Predicting N:K ratio 
  NK.sem.acid <- lme(NK ~ N + K,
                  random = ~1|cell,
                data = data.sem.model.acid),
  
  #Predicting N:P ratio 
  NP.sem.acid <- lme(NP ~ N + P,
                  random = ~1|cell,
                data = data.sem.model.acid),
  
  #Predicting K:P ratio 
  KP.sem.acid <- lme(KP ~ K + P,
                  random = ~1|cell,
                data = data.sem.model.acid),
  
  #Predicting biomass 
  bio.sem.acid <- lme(log(biomass) ~ MAT + MAP + PET + Nava + Pava + Kava + NK + NP + species_richness_quad,
                  random = ~1|cell,
                 data = data.sem.model.acid),
  
  #Predicting species richness 
  richness.sem.acid <- lme(log(species_richness) ~ biomass + Nava + NP + shan_evenness,
                        random = ~1|cell,
                      data = data.sem.model.acid),
  
  #Predicting evenness 
  evenness.sem.acid <- lme(pielou_evenness ~ NP + log(species_richness),
                        random = ~1|cell,
                      data = data.sem.model.acid))

#Run psem with new syntax
nutrients.psem.acid <- as.psem(SEM.list.nutrients.acid)

summ.nutrients.acid <- summary(nutrients.psem.acid)

coef.nutrients.acid <- summ.nutrients.acid$coefficients

AIC(nutrients.psem.acid, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.acid)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.acid, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_acid-soils-CORRECT.csv")
####

#### plotting SEM
plot(nutrients.psem.acid, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
### 6.1. creating dataset to run the SEM models (intermediate soils)----
data.sem.model.intermediate <- subset(data.sem.model, lith_simp=="intermediate")
data.sem.model.intermediate <- na.omit(data.sem.model.intermediate)

### 6.2. - piecewiseSEM having temperature, precipitation, N-deposition, soil age, N, P, K & N:P, K:P and N:K ratios as predictors----

#Fit piecewise model with random effects (using package nlme)

# Create component models and store in list
SEM.list.nutrients.intermediate <- list(
  
  #Predicting N values 
  nitro.sem.intermediate <- lme(Nava ~ MAT_quad + MAP_quad + PET + log(ndep) + soil_age,
                        random = ~1|cell,
                        data = data.sem.model.intermediate),
  
  #Predicting K values 
  potas.sem.intermediate <- lme(Kava ~ MAT + MAP + PET + soil_age,
                        random = ~1|cell,
                        data = data.sem.model.intermediate),
  
  #Predicting P values 
  phosp.sem.intermediate <- lme(Pava ~ MAT + MAP + PET + soil_age,
                        random = ~1|cell,
                        data = data.sem.model.intermediate),
  
  #Predicting N:K ratio 
  NK.sem.intermediate <- lme(NK ~ N + K,
                     random = ~1|cell,
                     data = data.sem.model.intermediate),
  
  #Predicting N:P ratio 
  NP.sem.intermediate <- lme(NP ~ N + P,
                     random = ~1|cell,
                     data = data.sem.model.intermediate),
  
  #Predicting K:P ratio 
  KP.sem.intermediate <- lme(KP ~ K + P,
                     random = ~1|cell,
                     data = data.sem.model.intermediate),
  
  #Predicting biomass 
  bio.sem.intermediate <- lme(log(biomass) ~ MAT + MAP + PET + Nava + Pava + Kava + NK + NP + species_richness_quad,
                      random = ~1|cell,
                      data = data.sem.model.intermediate),
  
  #Predicting species richness 
  richness.sem.intermediate <- lme(log(species_richness) ~ biomass + Nava + NP + shan_evenness,
                           random = ~1|cell,
                           data = data.sem.model.intermediate),
  
  #Predicting evenness 
  evenness.sem.intermediate <- lme(pielou_evenness ~ NP + log(species_richness),
                           random = ~1|cell,
                           data = data.sem.model.intermediate))

#Run psem with new syntax
nutrients.psem.intermediate <- as.psem(SEM.list.nutrients.intermediate)

summ.nutrients.intermediate <- summary(nutrients.psem.intermediate)

coef.nutrients.intermediate <- summ.nutrients.intermediate$coefficients

AIC(nutrients.psem.intermediate, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.intermediate)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.intermediate, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_intermediate-soils-CORRECT.csv")
####

#### plotting SEM
plot(nutrients.psem.intermediate, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

### testing some functions to assess predictors correlations & error terms
#cerror(log(biomass) %~~% log(species_richness), nutrients.psem.intermediate, data.sem.model.intermediate)

#with(data.sem.model.intermediate, cor.test(log(biomass), log(species_richness)))

####
##

#
### 7.1. creating dataset to run the SEM models (well-buffered soils)----
data.sem.model.buffer <- subset(data.sem.model, lith_simp=="well-buffered")

### 7.2. - piecewiseSEM having temperature, precipitation, N-deposition, soil age, N, P, K & N:P, K:P and N:K ratios as predictors----

#Fit piecewise model with random effects (using package nlme)

# Create component models and store in list
SEM.list.nutrients.buffer <- list(
  
  #Predicting N values 
  nitro.sem.buffer <- lme(Nava ~ MAT_quad + MAP_quad + PET + log(ndep) + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.buffer),
  
  #Predicting K values 
  potas.sem.buffer <- lme(Kava ~ MAT + MAP + PET + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.buffer),
  
  #Predicting P values 
  phosp.sem.buffer <- lme(Pava ~ MAT + MAP + PET + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.buffer),
  
  #Predicting N:K ratio 
  NK.sem.buffer <- lme(NK ~ N + K,
                           random = ~1|cell,
                           data = data.sem.model.buffer),
  
  #Predicting N:P ratio 
  NP.sem.buffer <- lme(NP ~ N + P,
                           random = ~1|cell,
                           data = data.sem.model.buffer),
  
  #Predicting K:P ratio 
  KP.sem.buffer <- lme(KP ~ K + P,
                           random = ~1|cell,
                           data = data.sem.model.buffer),
  
  #Predicting biomass 
  bio.sem.buffer <- lme(log(biomass) ~ MAT + MAP + PET + Nava + Pava + Kava + NK + NP + species_richness_quad,
                            random = ~1|cell,
                            data = data.sem.model.buffer),
  
  #Predicting species richness 
  richness.sem.buffer <- lme(log(species_richness) ~ biomass + Nava + NP + shan_evenness,
                                 random = ~1|cell,
                                 data = data.sem.model.buffer),
  
  #Predicting evenness 
  evenness.sem.buffer <- lme(pielou_evenness ~ NP + log(species_richness),
                                 random = ~1|cell,
                                 data = data.sem.model.buffer))

#Run psem with new syntax
nutrients.psem.buffer <- as.psem(SEM.list.nutrients.buffer)

summ.nutrients.buffer <- summary(nutrients.psem.buffer)

coef.nutrients.buffer <- summ.nutrients.buffer$coefficients

AIC(nutrients.psem.buffer, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.buffer)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.buffer, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_well-buffered-soils-CORRECT.csv")
####

#### plotting SEM
plot(nutrients.psem.buffer, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##


#
### 8.1. testing the piecewiseSEM with reversed nutrient ratios (e.g., from NP to PN) ----

###removing NA values prior to fitting the model
data.sem.model <- na.omit(data.sem.model)

### log-transforming evenness
#data.sem.model$log_evenness <- log(data.sem.model$shannon_evenness+1) ### not working in the SEM

# Create component models and store in list
SEM.list.nutrients_test <- list(
  
  #Predicting N values 
  nitro.sem_test <- lme(N ~ MAT_quad + MAP_quad + mean_10yr + soil_age,
                   random = ~1|plot_ID,
                   data = data.sem.model),
                   #control = lmeControl(maxIter = 1000)),
  
  #Predicting K values 
  potas.sem_test <- lme(K ~ MAT + MAP + soil_age,
                   random = ~1|plot_ID,
                   data = data.sem.model),
  
  #Predicting P values 
  phosp.sem_test <- lme(P ~ MAT + MAP + soil_age,
                   random = ~1|plot_ID,
                   data = data.sem.model),
  
  #Predicting N:K ratio 
  KN.sem_test <- lme(KN ~ N + K,
                random = ~1|plot_ID,
                data = data.sem.model),
  
  #Predicting N:P ratio 
  PN.sem_test <- lme(log(PN) ~ N + P,
                random = ~1|plot_ID,
                data = data.sem.model),
                #control = lmeControl(maxIter = 1000)),
  
  #Predicting K:P ratio 
  PK.sem_test <- lme(PK ~ K + P,
                random = ~1|plot_ID,
                data = data.sem.model),
  
  #Predicting biomass 
  bio.sem_test <- lme(log(biomass) ~ MAT + MAP + N + K + P + NK + NP + species_richness_quad,
                 random = ~1|plot_ID,
                 data = data.sem.model),
  
  #Predicting species richness 
  richness.sem_test <- lme(log(species_richness) ~ biomass + N + NP + shan_evenness,
                      random = ~1|plot_ID,
                      data = data.sem.model),
  
  #Predicting evenness 
  evenness.sem_test <- lme(pielou_evenness ~ NP + log(species_richness),
                      random = ~1|plot_ID,
                      data = data.sem.model))

#Run psem with new syntax
nutrients.psem_test <- as.psem(SEM.list.nutrients_test)
#fit <- lavaan::sem(SEM.list.nutrients, data = data.sem.model, cor = "spearman")

summ.nutrients_test <- summary(nutrients.psem_test)

coef.nutrients_test <- summ.nutrients_test$coefficients

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect

write.csv(coef.nutrients_test, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_test_nutrient-ratios.csv")
####

#### plotting SEM
plot(nutrients.psem_test, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

##
#### Don't run the analysis for calcareous soils
##

#
### 9.1. creating dataset to run the SEM models (calcareous soils)----
data.sem.model.calcareous <- subset(data.sem.model, lith_simp=="calcareous")
# subset(mydata, age > 30 & income > 50000)

### 9.2. - piecewiseSEM having temperature, precipitation, N-deposition, soil age, N, P, K & N:P, K:P and N:K ratios as predictors----

#Fit piecewise model with random effects (using package nlme)

# Create component models and store in list
SEM.list.nutrients.calcareous <- list(
  
  #Predicting N values 
  nitro.sem.calcareous <- lme(Nava ~ MAT_quad + MAP_quad + PET + log(ndep) + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.calcareous),
  
  #Predicting K values 
  potas.sem.calcareous <- lme(Kava ~ MAT + MAP + PET + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.calcareous),
  
  #Predicting P values 
  phosp.sem.calcareous <- lme(Pava ~ MAT + MAP + PET + soil_age,
                              random = ~1|cell,
                              data = data.sem.model.calcareous),
  
  #Predicting N:K ratio 
  NK.sem.calcareous <- lme(NK ~ N + K,
                           random = ~1|cell,
                           data = data.sem.model.calcareous),
  
  #Predicting N:P ratio 
  NP.sem.calcareous <- lme(NP ~ N + P,
                           random = ~1|cell,
                           data = data.sem.model.calcareous),
  
  #Predicting K:P ratio 
  KP.sem.calcareous <- lme(KP ~ K + P,
                           random = ~1|cell,
                           data = data.sem.model.calcareous),
  
  #Predicting biomass 
  bio.sem.calcareous <- lme(log(biomass) ~ MAT + MAP + PET + Nava + Pava + Kava + NK + NP + species_richness_quad,
                            random = ~1|cell,
                            data = data.sem.model.calcareous),
  
  #Predicting species richness 
  richness.sem.calcareous <- lme(log(species_richness) ~ biomass + Nava + NP + shan_evenness,
                                 random = ~1|cell,
                                 data = data.sem.model.calcareous),
  
  #Predicting evenness 
  evenness.sem.calcareous <- lme(pielou_evenness ~ NP + log(species_richness),
                                 random = ~1|cell,
                                 data = data.sem.model.calcareous))

#Run psem with new syntax
nutrients.psem.calcareous <- as.psem(SEM.list.nutrients.calcareous)

summ.nutrients.calcareous <- summary(nutrients.psem.calcareous)

coef.nutrients.calcareous <- summ.nutrients.calcareous$coefficients

AIC(nutrients.psem.calcareous, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.calcareous)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.calcareous, "Outputs/Tables/results_piecewiseSEM_GlobNut_1.0_calcareous-soils-ava.csv")
####

#### plotting SEM
plot(nutrients.psem.calcareous, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

####
####
#### END OF CODE
####
####