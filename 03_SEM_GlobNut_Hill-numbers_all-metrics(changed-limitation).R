## ---------------------------
##
## Script name: 04_SEM_GlobNut_1.0_v11
##
## Purpose of script: Calculating piecewise structural equation models
##
## Author: Leonardo H. Teixeira
##
## Date Created: 24-07-2024
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
data.nut <- readRDS("Data/Raw/02_GlobNut_subsampled.rds")
#data.nut <- read_csv("Data/Raw/globnut_SEM_bio.csv", col_names = T, na = "NA")
str(data.nut)
boxplot(data.nut$spec_ric ~ data.nut$lith_simp)


# 1.2. creating data.sem (lithology as integer)----
#write.csv(data.nut, "Data/Raw/data_nut_lithology_integer.csv")
#data.sem <- read_csv("Data/Raw/data_nut_lithology_integer.csv", col_names = T, na = "NA")
data.sem <- data.nut
data.sem$lith_integer <- data.sem$lith_simp
data.sem$lith_integer <- as.integer(as.factor(data.sem$lith_integer))
boxplot(data.sem$spec_ric ~ data.sem$lith_integer)

#
### 2. Checking data normality----
# checking data distribution with histograms
# Source functions to automate plots:
source("R/zzz_functions.R")

#Check the distribution of multiple continuous variables------
normality_check(data.sem, variables=c("biomass","N","P","K","NP","spec_ric",
                                      "pilou_eve","MAT","MAP","PET","ndep",
                                      "soil_age"))

#
### 2.1. Checking collinearity among variables----

### response variables
cor_num <- cor(data.sem[,c(9:11,13,15:20,22,24:25)],use="pairwise.complete.obs") # Correlation matrix
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
pairs(data.sem[,c(9:11,13,15:20,22,24:25)], lower.panel = NULL, use="pairwise.complete.obs")

# Option 3
library(psych)
pairs.panels(data.sem[,c(9:11,13,15:20,22,24:25)],
             density = FALSE, 
             ellipses = FALSE)

#
### 3. Structural equation modeling----
# SEM with soil age and rock acidity (lithology) as integers

# preparing variables
str(data.sem)

## create species richness (with correct name)
data.sem$species_richness <- data.sem$spec_ric

## create shannon diversity (with correct name)
data.sem$shannon_diversity <- data.sem$q1

## create simpson diversity (with correct name)
data.sem$simpson_diversity <- data.sem$q2

## create pielou evenness (with correct name)
data.sem$pielou_evenness <- data.sem$pilou_eve

## create N:K variable
#data.sem$NK <- data.sem$N/data.sem$K

## create K:P variable
#data.sem$KP <- data.sem$K/data.sem$P

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

### 3.1. creating dataset to run the SEM models (soil age as integer)----
data.sem.model <- select(data.sem, plot_ID, country, cell, plot_size, lat, lon, sample_year, biomass, N, P, NP, species_richness, pielou_evenness, shannon_diversity, simpson_diversity, MAT, MAP, PET, ndep, soil_age, lith_simp, lith_integer, lim)

str(data.sem.model)
data.sem.model <- as.data.frame(data.sem.model)
#data.sem.model$lith_merge <- as.factor(data.sem.model$lith_merge)
str(data.sem.model)

### 3.2. creating quadratic variables to use in the SEM----
# this needs to be done prior to running the SEM, otherwise variables won't be recognized
data.sem.model$MAT_quad <- data.sem.model$MAT^2
data.sem.model$MAP_quad <- data.sem.model$MAP^2
data.sem.model$species_richness_quad <- data.sem.model$species_richness^2
data.sem.model$biomass_quad <- data.sem.model$biomass^2
data.sem.model$NP_quad <- data.sem.model$NP^2

### 3.3. creating inverted ratios to use in the SEM----
# this is to check if partial correlations will remain the same with reversed nutrients ratio (i.e. from NP to PN)
#data.sem.model$PN <- (data.sem.model$P)/(data.sem.model$N)
#data.sem.model$KN <- (data.sem.model$K)/(data.sem.model$N)
#data.sem.model$PK <- (data.sem.model$P)/(data.sem.model$K)

### 3.4. creating proxy for nutrient availability via multiplying by biomass-----
data.sem.model$Nava <- (data.sem.model$N)*(data.sem.model$biomass)
data.sem.model$Pava <- (data.sem.model$P)*(data.sem.model$biomass)
#data.sem.model$Kava <- (data.sem.model$K)*(data.sem.model$biomass)

### 3.5. removing NA values prior to fitting the model----
data.sem.model <- na.omit(data.sem.model)

#
### 4.1. - creating dataset to run the SEM models
# (excluding plots with No limitation by N, P, K)----
#data.sem.model <- subset(data.sem.model, lim!="No limitation by N, P, K")
### this removed 143 plots!!

#hist(data.sem.model$Pava, breaks=100)
#quantile(data.sem.model$Pava)

#To avoid issues in piecewiseSEM, create column with transformed response variable or transform variable directly in dataframe
#example:
#data.no.contr$biomass.exoT <- (-1/(data.no.contr$biomass.exo+1))
###
##

#
### 4.2a. Checking normality of data.sem.model----
# checking data distribution with histograms
# Source functions to automate plots:
source("R/zzz_functions.R")

#Check the distribution of multiple continuous variables------
normality_check(data.sem.model, variables=c("biomass","NP","species_richness",
                                      "shannon_diversity","simpson_diversity",
                                      "MAT","MAP","ndep","soil_age",
                                      "lith_integer"))

#
### 4.2b. Checking collinearity among variables----

### response variables
cor_num_final <- cor(data.sem.model[,c(8:12,14:20,22,24)],use="pairwise.complete.obs") # Correlation matrix
cor_num_final

# Plot correlation matrix
# Option 1
library(corrplot)
corrplot(cor_num_final, 
         method="color", 
         type = "upper",
         addCoef.col = "black",  # Color of the correlation coefficients
         tl.col = "gray25",        # Color of the text labels
         tl.cex = 0.9,
         number.cex = 0.8,
         mar = c(0,0,0,8))
###
##

#
### 5.1. - piecewiseSEM for species diversity metrics (Harry's suggestion) -----

#Fit piecewise model with random effects (using package nlme)

# MAT, MAP, PET, N-deposition, soil age, rock acidity, N-avail, P-avail, N:P as predictors

# Create component models and store in list
SEM.list.nutrients.div <- list(
  
  #Predicting N:P ratio 
  NP.sem.div <- lme(NP ~ MAT_quad + MAP + PET + log(ndep) + soil_age + 
                       lith_integer,
                     random = ~1|cell,
                     data = data.sem.model),
  
  #Predicting biomass 
  bio.sem.div <- lme(log(biomass) ~ MAT + MAP + PET + log(ndep) + soil_age + 
                        lith_integer + log(NP) + species_richness + 
                        shannon_diversity + simpson_diversity,
                      random = ~1|cell,
                      data = data.sem.model),
  
  #Predicting species richness (log-transformed to get a normal distribution)
  richness.sem.div <- lme(log(species_richness) ~ log(biomass) + log(NP),
                           random = ~1|cell,
                           data = data.sem.model),

  #Predicting shannon diversity (log-transformed to get a normal distribution)
  shannon.sem.div <- lme(log(shannon_diversity) ~ log(biomass) + log(NP),
                          random = ~1|cell,
                          data = data.sem.model),

  #Predicting simpson diversity (log-transformed to get a normal distribution)
  simpson.sem.div <- lme(log(simpson_diversity) ~ log(biomass) + log(NP),
                        random = ~1|cell,
                        data = data.sem.model))

#Run psem with new syntax
nutrients.psem.div <- as.psem(SEM.list.nutrients.div)

summ.nutrients.div <- summary(nutrients.psem.div)

coef.nutrients.div <- summ.nutrients.div$coefficients

AIC(nutrients.psem.div, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.div)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.div, "Outputs/Tables/results_piecewiseSEM_GlobNut_hill-numbers_spp-diversity-metrics_v11_July24.csv")
####

#### plotting SEM
plot(nutrients.psem.div, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
### 6.1. - piecewiseSEM for species richness -----

#Fit piecewise model with random effects (using package nlme)

# MAT, MAP, PET, N-deposition, soil age, rock acidity, N-avail, P-avail, N:P as predictors

# Create component models and store in list
SEM.list.nutrients.rich <- list(
  
  #Predicting N:P ratio 
  NP.sem.rich <- lme(NP ~ MAT_quad + MAP + PET + log(ndep) + soil_age + 
                       lith_integer,
                  random = ~1|cell,
                data = data.sem.model),
  
  #Predicting biomass 
  bio.sem.rich <- lme(log(biomass) ~ MAT + MAP + PET + log(ndep) + soil_age + 
                        lith_integer + log(NP) + species_richness, #+
                        #shannon_diversity + simpson_diversity,
                  random = ~1|cell,
                 data = data.sem.model),
  
  #Predicting species richness
  richness.sem.rich <- lme(log(species_richness) ~ log(biomass) + log(NP),
                        random = ~1|cell,
                      data = data.sem.model))

#Run psem with new syntax
nutrients.psem.rich <- as.psem(SEM.list.nutrients.rich)

summ.nutrients.rich <- summary(nutrients.psem.rich)

coef.nutrients.rich <- summ.nutrients.rich$coefficients

AIC(nutrients.psem.rich, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.rich)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.rich, "Outputs/Tables/results_piecewiseSEM_GlobNut_hill-numbers_richness_v11_July24.csv")
####

#### plotting SEM
plot(nutrients.psem.rich, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
### 7.1. - piecewiseSEM for shannon diversity -----

#Fit piecewise model with random effects (using package nlme)

# MAT, MAP, PET, N-deposition, soil age, rock acidity, N-avail, P-avail, N:P as predictors

# Create component models and store in list
SEM.list.nutrients.shan <- list(
  
  #Predicting N:P ratio 
  NP.sem.shan <- lme(NP ~ MAT_quad + MAP + PET + log(ndep) + soil_age + 
                       lith_integer,
                     random = ~1|cell,
                     data = data.sem.model),
  
  #Predicting biomass 
  bio.sem.shan <- lme(log(biomass) ~ MAT + MAP + PET + log(ndep) + soil_age + 
                        lith_integer + log(NP) + shannon_diversity,
                      random = ~1|cell,
                      data = data.sem.model),
  
  #Predicting shannon diversity
  shannon.sem.shan <- lme(log(shannon_diversity) ~ log(biomass) + log(NP),
                          random = ~1|cell,
                          data = data.sem.model))

#Run psem with new syntax
nutrients.psem.shan <- as.psem(SEM.list.nutrients.shan)

summ.nutrients.shan <- summary(nutrients.psem.shan)

coef.nutrients.shan <- summ.nutrients.shan$coefficients

AIC(nutrients.psem.shan, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.shan)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.shan, "Outputs/Tables/results_piecewiseSEM_GlobNut_hill-numbers_shannon_v11_July24.csv")
####

#### plotting SEM
plot(nutrients.psem.shan, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

#
### 8.1. - piecewiseSEM for simpson diversity -----

#Fit piecewise model with random effects (using package nlme)

# MAT, MAP, PET, N-deposition, soil age, rock acidity, N-avail, P-avail, N:P as predictors

# Create component models and store in list
SEM.list.nutrients.simp <- list(
  
  #Predicting N:P ratio 
  NP.sem.simp <- lme(NP ~ MAT_quad + MAP + PET + log(ndep) + soil_age + 
                       lith_integer,
                     random = ~1|cell,
                     data = data.sem.model),
  
  #Predicting biomass 
  bio.sem.simp <- lme(log(biomass) ~ MAT + MAP + PET + log(ndep) + soil_age + 
                        lith_integer + log(NP) + simpson_diversity,
                      random = ~1|cell,
                      data = data.sem.model),

  #Predicting simpson diversity
  simpson.sem.simp <- lme(log(simpson_diversity) ~ log(biomass) + log(NP),
                          random = ~1|cell,
                          data = data.sem.model))

#Run psem with new syntax
nutrients.psem.simp <- as.psem(SEM.list.nutrients.simp)

summ.nutrients.simp <- summary(nutrients.psem.simp)

coef.nutrients.simp <- summ.nutrients.simp$coefficients

AIC(nutrients.psem.simp, AIC.type = "dsep", aicc = T)

LLchisq(nutrients.psem.simp)

# write.csv(coef.nutrients, "Outputs/Tables/results_piecewiseSEM_GlobNut.csv") ## results having only country as random effect
write.csv(coef.nutrients.simp, "Outputs/Tables/results_piecewiseSEM_GlobNut_hill-numbers_simpson_v11_July24.csv")
####

#### plotting SEM
plot(nutrients.psem.simp, node_attrs = list(shape = "rectangle", color = "black", fillcolor = "grey"))
####
##

####
####
#### END OF CODE
####
####
