## ---------------------------
##
## Script name: 03_SEM
##
## Purpose of script: Calculating piecewise structural equation models
##
## Author: Leonardo H. Teixeira & Annegreet veeken
##
## Date Created: 27-10-2025
##
## Email: leonardo.htp@gmail.com, veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes: Structural equation models are calculated using the R-package piecewiseSEM
##  
## References: 
## Lefcheck, J.S. (2016). PiecewiseSEM: Piecewise structural equation modeling in R for ecology, evolution, and systematics. Methods in Ecology and Evolution 7: 573–579. DOI:10.1111/2041-210X.12512
## https://jslefche.github.io/sem_book/
## Supplementary of doi/10.1126/science.1256330
##
## ---------------------------

# Loading Packages
library(tidyverse)
library(car)
library(piecewiseSEM)
library(easystats)
library(lme4)
library(corrplot)
library(semEff)

## opening & filtering data----
data.nut <- readRDS("outputs/02_GlobNut.rds")
str(data.nut)

# subsetting dataset - only plots with soil pH measurements (pH_field)
data.sem <- subset(data.nut, pH_field!= "NA")

## Preparing data for SEM
data.sem$habitat <- as.factor(data.sem$habitat) # transforming habitat type to factor

## Conceptual model
dagitty::dagitty("dag{
  ndep -> NP
  ndep -> biomass
  ndep -> pH_field
  NP -> spec_ric
  NP -> biomass
  plot_size -> spec_ric
  spec_pool -> spec_ric
  biomass -> spec_ric
  MAP -> NP
  MAP -> biomass
  MAP -> pH_field
  MAT -> NP
  MAT -> biomass
  MAT -> pH_field
  PET -> NP
  PET -> biomass
  PET -> pH_field
  pH_field -> NP
  pH_field -> spec_ric
}") %>% plot


# variable selection and transformations
vars_to_center <- c("ndep", "MAT", "MAP", "PET","spec_ric", "biomass", "NP", "pH_field", "spec_pool", "plot_size")
vars_to_log <- c("ndep", "biomass", "MAP", "NP", "spec_ric", "q1", "q2", "plot_size")
vars_to_quad <- c("ndep", "MAT", "MAP", "biomass", "NP", "pH_field")

#  creating dataset to run the SEM models 
data.sem.model <- data.sem %>%
  dplyr::select(plot_ID, country, lon, lat, dg_cell, plot_size, sample_year, biomass,  NP, spec_ric,
                q1, q2, MAT, MAP, PET, ndep, pH_field, habitat, spec_pool) %>% 
  # removing NA values prior to fitting the model
  drop_na() %>% 
  # log transform to approach normality and improve model fit
  dplyr::mutate(across(all_of(vars_to_log), log)) %>% 
  # center to reduce co-linearity between linear and quadratic terms
  dplyr::mutate(across(all_of(vars_to_center),
           ~ as.numeric(scale(.x, center = TRUE, scale = FALSE)))   # ensure vector, not matrix
  ) %>%
  # creating quadratic variables for the the SEM this needs to be done prior to running the SEM, otherwise variables won't be recognized
  dplyr::mutate(across(all_of(vars_to_quad),
                       ~ .^2, .names = "{.col}_quad"   # ensure vector, not matrix
  ))

#  Checking correlation among variables
cor_mat <- cor(data.sem.model[,c("biomass","NP","spec_ric","MAT","MAP","PET","ndep",
                                      "pH_field","plot_size", "spec_pool")], use="pairwise.complete.obs") # Correlation matrix
colnames(cor_mat) <- c("Biomass", "N/P", "Species richness", "MAT", "MAP", "PET", "N deposition", "pH", "Plot size", "Species pool") # nicer labels
rownames(cor_mat) <- c("Biomass", "N/P", "Species richness", "MAT", "MAP", "PET", "N deposition", "pH", "Plot size", "Species pool") # nicer labels

# Plot & save correlation matrix
png("figures/202511-Figures-resubmission/Supp-correlation-matrix.png", res = 300,
    width = 7, height = 7, units = "in")
corrplot(cor_mat, 
         method="color", 
         type = "upper",
         addCoef.col = "black",  
         tl.col = "gray25",       
         tl.cex = 0.9,
         number.cex = 0.8,
         mar = c(0,0,0,8))
dev.off()

## Sub models for piecewiseSEM for species richness -----
#Fit piecewise model with random effects (using package lme4)
## Predicting N:P ratio 
NP.lmm.rich <- lmer(NP ~ MAT +  MAP + PET + ndep + pH_field + I(pH_field^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                        data = data.sem.model)
#model checks
check_model(NP.lmm.rich)
model_parameters(NP.lmm.rich)

# calculate composites
# Temperature hypothesized to be quadratic but linear & quadratic is not significant -> don't use composite
np.ph.lm <- lmer(NP ~ pH_field +  I(pH_field^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) , data = data.sem.model)
check_model(np.ph.lm)
data.sem.model$comp_np_ph <- summary(np.ph.lm)$coefficients[2,1] * data.sem.model$pH_field + summary(np.ph.lm)$coefficients[3,1] * data.sem.model$pH_field_quad

## Predicting soil pH 
ph.lmm.rich <- lmer(pH_field ~  ndep + MAT + MAP +  PET + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                   data = data.sem.model)
# model checks
check_model(ph.lmm.rich)
model_parameters(ph.lmm.rich)

## Predicting biomass 
bio.lmm.rich <- lmer(biomass ~ MAT + I(MAT^2) + MAP + PET + ndep + NP + pH_field +  I(pH_field^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                    data = data.sem.model)
# model checks
check_model(bio.lmm.rich)
model_parameters(bio.lmm.rich)
# Calculate composite
bio.ph.lm <- lmer(biomass ~ pH_field +  I(pH_field^2) + (1|dg_cell) + (1|habitat) + (1|sample_year), data = data.sem.model)
check_model(bio.ph.lm)
data.sem.model$comp_bio_ph <- summary(bio.ph.lm)$coefficients[2,1] * data.sem.model$pH_field + summary(bio.ph.lm)$coefficients[3,1] * data.sem.model$pH_field_quad
bio.mat.lm <- lmer(biomass ~ MAT +  I(MAT^2) + (1|dg_cell) + (1|habitat) + (1|sample_year), data = data.sem.model)
check_model(bio.mat.lm)
data.sem.model$comp_bio_mat <- summary(bio.mat.lm)$coefficients[2,1] * data.sem.model$MAT + summary(bio.mat.lm)$coefficients[3,1] * data.sem.model$MAT_quad

#Predicting species richness 
richness.lmm.rich <- lmer(spec_ric ~ MAT + MAP + PET +
                           biomass + I(biomass^2) + NP + I(NP^2) + 
                           ndep + I(ndep^2) + pH_field + I(pH_field^2) + plot_size + spec_pool + (1|dg_cell) + (1|habitat) + (1|sample_year),
                         data = data.sem.model)
# model checks
check_model(richness.lmm.rich)
model_parameters(richness.lmm.rich)
# calculate composites
rich.bio.lm <- lmer(spec_ric ~ biomass +  I(biomass^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) , data = data.sem.model)
check_model(rich.bio.lm)
data.sem.model$comp_rich_bio <- summary(rich.bio.lm)$coefficients[2,1] * data.sem.model$biomass + summary(rich.bio.lm)$coefficients[3,1] * data.sem.model$biomass_quad
rich.np.lm <- lmer(spec_ric ~ NP +  I(NP^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) , data = data.sem.model)
check_model(rich.bio.lm)
data.sem.model$comp_rich_np <- summary(rich.np.lm)$coefficients[2,1] * data.sem.model$NP + summary(rich.np.lm)$coefficients[3,1] * data.sem.model$NP_quad
rich.ph.lm <- lmer(spec_ric ~ pH_field +  I(pH_field^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) , data = data.sem.model)
check_model(rich.ph.lm)
data.sem.model$comp_rich_ph <- summary(rich.ph.lm)$coefficients[2,1] * data.sem.model$pH_field + summary(rich.ph.lm)$coefficients[3,1] * data.sem.model$pH_field_quad
rich.ndep.lm <- lmer(spec_ric ~ ndep +  I(ndep^2) + (1|dg_cell) + (1|habitat) + (1|sample_year) , data = data.sem.model)
check_model(rich.ndep.lm)
data.sem.model$comp_rich_ndep <- summary(rich.ndep.lm)$coefficients[2,1] * data.sem.model$ndep + summary(rich.ndep.lm)$coefficients[3,1] * data.sem.model$ndep_quad


# SEM Without composites to understand basis-set
NP.sem.rich <- lmer(NP ~ MAT + MAP + PET + ndep + pH_field + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                     data = data.sem.model)
ph.sem.rich <- lmer(pH_field ~  ndep + MAP + MAT + PET + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                     data = data.sem.model)
bio.sem.rich <- lmer(biomass ~ MAT  + MAP + PET + ndep + NP  + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                      data = data.sem.model)
richness.sem.rich <- lmer(spec_ric ~
                             biomass + NP + spec_pool +
                             ndep + pH_field + plot_size + (1|dg_cell) + (1|habitat) + (1|sample_year),
                           data = data.sem.model)
nutrients.psem.rich <- psem(NP.sem.rich, ph.sem.rich, bio.sem.rich, richness.sem.rich)
basis_set <-  basisSet(nutrients.psem.rich)

## SEM specifications -----
# SEM 1 - based on conceptual model
NP.sem.rich1 <- lmer(NP ~ MAT + MAP + PET + ndep + comp_np_ph + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                    data = data.sem.model)
ph.sem.rich1 <- lmer(pH_field ~  ndep + MAP + MAT + PET + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                    data = data.sem.model)
bio.sem.rich1 <- lmer(biomass ~ comp_bio_mat  + MAP + PET + ndep + NP  + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                     data = data.sem.model)
richness.sem.rich1 <- lmer(spec_ric ~
                            comp_rich_bio + comp_rich_np + spec_pool +
                            comp_rich_ndep + comp_rich_ph + plot_size + (1|dg_cell) + (1|habitat) + (1|sample_year),
                          data = data.sem.model)
nutrients.psem.rich1 <- psem(NP.sem.rich1, ph.sem.rich1, bio.sem.rich1, richness.sem.rich1)
basis_set1 <- basisSet(nutrients.psem.rich1)
indep_test1 <- basis_set1[c(2:4, # Climate & species richness
                            18:20, # species pool & endogenous variables
                            27:29, # plot size & endogenous variables 
                            32)] # pH - biomass
sem1 <- summary(nutrients.psem.rich1, basis.set = indep_test1)
aic1 <- AIC_psem(nutrients.psem.rich1, basis.set = indep_test1)

# SEM 2 - based on conceptual model + climate and species richness relationship
NP.sem.rich2 <- lmer(NP ~ MAT + MAP + PET + ndep + comp_np_ph + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                     data = data.sem.model)
ph.sem.rich2 <- lmer(pH_field ~  ndep + MAP + MAT + PET + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                     data = data.sem.model)
bio.sem.rich2 <- lmer(biomass ~ comp_bio_mat + MAP + PET + ndep + NP  + (1|dg_cell) + (1|habitat) + (1|sample_year) ,
                      data = data.sem.model)
richness.sem.rich2 <- lmer(spec_ric ~ MAT + MAP + PET +
                             comp_rich_bio + comp_rich_np + spec_pool +
                             comp_rich_ndep + comp_rich_ph + plot_size + (1|dg_cell) + (1|habitat) + (1|sample_year),
                           data = data.sem.model)
nutrients.psem.rich2 <- psem(NP.sem.rich2, ph.sem.rich2, bio.sem.rich2, richness.sem.rich2)
basis_set2 <- basisSet(nutrients.psem.rich2)
indep_test2 <- basis_set2[c(15:17, # species pool & endogenous variables
                            24:26, # plot size & endogenous variables 
                            29)] # pH - biomass
sem2 <- summary(nutrients.psem.rich2, basis.set = indep_test2)
aic2 <- AIC_psem(nutrients.psem.rich2, basis.set = indep_test2)

# Table for supplementary
sem2$coefficients %>%
  rename(`SE` = Std.Error,`Critical value` = Crit.Value, `P-value` = P.Value, `Standarized estimate` = Std.Estimate, `-` = last_col()) %>% 
  mutate(DF = round(DF,1), 
         `Critical value` = round(`Critical value`, 2),
         Response = case_match(Response,
                            "NP" ~ "N/P", 
                           "spec_ric" ~ "Species richness",
                           "ndep" ~ "N deposition", 
                           "pH_field" ~ "pH", 
                           "plot_size" ~ "Plot size", 
                           "spec_pool" ~ "Species pool",
                           "biomass" ~ "Biomass",
                           .default = Response),
         Predictor = case_match(Predictor,
                               "NP" ~ "N/P", 
                               "spec_ric" ~ "Species richness",
                               "ndep" ~ "N deposition", 
                               "pH_field" ~ "pH", 
                               "plot_size" ~ "Plot size", 
                               "spec_pool" ~ "Species pool",
                               "biomass" ~ "Biomass",
                               "comp_np_ph" ~ "pH (composite)",
                               "comp_bio_mat" ~ "MAT (composite)",
                               "comp_rich_bio" ~ "Biomass (composite)",
                               "comp_rich_ph" ~ "pH (composite)",
                               "comp_rich_np" ~ "N/P (composite)",
                               "comp_rich_ndep" ~ "N deposition (composite)",
                               .default = Predictor)) %>% 
  knitr::kable()

# does adding species richness - climate improve model?
aic1$AIC - aic2$AIC
anova(nutrients.psem.rich1, nutrients.psem.rich2)

# For figure (made in illustrator)
# line size for figure
linesize <- sem1$coefficients %>% 
  select(-last_col()) %>% 
  mutate(linesize = case_when(P.Value > 0.05 ~ NA,
                              between(abs(Std.Estimate), 0, 0.1) ~ 1,
                              between(abs(Std.Estimate), 0.1, 0.2) ~ 2,
                              between(abs(Std.Estimate), 0.2, 0.3) ~ 3,
                              between(abs(Std.Estimate), 0.3, 0.4) ~ 4,
                              between(abs(Std.Estimate), 0.4, 0.5) ~ 5,
  ))

linesize <- sem2$coefficients %>% 
  select(-last_col()) %>% 
  mutate(linesize = case_when(P.Value > 0.05 ~ NA,
                              between(abs(Std.Estimate), 0, 0.1) ~ 1,
                              between(abs(Std.Estimate), 0.1, 0.2) ~ 2,
                              between(abs(Std.Estimate), 0.2, 0.3) ~ 3,
                              between(abs(Std.Estimate), 0.3, 0.4) ~ 4,
                              between(abs(Std.Estimate), 0.4, 0.5) ~ 5,
                              between(abs(Std.Estimate), 0.5, 0.6) ~ 6,
  ))

