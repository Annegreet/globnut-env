## ---------------------------
##
## Script name: 01_Plot_selection.R
##
## Purpose of script: Plot selection for further analysis 
##
## Author: Annegreet Veeken
##
## Date Created: 2023-11-01
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##  
## References:
##
## ---------------------------

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")

## Load data
data_dir <- "Z:/_GlobNut1.0/"
# Predictors
ndep <- read.csv(paste0(data_dir, "Ndeposition.csv"))
clim <- read.csv(paste0(data_dir, "ERA5_climate.csv"))
age <- read.csv(paste0(data_dir, "soilage.csv")) %>% 
  dplyr::select(-lat, -lon)
bed <-  read.csv(paste0(data_dir, "lithology.csv")) %>% 
  dplyr::select(plot_ID, lith = value_chr)
# response
npk <- read.csv(paste0(data_dir, "GlobNut1.0_nutrients.csv")) %>% 
  # add column with nutrient limitation
  mutate(lim = case_when(NP < 10 & NK < 2.1 & N < 2.0 ~ "N-limitation", 
                         NP > 16 & KP > 3.4 & P < 0.11 ~ "P-limitation", 
                         NK >= 2.1 & KP < 3.4 & K < 0.8 ~  "K(co)-limitation", 
                         NP >= 10 & NP < 16 & KP >= 3.4 ~ "Co-limitation N-P",
                         P >= 0.11 & N >= 2.0 & K >= 0.8 ~ "No limitation by N, P, K",
                         TRUE ~ "Unclear"),
         lim_ratio = case_when(NP < 10 & NK < 2.1 ~ "N-limitation", 
                               NP > 16 & KP > 3.4 ~ "P-limitation", 
                               NK >= 2.1 & KP < 3.4 ~  "K(co)-limitation", 
                               NP >= 10 & NP < 16 & KP > 3.4 ~ "Co-limitation N-P",
                               TRUE ~ "Unclear")) 
meta <- read.csv(paste0(data_dir, "GlobNut1.0_metadata.csv"))
grid <- readRDS("outputs/01_Globnut_grid_res15.rds") %>% 
  rename(plot_ID = globnut.plot_ID) %>% 
  dplyr::select(-lat, -lon)
spec <- readRDS("outputs/01_Species_indices.rds")

# dataframe for analysis
globnut <- grid %>%
  # Join datasets
  left_join(spec, by = c("plot_ID")) %>% 
  left_join(npk, by = "plot_ID") %>% 
  left_join(clim, by = "plot_ID") %>% 
  left_join(ndep, by = "plot_ID") %>% 
  left_join(age, by = "plot_ID") %>% 
  left_join(meta, by = "plot_ID") %>% 
  left_join(bed, by = "plot_ID") %>% 
  # filter plots that have been fertilized
  filter(!harm_fert_appl %in% c(1,2))  %>% 
  # select relevant columns
  dplyr::select(plot_ID, country, cell, plot_size, sample_year = year, lat, lon, spec_ric, grass_cover, 
                pilou_eve, soil_age, lith, ndep = sum_5yr, MAT, MAP, N, P, K, NP, lim, biomass) %>% 
  # keep only plots with complete data
  drop_na(-plot_size, -sample_year) %>% 
  # filter K-limited and unclear limitation plots
  filter(!lim %in% c("Unclear",  "K(co)-limitation")) %>% 
  # filter out outliers with high biomass (potential filter)
  mutate(z_biomass = (biomass-mean(biomass))/sd(biomass)) 
saveRDS(globnut, "outputs/01_GlobNut.rds")

# sub sampled data set to check the effect of oversample areas - 5 per cell
set.seed(123)
globnut_samp <- globnut %>% 
  group_by(cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut, "outputs/01_GlobNut_subsampled.rds")