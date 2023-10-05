## ---------------------------
##
## Script name: 01_Species_diversity_indices
##
## Purpose of script: Calculating species diversity indices
##
## Author: Annegreet Veeken
##
## Date Created: 2023-09-15
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

options(scipen = 6, digits = 4) # avoid scientific notation

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(vegan)) install.packages("vegan")

## Load data
data_dir <- "Z:/_GLOBNUT1.0/" # directory with raw Globnut 1.0 data
spec_raw <- read.csv(paste0(data_dir, "GlobNut1.0_species.csv"))

## Data wrangling
spec <- spec_raw %>% 
  # only vascular plants
  dplyr::filter(vascular_plant == 1) %>% 
  # exclude unknowns
  dplyr::filter(!str_detect(species_new, pattern = "Unknown|unknown|unkown")) %>% 
  # convert subspecies to species
  mutate(species_new = str_remove(species_new, pattern = " subsp.*| var\\.\\b.")) %>%
  group_by(plot_ID, species_new) %>% # aggregate by species
  summarise(cover = sum(cover, na.rm = TRUE)) %>% 
  # remove zero abundances
  dplyr::filter(!cover == 0|!is.na(cover))

## Species richness and evenness
spec_ric <- spec %>% 
  # calculate species richness by plot 
  group_by(plot_ID) %>% 
  summarise(spec_ric = n_distinct(species_new)) %>% 
  na.omit()

spec_wide <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  dplyr::select(-covertot) %>% 
  # convert to wide to calculate diversity measures
  pivot_wider(names_from = species_new, values_from = cover, values_fill = 0) 
spec_ric$shan <- diversity(spec_wide[,-1]) # Shannon entropy
spec_ric$shan_eve <- exp(spec_ric$shan) / spec_ric$spec_ric # shannon evenness
spec_ric$pilou_eve <- spec_ric$shan / log(spec_ric)

saveRDS(spec_ric, "Outputs/01_Species_diversity.rds")

