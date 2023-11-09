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

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(vegan)) install.packages("vegan")
if(!require(readxl)) install.packages("readxl")

## Load data
data_dir <- "Z:/_GLOBNUT1.0/" # directory with Globnut 1.0 data
spec_raw <- read.csv(paste0(data_dir, "GlobNut1.0_species.csv"))

## Data wrangling
spec <- spec_raw %>% 
  # only vascular plants
  dplyr::filter(vascular_plant == 1) %>% 
  # convert subspecies to species
  mutate(species_new = str_remove(species_new, pattern = " subsp.*| var\\.")) %>%
  group_by(plot_ID, family_new, species_new) %>% # aggregate by species
  summarise(cover = sum(cover, na.rm = TRUE)) %>% 
  # remove zero abundances
  dplyr::filter(cover != 0)  %>% 
  # add genus column 
  mutate(genus_new = word(species_new, 1))

## Species richness and evenness
spec_ric <- spec %>% 
  # calculate species richness by plot 
  group_by(plot_ID) %>% 
  summarise(spec_ric = n_distinct(species_new)) 

spec_wide <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  dplyr::select(species_new, cover) %>% 
  # convert to wide to calculate diversity measures
  pivot_wider(names_from = species_new, values_from = cover, values_fill = 0) 
spec_ric$shan <- diversity(spec_wide[,-1]) # Shannon entropy
spec_ric$shan_eve <- exp(spec_ric$shan) / spec_ric$spec_ric # shannon evenness
spec_ric$pilou_eve <- spec_ric$shan / log(spec_ric$spec_ric)

## Percentage grass
grass_globnut <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  # create grouping column for grass
  mutate(grass = ifelse(family_new == "Poaceae", "Grass", "Not grass") %>% as.factor(.)) %>% 
  group_by(plot_ID, grass, .drop = FALSE) %>% 
  summarise(grass_cover = sum(cover)) %>% 
  filter(grass == "Grass") %>% 
  dplyr::select(plot_ID, grass_cover)

## Percentage forbs
forb_globnut <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  # create grouping column for forb
  mutate(forb = ifelse(family_new %in% 
                         c("Poaceae", "Juncaceae","Cyperaceae"), 
                       "Not forb", "Forb") %>% as.factor(.)) %>% 
  group_by(plot_ID, forb, .drop = FALSE) %>% 
  summarise(forb_cover = sum(cover)) %>% 
  filter(forb == "Forb") %>% 
  dplyr::select(plot_ID, forb_cover)
spec_ric <- spec_ric %>% 
  left_join(grass_globnut, by = "plot_ID") %>% 
  left_join(forb_globnut, by = "plot_ID")

saveRDS(spec_ric, "Outputs/01_Species_indices.rds")