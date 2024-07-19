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
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(vegan)) install.packages("vegan")
if (!require(readxl)) install.packages("readxl")

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
  mutate(genus_new = word(species_new, 1)) %>% 
  filter(!species_new == "")

## Species richness and evenness
spec_ric <- spec %>% 
  # calculate species richness by plot 
  group_by(plot_ID) %>% 
  summarise(spec_ric = n_distinct(species_new)) 

spec_wide <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  dplyr::select(plot_ID,species_new, cover) %>% 
  # convert to wide to calculate diversity measures
  pivot_wider(names_from = species_new, values_from = cover, values_fill = 0) 
spec_ric$shan <- diversity(spec_wide[,-1], index = "shannon") # Shannon entropy
spec_ric$shan_eve <- exp(spec_ric$shan) / spec_ric$spec_ric # shannon evenness
spec_ric$pilou_eve <- spec_ric$shan / log(spec_ric$spec_ric)
spec_ric$q1 <- exp(spec_ric$shan)
spec_ric$q2 <- diversity(spec_wide[,-1], index = "invsimpson")

saveRDS(spec_ric, "Outputs/01_Species_indices.rds")

## Species list for supplementary ----
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

# select plot ids for species list
ids <- meta  %>% 
  # filter for Eurasian sites
  dplyr::filter(between(lat, 35, 72)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  # filter plots that have been fertilized
  filter(!harm_fert_appl %in% c(1,2)) %>% 
  # exclude unclear limitation %>% 
  left_join(npk, by = "plot_ID") %>% 
  # filter K-limited and unclear limitation plots
  filter(!lim %in% c("Unclear",  "No limitation by N, P, K")) %>% 
  # filter out outliers with high biomass (potential filter)
  mutate(z_biomass = (biomass-mean(biomass, na.rm = TRUE))/sd(biomass, na.rm = TRUE)) %>% 
  filter(z_biomass < 4) %>% 
  pull(plot_ID)

spec_list <- spec %>% 
  left_join(npk[,c("plot_ID","lim")], by = "plot_ID") %>% 
  filter(plot_ID %in% ids) %>% 
  group_by(species_new) %>% 
  mutate(n_plots = n_distinct(plot_ID),
         percentage_plots = n_plots/length(ids)*100) %>% 
  group_by(species_new, lim) %>% 
  mutate(n_lim = n_distinct(plot_ID),
         percentage_lim = n_lim/n_plots*100) %>% 
  dplyr::select(species_new, n_plots, percentage_plots, lim,  percentage_lim) %>%
  distinct() %>% 
  pivot_wider( names_from = lim, values_from = percentage_lim) %>% 
  # only species
  filter(str_detect(species_new, pattern = " ")) %>% 
  mutate(across(where(~is.numeric(.)), ~round(., digits = 2)))
write.csv(spec_list, "outputs/Species_list.csv", row.names = FALSE, na = "")


