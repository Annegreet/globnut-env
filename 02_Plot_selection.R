## ---------------------------
##
## Script name: 02_Plot_selection.R
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
## Augusto, L., Achat, D. L., Jonard, M., Vidal, D. & Ringeval, B. (2017). Soil parent material—A major driver of plant nutrient limitations in terrestrial ecosystems. Global Change Biology, 23(9), 3808–3824. https://doi.org/10.1111/gcb.13691
## Palpurina, S., Chytrý, M., Hölzel, N., Tichý, L., Wagner, V., Horsák, M., Axmanová, I., Hájek, M., Hájková, P., Freitag, M., Lososová, Z., Mathar, W., Tzonev, R., Danihelka, J. & Dřevojan, P. (2019). The type of nutrient limitation affects the plant species richness–productivity relationship: Evidence from dry grasslands across Eurasia. Journal of Ecology, 107(3), 1038–1050. https://doi.org/10.1111/1365-2745.13084
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
  dplyr::select(plot_ID, lith = value_chr) %>% 
  # simplify lithology category same way as Augusto et al (2017)
  mutate(lith_simp = case_when(lith %in% c("pa", "ss", "va", "su") ~ "acid",
                               lith %in% c("mt", "vi") ~ "intermediate",
                               lith %in% c("sc", "sm","pb", "vb") ~ "well-buffered")) 
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
globnut_raw <- grid %>%
  # Join datasets
  left_join(spec, by = c("plot_ID")) %>% 
  left_join(npk, by = "plot_ID") %>% 
  left_join(clim, by = "plot_ID") %>% 
  left_join(ndep, by = "plot_ID") %>% 
  left_join(age, by = "plot_ID") %>% 
  left_join(meta, by = "plot_ID") %>% 
  left_join(bed, by = "plot_ID")  %>% 
  # globnut plots with complete data
  drop_na(lat, spec_ric, biomass, N, P)

globnut <- globnut_raw %>% 
  # filter plots that have been fertilized
  filter(!harm_fert_appl %in% c(1,2)) %>% 
  # select relevant columns
  dplyr::select(plot_ID, country, cell, plot_size, sample_year = year, lat, lon,
                spec_ric, grass_cover, forb_cover, pilou_eve, soil_age, lith_simp, 
                ndep = sum_5yr, MAT, MAP, PET, AI, N, P, K, NP, lim, biomass) %>% 
  # keep only plots with complete data
  drop_na(-plot_size, -sample_year) %>% 
  # filter K-limited and unclear limitation plots
  filter(!lim %in% c("Unclear",  "K(co)-limitation")) %>% 
  # filter out outliers with high biomass (potential filter)
  mutate(z_biomass = (biomass-mean(biomass))/sd(biomass)) %>% 
  filter(z_biomass < 4) %>% 
  dplyr::select(-z_biomass)
saveRDS(globnut, "outputs/02_GlobNut.rds")

# sub sampled data set to check the effect of oversample areas - 5 per cell
set.seed(123)
globnut_samp <- globnut %>% 
  group_by(cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut_samp, "outputs/02_GlobNut_subsampled.rds")

# check for effect of grid size
grid14 <- readRDS("outputs/01_Globnut_grid_res14.rds")
globnut_samp <- globnut %>%
  dplyr::select(-cell) %>% 
  left_join(grid14, by = c("plot_ID" = "globnut.plot_ID")) %>% 
  group_by(cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut_samp, "outputs/02_GlobNut_subsampled14.rds")

grid16 <- readRDS("outputs/01_Globnut_grid_res16.rds")
globnut_samp <- globnut %>%
  dplyr::select(-cell) %>% 
  left_join(grid16, by = c("plot_ID" = "globnut.plot_ID")) %>% 
  group_by(cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut_samp, "outputs/02_GlobNut_subsampled16.rds")

# check effect of ndeposition cummulative 10
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
  dplyr::select(plot_ID, country, cell, plot_size, sample_year = year,
                lat, lon,spec_ric, grass_cover, forb_cover, pilou_eve, soil_age, lith_simp, 
                ndep = sum_10yr, MAT, MAP, N, P, K, NP, lim, biomass) %>% 
  # keep only plots with complete data
  drop_na(-plot_size, -sample_year) %>% 
  # filter K-limited and unclear limitation plots
  filter(!lim %in% c("Unclear",  "K(co)-limitation")) %>% 
  # filter out outliers with high biomass (potential filter)
  mutate(z_biomass = (biomass-mean(biomass))/sd(biomass)) %>% 
  filter(z_biomass < 4) %>% 
  dplyr::select(-z_biomass) %>% 
  group_by(cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut, "outputs/02_GlobNut_subsampled_ndep10yr.rds")
