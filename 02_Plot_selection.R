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
## Palpurina, S., Chytrý, M., Hölzel, N., Tichý, L., Wagner, V., Horsák, M., Axmanová, I., Hájek, M., Hájková, P., Freitag, M., Lososová, Z., Mathar, W., Tzonev, R., Danihelka, J. & Dřevojan, P. (2019). The type of nutrient limitation affects the plant species richness–productivity relationship: Evidence from dry grasslands across Eurasia. Journal of Ecology, 107(3), 1038–1050. https://doi.org/10.1111/1365-2745.13084
## ---------------------------

## Load packages
library(tidyverse)
library(raster)
library(terra)

## Load data
data_dir <- "Z:/_GLOBNUT1.0/"
ndepos <- read.csv(paste0(data_dir, "ndeposition_EMEP_zhu.csv"))
clim <- read.csv(paste0(data_dir, "ERA5_climate.csv"))
elev <- readRDS("env_data/outputs/elevation.rds")
pH <- readRDS("outputs/pH-data.rds")
eunis <- read.csv(paste0(data_dir, "Globnut1.0_EUNIS.csv"))
npk <- read.csv(paste0(data_dir, "GlobNut1.0_nutrients.csv")) %>% 
  # add column with nutrient limitation
  mutate(lim = case_when(NP <= 13.5 & NK <= 2.1 ~ "N-limitation",
                         NP > 16 & KP > 3.4 ~ "P-limitation",
                         NK > 2.1 & KP <= 3.4 ~  "K(co)-limitation",
                         NP >= 13.5 & NP <= 16  ~ "Co-limitation N-P",
                         P >= 0.11 & N >= 2 & K >= 0.8 ~ "No limitation by N, P, K",
                         TRUE ~ "No limitation by N, P, K")) %>%
  mutate(lim2 = case_when(NP <= 10 & NK <= 2.1 & N < 2 ~ "N-limitation", 
                         NP > 16 & KP > 3.4 & P < 0.11 ~ "P-limitation", 
                         NK > 2.1 & KP <= 3.4 & K < 0.8 ~  "K(co)-limitation", 
                         NP >= 10 & NP <= 16 & P < 0.11 & N < 2 ~ "Co-limitation N-P",
                         P >= 0.11 & N >= 2 & K >= 0.8 ~ "No limitation by N, P, K",
                         TRUE ~ "No limitation by N, P, K"))
meta <- read.csv(paste0(data_dir, "GlobNut1.0_metadata.csv"))
grid <- readRDS("outputs/01_Globnut_grid_res15.rds") %>%
  rename(plot_ID = globnut.plot_ID) %>%
  dplyr::select(-lat, -lon, dg_cell = cell)
spec_pool <- readRDS("env_data/outputs/Cai_spec_pool.rds")
# Use EMEP raster for site effects
masras <- raster("Z:/Organized-globnut/Geo-data/EMEP/Data_EMEP_report2023/EMEP01_rv5.0_year.2022met_2021emis.nc")
masras <- extend(masras, extent(-25,150,35,75))
values(masras) <- 1:ncell(masras)
# Get the corresponding grid cells for each plot
meta$cell <- terra::extract(x = masras, y = as.matrix(meta[,c("lon","lat")]))
meta$cell_lon <- xFromCell(masras, meta$cell) 
meta$cell_lat <- yFromCell(masras, meta$cell) 

spec <- readRDS("outputs/01_Species_indices.rds")

# all available data
globnut_raw <- meta %>%
  # Join datasets
  left_join(spec, by = c("plot_ID")) %>% 
  left_join(npk, by = "plot_ID") %>% 
  left_join(pH, by = "plot_ID") %>% 
  left_join(clim, by = "plot_ID") %>% 
  left_join(ndepos, by = "plot_ID") %>% 
  left_join(elev, by = "plot_ID") %>%
  left_join(eunis, by = "plot_ID") %>% 
  left_join(grid, by = "plot_ID") %>% 
  left_join(spec_pool, by = "plot_ID") %>% 
  drop_na(lon, lat, spec_ric, q1, q2, biomass, N, P)
saveRDS(globnut_raw, "outputs/02_Globnut_raw.rds")

# selection of plots for analysis
globnut <- globnut_raw %>% 
  # filter K-limited and no limitation plots
  filter(!lim %in% c("K(co)-limitation", "No limitation by N, P, K")) %>% 
  # filter out outliers with high biomass 
  drop_na(biomass) %>% 
  mutate(z_biomass = (biomass - mean(biomass))/sd(biomass)) %>% 
  filter(z_biomass < 4) %>% 
  # filter out plots that are fertilized
  filter(harm_fert_appl == 0) %>% 
  # select relevant columns
  dplyr::select(cont_ID, plot_ID, country, cell, dg_cell, cell_lon, cell_lat, 
                plot_size, sample_year = year, 
                lat, lon, elev, spec_ric, q1, q2, biomass, spec_pool,
                ndep = sum_5yr, MAT, MAP, PET, pH, pH_field, pH_soilgrids, N, P,
                K, NP, lim, habitat = sub_class, data_source) %>% 
  # remove globnut plots with incomplete data
  drop_na(MAT, MAP, pH, ndep,  plot_size, habitat, spec_pool) 

saveRDS(globnut, "outputs/02_GlobNut.rds")

# selection based on Palpurina criteria
globnut <- globnut_raw %>% 
  # filter K-limited and no limitation plots
  filter(!lim2 %in% c("K(co)-limitation", "No limitation by N, P, K")) %>%  # has higher number of "No limitation" 
  # filter out outliers with high biomass 
  drop_na(biomass) %>% 
  mutate(z_biomass = (biomass - mean(biomass))/sd(biomass)) %>% 
  filter(z_biomass < 4) %>% 
  # filter out plots that are fertilized
  filter(harm_fert_appl == 0) %>% 
  # select relevant columns
  dplyr::select(cont_ID, plot_ID, country, cell, dg_cell, cell_lon, cell_lat, 
                plot_size, sample_year = year, 
                lat, lon, elev, spec_ric, q1, q2, biomass, spec_pool,
                ndep = sum_5yr, MAT, MAP, PET, pH, pH_field, pH_soilgrids, N, P,
                K, NP, lim2, habitat = sub_class, data_source) %>% 
  # remove globnut plots with incomplete data
  drop_na(lon, lat, spec_ric, q1, q2, biomass, N, P, MAT, MAP, pH,
          ndep,  plot_size, habitat, spec_pool) 
saveRDS(globnut, "outputs/02_GlobNut_critical_ratios_palpurina.rds")