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
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(raster)) install.packages("raster")

## Load data
data_dir <- "~/Data/Globnut_offline/_GLOBNUT1.0/"
ndepos <- read.csv(paste0(data_dir, "ndeposition.csv"))
clim <- read.csv(paste0(data_dir, "ERA5_climate.csv"))
elev <- readRDS("env_data/outputs/elevation.rds")
pH <- readRDS("outputs/pH-data.rds")
eunis <- read.csv(paste0(data_dir, "Globnut1.0_EUNIS.csv"))
npk <- read.csv(paste0(data_dir, "GlobNut1.0_nutrients.csv")) %>% 
  # add column with nutrient limitation
  mutate(lim = case_when(NP <= 13.5 & NK <= 2.1 ~ "N-limitation", 
                         NP > 16 & KP > 3.4 ~ "P-limitation", 
                         NK > 2.1 & KP <= 3.4 ~  "K(co)-limitation", 
                         NP >= 13.5 & NP <= 16 ~ "Co-limitation N-P",
                         P >= 0.11 & N >= 2 & K >= 0.8 ~ "No limitation by N, P, K",
                         TRUE ~ "No limitation by N, P, K"))
meta <- read.csv(paste0(data_dir, "GlobNut1.0_metadata.csv"))
grid <- readRDS("outputs/01_Globnut_grid_res15.rds") %>%
  rename(plot_ID = globnut.plot_ID) %>%
  dplyr::select(-lat, -lon, dg_cell = cell)
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
  left_join(grid, by = "plot_ID")
  
saveRDS(globnut_raw, "outputs/02_Globnut_raw.rds")

# selection of plots for analysis
globnut <- globnut_raw %>% 
  # select plots in EUrasia
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  # filter out plots that have been fertilized
  filter(!harm_fert_appl %in% c(1,2)) %>% 
  # filter K-limited and no limitation plots
  filter(!lim %in% c("K(co)-limitation", "No limitation by N, P, K")) %>% 
  # filter out outliers with high biomass 
  drop_na(biomass) %>% 
  mutate(z_biomass = (biomass - mean(biomass))/sd(biomass)) %>% 
  filter(z_biomass < 4) %>% 
  # select relevant columns
  dplyr::select(cont_ID, plot_ID, country, cell, dg_cell, cell_lon, cell_lat, 
                plot_size, sample_year = year, 
                lat, lon, elev, spec_ric, q1, q2, biomass,
                ndep = sum_5yr, MAT, MAP, PET, pH, pH_field, pH_soilgrids, N, P,
                K, NP, lim, habitat = sub_class) %>% 
  # remove globnut plots with incomplete data
  drop_na(lon, lat, spec_ric, q1, q2, biomass, N, P, MAT, MAP, pH,
          ndep,  plot_size, habitat) 
saveRDS(globnut, "outputs/02_GlobNut.rds")


# sub sampled data set to account for oversampled areas - max 5 plots per 3.6km2 cell
set.seed(123)
globnut_samp <- globnut %>% 
  group_by(dg_cell) %>%
  slice_sample(n = 5, replace = FALSE)
saveRDS(globnut_samp, "outputs/02_GlobNut_subsampled.rds")

