## ---------------------------
##
## Script name: 01_Gridding_globnut.R
##
## Purpose of script: Gridding GlobNut database, reducing sampling bias
##
## Author: Annegreet Veeken
##
## Date Created: 2023-10-19
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##  
## References:
## - script based on dggridR vignette: https://cran.r-project.org/web/packages/dggridR/vignettes/dggridR.html
## ---------------------------

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(sf)) install.packages("sf")
if (!require(dggridR)) install.packages("dggridR")
if (!require(rnaturalearth)) install.packages("rnaturalearth")
if (!require(rnaturalearthdata)) install.packages("rnaturalearthdata")

## Load data
data_dir <- "Z:/_GLOBNUT1.0/" # directory with raw Globnut 1.0 data
meta_env <- read.csv(paste0(data_dir, "Globnut1.0_env_variable_meta.csv"))
globnut <- read.csv(paste0(data_dir, "Globnut1.0_metadata.csv")) %>% 
  filter(!is.na(lat) | !is.na(lon)) 

# Construct a global grid with cells approximately 3.6 km2 (res 15), meaning the distance between center of adjacent cells is ~1.9km
dggs <- dgconstruct(res = 15, metric = TRUE) # see table for different resolutions https://github.com/r-barnes/dggridR or dggetres()

# Get the corresponding grid cells for each plot
globnut$cell <- dgGEO_to_SEQNUM(dggs,globnut$lon,globnut$lat)$seqnum

# check the plot count by cell
plotcounts  <- globnut %>% group_by(cell) %>% summarise(count = n())
quantile(plotcounts$count)

# get center coordinates of the cells
gn_grid <- dgSEQNUM_to_GEO(dggs,globnut$cell) %>% 
  data.frame(globnut$plot_ID, cell = globnut$cell, lon = .$lon_deg, lat = .$lat_deg) 
# saveRDS(gn_grid, "outputs/01_Globnut_grid_res15.rds")





