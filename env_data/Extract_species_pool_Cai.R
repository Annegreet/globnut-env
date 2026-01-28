## ---------------------------
##
## Script name: Extract species pool Cai et al    
##
## Purpose of script: Extract species richness information from Cai et al paper, to approximate regional species pool
##
## Author: Annegreet Veeken
##
## Date Created: 2025-07-16
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes: downloaded data is from the XGboost model with cubic polynomial trend surface https://gift.uni-goettingen.de/shiny/predictions/ 
##  
## References: https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.18533
##
## ---------------------------

## Load packages
if(!require(tidyverse)) install.packages(tidyverse)
if(!require(dggridR)) install.packages(dggridR)

## Load data
load("~/Data/Cai-etal/sr_Random-Forest-trend-surface_Prediction_7774_LongLat.RData")
globnut <- read_rds("outputs/02_GlobNut.rds")

# assign globnut plots to same grid as cai
dggs <- dgconstruct(res = 8, metric = TRUE) # see table for different resolutions https://github.com/r-barnes/dggridR or dggetres()
globnut$gr_8_ID <- dgGEO_to_SEQNUM(dggs,globnut$lon,globnut$lat)$seqnum

# join datasets
spec_pool <- left_join(globnut, predictions_grid, by = c( "gr_8_ID")) %>% 
  select(plot_ID, 
         ndep ,
         spec_pool = value)

saveRDS(spec_pool, "env_data/outputs/Cai_spec_pool.rds")

