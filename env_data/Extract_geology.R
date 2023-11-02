## ---------------------------
##
## Script name: Extract_geology.R
##
## Purpose of script: 
##
## Author: Annegreet Veeken
##
## Date Created: 2023-07-20
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes: sf package was used because it is supposedly faster than terra
##  
## References:
##
## ---------------------------


## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(terra)) install.packages("terra")
if(!require(sf)) install.packages("sf")

## Load data
globnut <- read.csv("Z:/_Globnut1.0/GlobNut1.0_metadata.csv")
glim <- rast("Z:/geo_data/GLiM/glim_wgs84_0point5deg.txt.asc")# global lithological map
glim_names <- read.table("Z:/geo_data/GLiM/Classnames.txt", sep = ";",
                         header = TRUE) %>% 
  rename(lith_code = xx, lith_ID = Value_) %>% 
  mutate(description = case_when(lith_code == "su" ~ "Unconsolidated Sediments",
                                 lith_code == "vb" ~ "Basic volcanic rocks",
                                 lith_code == "ss" ~ "Siliciclastic sedimentary rocks",
                                 lith_code == "pb" ~ "Basic plutonic rocks",
                                 lith_code == "sm" ~ "Mixed sedimentary rocks",
                                 lith_code == "sc" ~ "Carbonate sedimentary rocks",
                                 lith_code == "va" ~ "Acid volcanic rocks",
                                 lith_code == "mt" ~ "Metamorphics",
                                 lith_code == "pa" ~ "Acid plutonic rocks",
                                 lith_code == "vi" ~ "Intermediate volcanic rocks",
                                 lith_code == "wb" ~ "Water Bodies",
                                 lith_code == "py" ~ "Pyroclastics mentioned",
                                 lith_code == "pi" ~ "Intermediate plutonic rocks",
                                 lith_code == "ev" ~ "Evaporites",
                                 lith_code == "nd" ~ "No Data",
                                 lith_code == "ig" ~ "Ice and Glaciers"))

# GliM - lithology
crs(glim)
lith_globnut <- cbind(globnut[, c("plot_ID","lon","lat")], 
                      terra::extract(x = glim, 
                                     y = data.frame(globnut[, c("lon","lat")]))) %>% 
  dplyr::select(-ID, lith_ID = glim_wgs84_0point5deg.txt) %>% 
  left_join(glim_names, by = "lith_ID") %>% 
  mutate(var_name = "lith",
         unit = NA,
         data_citation = "Hartmann, Jörg; Moosdorf, Nils (2012): Global Lithological Map Database v1.0 (gridded to 0.5° spatial resolution). PANGAEA, https://doi.org/10.1594/PANGAEA.788537",
         data_url = "https://doi.pangaea.de/10.1594/PANGAEA.788537",
         data_source = "GliM",
         orig_res = "0.5x0.5 degree",
         obs_year = NA) %>% 
  dplyr::select(plot_ID, lat, lon, obs_year, var_name, value_chr = lith_code, unit, description, data_source, 
                orig_res, data_url, data_citation)
saveRDS(lith_globnut, file = "env_data/outputs/GLiM_lithology.rds")
write.csv(lith_globnut, file = "Z:/_GLOBNUT1.0/lithology.csv")
# plot
lith_globnut %>%
  group_by(var_name) %>% 
  count(description) %>% 
  ggplot(aes(x = description, y = n)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~var_name, scales = "free")
