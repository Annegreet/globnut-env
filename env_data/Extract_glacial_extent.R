## ---------------------------
##
## Script name: Extract_glacial_extent.R
##
## Purpose of script: 
##
## Author: Annegreet Veeken
##
## Date Created: 2023-07-12
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##  
## References:
## - https://envidat.ch/#/metadata/chelsa_trace
## - https://osf.io/7jen3/?view_only=
## - RABASSA, J., CORONATO, A. & MARTÍNEZ, O. (2011). Late Cenozoic glaciations in Patagonia and Tierra del Fuego: an updated review. Biological Journal of the Linnean Society, 103(2), 316–335. https://doi.org/10.1111/j.1095-8312.2011.01681.x
## ---------------------------

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(raster)) install.packages("raster")
if(!require(terra)) install.packages("terra")

## Load data
# globnut coordinates
data_dir <- "Z:/_GLOBNUT1.0/"
meta <- read.csv(paste0(data_dir, "GlobNut1.0_metadata.csv"))

## CHELSA data (LGM) ----
if(1){
# Glacial extent from CHELSA-TraCE21k - directories
f <- list.files(
  "Z:/geo_data/CHELSA_glacial_extent",
  pattern = "*.tif",
  full.names = T) %>% 
  str_subset(pattern = "gle") # filter .tif files for glacial extent (glz files are elevation glacial elevation m.a.s.l)
# load, stack and convert to SpatRaster
glacial <- stack(f) %>% rast()

# Rename layers
# look-up table time BP (see documentation at download website https://envidat.ch/#/metadata/chelsa_trace)
time_lookup <- data.frame(time_ID = 20:-200 %>% as.character() %>% 
                            str_replace(pattern = "-", replacement = "."), # - in file names gets replaced by . when loading and stacking 
                          time_BP = seq(0, 22000, by = 100))

glacial_globnut <- cbind(meta[, c("plot_ID","lon","lat")],
                         terra::extract(x = glacial, y = data.frame(meta[, c("lon","lat")]))) %>% 
  pivot_longer(cols = starts_with("CHELSA_TraCE21k_gle_"), 
               names_to = "time_ID", values_to = "glacial") %>% 
  mutate(time_ID = str_remove(time_ID, "CHELSA_TraCE21k_gle_") %>% 
           str_remove("_V1.0")) %>% 
  left_join(time_lookup, by = "time_ID") 

chelsa_lgm <- glacial_globnut %>% 
  # extract time point of deglaciation for every site
  group_by(plot_ID) %>% 
  arrange(time_BP, .by_group = TRUE) %>% 
  mutate(deglac = case_when(all(is.na(glacial)) ~ ">21000",# glaciation is longer than 21000 years ago
                            all(glacial == 1) ~ "current glaciation", # current glaciation
                            any(is.na(glacial)) ~ "deglaciation")) %>% 
  slice(which(first(deglac == "current glaciation")), first(which(!is.na(glacial))),
        which(first(deglac == ">21000"))) %>% 
  mutate(value = case_when(deglac == "current glaciation" ~ 0,
                               deglac == ">21000" ~ Inf,
                               deglac == "deglaciation" ~ time_BP),
          obs_year = NA, 
           data_source = "CHELSA_TraCE21k",
           var_name = "glacial_extent",
           description = "Years since last glaciation - Inf = no glaciation for more than 21000 years",
           ID = NULL, # remove column
           unit = "years ago",
           orig_res = "30 arcsec",
           data_url = "https://envidat.ch/#/metadata/chelsa_trace",
           data_citation = "Karger, D. N., Nobis, M. P., Normand, S., Graham, C. H., Zimmermann, N. E. (2020). CHELSA-TraCE21k: Downscaled transient temperature and precipitation data since the last glacial maximum. EnviDat. https://www.doi.org/10.16904/envidat.211.") %>%
  # order columns
    dplyr::select(plot_ID, obs_year, var_name, value, unit, description, data_source,
                  orig_res, data_url, data_citation)
 
}

## Before LGM (Northern hemisphere)----
if(1){
# get file paths of shapefiles by time slice
file_paths <- list.files("Z:/geo_data/Glacial_extent/", recursive = TRUE, include.dirs = TRUE, full.names = TRUE) %>% 
  str_subset("hypothesised ice-sheet reconstructions") %>% 
  str_subset("best_estimate.shp") %>% 
  str_subset(".xml", negate = TRUE)

# read in a shape file to check
glac <- vect(file_paths[7]) %>% terra::union()
# reproject plot coordinates to crs of shape files
plot_coord <- meta %>%
  select(plot_ID, lon, lat) %>%
  na.omit() %>%
  # convert to vector
  vect(., geom = c("lon", "lat"), crs = "epsg:4326") %>%
  # reproject
  project(crs(glac))
# check
ggplot(st_as_sf(glac)) +
  geom_sf() +
  geom_sf(data = st_as_sf(plot_coord))


# find intersection for all shape files
glac_time <- list()
for (i in 1:length(file_paths)) {
  # load vector
  glac_shp <- vect(file_paths[i])
  # extract time period
  time_period <- str_extract(file_paths[i], pattern = "\\/([^\\/]*)(?=\\_best_estimate\\.shp)") %>% 
    str_remove("/")
  print(time_period)
  
  # reproject plot coordinates to crs of shape file (LGM is in another projection than others - so do it in loop)
  plot_coord <- meta %>% 
    dplyr::select(plot_ID, lon, lat) %>% 
    na.omit() %>% 
    # convert to vector
    vect(., geom = c("lon", "lat"), crs = "epsg:4326") %>% 
    # reproject to crs of shp file
    project(crs(glac_shp))
  
  # extract glacial extent by globnut coordinate
  glac_point <- cbind(time_period = time_period, 
                      # data.frame(plot_coord), 
                      terra::extract(glac_shp, plot_coord)) %>% 
    # remove empty
    filter(!if_all(3:ncol(.), is.na))
  
  # append time slice to list
  glac_time[[length(glac_time) + 1]] <- glac_point
}

## Create df with time since last glaciation by plot
plot_ids <- data.frame(plot_coord) %>% rowid_to_column()
north_glac <- bind_rows(glac_time) %>% 
  # add plot_ID
  left_join(plot_ids, by = c("id.y" = "rowid")) %>% 
  dplyr::select(plot_ID, time_period) %>% 
  # set order to time period variable
  mutate(time_period = factor(time_period, levels = c("LGM",
                                                      "30ka",
                                                      "35ka",
                                                      "40ka",
                                                      "45ka", "MIS4","MIS5a","MIS5b", "MIS5c",
                                                      "MIS5d","MIS6","MIS8","MIS10",
                                                      "MIS12", "MIS16","MIS20-24",
                                                      "EarlyMatuyama",
                                                      "LateGauss"), ordered = TRUE)) %>% 
  group_by(plot_ID) %>% 
  summarise(last_glac = min(time_period), # time period of last glaciation on the site
            n_glac = n_distinct(time_period)) 

}

## Compile soil age variable per plot ----
age <- chelsa_lgm %>% 
  left_join(meta, by = "plot_ID") %>% 
  left_join(north_glac, by = "plot_ID") %>% 
  dplyr::select(plot_ID, lat, lon, lgm_glacial_extent = value, last_glac)  %>% 
  dplyr::filter(!is.na(lat)) %>% 
  mutate(soil_age = case_when(lgm_glacial_extent < 21000 ~ "1", # 21 000 years or less
                              last_glac %in% c("30ka","35ka","40ka","MIS4","Mis5a","Mis5b","MIS5c","MIS5d", "MIS6","MIS8", "MIS10", "MIS12") ~ "2", #30ka and 4770000 
                              is.na(last_glac) ~ "3")) # older than 477000
write.csv(age, "Z:/_GLOBNUT1.0/soilage.csv")
saveRDS(age, "env_data/outputs/Soilage.rds")
