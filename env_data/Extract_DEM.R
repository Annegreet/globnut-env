## ---------------------------
##
## Script name: 01_extract_DEM
##
## Purpose of script: Extract height values for plots in globnut 
##
## Author: Annegreet Veeken
##
## Date Created: 2023-06-05
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
## - download NASA data script by Julian Schrader (Globnut/geodata/DEM_download_globnut.R)   
## - how to cite: https://lpdaac.usgs.gov/
## References:
##
## ---------------------------


## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(terra)) install.packages("terra")
if(!require(httr)) install.packages("httr")

## Load data
# globnut coordinates
globnut <- read.csv("Z:/_GLOBNUT1.0/GlobNut1.0_metadata.csv")

# NASA DEM ----
## create df with coordinates and filenames for the .hgt files
plot_coo <-  globnut[,c("plot_ID", "lat", "lon")]

plot_coo <-  as.data.frame(sapply(plot_coo, as.numeric))

#deletes all erroneous coordinates (should ideally be none)
plot_coo <-  plot_coo[which(plot_coo$lat > -90 & plot_coo$lat < 90),]
plot_coo <-  plot_coo[which(plot_coo$lon > -180 & plot_coo$lon < 180),]

#First, create two empty cols
plot_coo$lat_DEM <-  NA
plot_coo$lon_DEM <-  NA


for(i in 1:nrow(plot_coo)){
  
  #For Lat: Rounds coordinates, deletes "-" in case its in southern hemisphere and adds "N" and "S"
  if(plot_coo$lat[i] >= 0) {plot_coo$lat_DEM[i] = paste0("N", "0", abs(floor(plot_coo$lat[i])), sep = "")
  } else {plot_coo$lat_DEM[i] = paste("S", "0", abs(floor(plot_coo$lat[i])), sep = "")}
  
  #For Long: Rounds coordinates, deletes "-" in case its in west  and adss "E" and "W"
  if(plot_coo$lon[i] >= 0) {plot_coo$lon_DEM[i] = paste("E", "00", abs(floor(plot_coo$lon[i])), sep = "")
  } else {plot_coo$lon_DEM[i] = paste("W", "00", abs(floor(plot_coo$lon[i])), sep = "")}                                                   
  
  #Makes sure Lat names have 3 characters 
  if(nchar(plot_coo$lat_DEM[i]) > 3) {str_sub(plot_coo$lat_DEM[i], 2, -3) = ""; plot_coo$lat_DEM[i]}
  
  #Makes sure Long names have 4 characters
  if(nchar(plot_coo$lon_DEM[i]) > 4) {str_sub(plot_coo$lon_DEM[i], 2, -4) = ""; plot_coo$lon_DEM[i]}
  if(nchar(plot_coo$lon_DEM[i]) > 4) {str_sub(plot_coo$lon_DEM[i], 2, -4) = ""; plot_coo$lon_DEM[i]}
  
}

plot_coo_unique <- unique(plot_coo[,c("lat_DEM", "lon_DEM")])
# add file name column
plot_coo$dem_file <- paste0(plot_coo$lat_DEM, plot_coo$lon_DEM)

## Download data from NASA 
#Specify path here where download data should be stored
path <- "Z:/Organized-globnut/Geo-data/Elevation_DEM/20230605-DEM30_download/"

#This loop is for downloading the DEM30 data from NASA
# for(i in 1:nrow(plot_coo_unique)){
#   print(i)
#   
#   destfile = paste(path, plot_coo_unique$lat_DEM[i], plot_coo_unique$lon_DEM[i], ".zip", sep = "")
#   
#   if(paste(plot_coo_unique$lat_DEM[i], plot_coo_unique$lon_DEM[i], ".zip", sep = "") %in% list.files(path)){next} #only downloads new tiles
#   
#   
#   path_download = paste("https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/", 
#                         plot_coo_unique$lat_DEM[i], plot_coo_unique$lon_DEM[i], ".SRTMGL1.hgt.zip", sep = "")
#   
#   GET(path_download, 
#       authenticate("USERNAME", "PASSWORD"),#add username and PW here
#       write_disk(destfile)) # destfile is an empty file that corresponds to the file of the download including the path (similar to the destfile in the download.file function)
#   
# }
# 

# File path and name of the .hgt file
hgt_files <- list.files(path, full.names = T)

# Unzip the .hgt files
# purrr::map(hgt_files, ~unzip(., exdir = path, overwrite = FALSE))

# Get the list of unzipped files
unzipped_files <- list.files(path, pattern = ".hgt$", full.names = TRUE)
dem_names <- list.files(path, pattern = ".hgt$") %>% 
  str_remove(., pattern = ".hgt")

hgt_list <- list()
for(i in 1:length(dem_names)){
  hgt <- unzipped_files %>% 
    # filter and relevant raster
    str_subset(., pattern = dem_names[i]) %>% 
    terra::rast(.) %>% 
    # extract height from relevant raster
    terra::extract(x = ., y = plot_coo[plot_coo$dem_file == dem_names[i], c("lon","lat")]) %>% 
    dplyr::select(-1)
  # add plot_id, lat and lon
  hgt$plot_ID <- plot_coo$plot_ID[plot_coo$dem_file == dem_names[i]]
  hgt$lat <- plot_coo$lat[plot_coo$dem_file == dem_names[i]]
  hgt$lon <- plot_coo$lon[plot_coo$dem_file == dem_names[i]]
  colnames(hgt) <- c("value", "plot_ID", "lat", "lon")
  # append to list
  hgt_list[[length(hgt_list)+1]] <- hgt
  }

# create 1 df out of list
hgt <- bind_rows(hgt_list) %>%
  # NA values for missing observations
  right_join(globnut[,c("plot_ID","lat", "lon")], by = c("plot_ID","lat","lon")) %>% 
  mutate(data_source = "NASA",
         var_name = "elev",
         unit = "m.a.s.l.", 
         description = "NASA DEM extracted meters above sea level",
         data_url = "https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/",
         data_citation = NA,
         orig_res = NA,
         obs_year = NA) %>% 
  # order columns
  dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
                orig_res, data_url, data_citation)
saveRDS(hgt, "GlobNut_Env_var/Outputs/NASA_elevation.rds")

# # Terrain variables
# ter_list <- list()
# for(i in 1:length(dem_names)){
#   ter <- unzipped_files %>% 
#     # filter and relevant raster
#     str_subset(., pattern = dem_names[i]) %>% 
#     rast(.) %>% 
#     # calculate terrain variables in spatraster format
#     terra::terrain(., v = c("slope", "aspect", "TRI"), unit = "degrees",
#                    neighbors = 8 ) %>% 
#     # calculate terrain variable
#     terra::extract(x = ., y = plot_coo[plot_coo$dem_file == dem_names[i], c("lon","lat")]) %>% 
#     dplyr::select(-1)
#   # add plot_id, lat and lon
#   ter$plot_ID <- plot_coo$plot_ID[plot_coo$dem_file == dem_names[i]]
#   ter$lat <- plot_coo$lat[plot_coo$dem_file == dem_names[i]]
#   ter$lon <- plot_coo$lon[plot_coo$dem_file == dem_names[i]]
#   # append to list
#   ter_list[[length(ter_list)+1]] <- ter
# }
# 
# var_desc <- data.frame(var_name = c("slope", "aspect", "TRI"),
#                        unit = c("degrees", "degrees", NA),
#                        description = c("Slope based on 8 neighboring cells, terra terrain function",
#                                        "Aspect based on 8 neighboring cells, terra terrain function",
#                                        "Terrain Ruggedness Index, mean of differences between 8 neighboring cells, terra terrain function"))
# ter <- bind_rows(ter_list) %>% 
#   # add NA values for missing observations
#   right_join(globnut[,c("plot_ID","lat", "lon")], by = c("plot_ID","lat","lon")) %>% 
#   pivot_longer(cols = slope:TRI, names_to = "var_name", values_to = "value") %>% 
#   left_join(var_desc, by = "var_name") %>% 
#   mutate(data_source = "NASA",
#          data_url = "https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/",
#          data_citation = NA,
#          orig_res = NA,
#          obs_year = NA) %>% 
#   # order columns
#   dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
#                 orig_res, data_url, data_citation) %>% 
#   bind_rows(hgt)
# 
# base::saveRDS(ter, file = "GlobNut_Env_var/Outputs/NASA_topography.rds")

# Arctic DEM ----
# NASA DEM doesn't cover the polar region, use arcticDEM instead https://www.pgc.umn.edu/data/arcticdem/
# ArcticDEM is tile 
# identify tile names to with Globnut plots
arctic_mos <- vect("Z:/Organized-globnut/Geo-data/Elevation_DEM/ArcticDEM_Mosaic_Index_latest_shp/ArcticDEM_Mosaic_Index_v4_1_10m.shp")

plot_coo$id.y <- 1:nrow(plot_coo) 
mosaic_names <- arctic_mos %>% 
  project("epsg:4326") %>% 
  terra::extract(x = ., y = plot_coo[, c("lon","lat")]) %>% 
  left_join(plot_coo[,c("id.y","plot_ID", "lon","lat")], by = "id.y") %>% 
  drop_na() 

mosaic_ID <- mosaic_names %>% 
  pull(tile) %>% 
  unique()

# path_download <- paste0("https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/10m/", 
#                       mosaic_ID,"/",mosaic_ID, "_10m_v4.1.tar.gz")
# path_save <- paste0("Z:/geo_data/Elevation_DEM/ArcticDEM_download/",mosaic_ID, ".tar.gz")
path <- "Z:/Organized-globnut/Geo-data/Elevation_DEM/ArcticDEM_download/"

# for(i in 1:length(mosaic_ID)){
#   print(i)
#   
#   destfile = paste(path, mosaic_ID[i], "_10m_v4.1.tar.gz", sep = "")
#   
#   if(paste(mosaic_ID[i], "_10m_v4.1.tar.gz", sep = "") %in% list.files(path)){next} #only downloads new tiles
#   
#   path_download = paste("https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/10m/", 
#                         mosaic_ID[i],"/",mosaic_ID[i], "_10m_v4.1.tar.gz", sep = "")
#   
#   GET(path_download, 
#       write_disk(destfile)) # destfile is an empty file that corresponds to the file of the download including the path (similar to the destfile in the download.file function)
#   untar(destfile, files = "*_dem.tif", exdir = path)
# }

unzipped_files <- list.files(path, pattern = "dem.tif", full.names = TRUE)

hgt_list_arc <- list()

for(i in 1:length(mosaic_ID)){
  print(mosaic_ID[i])
  
  # convert globnut coordinate to espg:3413 
  globnut_vec <- mosaic_names[mosaic_names$tile == mosaic_ID[i], c("lon","lat")] %>%  
    vect(., crs = "EPSG:4326", geom=c("lon", "lat")) %>% 
    # project to EPSG:3413 from EPSG:4326
    terra::project(., "EPSG:3413") 
  
  hgt <- unzipped_files %>% 
    # filter and relevant raster
    str_subset(., pattern = paste0("/", mosaic_ID[i])) %>% 
    rast(.) %>% 
    # extract height from relevant raster
    terra::extract(x = ., y = globnut_vec) %>% 
    dplyr::select(-1)
  # add plot_id, lat and lon
  hgt$plot_ID <- mosaic_names$plot_ID[mosaic_names$tile == mosaic_ID[i]]
  hgt$lat <- mosaic_names$lat[mosaic_names$tile == mosaic_ID[i]]
  hgt$lon <- mosaic_names$lon[mosaic_names$tile == mosaic_ID[i]]
  colnames(hgt) <- c("value", "plot_ID", "lat", "lon")
  # append to list
  hgt_list[[length(hgt_list)+1]] <- hgt
}

hgt <- bind_rows(hgt_list) %>%
  filter(!is.na(value)) %>% 
  mutate(data_source = "ArcticDEM",
         var_name = "elev",
         unit = "m.a.s.l.", 
         description = "ArcticDEM extracted meters above sea level",
         data_url = "https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/10m/",
         data_citation = "Porter, Claire; Morin, Paul; Howat, Ian; Noh, Myoung-Jon; Bates, Brian; Peterman, Kenneth; Keesey, Scott; Schlenk, Matthew; Gardiner, Judith; Tomko, Karen; Willis, Michael; Kelleher, Cole; Cloutier, Michael; Husby, Eric; Foga, Steven; Nakamura, Hitomi; Platson, Melisa; Wethington, Michael, Jr.; Williamson, Cathleen; Bauer, Gregory; Enos, Jeremy; Arnold, Galen; Kramer, William; Becker, Peter; Doshi, Abhijit; D’Souza, Cristelle; Cummens, Pat; Laurier, Fabien; Bojesen, Mikkel, 2018, “ArcticDEM”, https://doi.org/10.7910/DVN/OHHUKH, Harvard Dataverse, V1, [Date Accessed: 17-8-2023]",
         orig_res = "10 m",
         obs_year = NA) %>% 
  # order columns
  dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
                orig_res, data_url, data_citation)
base::saveRDS(hgt, file = "GlobNut_Env_var/Outputs/ArcticDEM_elevation.rds")

# Terrain for arctic dem not currently working properly - only retrieving 
# ter_list <- list()
# for(i in 1:length(mosaic_ID)){
#   print(mosaic_ID[i])
#   
#   # convert globnut coordinate to espg:3413 
#   globnut_vec <- mosaic_names[mosaic_names$tile == mosaic_ID[i], c("lon","lat")] %>%  
#     vect(., crs = "EPSG:4326", geom=c("lon", "lat")) %>% 
#     # project to EPSG:3413 from EPSG:4326
#     terra::project(., "EPSG:3413") 
# 
#   ter <- unzipped_files %>% 
#     # filter and relevant raster
#     str_subset(., pattern = paste0("/", mosaic_ID[5])) %>% 
#     rast(.) %>% 
#     # calculate terrain variables in spatraster format
#     terra::terrain(., v = c("slope", "aspect", "TRI"), unit = "degrees",
#                    neighbors = 8 ) %>% 
#     # calculate terrain variable
#     terra::extract(x = ., y = globnut_vec) %>% 
#     dplyr::select(-1)
#   # add plot_id, lat and lon
#   ter$plot_ID <- mosaic_names$plot_ID[mosaic_names$tile == mosaic_ID[i]]
#   ter$lat <- mosaic_names$lat[mosaic_names$tile == mosaic_ID[i]]
#   ter$lon <- mosaic_names$lon[mosaic_names$tile == mosaic_ID[i]]
#   # append to list
#   ter_list[[length(ter_list)+1]] <- ter
# }
# 
# var_desc <- data.frame(var_name = c("slope", "aspect", "TRI"),
#                        unit = c("degrees", "degrees", NA),
#                        description = c("Slope based on 8 neighboring cells, terra terrain function",
#                                        "Aspect based on 8 neighboring cells, terra terrain function",
#                                        "Terrain Ruggedness Index, mean of differences between 8 neighboring cells, terra terrain function"))
# ter <- bind_rows(ter_list) %>% 
#   pivot_longer(cols = slope:TRI, names_to = "var_name", values_to = "value") %>% 
#   filter(!is.na(value)) %>% 
#   left_join(var_desc, by = "var_name") %>% 
#   mutate(data_source = "ArcticDEM",
#          var_name = "elev",
#          unit = "m.a.s.l.", 
#          description = "ArcticDEM extracted meters above sea level",
#          data_url = "https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/latest/10m/",
#          data_citation = "Porter, Claire; Morin, Paul; Howat, Ian; Noh, Myoung-Jon; Bates, Brian; Peterman, Kenneth; Keesey, Scott; Schlenk, Matthew; Gardiner, Judith; Tomko, Karen; Willis, Michael; Kelleher, Cole; Cloutier, Michael; Husby, Eric; Foga, Steven; Nakamura, Hitomi; Platson, Melisa; Wethington, Michael, Jr.; Williamson, Cathleen; Bauer, Gregory; Enos, Jeremy; Arnold, Galen; Kramer, William; Becker, Peter; Doshi, Abhijit; D’Souza, Cristelle; Cummens, Pat; Laurier, Fabien; Bojesen, Mikkel, 2018, “ArcticDEM”, https://doi.org/10.7910/DVN/OHHUKH, Harvard Dataverse, V1, [Date Accessed: 17-8-2023]",
#          orig_res = "10 m",
#          obs_year = NA) %>% 
#   # order columns
#   dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
#                 orig_res, data_url, data_citation) %>% 
#   bind_rows(hgt)


## Compile results NASA and arctic DEM
nasa_elev <- readRDS("GlobNut_Env_var/Outputs/NASA_elevation.rds")
arctic_elev <- readRDS("GlobNut_Env_var/Outputs/ArcticDEM_elevation.rds")

elev <- nasa_elev %>% 
  left_join(arctic_elev %>% 
              dplyr::select(plot_ID, value), by = "plot_ID", suffix = c("_nasa", "_arctic")) %>% 
  mutate(elev = ifelse(is.na(value_nasa), value_arctic, value_nasa)) %>% 
  select(plot_ID, elev) %>% 
  distinct()
saveRDS(elev, "~/Projects/Ongoing/globnut-env/env_data/outputs/elevation.rds")

