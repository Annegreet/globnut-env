## ---------------------------
##
## Script name: Extract_N_deposition.R
##
## Purpose of script: Extracting Nitrogen deposition data from the EMEP and Ackerman dataset
##
## Author: Annegreet Veeken
##
## Date Created: 2023-06-22
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##  
## References:
## https://www.emep.int/mscw/mscw_moddata.html
## https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018GB005990
## ---------------------------


## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ncdf4)) install.packages("ncdf4")
if (!require(ncdump)) install.packages("ncdump")
if (!require(raster)) install.packages("raster")
if (!require(terra)) install.packages("terra")
if (!require(sf)) install.packages("sf")
if (!require(rnaturalearth)) install.packages("rnaturalearth")
if (!require(rnaturalearthdata)) install.packages("rnaturalearthdata")

## EMEP data (europe/russia) ----
if (1) {
## Download EMEP data (only need to do this once!) https://www.emep.int/mscw/mscw_moddata.html
  # The NetCDF files have the following convention:
  #   {GRID}_{MODEL VERSION}_{TIME RESOLUTION}.{YEAR}met_{EMISSION YEAR}emis_({REPORTING YEAR}).nc
globnut_dir <- "Z:/geo_data/EMEP/Data_EMEP_report2023/"
years <- 1990:2020
emep_urls <- paste0("https://thredds.met.no/thredds/fileServer/data/EMEP/2023_Reporting/EMEP01_rv5.0_year.",
                    years, "met_", years,"emis_rep2023.nc") # EMEP is available from 1990-2020 on a yearly, monthly, daily and hourly resolution. Here yearly is dowloaded for last 10 year
# 2021 and 2022 has another url
emep_urls <- c(emep_urls, "https://thredds.met.no/thredds/fileServer/data/EMEP/2023_Reporting/EMEP01_rv5.0_year.2021met_2021emis.nc", "https://thredds.met.no/thredds/fileServer/data/EMEP/2023_Reporting/EMEP01_rv5.0_year.2022met_2021emis.nc")
file_names <- c(paste0(globnut_dir,"EMEP01_rv5.0_year.",
                     years, "met_",years,"emis_rep2023.nc"), 
                paste0(globnut_dir,c("EMEP01_rv5.0_year.2021met_2021emis.nc", "EMEP01_rv5.0_year.2022met_2021emis.nc")))
names(file_names) <- c(years, 2021, 2022)
# download
# map2(emep_urls, file_names, ~download.file(url = .x, destfile = .y, mode = 'wb'))


## Load data
# globnut coordinates
globnut <- read.csv("Z:/_GLOBNUT1.0/GlobNut1.0_metadata.csv") %>% 
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  dplyr::select(plot_ID, lon, lat) %>%
  na.omit()

## Extract data from EMEP grid
# EMEP data 2020
emep_2020 <- nc_open(file_names["2020"])
emep_meta <- NetCDF(file_names["2020"])
# plot(emep_meta)
print(emep_2020) # check layers in file
# latitude and longitude of the EMEP grid
lat <- ncvar_get(emep_2020, "lat") 
lon <- ncvar_get(emep_2020, "lon")

# europe shapefiles
country <- c("Albania","Aland","Andorra","Austria","Belgium","Bulgaria",
             "Bosnia and Herzegovina","Belarus","Switzerland","Czech Republic",
             "Germany","Denmark","Spain","Estonia","Finland","France","Faroe Islands",
             "United Kingdom","Guernsey","Greece","Croatia","Hungary","Isle of Man",
             "Ireland","Iceland","Italy","Jersey","Kosovo","Liechtenstein","Lithuania",
             "Luxembourg","Latvia","Monaco","Moldova","Macedonia","Malta","Montenegro",
             "Netherlands", "The Netherlands", "Norway","Poland","Portugal","Romania","Russia","San Marino",
             "Republic of Serbia","Slovakia","Slovenia","Sweden","Ukraine","Vatican")
# fix geometry of Russia (sf_use_s2() must be FALSE for this)
sf::sf_use_s2(FALSE)
europe <- ne_countries(scale = "medium", continent = "Europe", 
                       country = country,
                       returnclass = "sf") %>% 
  st_make_valid() %>% 
  # crop to extent of Europe
  st_crop(., c(xmin = min(lon), ymin = min(lat), xmax = max(lon), ymax = max(lat)))
sf::sf_use_s2(TRUE) # back to default

## Read and stack N deposition from EMEP grid
# nitrogen oxides NOx(NO2+NO) (dry)
nox_raster <- stack(as.list(file_names),
                     varname = "DDEP_OXN_m2Grid", lvar = 3, level = 1, ncd = TRUE) %>% 
  # convert to terra spatraster
  rast()
# Reduced Nitrogen NHx (NH3+NH4+)
nhx_raster <- stack(as.list(file_names),
                    varname = "DDEP_RDN_m2Grid", lvar = 3, level = 1, ncd = TRUE) %>% 
  # convert to terra spatraster
  rast()
# wet
wnox_raster <- stack(as.list(file_names),
                    varname = "WDEP_OXN", lvar = 3, level = 1, ncd = TRUE) %>% 
  # convert to terra spatraster
  rast()
wnhx_raster <- stack(as.list(file_names),
                     varname = "WDEP_RDN", lvar = 3, level = 1, ncd = TRUE) %>% 
  # convert to terra spatraster
  rast()


# plot raster
nox_raster %>% 
  as.data.frame(., xy = T) %>% 
  pivot_longer(cols = X1990:X2022, names_to = "year", values_to = "DDEP_OXN_m2Grid") %>% 
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = DDEP_OXN_m2Grid)) +
  geom_sf(data = europe, alpha = 0.5) +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~year)

nhx_raster %>% 
  as.data.frame(., xy = T) %>% 
  pivot_longer(cols = X1990:X2022, names_to = "year", values_to = "DDEP_RDN_m2Grid") %>% 
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = DDEP_RDN_m2Grid)) +
  geom_sf(data = europe, alpha = 0.5) +
  scale_fill_viridis_c(trans = "log") +
  facet_wrap(~year)

## Extract N-deposition by globnut coordinate
#NOX
# dry deposition NOX
nox_globnut <- cbind(globnut[, c("plot_ID","lon","lat")], 
                    terra::extract(x = nox_raster, y = data.frame(globnut[, c("lon","lat")]))) 
# wet deposition NOX
wnox_globnut <- cbind(globnut[, c("plot_ID","lon","lat")], 
                     terra::extract(x = wnox_raster, y = data.frame(globnut[, c("lon","lat")]))) 
tot_nox <- bind_rows(dry = nox_globnut, wet = wnox_globnut, .id = "type") %>% 
  pivot_longer(cols = X1990:X2022, names_to = "obs_year", values_to = "value") %>% 
  group_by(plot_ID,lat,lon, obs_year) %>% 
  # sum wet and dry deposition
  summarise(value = sum(value)) %>% 
  mutate(obs_year = str_remove(obs_year, "X") %>% as.numeric,
         data_source = "EMEP",
         var_name = "NOx_dep",
         description = "wet + dry NOx deposition from EMEP model",
         ID = NULL,
         unit = "mg N/m2",
         orig_res = "0.1x0.1 degree",
         data_url = "https://www.emep.int/mscw/mscw_moddata.html",
         data_citation = "Simpson, D., Benedictow, A., Berge, H., Bergström, R., Emberson, L. D., Fagerli, H., Flechard, C. R., Hayman, G. D., Gauss, M., Jonson, J. E., Jenkin, M. E., Nyíri, A., Richter, C., Semeena, V. S., Tsyro, S., Tuovinen, J.-P., Valdebenito, Á., and Wind, P.: The EMEP MSC-W chemical transport model – technical description, Atmos. Chem. Phys., 12, 7825–7865, https://doi.org/10.5194/acp-12-7825-2012, 2012
         Data produced by EMEP/MSC-W are licensed under Norwegian license for public data (NLOD) and Creative Commons 4.0 BY International. Credit should be given to The Norwegian Meteorological institute, shortened `MET Norway`, as the source of data. Some suggestions: `Data from The Norwegian Meteorological Institute`, `Based on data from MET Norway`") %>% 
  # order columns
  dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
         orig_res, data_url, data_citation)
saveRDS(tot_nox, file = "env_data/outputs/EMEP_NOxdeposition.rds")
saveRDS(tot_nox, file = "env_data/outputs/EMEP_NOxdeposition.rds")
#NHx
# dry deposition NHX
nhx_globnut <- cbind(globnut[, c("plot_ID","lon","lat")], 
                     terra::extract(x = nhx_raster, y = data.frame(globnut[, c("lon","lat")])))
# wet deposition NHX
wnhx_globnut <- cbind(globnut[, c("plot_ID","lon","lat")], 
                      terra::extract(x = wnhx_raster, y = data.frame(globnut[, c("lon","lat")]))) 
tot_nhx <- bind_rows(dry = nhx_globnut, wet = wnhx_globnut, .id = "type") %>% 
  pivot_longer(cols = X1990:X2022, names_to = "obs_year", values_to = "value") %>% 
  group_by(plot_ID,lat,lon, obs_year) %>% 
  # sum wet and dry deposition
  summarise(value = sum(value)) %>% 
  mutate(obs_year = str_remove(obs_year, "X") %>% as.numeric,
         var_name = "NHx_dep",
         data_source = "EMEP",
         description = "wet and dry NHx deposition from EMEP model",
         ID = NULL,
         unit = "mg N/m2",
         orig_res = "0.1x0.1 degree",
         data_url = "https://www.emep.int/mscw/mscw_moddata.html",
         data_citation = "Simpson, D., Benedictow, A., Berge, H., Bergström, R., Emberson, L. D., Fagerli, H., Flechard, C. R., Hayman, G. D., Gauss, M., Jonson, J. E., Jenkin, M. E., Nyíri, A., Richter, C., Semeena, V. S., Tsyro, S., Tuovinen, J.-P., Valdebenito, Á., and Wind, P.: The EMEP MSC-W chemical transport model – technical description, Atmos. Chem. Phys., 12, 7825–7865, https://doi.org/10.5194/acp-12-7825-2012, 2012
         Data produced by EMEP/MSC-W are licensed under Norwegian license for public data (NLOD) and Creative Commons 4.0 BY International. Credit should be given to The Norwegian Meteorological institute, shortened `MET Norway`, as the source of data. Some suggestions: `Data from The Norwegian Meteorological Institute`, `Based on data from MET Norway`") %>% 
  # order columns
  dplyr::select(plot_ID, lat, lon, obs_year, var_name, value, unit, description, data_source, 
         orig_res, data_url, data_citation)

saveRDS(tot_nhx, file = "env_data/outputs/EMEP_NHxdeposition.rds")
saveRDS(tot_nhx, file = "env_data/outputs/EMEP_NHxdeposition.rds")
}

# Supplement missing values with global nitrogen data from Ackerman etal. (2019) ----
if (1) {
# load in data
tot_nhx <- readRDS("env_data/outputs/EMEP_NHxdeposition.rds")
tot_nox <- readRDS("env_data/outputs/EMEP_NOxdeposition.rds")
glob_ndep <- read.csv("Z:/geo_data/Global_n_deposition/oxidized_reduced_N_deposition(1).csv") 

# convert df to vector, then rasterize
glob_ndep_vect <- glob_ndep %>% 
  vect(geom = c("longitude","latitude"), crs = st_crs(4326)) 
# create empty raster with extent and resolution of the data
r <- rast(ext(glob_ndep_vect), ncol = n_distinct(glob_ndep$longitude), 
          nrow = n_distinct(glob_ndep$latitude), crs = st_crs(4326)) 
# oxidized N
glob_nox <- terra::rasterize(glob_ndep_vect, r, "total_oxidized_2016") 
plot(glob_nox)
# reduced N
glob_nhx <- terra::rasterize(glob_ndep_vect, r, "total_reduced_2016") 
plot(glob_nhx)

# # Get plot IDs for plots without data
# no_data <- tot_nhx %>%
#   filter(is.na(value)) %>%
#   pull(plot_ID) %>%
#   unique()
# no_data <- meta %>%
#   filter(plot_ID %in% no_data) %>%
#   dplyr::select(plot_ID, lat, lon) %>%
#   mutate(across(.cols = everything(), ~as.numeric(.))) %>%
#   na.omit()

# extract 
glob_nox_plot <- cbind(globnut[, c("plot_ID","lon","lat")], terra::extract(x = glob_nox, y = globnut[, c("lon","lat")])) %>% 
  rename(value = last) %>% 
  mutate(obs_year = 2016,
         var_name = "NOx_dep_glob",
         data_source = "Ackerman et al 2019",
         description = "wet and dry NOx deposition",
         ID = NULL,
         unit = "mg N/m2",
         orig_res = "2x2.5 degree",
         data_url = "https://conservancy.umn.edu/handle/11299/197613",
         data_citation = "Ackerman, D., Millet, D. B., & Chen, X. (2019). Global estimates of inorganic nitrogen deposition across four decades. Global Biogeochemical Cycles, 33, 100–107. https://doi.org/10.1029/2018GB005990")
glob_nhx_plot <- cbind(globnut[, c("plot_ID","lon","lat")], terra::extract(x = glob_nhx, y = globnut[, c("lon","lat")])) %>% 
  rename(value = last) %>% 
  mutate(obs_year = 2016,
         var_name = "NHx_dep",
         data_source = "Ackerman et al 2019",
         description = "wet and dry NHx deposition",
         ID = NULL,
         unit = "mg N/m2",
         orig_res = "2x2.5 degree",
         data_url = "https://conservancy.umn.edu/handle/11299/197613",
         data_citation = "Ackerman, D., Millet, D. B., & Chen, X. (2019). Global estimates of inorganic nitrogen deposition across four decades. Global Biogeochemical Cycles, 33, 100–107. https://doi.org/10.1029/2018GB005990")
glob_ndep_plot <- bind_rows(glob_nox_plot, glob_nhx_plot)
saveRDS(glob_ndep_plot, "env_data/outputs/Ackerman_ndeposition.rds")
}

## Calculate cummulative N-deposition per plot -----
# load data
nox <- readRDS("env_data/outputs/EMEP_NOxdeposition.rds")
nhx <- readRDS("env_data/outputs/EMEP_NHxdeposition.rds")
missing_emep <- nox %>% 
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  filter(is.na(value)) %>% 
  pull(plot_ID) %>% 
  unique
ack <- readRDS("env_data/outputs/Ackerman_ndeposition.rds") %>% 
  # only plots not in emep data
  filter(plot_ID %in% missing_emep) %>% 
  mutate(obs_year = NA)
meta <- read.csv("Z:/_GLOBNUT1.0/GlobNut1.0_metadata.csv") %>% 
  # for the plots sampled in 1989, take the ndep observation from 1990
  mutate(year = ifelse(year == 1989, 1990, year))

## 5-year sum
## imputing missing years with last observation (i.e plots sampled between 1990 - 1994) or only observation of 2016 in case of ackerman plots
# create dataframe with years needed for the calculation of n depositoin
yr5 <- meta[, c("plot_ID", "year")] %>% 
  # calculate mean for 5 years prior sampling
  mutate(start_year = year - 4, end_year = year, year = NULL) %>% 
  na.omit()
l <- list() # surely this can be done nicer than a loop?
for (i in 1:nrow(yr5)) {
  t <- expand.grid(plot_ID = yr5[i,]$plot_ID, obs_year = yr5[i,]$start_year:yr5[i,]$end_year)
  l[[i]] <- t
}
yr5 <- bind_rows(l)

## combining datasets of emep and ackerman
ndep <- bind_rows(nhx, nox) %>% 
  filter(!is.na(value)) %>% 
  bind_rows(ack) %>% #add ackerman data
  left_join(meta[,c("plot_ID","year")], by = "plot_ID") %>% # join with meta data for sample year 
  mutate(obs_year = ifelse(is.na(obs_year), year, obs_year)) %>% # set to ackerman obs year to year of sampling
  # calculate total Ndep (Nox+nhx) by plot
  group_by(plot_ID, obs_year) %>% 
  summarise(ndep = sum(value, na.rm = TRUE) %>% na_if(0)) 

ndep_5yr_imp <- yr5 %>% # years needed
  # available data
  left_join(ndep, by = c("plot_ID", "obs_year")) %>%  
  # impute missing 
  group_by(plot_ID) %>% 
  mutate(n_obs = sum(!is.na(ndep))) %>% # number of observations
  fill(ndep, .direction = "updown") %>% 
  # calculate mean and sum of ndep
  summarise(mean_5yr = mean(ndep),
            sd_5yr = sd(ndep),
            sum_5yr = sum(ndep),
            n_obs_5yr = unique(n_obs)) 

## 10 year sum
yr10 <- meta[, c("plot_ID", "year")] %>% 
  # calculate mean for 10 years prior sampling
  mutate(start_year = year - 9, end_year = year, year = NULL) %>% 
  na.omit()

l <- list() # surely this can be done nicer than a loop?
for (i in 1:nrow(yr10)) {
  t <- expand.grid(plot_ID = yr10[i,]$plot_ID, obs_year = yr10[i,]$start_year:yr10[i,]$end_year)
  l[[i]] <- t
}
yr10 <- bind_rows(l)

# 10 year sum and average
ndep_10yr_imp <- yr10 %>% # years needed
  # available data
  left_join(ndep, by = c("plot_ID", "obs_year")) %>%  
  # impute missing 
  group_by(plot_ID) %>% 
  mutate(n_obs = sum(!is.na(ndep))) %>% # number of observations
  fill(ndep, .direction = "updown") %>% 
  summarise(mean_10yr = mean(ndep),
            sd_10yr = sd(ndep),
            sum_10yr = sum(ndep),
            n_obs_10yr = unique(n_obs)) 

ndep <- left_join(ndep_5yr_imp, ndep_10yr_imp, by = "plot_ID")

# write to csv for globnut 1.0
write.csv(ndep, "Z:/_GLOBNUT1.0/ndeposition.csv", row.names = FALSE,  fileEncoding = "UTF-8")


## Base map for fig 1
# EMEP
dnox_2020 <- raster(nox_raster$X2020)
dnhx_2020 <- raster(nhx_raster$X2020)
wnox_2020 <- raster(wnox_raster$X2020)
wnhx_2020 <- raster(wnhx_raster$X2020)
# cumulative
ndep_emep <- calc(stack(dnox_2020, dnhx_2020,wnox_2020,wnhx_2020), sum) * 0.01  # tp g/ha
# ndep_emep <- raster::aggregate(ndep_emep, fact =5, fun = mean)
plot(ndep_emep)
# global
ndep_global <- calc(stack(raster(glob_nhx), raster(glob_nox)), sum ) * 0.01 # to g/ha
plot(ndep_global)

# combine emep and ackerman

# crop and resample ackerman to resolution of emep
ras_extent <- raster(xmn = 90, xmx = 150, ymn = 30, ymx = 82, res = raster::res(ndep_emep))
resample_ndep <- raster::resample(ndep_global, ras_extent, method = "bilinear") 
ndep_comb <- raster::merge(ndep_emep, resample_ndep)
raster::writeRaster(ndep_comb, "figures/Files_Ton/ndeposition.tif",format = "GTiff", overwrite = T)


# Euroasia shape file
sf::sf_use_s2(FALSE)
europe <- ne_countries(scale = "medium",
                       # country = country,
                       returnclass = "sf") %>%
  # fix geometry of russia (sf_use_s2() must be FALSE for this)
  st_make_valid() %>%
  st_crop(., c(xmin = -25, ymin = 35, xmax = 150, ymax = 72))

ndep  <- mask(ndep_comb, europe)

r_df <- as.data.frame(ndep, xy = TRUE)
colnames(r_df) <- c("x", "y", "value")

r_df <- r_df %>%
  mutate(category = cut(value, breaks = c(-Inf, 5, 10, 20, 40, Inf),
                        labels = c("< 5", "5-10", "10-20", "20-40", "> 40")))

ggplot(r_df, aes(x = x, y = y, fill = category)) +
  geom_raster() +
  scale_fill_manual(values = c("< 5" = "grey90",
                               "5-10" = "grey70",
                               "10-20" = "grey50",
                               "20-40" = "grey30",
                               "> 40" = "grey10"),na.translate = FALSE) +
  labs(fill = "kg N/ha/yr") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave("figures/Files_Ton/Base_map_ndep.pdf", dpi = 300)
