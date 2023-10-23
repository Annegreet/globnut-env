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
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(vegan)) install.packages("vegan")
if(!require(readxl)) install.packages("readxl")
if(!require(TNRS)) install.packages("TNRS")

## Load data
data_dir <- "Z:/_GLOBNUT1.0/" # directory with Globnut 1.0 data
spec_raw <- read.csv(paste0(data_dir, "GlobNut1.0_species.csv"))
nfix <- read_xlsx("Nodb.xlsx", skip = 1, sheet = 1) %>% 
  drop_na(genus)

## Data wrangling
spec <- spec_raw %>% 
  # only vascular plants
  dplyr::filter(vascular_plant == 1) %>% 
  # exclude unknowns
  # dplyr::filter(!str_detect(species_new, pattern = "Unknown|unknown|unkown")) %>% 
  # convert subspecies to species
  mutate(species_new = str_remove(species_new, pattern = " subsp.*| var\\.\\b.")) %>%
  group_by(plot_ID, family_new, species_new) %>% # aggregate by species
  summarise(cover = sum(cover, na.rm = TRUE)) %>% 
  # remove zero abundances
  dplyr::filter(cover != 0)  %>% 
  # add genus column 
  mutate(genus = word(species_new, 1))

## Species richness and evenness
spec_ric <- spec %>% 
  # calculate species richness by plot 
  group_by(plot_ID) %>% 
  summarise(spec_ric = n_distinct(species_new)) 

spec_wide <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  dplyr::select(species_new, cover) %>% 
  # convert to wide to calculate diversity measures
  pivot_wider(names_from = species_new, values_from = cover, values_fill = 0) 
spec_ric$shan <- diversity(spec_wide[,-1]) # Shannon entropy
spec_ric$shan_eve <- exp(spec_ric$shan) / spec_ric$spec_ric # shannon evenness
spec_ric$pilou_eve <- spec_ric$shan / log(spec_ric$spec_ric)

## Percentage grass
grass_globnut <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  # create grouping column for grass
  mutate(grass = ifelse(family_new == "Poaceae", "Grass", "Not grass") %>% as.factor(.)) %>% 
  group_by(plot_ID, grass, .drop = FALSE) %>% 
  summarise(grass_cover = sum(cover)) %>% 
  filter(grass == "Grass") %>% 
  select(plot_ID, grass_cover)

## Percentage Nfix
# harmonize genus names the same way as globnut
tax_check <- TNRS(nfix$genus)
nfix_harm <- nfix %>% 
  left_join(tax_check[,c("Name_submitted", "Accepted_name")], by = c("genus" = "Name_submitted")) %>% 
  # harmonisation merges genera, creating duplicates for consensus column, check if conflicts arise
  group_by(Accepted_name) %>% 
  summarise(consensus = paste(unique(`Consensus estimate`), collapse = ", "),
            n = n_distinct(`Consensus estimate`)) %>%  #used to check for conflicts
  # Zygia has conflict between consensus estimates, chose likely_Rhizobia instead of None, n studies is highest for this
  mutate(consensus = case_match(consensus, "Rhizobia, None"~ "likely_Rhizobia",
                                "Rhizobia, likely_Rhizobia" ~ "likely_Rhizobia",
                                "None, unlikely_Rhizobia" ~ "unlikely_Rhizobia",
                                "likely_present, Present" ~ "likely_present",
                                "likely_Rhizobia, Rhizobia" ~ "likely_Rhizobia",
                                .default = consensus)) %>% 
  select(Accepted_name, consensus)

nfix_globnut <- spec %>% 
  #  to 100 %
  group_by(plot_ID) %>% 
  mutate(covertot = sum(cover), cover = cover/covertot) %>% 
  # join Nfix data
  left_join(nfix_harm, by = c("genus" = "Accepted_name")) %>%
  # create one column assigning genera to nodulated or not nodulated
  mutate(nfix = case_when(consensus %in% c("None", "unlikely_Frankia", "unlikely_Rhizobia") ~ "Non-nodulated",
                          is.na(consensus) ~ "Non-nodulated", # genera not in NoddB are likely not nfixers
                          consensus %in% c("Frankia", "likely_present", "likely_Rhizobia",
                                                      "Nostocaceae", "Present", "Rhizobia") ~ "Nodulated") %>% 
           as.factor(.)) %>% 
  group_by(plot_ID, nfix, .drop = FALSE) %>% 
  summarise(nfix_cover = sum(cover)) %>%  
  filter(nfix == "Nodulated") %>% 
  select(plot_ID, nfix_cover)

spec_ric <- spec_ric %>% 
  left_join(grass_globnut, by = "plot_ID") %>%
  left_join(nfix_globnut, by = "plot_ID")

saveRDS(spec_ric, "Outputs/01_Species_indices.rds")

