## ---------------------------
##
## Script name: 03_geography_limitation 
##
## Purpose of script: plot the distirbution of sites and their limitation (fig1)
## Plot species richness - biomass by limitation (Supplementary)
##
## Author: Annegreet Veeken
##
## Date Created: 2024-02-19
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
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(easystats)) install.packages("easystats")
if (!require(rnaturalearth)) install.packages("rnaturalearth")
if (!require(rnaturalearthdata)) install.packages("rnaturalearthdata")
if (!require(patchwork)) install.packages("patchwork")
if (!require(sf)) install.packages("sf")
if (!require(ggExtra)) install.packages("ggExtra")
if (!require(gridExtra)) install.packages("gridEXtra")
if (!require(quantreg)) install.packages("quantreg")
if (!require(mgcv)) install.packages("mgcv")

## Load data
globnut <- readRDS("outputs/02_GlobNut.rds") %>% 
  mutate(biomass_cent = biomass - mean(biomass)) 

# Euroasia shape file
world <- ne_countries(scale = "medium", returnclass = "sf")
eurasia <- world[world$continent %in% c("Asia", "Europe"),]

ecoreg_globnut <- readRDS("outputs/Ecoregions.rds")

# add to globnut
globnut <- globnut %>% 
  left_join(ecoreg_globnut[,c("plot_ID", "Ecoreg.ECO_NAME")], by = "plot_ID") %>% 
  rename(ecoregion = Ecoreg.ECO_NAME)

region_alias <- data.frame(ecoregion = c("Alps conifer and mixed forests",
                                      "Central European mixed forests",
                                      "East Siberian taiga","European Atlantic mixed forests",
                                      "Kazakh steppe","Pannonian mixed forests",
                                      "Western European broadleaf forests"),
                           region = c("alps", "north_eu", "jak", "west_eu", "kaz", 
                                     "east_eu", "central_eu"),
                           label = c("Alps", "Northeastern Europe", "Yakutsk", "Western Europe",
                                     "Kazakhstan", "Eastern Europe", "Central Europe"))

# select plots for vizualization
globnut <- globnut %>% 
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  left_join(region_alias, by = "ecoregion")
  
# subset for TOn 
sub_ton <- globnut %>% 
  dplyr::select(plot_ID, lat,lon, limitation = lim, region)
write.csv(sub_ton, "figures/Files_Ton/Latlon_studysites.csv", row.names = F)

# labels
region_labels <- globnut %>% 
  group_by(region) %>% 
  summarise(totplots = n()) %>% 
  left_join(region_alias, by = "region") %>% 
  mutate(plot_lab = paste0(label, " (N = ", totplots, ")"))
xlabs <- c(0, 50, 100, 150)
ylabs <- seq(35, 70, by = 10)

# colors
colours_lim <- c("Co-limitation N-P" = '#1b9e77',"N-limitation" = '#d95f02',
                 "P-limitation" = '#7570b3')

region_color <- c(jak = '#a65628',
                  north_eu = '#377eb8',
                  kaz = '#4daf4a',
                  west_eu = '#984ea3',
                  central_eu = '#ff7f00',
                  east_eu = '#e41a1c',
                  alps = '#f781bf')

# barplot per region as annotation
lim_region <- globnut %>% 
  group_by(region, lim) %>% 
  summarise(n_lim = n()) %>% 
  group_by(region) %>% 
  mutate(percent_lim = n_lim/sum(n_lim)*100)

p_west <- lim_region %>%
  filter(region == "west_eu") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "west_eu"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_cent <- lim_region %>%
  filter(region == "central_eu") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "central_eu"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_north <- lim_region %>%
  filter(region == "north_eu") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "north_eu"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_east <- lim_region %>%
  filter(region == "east_eu") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "east_eu"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_jak <- lim_region %>%
  filter(region == "jak") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "jak"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_kaz <- lim_region %>%
  filter(region == "kaz") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "kaz"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

p_alps <- lim_region %>%
  filter(region == "alps") %>%
  ggplot(aes(x = percent_lim, fill = lim, y = region)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    alpha = 0.5
  ) +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title =  region_labels$plot_lab[region_labels$region == "alps"]) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.ticks.length = unit(0, "pt"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 6)
  )

# convert points and annotation to sf object
regions <- data.frame(region = c("west_eu", "north_eu","central_eu",
                                 "alps", "east_eu", "jak", "kaz"),
                      xmin = c(2, 10, 8, 6, 16, 126, 65),
                      ymin = c(50, 49, 47, 45, 46.5, 60, 48),
                      xmax = c(8, 26, 17, 14, 24, 133, 80),
                      ymax = c(54, 54.5, 50, 47, 50, 64, 56)) 
write_csv(regions,"figures/Files_Ton/Regions_Globnut_map.csv")
coords_sf <- st_as_sf(globnut, coords = c("lon", "lat"), crs = 4326)
rect_sf <- regions %>% 
  rowwise() %>%
  mutate(geometry = list(st_polygon(list(rbind(
    c(lon = xmin, lat = ymin), 
    c(lon = xmin, lat = ymax), 
    c(lon = xmax, lat = ymax), 
    c(lon = xmax, lat = ymin), 
    c(lon = xmin, lat = ymin)
  ))))) %>% 
  st_as_sf(crs = 4326)

annon_arrow <- data.frame(x = c(0.15,0.31),
                          xend = c(0.235,0.31),
                          y = c(0.32,0.22),
                          yend = c(0.32,0.5))

p <- ggplot() +
  geom_sf(data = eurasia, fill = "grey90", col = "grey90") +
  geom_sf(data = coords_sf, aes(col = lim), alpha = 0.25, size = 0.5) +
  geom_sf(data = rect_sf, fill = NA) +
  scale_color_manual(values = colours_lim) +
  coord_sf(crs = "+proj=lcc +lat_1=45 +lat_2=55 +lat_0=50 +lon_0=35",
           xlim = c(-4000000, 4500000), ylim = c(-1000000, 4500000)) +
  theme_bw() +
  theme(legend.position = c(0.15,0.8),
        legend.margin = margin(t = 0, r = 3, b = 0, l = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
       legend.spacing.y = unit(0.1, 'cm'), # brings legend items closer vertically
       legend.key = element_blank(), # removes the white boxes in the legend keys
       legend.background = element_blank(),
       axis.title = 
       ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   size = 3)))
p <- p +
  annotation_custom(grob = grid::segmentsGrob(
                    x0 = unit(annon_arrow$x, "npc"), 
                    y0 = unit(annon_arrow$y, "npc"), 
                    x1 = unit(annon_arrow$xend, "npc"), 
                    y1 = unit(annon_arrow$yend, "npc"))) 

(p_map <- p + 
    inset_element(p_cent, left = 0.25, bottom = 0.5, right = 0.35, top = 0.57, ignore_tag = TRUE) +
    inset_element(p_west, left = 0.05, bottom = 0.3, right = 0.15, top = 0.37, ignore_tag = TRUE) +
    inset_element(p_east, left = 0.39, bottom = 0.14, right = 0.49, top = 0.21, ignore_tag = TRUE) +
    inset_element(p_north, left = 0.35, bottom = 0.33, right = 0.45, top = 0.4, ignore_tag = TRUE) +
    inset_element(p_kaz, left = 0.7, bottom = 0.2, right = 0.8, top = 0.27, ignore_tag = TRUE) +
    inset_element(p_alps, left = 0.15, bottom = 0.14, right = 0.25, top = 0.21, ignore_tag = TRUE) +
    inset_element(p_jak, left = 0.8, bottom = 0.8, right = 0.9, top = 0.87, ignore_tag = TRUE))

## Proportion Limitation by long/latitude bin ----
prop_lon <- globnut %>% 
  mutate(lon_bin = cut(lon, breaks = seq(-10, 150, by = 1),labels = seq(-9.5, 149.5, by = 1))) %>% 
  group_by(lon_bin,lim) %>% 
  summarise(n = n()) %>% 
  group_by(lon_bin) %>% 
  mutate(prop = n/sum(n),
         tot_plots = sum(n)) %>%  
  filter(tot_plots >= 10) %>% 
  mutate(lon_bin = as.numeric(as.character(lon_bin)),
         totprop = sum(prop),
         lim = case_match(lim,"N-limitation"~ "nlim",
                          "P-limitation" ~ "plim",
                          "Co-limitation N-P" ~ "colim")) 

lon_nlim <- gam(prop ~ s(lon_bin), data = prop_lon[prop_lon$lim == "nlim",], family = betar(link = "logit"))
lon_plim <- gam(prop ~ s(lon_bin), data = prop_lon[prop_lon$lim == "plim",], family = betar(link = "logit"))
lon_colim <- gam(prop ~ s(lon_bin), data = prop_lon[prop_lon$lim == "colim",], family = betar(link = "logit"))

pred <- data.frame(lon_bin =  seq(-9.5, 149.5, by = 1))
pred_lon_nlim <- predict.gam(lon_nlim, newdata = pred, se.fit = T)
pred_lon_plim <- predict.gam(lon_plim, newdata = pred, se.fit = T)
pred_lon_colim <- predict.gam(lon_colim, newdata = pred, se.fit = T)

pred$nlim_mean <- plogis(pred_lon_nlim$fit)
pred$nlim_up <- plogis(pred_lon_nlim$fit + qnorm(0.975) * pred_lon_nlim$se.fit)
pred$nlim_low <- plogis(pred_lon_nlim$fit - qnorm(0.975) * pred_lon_nlim$se.fit)
pred$plim_mean <- plogis(pred_lon_plim$fit)
pred$plim_up <- plogis(pred_lon_plim$fit + qnorm(0.975) * pred_lon_plim$se.fit)
pred$plim_low <- plogis(pred_lon_plim$fit - qnorm(0.975) * pred_lon_plim$se.fit)
pred$colim_mean <- plogis(pred_lon_colim$fit)
pred$colim_up <- plogis(pred_lon_colim$fit + qnorm(0.975) * pred_lon_colim$se.fit)
pred$colim_low <- plogis(pred_lon_colim$fit - qnorm(0.975) * pred_lon_colim$se.fit)

p_lon <- pred %>% 
  pivot_longer(cols = -lon_bin, names_to = c("lim","est"), names_sep = "_", values_to =  "value") %>% 
  pivot_wider(names_from = "est", values_from = "value") %>% 
  ggplot(aes(x = lon_bin, y = mean, fill = lim, color = lim)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.5, color = NA) +
  # geom_point(data = prop_lon, aes(x = lon_bin, y = c(mean = prop))) +
  scale_color_manual(values =  c(colim = '#1b9e77', nlim = '#d95f02', plim = '#7570b3')) +
  scale_fill_manual(values =  c(colim = '#1b9e77' , nlim = '#d95f02', plim = '#7570b3')) +
  scale_x_continuous("Longitude", breaks = xlabs, labels = paste0(xlabs,'°E')) +
  scale_y_continuous("% plots with limitation", labels = scales::percent, limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none")

# latitude
prop_lat <- globnut %>% 
  mutate(lat_bin = cut(lat, breaks = seq(35, 70, by = 1),
                       labels = seq(35.5,69.5, by = 1))) %>% 
  group_by(lat_bin,lim) %>% 
  summarise(n = n()) %>% 
  group_by(lat_bin) %>% 
  mutate(prop = n/sum(n),
         tot_plots = sum(n)) %>%  
  filter(tot_plots >= 10) %>% 
  mutate(lat_bin = as.numeric(as.character(lat_bin)),
         totprop = sum(prop)) 

lat_nlim <- gam(prop ~ s(lat_bin), data = prop_lat[prop_lat$lim == "N-limitation",], family = betar(link = "logit"))
lat_plim <- gam(prop ~ s(lat_bin), data = prop_lat[prop_lat$lim == "P-limitation",], family = betar(link = "logit"))
lat_colim <- gam(prop ~ s(lat_bin), data = prop_lat[prop_lat$lim == "Co-limitation N-P",], family = betar(link = "logit"))

pred_lat <- data.frame(lat_bin =  seq(35.5, 69.5, by = 1))
pred_lat_nlim <- predict.gam(lat_nlim, newdata = pred_lat, se.fit = T)
pred_lat_plim <- predict.gam(lat_plim, newdata = pred_lat, se.fit = T)
pred_lat_colim <- predict.gam(lat_colim, newdata = pred_lat, se.fit = T)

pred_lat$nlim_mean <- plogis(pred_lat_nlim$fit)
pred_lat$nlim_up <- plogis(pred_lat_nlim$fit + qnorm(0.975) * pred_lat_nlim$se.fit)
pred_lat$nlim_low <- plogis(pred_lat_nlim$fit - qnorm(0.975) * pred_lat_nlim$se.fit)
pred_lat$plim_mean <- plogis(pred_lat_plim$fit)
pred_lat$plim_up <- plogis(pred_lat_plim$fit + qnorm(0.975) * pred_lat_plim$se.fit)
pred_lat$plim_low <- plogis(pred_lat_plim$fit - qnorm(0.975) * pred_lat_plim$se.fit)
pred_lat$colim_mean <- plogis(pred_lat_colim$fit)
pred_lat$colim_up <- plogis(pred_lat_colim$fit + qnorm(0.975) * pred_lat_colim$se.fit)
pred_lat$colim_low <- plogis(pred_lat_colim$fit - qnorm(0.975) * pred_lat_colim$se.fit)

p_lat <- pred_lat %>% 
  pivot_longer(cols = -lat_bin, names_to = c("lim","est"), names_sep = "_", values_to =  "value") %>% 
  pivot_wider(names_from = "est", values_from = "value") %>% 
  ggplot(aes(x = lat_bin, y = mean, fill = lim, color = lim)) +
  geom_line() +
  geom_ribbon(aes(ymin = low, ymax = up), color = NA, alpha = 0.5) +
  scale_color_manual(values =  c(colim = '#1b9e77', nlim = '#d95f02', plim = '#7570b3')) +
  scale_fill_manual(values =  c(colim = '#1b9e77', nlim = '#d95f02', plim = '#7570b3')) +
  scale_x_continuous("Latitude", breaks = ylabs, labels = paste0(ylabs,'°N')) +
  scale_y_continuous("", labels = scales::percent, limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none")

## combine plots fig 1 ----
design <- "AAAAAA
           AAAAAA
           BBBCCC"
p_all <- p_map + free(p_lon) + free(p_lat) +
  plot_layout(design = design) + plot_annotation(tag_levels = "a")
p_all

ggsave("figures/fig1_geography_lim_regions_Nlim14.png", p_all, dpi = 300, height = 7, width = 7)
ggsave("figures/Files_Ton/fig1_geography_lim_regions.pdf", p_all, dpi = 300, height = 7, width = 7)

# Biomass - species plot () ----
globnut <- readRDS("outputs/02_GlobNut.rds") %>% 
  ungroup() %>% 
  filter(!lim == "No limitation by N, P, K") %>% 
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  drop_na %>% 
  mutate(biomass_cent = biomass - mean(biomass))

n_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), 
           data = globnut[globnut$lim == "N-limitation",])
co_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size),
            data = globnut[globnut$lim == "Co-limitation N-P",])
p_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), 
           data = globnut[globnut$lim == "P-limitation",])

n_qr <- rq(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "N-limitation",],
            tau = 0.9)
co_qr <- rq(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "Co-limitation N-P",],
           tau = 0.9)
p_qr <- rq(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "P-limitation",],
            tau = 0.9)

df <- data.frame(biomass_cent = seq(-300,1200), plot_size = 4)
pred_qr <- rbind(predict(n_qr, df) %>% 
                   as.data.frame() %>% 
                   mutate(spec_ric = exp(.),  
                          biomass_cent = df$biomass_cent,
                          lim = "N-limitation",
                          biomass = biomass_cent + mean(globnut$biomass)),
                 predict(co_qr, df) %>% 
                   as.data.frame() %>% 
                   mutate(spec_ric = exp(.),  
                          biomass_cent = df$biomass_cent,
                          lim = "Co-limitation N-P",
                          biomass = biomass_cent + mean(globnut$biomass)),
                 predict(p_qr, df) %>% 
                   as.data.frame() %>% 
                   mutate(spec_ric = exp(.),  
                          biomass_cent = df$biomass_cent,
                          lim = "P-limitation",
                          biomass = biomass_cent + mean(globnut$biomass)))
                 

pred_lm <- rbind(ggeffects::ggpredict(n_lm, terms = "biomass_cent [-300:1200]") %>% 
                   as.data.frame() %>% 
        mutate(lim = "N-limitation",
               biomass = x + mean(globnut$biomass),
               spec_ric = predicted),
     ggeffects::ggpredict(co_lm, terms = "biomass_cent [-300:1200]") %>% 
        as.data.frame() %>% 
        mutate(lim = "Co-limitation N-P",
               biomass = x + mean(globnut$biomass),
               spec_ric = predicted),
     ggeffects::ggpredict(p_lm, terms = "biomass_cent [-300:1200]") %>% 
        as.data.frame() %>% 
        mutate(lim = "P-limitation",
               biomass = x + mean(globnut$biomass),
               spec_ric = predicted))
# significance between limitation categories
m <- aov(log(biomass) ~ lim, data = globnut)
plot(m)
TukeyHSD(m)

m <- aov(log(spec_ric) ~ lim, data = globnut)
plot(m)

TukeyHSD(m)

colours_lim <- c('#1b9e77','#d95f02','#7570b3')
g <- ggplot() +
  geom_point(data = globnut, aes(x = biomass, y = spec_ric, col = lim), alpha = 0.25, size = 0.5) +
  geom_line(data = pred_lm, aes(x = biomass, y = spec_ric, col = lim)) + 
  geom_line(data = pred_qr, aes(x = biomass, y = spec_ric, col = lim), linetype = "dashed") +
  scale_color_manual(values = colours_lim) +
  scale_x_continuous(Biomass~(g/m^2), limits = c(0,1600)) +
  scale_y_continuous("Species richness", limits = c(0, 150)) +
  theme_bw() +
  theme(legend.position = c(0.8,0.9),
        legend.title = element_blank(),
        legend.background = element_blank())
g
g_mar <- ggMarginal(g, type = "boxplot", size = 5, groupFill = TRUE, groupColour = TRUE)

ggsave("figures/Files_Ton/Sup_bio_specRic_lim.pdf", g_mar, 
       dpi = 300, height = 5, width = 5)

ggsave("figures/Sup_bio_specRic_limN14.png",g_mar, dpi = 300, height = 5, width = 5)
