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
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rnaturalearth)) install.packages("rnaturalearth")
if(!require(rnaturalearthdata)) install.packages("rnaturalearthdata")
if(!require(sf)) install.packages("sf")
if(!require(ggExtra)) install.packages("ggExtra")
if(!require(gridExtra)) install.packages("gridEXtra")
if(!require(quantreg)) install.packages("quantreg")

## Load data
globnut <- readRDS("outputs/02_GlobNut_subsampled.rds") %>% 
  mutate(biomass_cent = biomass - mean(biomass))

# colors
colours_lim <- c('#1b9e77','#d95f02','#7570b3')

# Europe shape fi;e
sf::sf_use_s2(FALSE)
europe <- ne_countries(scale = "medium", 
                       # country = country,
                       returnclass = "sf") %>% 
  # fix geometry of russia (sf_use_s2() must be FALSE for this)
  st_make_valid() %>% 
  # crop to extent of Europe
  st_crop(., c(xmin = -30, ymin = 35, xmax = 150, ymax = 72))

globnut <- globnut %>% 
  dplyr::filter(between(lat, 35, 75)) %>% 
  dplyr::filter(between(lon, -25, 150)) %>% 
  filter(!lim %in% c("No limitation by N, P, K")) 
sf_use_s2(TRUE)
xlabs <- c(0,50,100,150)
ylabs <- seq(35,70, by = 10)
p <- ggplot() +
  geom_sf(data = europe, fill = "grey80", col = "white") +
  geom_point(data = globnut, aes(x = lon, y = lat, fill = lim, col = lim), alpha = 0.25, size =1) +
  scale_color_manual(values = colours_lim) +
  coord_sf(datum = st_crs(4326)) +
  scale_x_continuous("Longitude", breaks = xlabs, labels = paste0(xlabs,'°E')) +
  scale_y_continuous("Latitude", breaks = ylabs, labels = paste0(ylabs,'°N')) +
  theme_bw() +
  theme(legend.position = c(0.1,0.8),
        legend.margin = margin(t = 0, r = 3, b = 0, l = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
       legend.spacing.y = unit(0.1, 'cm'), # brings legend items closer vertically
       legend.key = element_blank(), # removes the white boxes in the legend keys
       legend.background = element_blank()
       ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_mar <-  ggMarginal(p, type = "histogram", size = 10, binwidth = 5, groupFill = TRUE, groupColour = TRUE)
p_mar


ggsave("figures/fig1_geography_lim_subsamples.png", p_mar, dpi = 300, height = 7, width = 7)


# Biomass - species plot ()
n_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "N-limitation",])
co_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "Co-limitation N-P",])
p_lm <- lm(log(spec_ric) ~ biomass_cent + I(biomass_cent^2) + log(plot_size), data = globnut[globnut$lim == "P-limitation",])

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
                 

pred_lm <- rbind(ggeffect(n_lm, terms = "biomass_cent [-300:1200]") %>% 
                   as.data.frame() %>% 
        mutate(lim = "N-limitation",
               biomass = x + mean(globnut$biomass),
               spec_ric = exp(predicted)),
      ggeffect(co_lm, terms = "biomass_cent [-300:1200]") %>% 
        as.data.frame() %>% 
        mutate(lim = "Co-limitation N-P",
               biomass = x + mean(globnut$biomass),
               spec_ric = exp(predicted)),
      ggeffect(p_lm, terms = "biomass_cent [-300:1200]") %>% 
        as.data.frame() %>% 
        mutate(lim = "P-limitation",
               biomass = x + mean(globnut$biomass),
               spec_ric = exp(predicted)))
# significance between limitation categories
m <- aov(log(biomass) ~ lim, data = globnut)
plot(m)
TukeyHSD(m)

m <- aov(log(spec_ric) ~ lim, data = globnut)
plot(m)
TukeyHSD(m)


g <- ggplot() +
  geom_point(data = globnut, aes(x = biomass, y = spec_ric, col = lim), alpha = 0.25, size = 0.5) +
  geom_line(data = pred_lm, aes(x = biomass, y = spec_ric, col = lim)) + 
  geom_line(data = pred_qr, aes(x = biomass, y = spec_ric, col = lim), linetype = "dashed") +
  scale_color_manual(values = colours_lim) +
  scale_x_continuous(Biomass~(g/m^2), limits = c(0,1600)) +
  scale_y_continuous("Species richness", limits = c(0, 150)) +
  geom_segment(aes(ymin = ))
  theme_bw() +
  theme(legend.position = c(0.8,0.9),
        legend.title = element_blank(),
        legend.background = element_blank())
g
g_mar <- ggMarginal(g, type = "boxplot", size = 5, groupFill = TRUE, groupColour = TRUE)


ggsave("figures/Sup_bio_specRic_lim.png",g_mar, dpi = 300, height = 5, width = 5)

p2 <- globnut %>% 
  mutate(lat_bin = cut(lat, breaks = seq(30,80,by = 5)),
         long_bin = cut(lon, breaks = seq(0,150, by = 10))) %>% 
  group_by(lat_bin, lim) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = lat_bin, y = n, fill = lim)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours_lim) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")
p3 <- globnut %>% 
  mutate(lat_bin = cut(lat, breaks = seq(30,80,by = 5)),
         long_bin = cut(lon, breaks = seq(0,150, by = 10))) %>% 
  group_by(long_bin, lim) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = long_bin, y = n, fill = lim)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours_lim) +
  theme_void() +
  theme(legend.position = "none")

grid.arrange(p,p3, p2, layout_matrix = rbind(c(NA,2,2,NA),
                                             c(NA,1,1,NA),
                                              c(NA,1,1,3)),
             heights = c(1,1,3))
