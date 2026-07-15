## ---------------------------
##
## Script name: 03_geography_limitation 
##
## Purpose of script: Produces figure 1a and b, figure S1a and b, and supplementary figures S3 and S5
##
## Author: Annegreet Veeken
##
## Date Created: 2024-02-19
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes: Figures in the manuscript have been post-processed in Adobe Illustrator 
##  so this script does not produce figures exactly as they appear in the manuscript
##
## ---------------------------


## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rnaturalearth)) install.packages("rnaturalearth")
if (!require(rnaturalearthdata)) install.packages("rnaturalearthdata")
if (!require(patchwork)) install.packages("patchwork")
if (!require(sf)) install.packages("sf")
if (!require(ggExtra)) install.packages("ggExtra")
if (!require(terra)) install.packages("terra")
if (!require(tidyterra)) install.packages("tidyterra")
if (!require(RColorBrewer)) install.packages("RcolorBrewer")

## Load data
lim_order <- c("N-limitation" ,
               "Co-limitation N-P","P-limitation" )
globnut <- readRDS("outputs/02_GlobNut.rds") %>% 
  mutate(lim = factor(lim, levels = lim_order),
         ndep = log(ndep/5))  %>% # to g/ha/yr 
  arrange(lim)  
  
# ndep
world_vec <- ne_countries(scale = "medium", returnclass = "sv")
eurasia_vec <- world_vec[world_vec$continent %in% c("Asia", "Europe"),]
ndep <- rast("figures/Files_Ton/ndeposition.tif") %>% 
  mutate(category = cut(ndeposition, breaks = c(-Inf,5,10,15,20,30,Inf),
                        labels = c("< 5","5-10", "10-15", "15-20","20-30",">30"))) 
# remove values in the sea
ndep <- mask(ndep, eurasia_vec)

# Euroasia shape file
world <- ne_countries(scale = "medium", returnclass = "sf")
eurasia <- world[world$continent %in% c("Asia", "Europe"),]

xlabs <- c(0, 50, 100, 150)
ylabs <- seq(35, 70, by = 10)

# colors
colours_lim <- c("N-limitation" = '#d95f02',"Co-limitation N-P" = '#1b9e77',
                 "P-limitation" = '#7570b3')
coords_sf <- st_as_sf(globnut, coords = c("lon", "lat"), crs = 4326)

# Figure 1 with critical ratio of 13.5 for co-limitation ----
# fig 1a - Map with marginal histograms
p <- ggplot(ndep, aes(x = x, y = y, fill = category)) +
  geom_raster() +
  scale_fill_grey(Annual~N~deposition~(kg~ha^-1)(2020), start = 0.9, end = 0.3, 
                  na.value = "transparent", na.translate = F) +
  geom_point(data = globnut, aes(x = lon, y = lat, color = lim, shape = lim), 
             size = 0.05, alpha = 0.5, inherit.aes = FALSE) +
  scale_shape_manual(name = " ", values = c(3, 5, 17)) +
  scale_color_manual(name = "", values = colours_lim) +
  scale_x_continuous(name = "", limits = c(-10, 150), expand = c(0, 0), 
                     breaks = c(0,50,100,150)) +
  scale_y_continuous(name = "", limits = c(35, 75), expand = c(0, 0),
                     breaks = c(40,50,60,70)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title.position = "bottom",
        plot.margin = margin(0, 0, 0, 0)) +
  guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1, color = colours_lim), 
                              direction = "vertical", order = 1),
         color = guide_legend(override.aes = list(shape = 15, size = 4, alpha = 1), 
                              direction = "vertical", order = 2))
p
# Top marginal - with y-axis
p_top <- ggplot(globnut, aes(x = lon, fill = lim)) +
  geom_histogram(position = "stack", binwidth = 10, boundary = -10, closed = "left") +
  scale_fill_manual("",values = colours_lim) +
  scale_x_continuous(limits = c(-10, 150), expand = c(0, 0),
                     breaks = seq(-10, 150, by = 10)) +
  scale_y_continuous("",expand = expansion(mult = c(0, 0.1)), breaks = seq(0,450, by =150)) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_line())

# Right marginal - with y-axis (which becomes x-axis after coord_flip)
p_right <- ggplot(globnut, aes(x = lat, fill = lim)) +
  geom_histogram(position = "stack", binwidth = 5, boundary = 40, closed = "left") +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(limits = c(35, 75), expand = c(0, 0),
                     breaks = seq(40, 70, by = 5)) +
  scale_y_continuous("",expand = expansion(mult = c(0, 0.1)), breaks =  seq(0,900, by =150)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = -90),
        axis.ticks.x = element_line()) 

# Combine
p_mag <- p_top + plot_spacer() + p + p_right + 
  plot_layout(ncol = 2, nrow = 2, 
              widths = c(6, 1), 
              heights = c(1, 6))
ggsave("figures/202511-Figures-resubmission/PDFs-for-Ton/fig1a-NP13.5.pdf", p_mag, dpi = 300, height = 5, width = 7)

## fig 1b - Proportion Limitation type by n-deposition category ----
lim_ndep <- globnut %>% 
  mutate(ndep_cat = cut(exp(ndep), breaks = c(-Inf,5,10,15,20,30,Inf),
                        labels = c("< 5","5-10", "10-15", "15-20","20-30",">30")),
         lim = factor(lim, levels = c("N-limitation", "Co-limitation N-P", "P-limitation"))) %>% 
  group_by(ndep_cat, lim) %>% 
  summarise(n = n())

ggplot(lim_ndep, aes(x = ndep_cat, y = n, fill = lim)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colours_lim) +
  theme_bw() +
  xlab(Average~annual~N~deposition~(kg~ha^-1)) +
  scale_y_continuous("Number of plots with limitation") +
  theme(legend.title = element_blank())
ggsave("figures/202511-Figures-resubmission/PDFs-for-Ton/fig1b-13.5.pdf", dpi = 300, height = 4, width = 6.5)

## Figure 1 but with critical ratios and concentrations of Palpurina ----
rm("globnut", "p_mag", "p_right", "p_top","lim_ndep")

globnut2 <- readRDS("outputs/02_GlobNut_critical_ratios_palpurina.rds") %>% 
  mutate(lim = factor(lim2, levels = lim_order),
         ndep = log(ndep/5))  %>% # to g/ha/yr 
  arrange(lim)  

p2 <- ggplot(ndep, aes(x = x, y = y, fill = category)) +
  geom_raster() +
  scale_fill_grey(Annual~N~deposition~(kg~ha^-1)(2020), start = 0.9, end = 0.3, 
                  na.value = "transparent",na.translate = F) +
  geom_point(data = globnut2, aes(x = lon, y = lat, color = lim, shape = lim), 
             size = 0.05, alpha = 0.5, inherit.aes = FALSE) +
  scale_shape_manual(name = " ", values = c(3, 5, 17)) +
  scale_color_manual(name = "", values = colours_lim) +
  scale_x_continuous(name = "", limits = c(-10, 150), expand = c(0, 0), 
                     breaks = c(0,50,100,150)) +
  scale_y_continuous(name = "", limits = c(35, 75), expand = c(0, 0),
                     breaks = c(40,50,60,70)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title.position = "bottom",
        plot.margin = margin(0, 0, 0, 0)) +
  guides(shape = guide_legend(override.aes = list(size = 4, alpha = 1, color = colours_lim), 
                              direction = "vertical", order = 1),
         color = guide_legend(override.aes = list(shape = 15, size = 4, alpha = 1), 
                              direction = "vertical", order = 2))
p2
# Top marginal - with y-axis
p_top2 <- ggplot(globnut2, aes(x = lon, fill = lim)) +
  geom_histogram(position = "stack", binwidth = 10, boundary = -10, closed = "left") +
  scale_fill_manual("",values = colours_lim) +
  scale_x_continuous(limits = c(-10, 150), expand = c(0, 0),
                     breaks = seq(-10, 150, by = 10)) +
  scale_y_continuous("",expand = expansion(mult = c(0, 0.1)),breaks =  seq(0,450, by =150)) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_line())

# Right marginal - with y-axis (which becomes x-axis after coord_flip)
p_right2 <- ggplot(globnut2, aes(x = lat, fill = lim)) +
  geom_histogram(position = "stack", binwidth = 5, boundary = 40, closed = "left") +
  scale_fill_manual(values = colours_lim) +
  scale_x_continuous(limits = c(35, 75), expand = c(0, 0),
                     breaks = seq(40, 70, by = 5)) +
  scale_y_continuous("",expand = expansion(mult = c(0, 0.1)),breaks =  seq(0,900, by =150)) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = -90),
        axis.ticks.x = element_line()) 

# Combine
p_mag2 <- p_top2 + plot_spacer() + p2 + p_right2 + 
  plot_layout(ncol = 2, nrow = 2, 
              widths = c(6, 1), 
              heights = c(1, 6))
ggsave("figures/202511-Figures-resubmission/PDFs-for-Ton/fig1a-NP10.pdf", p_mag2, dpi = 300, height = 5, width = 7)

## fig 1b - Proportion Limitation type by n-deposition category ----
lim_ndep2 <- globnut2 %>% 
  mutate(ndep_cat = cut(exp(ndep), breaks = c(-Inf,5,10,15,20,30,Inf),
                                       labels = c("< 5","5-10", "10-15", "15-20","20-30",">30")),
         lim = factor(lim, levels = c("N-limitation", "Co-limitation N-P", "P-limitation"))) %>% 
  group_by(ndep_cat, lim) %>% 
  summarise(n = n())

P_hist2 <- ggplot(lim_ndep2, aes(x = ndep_cat, y = n, fill = lim)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colours_lim) +
  theme_bw() +
  xlab(Average~annual~N~deposition~(kg~ha^-1)) +
  scale_y_continuous("Number of plots with limitation") +
  theme(legend.title = element_blank())

ggsave("figures/202511-Figures-resubmission/PDFs-for-Ton/fig1b-NP10.pdf", P_hist2, dpi = 300, height = 4, width = 6.5)


# Supplementary figure N & P concentrations ----
if(0){
m1 <-   summary(lm(N ~ ndep, data = globnut))
m2 <- summary(lm(P ~ ndep, data = globnut))
p1 <- lm(N ~ ndep, data =globnut ) %>% 
  ggeffects::ggpredict(terms = "ndep[all]") %>%
  as.data.frame() %>% 
  ggplot(aes(x = exp(x), y = predicted)) +
  geom_line() +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high), fill = "grey70", alpha = 0.5) + 
  annotate("text", x = 5, y = 1.4,
           label = paste0("y ~ ", round(m1$coefficients[1,1],1), " + ", 
                          round(m1$coefficients[2,1],3), "x, p = ", round(m1$coefficients[2,4],2))) +
  scale_x_continuous(Average~annual~N~deposition~(kg~ha^-1),
                     trans = "log",
                     breaks = c(1,2,5,10,25,60),
                     labels = c(1,2,5,10,25,60),
                     limits = c(1,85)) +
  scale_y_continuous(N~`in`~dry~weight~(g~g^-1)) +
  theme_bw()
p1

p2 <- lm(P ~ ndep, data = globnut ) %>% 
  ggeffects::ggpredict(terms = "ndep[all]") %>%
  as.data.frame() %>% 
  ggplot(aes(x = exp(x), y = predicted)) +
  geom_line() +
  geom_ribbon( aes(ymin = conf.low, ymax = conf.high), fill = "grey70", alpha = 0.5) + 
  annotate("text", x = 5, y = 0.18,
           label = paste0("y ~ ", round(m2$coefficients[1,1],1), " + ", round(m2$coefficients[2,1],3), "x, p < 0.001")) +
  scale_x_continuous(Average~annual~N~deposition~(kg~ha^-1),
                     trans = "log",
                     breaks = c(1,2,5,10,25,60),
                     labels = c(1,2,5,10,25,60),
                     limits = c(1,85)) +
  scale_y_continuous("P in dry weight (g/g)") +
  theme_bw()
p1 + p2 + plot_layout(axes = "collect") + plot_annotation(tag_levels = "a")
ggsave("figures/202511-Figures-resubmission/Sup_ndep_N_P.png", width = 7.24, height = 4, dpi = 300)


## pH map supplementary
ggplot() +
  geom_sf(data = eurasia) +
  geom_point(data = globnut[is.na(globnut$pH_field),], aes(x = lon, y = lat), 
             shape = 21, fill = NA, size = 0.5) +
  geom_point(data = globnut[!is.na(globnut$pH_field),], aes(x = lon, y = lat, col = pH_field)) +
  scale_x_continuous(name = "", limits = c(-10, 150), expand = c(0, 0), 
                     breaks = c(0,50,100,150)) +
  scale_y_continuous(name = "", limits = c(35, 75), expand = c(0, 0),
                     breaks = c(40,50,60,70)) +
  scale_color_viridis_c("Field measured pH") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("figures/202511-Figures-resubmission/Sup-pH_map.png", dpi = 300, width = 7.24, height = 6)
}

## Table with % limitation by ndep category 
globnut_pal <- readRDS("outputs/02_GlobNut_critical_ratios_palpurina.rds")
perc_pal <- globnut_pal %>% 
  mutate(ndep_cat = cut(ndep/5, breaks = seq(0, 30, by = 5))) %>% 
  group_by(ndep_cat,lim2) %>% 
  summarise(n_plots = n()) %>% 
  group_by(ndep_cat) %>% 
  mutate(perc_lim = round(n_plots/sum(n_plots) * 100, digits = 1)) 
globnut %>% 
  mutate(ndep_cat = cut(exp(ndep), breaks = seq(0, 30, by = 5))) %>% 
  group_by(ndep_cat,lim) %>% 
  summarise(n_plots = n()) %>% 
  group_by(ndep_cat) %>% 
  mutate(perc_lim = round(n_plots/sum(n_plots) *100, digits = 1)) %>% 
  left_join(perc_pal, by = c("ndep_cat", "lim" = "lim2")) %>% 
  rename(`N deposition range` = ndep_cat ,
         `Limitation type` = lim,
         `Number of plots` = n_plots.x,
         `%`= perc_lim.x,
         `Number of plots (Palpurina method)` = n_plots.y,
         `% (Palpurian method)`= perc_lim.y
  ) %>% knitr::kable()

