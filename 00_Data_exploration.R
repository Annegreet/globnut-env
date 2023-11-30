## ---------------------------
##
## Script name: 00_Data_exploration
##
## Purpose of script: Producing scatters according to SEM proposal
##
## Author: Annegreet Veeken
##
## Date Created: 2023-09-21
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

options(scipen = 6, digits = 4) # avoid scientific notation

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")

## Load data
data_dir <- "Z:/_GLOBNUT1.0/"
cont <- read.csv(paste0(data_dir, "GlobNut1.0_contributors.csv"), sep = ";")
spec <- read.csv(paste0(data_dir, "GlobNut1.0_species.csv"))
npk <- read.csv(paste0(data_dir, "GlobNut1.0_nutrients.csv")) %>% 
  mutate(lim = case_when(NP < 13.5 ~ "N-limitation",
                         NP > 16 ~ "P-limitation",
                         TRUE ~ "Unclear")) 
env <- read.csv(paste0(data_dir, "GlobNut1.0_env_variable.csv"))
meta <- read.csv(paste0(data_dir, "GlobNut1.0_metadata.csv")) 
ndep <- read.csv(paste0(data_dir, "GlobNut1.0_ndeposition.csv"))
not_fert <- filter(meta, harm_fert_appl == 0) %>% pull(plot_ID)
spec_ric <- readRDS("outputs/01_Species_diversity.rds")

## Plots
# plot lm equation and r square
lm_eq <- function(x, y) {
  m <- lm(y ~ x)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


# Temperature -----
mat <- npk %>% 
  left_join(spec_ric, by = "plot_ID") %>% 
  left_join(env[,c("plot_ID", "worldclim_bio_1")], by = "plot_ID")

ggplot(mat, aes(x = worldclim_bio_1, y = N)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(mat$worldclim_bio_1, na.rm = TRUE), y = max(mat$N, na.rm = TRUE), label = lm_eq(x = mat$worldclim_bio_1, y = mat$N), hjust=1, parse = TRUE) +
  xlab("MAT") +
  theme_bw()
ggsave("figures/mat_n.png", dpi = 300, height = 5, width = 7)

ggplot(mat, aes(x = worldclim_bio_1, y = P)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(mat$worldclim_bio_1, na.rm = TRUE), y = max(mat$P, na.rm = TRUE), label = lm_eq(x = mat$worldclim_bio_1, y = mat$P), hjust=1, parse = TRUE) +
  xlab("MAT") +
  theme_bw()
ggsave("figures/mat_p.png", dpi = 300, height = 5, width = 7)

ggplot(mat, aes(x = worldclim_bio_1, y = K)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(mat$worldclim_bio_1, na.rm = TRUE), y = max(mat$K, na.rm = TRUE), label = lm_eq(x = mat$worldclim_bio_1, y = mat$K), hjust=1, parse = TRUE) +
  xlab("MAT") +
  theme_bw()
ggsave("figures/mat_k.png", dpi = 300, height = 5, width = 7)

ggplot(mat, aes(x = worldclim_bio_1, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  # geom_smooth(method = "lm", color = "grey50", formula = y~x + I(x^2)) +
  geom_smooth(method = "lm", color = "grey50") +
  stat_regline_equation()
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(mat$worldclim_bio_1, na.rm = TRUE), y = max(mat$biomass, na.rm = TRUE), label = lm_eq(x = mat$worldclim_bio_1, y = mat$biomass), hjust=1, parse = TRUE) +
  xlab("MAT") +
  theme_bw()
ggsave("figures/mat_biomass.png", dpi = 300, height = 5, width = 7)

# Precipitation ----
map <- npk %>% 
  left_join(spec_ric, by = "plot_ID") %>% 
  left_join(env[,c("plot_ID", "worldclim_bio_12")], by = "plot_ID")

ggplot(map, aes(x = worldclim_bio_12, y = N)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(map$worldclim_bio_12, na.rm = TRUE), y = max(map$N, na.rm = TRUE), label = lm_eq(x = map$worldclim_bio_12, y = map$N), hjust=1, parse = TRUE) +
  xlab("MAP") +
  theme_bw()
ggsave("figures/map_N.png", dpi = 300, height = 5, width = 7)

ggplot(map, aes(x = worldclim_bio_12, y = P)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(map$worldclim_bio_12, na.rm = TRUE), y = max(map$P, na.rm = TRUE), label = lm_eq(x = map$worldclim_bio_12, y = map$N), hjust=1, parse = TRUE) +
  xlab("MAP") +
  theme_bw()
ggsave("figures/map_p.png", dpi = 300, height = 5, width = 7)

ggplot(map, aes(x = worldclim_bio_12, y = K)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(map$worldclim_bio_12, na.rm = TRUE), y = max(map$K, na.rm = TRUE), label = lm_eq(x = map$worldclim_bio_12, y = map$K), hjust=1, parse = TRUE) +
  xlab("MAP") +
  theme_bw()
ggsave("figures/map_k.png", dpi = 300, height = 5, width = 7)

ggplot(map, aes(x = worldclim_bio_12, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(map$worldclim_bio_12, na.rm = TRUE), y = max(map$biomass, na.rm = TRUE), label = lm_eq(x = map$worldclim_bio_12, y = map$biomass), hjust=1, parse = TRUE) +
  xlab("MAP") +
  theme_bw()
ggsave("figures/map_biomass.png", dpi = 300, height = 5, width = 7)

# N deposition ----
npk_ndep <- npk %>% 
  left_join(spec_ric, by = "plot_ID") %>% 
  left_join(ndep, by = "plot_ID") %>% 
  filter(n_obs_5yr >= 5)

ggplot(npk_ndep, aes(x = sum_5yr, y = N)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_ndep$sum_5yr, na.rm = TRUE), y = max(npk_ndep$N, na.rm = TRUE), label = lm_eq(x = npk_ndep$sum_5yr, y = npk_ndep$N), hjust=1, parse = TRUE) +
  xlab("5 year sum N deposiotn") +
  theme_bw()
ggsave("figures/ndep_n.png", dpi = 300, height = 5, width = 7)

# unfertilized plots
npk_ndep <- npk %>% 
  left_join(spec_ric, by = "plot_ID") %>% 
  left_join(ndep, by = "plot_ID") %>% 
  filter(n_obs_5yr >= 5) %>% 
  filter(plot_ID %in% not_fert)

ggplot(npk_ndep, aes(x = sum_5yr, y = N)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_ndep$sum_5yr, na.rm = TRUE), y = max(npk_ndep$N, na.rm = TRUE), label = lm_eq(x = npk_ndep$sum_5yr, y = npk_ndep$N), hjust=1, parse = TRUE) +
  xlab("5 year sum N deposiotn") +
  ggtitle("Unfertilized plots") +
  theme_bw()
ggsave("figures/ndep_n_unfert.png", dpi = 300, height = 5, width = 7)


# Nutrients ----
npk_sr <- npk %>% 
  left_join(spec_ric, by = "plot_ID")
ggplot(npk_sr, aes(x = N, y = NP)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$N, na.rm = TRUE), y = max(npk_sr$NP, na.rm = TRUE), label = lm_eq(x = npk_sr$N, y = npk_sr$NP), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/n_np.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = N, y = NK)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$N, na.rm = TRUE), y = max(npk_sr$NK, na.rm = TRUE), label = lm_eq(x = npk_sr$N, y = npk_sr$NK), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/n_nk.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = P, y = NP)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  # geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  # annotate("text", x = max(npk_sr$P, na.rm = TRUE), y = max(npk_sr$NP, na.rm = TRUE), label = lm_eq(x = npk_sr$P, y = npk_sr$NP), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/p_np.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = P, y = KP)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  # geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  # annotate("text", x = max(npk_sr$P, na.rm = TRUE), y = max(npk_sr$KP, na.rm = TRUE), label = lm_eq(x = npk_sr$P, y = npk_sr$KP), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/p_kp.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = K, y = NK)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  # geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  # annotate("text", x = max(npk_sr$P, na.rm = TRUE), y = max(npk_sr$KP, na.rm = TRUE), label = lm_eq(x = npk_sr$P, y = npk_sr$KP), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/k_NK.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = K, y = KP)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$K, na.rm = TRUE), y = max(npk_sr$KP, na.rm = TRUE), label = lm_eq(x = npk_sr$K, y = npk_sr$KP), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/k_KP.png", dpi = 300, height = 5, width = 7)


ggplot(npk_sr, aes(x = N, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$N, na.rm = TRUE), y = max(npk_sr$biomass, na.rm = TRUE), label = lm_eq(x = npk_sr$N, y = npk_sr$biomass), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/N_biomass.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = P, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$P, na.rm = TRUE), y = max(npk_sr$biomass, na.rm = TRUE), label = lm_eq(x = npk_sr$P, y = npk_sr$biomass), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/P_biomass.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = K, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$K, na.rm = TRUE), y = max(npk_sr$biomass, na.rm = TRUE), label = lm_eq(x = npk_sr$K, y = npk_sr$biomass), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/k_biomass.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = N, y = spec_ric)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$N, na.rm = TRUE), y = max(npk_sr$spec_ric, na.rm = TRUE), label = lm_eq(x = npk_sr$N, y = npk_sr$spec_ric), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/N_spec_ric.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = NP, y = spec_ric)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$NP, na.rm = TRUE), y = max(npk_sr$spec_ric, na.rm = TRUE), label = lm_eq(x = npk_sr$NP, y = npk_sr$spec_ric), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/NP_spec_ric.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = NP, y = biomass)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$NP, na.rm = TRUE), y = max(npk_sr$biomass, na.rm = TRUE), label = lm_eq(x = npk_sr$NP, y = npk_sr$biomass), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/NP_biomass.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = NP, y = shan_eve)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$NP, na.rm = TRUE), y = max(npk_sr$shan_eve, na.rm = TRUE), label = lm_eq(x = npk_sr$NP, y = npk_sr$shan_eve), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/NP_shan_eve.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = biomass, y = spec_ric)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$biomass, na.rm = TRUE), y = max(npk_sr$spec_ric, na.rm = TRUE), label = lm_eq(x = npk_sr$biomass, y = npk_sr$spec_ric), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/biomass_specric.png", dpi = 300, height = 5, width = 7)

ggplot(npk_sr, aes(x = shan_eve, y = spec_ric)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  # geom_smooth(method = "lm", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  # annotate("text", x = max(npk_sr$shan_eve, na.rm = TRUE), y = max(npk_sr$spec_ric, na.rm = TRUE), label = lm_eq(x = npk_sr$shan_eve, y = npk_sr$spec_ric), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/speceve_specric.png", dpi = 300, height = 5, width = 7)

# Quadratic relationships ----
lm_eq <- function(x, y) {
  m <- lm(y ~ x + x^2)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
ggplot(npk_sr, aes(x = biomass, y = spec_ric)) +
  geom_point(aes(color = lim), alpha = 0.7) +
  geom_smooth(method = "gam", color = "grey50") +
  scale_color_manual("Limitation", values = c("N-limitation" = "darkorange", "P-limitation" = "purple", "Unclear" = "cyan4")) +
  annotate("text", x = max(npk_sr$biomass, na.rm = TRUE), y = max(npk_sr$spec_ric, na.rm = TRUE), 
           label = lm_eq(x = npk_sr$biomass, y = npk_sr$spec_ric), hjust=1, parse = TRUE) +
  theme_bw()
ggsave("figures/biomass_specric.png", dpi = 300, height = 5, width = 7)

