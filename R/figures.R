############################################################
# Supplementary Material – Figure Reproduction Script
# Manuscript: "How much sampling is enough in multivariate community ecology? A simulation and cost-benefit framework"
# Author: Arturo Sanchez-Porras, Edlin Guerra-Castro
# Notes:
# - This script assumes that all required data files are available in the /data/ directory.
# - Figures are generated using ggplot2 and base R graphics.
# - Each figure chunk is labeled according to its appearance in the manuscript.
############################################################

#Required packages:
library(tidyverse)
library(vegan)
library(SSP)
library(ecocbo)
library(ggplot2)
library(patchwork)


# Figure 2 ----------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
dataMacro <- read_csv("./data/dataMacroMultiBeach.csv") |>
  as.data.frame() |>
  filter(locality %in% c("El Cuyo", "Celestun", "Sisal")) |>
  mutate(dummy = 1)

MacroTransformed <- dataMacro |>
  select(-c("locality")) |>
  mutate(across(where(is.numeric), ~ log(.x + 1)))

nmdsMacro <- vegan::metaMDS(
  MacroTransformed[, -1],
  distance = "bray",
  k = 2,
  trymax = 30,
  autotransform = FALSE
)

scoresMacro <- as.data.frame(vegan::scores(nmdsMacro, display = "sites")) |>
  bind_cols(MacroTransformed[1]) |>
  bind_cols(dataMacro[1])
scoresMacro$zone[scoresMacro$zone == "Dry"] <- "Supratidal"
scoresMacro$zone[scoresMacro$zone == "Wet"] <- "Intertidal"

scoresMacro$transect <- rownames(scoresMacro)

fig2.a <- scoresMacro |>
  mutate(Zone = zone)|>
  ggplot(aes(x = NMDS1, y = NMDS2, color = zone)) +
  stat_ellipse(level = 0.8,
               alpha = 0.3,
               show.legend = FALSE)+
  geom_point(size = 3) +
  annotate("text",
           label = "A",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)+
  annotate("text",
           x = Inf,
           y = Inf ,
           vjust = 1.0,
           hjust = 1.0,
           size = 3,
           label = paste0("2D Stress = ", round(nmdsMacro$stress, 3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linewidth = 0.5))

fig2.a

#Mangrove root epibionts
dataEpi <- SSP::epibionts

epiTransformed <- dataEpi |>
  mutate(site = as.character(site)) |>
  mutate(across(where(is.numeric), ~ .x^(1 / 4))) |>
  as.data.frame()


nmdsEpi <- vegan::metaMDS(
  epiTransformed[, -c(1:2)],
  distance = "bray",
  k = 2,
  trymax = 30,
  autotransform = FALSE
)

scoresEpi <- as.data.frame(vegan::scores(nmdsEpi, display = "sites")) |>
  bind_cols(epiTransformed[c(1:2)])

fig2.b <- scoresEpi |>
  mutate(Sector = factor(sector,
                         levels = c("E", "M", "I"),
                         labels = c("External", "Intermediate", "Internal")))|>
  ggplot(aes(x = NMDS1, y = NMDS2, color = Sector)) +
  geom_point(size = 3) +
  annotate("text",
           label = "B",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)+
  annotate("text",
           x = Inf,
           y = Inf ,
           vjust = 1.0,
           hjust = 1.0,
           size = 3,
           label = paste0("2D Stress = ", round(nmdsEpi$stress, 3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        # axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linewidth = 0.5))

fig2.b

# Biomass of coral reef fish
dataFish <- read_csv("data/dataFish.csv") |>
  select(-`SAMPLE CODE`) |>
  as.data.frame()
load("./data/dataFish.rda")

FishTransformed <- dataFish |>
  mutate(LOCATION = factor(LOCATION)) |>
  mutate(across(where(is.numeric), ~ log(.x + 1)))


nmdsFish <- vegan::metaMDS(
  FishTransformed[, -c(1:2)],
  distance = "bray",
  k = 2,
  trymax = 30,
  autotransform = FALSE
)

scoresFish <- as.data.frame(vegan::scores(nmdsFish, display = "sites")) |>
  bind_cols(FishTransformed[c(1:2)])

scoresFish$transect <- rownames(scoresFish)

fig2.c <- scoresFish |>
  mutate(Region = factor(REGION,
                         levels = c("NORTH", "EAST", "SOUTH", "WEST"),
                         labels = c("North", "East", "South", "West")))|>
  ggplot(aes(x = NMDS1, y = NMDS2, color = Region)) +
  stat_ellipse(level = 0.8,
               alpha = 0.3,
               show.legend = FALSE)+
  geom_point(size = 3) +
  annotate("text",
           label = "C",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)+
  annotate("text",
           x = Inf,
           y = Inf ,
           vjust = 1.0,
           hjust = 1.0,
           size = 3,
           label = paste0("2D Stress = ", round(nmdsFish$stress, 3))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linewidth = 0.5))

fig2.c

fig_mds <- fig2.a / fig2.b / fig2.c

ggsave(
  "./figures/figure_2.pdf",
  plot = fig_mds,
  device = "pdf",
  width = 18,
  height = 24,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_2.png",
  plot = fig_mds,
  device = "png",
  width = 18,
  height = 24,
  units = "cm",
  dpi = 600
)

# Figure 3 -------------------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
dataMacro <- dataMacro[, names(dataMacro) != "dummy"]
simResMacro <- prep_data(
  as.data.frame(dataMacro[,-1]),
  type = "counts",
  Sest.method = "average",
  cases = 50,
  N = 100,
  M = 1,
  n = 30,
  m = 1,
  k = 30,
  transformation = "none",
  method = "bray",
  dummy = TRUE,
  useParallel = TRUE,
  model = "single.factor",
  jitter.base = 0.5
)

betaMacro <- sim_beta(simResMacro)

cboMacro <- sim_cbo(betaMacro, cn = 25, cm = 60)

fig3.a <- plot_power(
  betaMacro,
  method = "both",
  completePlot = TRUE,
  cbo = cboMacro
)

fig3.a

#Mangrove root epibionts
simResEpi <- prep_data(
  epiTransformed,
  type = "cover",
  Sest.method = "average",
  cases = 50,
  N = 100,
  M = 20,
  n = 15,
  m = 15,
  k = 30,
  transformation = "none",
  method = "bray",
  dummy = TRUE,
  useParallel = TRUE,
  model = "nested.symmetric",
  jitter.base = 0.5
)

betaEpi <- sim_beta(simResEpi)

cboEpi <- sim_cbo(betaEpi, cn = 45, cm = 20)

fig3.b <- plot_power(betaEpi, method = "both", cbo = cboEpi)

fig3.b

# Biomass of coral reef fish
simResFish <- prep_data(
  FishTransformed,
  type = "cover",
  Sest.method = "average",
  cases = 50,
  N = 100,
  M = 20,
  n = 15,
  m = 10,
  k = 30,
  transformation = "none",
  method = "bray",
  dummy = TRUE,
  useParallel = TRUE,
  model = "nested.symmetric",
  jitter.base = 0.5
)

betaFish <- sim_beta(simResFish)

cboFish <- sim_cbo(betaFish, cn = 35, cm = 90)

fig3.c <- plot_power(
  betaFish,
  method = "both",
  cbo = cboFish
)

fig3.c

fig_ecocbo <- fig3.a / fig3.b / fig3.c

ggsave(
  "./figures/figure_3.pdf",
  plot = fig_ecocbo,
  device = "pdf",
  width = 18,
  height = 24,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_3.png",
  plot = fig_ecocbo,
  device = "png",
  width = 18,
  height = 24,
  units = "cm",
  dpi = 600
)


# Figure S1.1 -------------------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
ES1 <- sim_ES(dataMacro[,-1],
              steps = 10,
              type = "counts",
              Sest.method = "average",
              cases = 5,
              N = 100,
              M = NULL,
              n = 30,
              m = NULL,
              k = 30,
              transformation = "none",
              method = "bray",
              dummy = TRUE,
              useParallel = TRUE,
              model = "single.factor",
              jitter.base = 0.5)

figS1_1.a <- plot(ES1, type = "stacked", inferential_var = "omega2")+
  annotate("text",
           label = "A",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)

#Mangroove root epibionts
ES2 <- sim_ES(SSP::epibionts, steps = 6, cases = 4, M = 20, n = 20, m = 15, k = 30, model = "nested.symmetric")

figS1_1.b <- plot.effect_size_data_nested(ES2, type = "stacked", inferential_var = "omega2")+
  annotate("text",
           label = "B",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)

fig_ES_simper <- figS1_1.a / figS1_1.b

ggsave(
  "./figures/figure_S1_1.pdf",
  plot = fig_ES_simper,
  device = "pdf",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_S1_1.png",
  plot = fig_ES_simper,
  device = "png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

# Figure S1.2 -------------------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
figS1_2.a <- plot(ES1, type = "scatter", inferential_var = "omega2", add_smooth = TRUE)+
  annotate("text",
           label = "A",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)

#Mangroove root epibionts
figS1_2.b <- plot.effect_size_data_nested(ES2, type = "scatter", inferential_var = "omega2", add_smooth = TRUE)+
  annotate("text",
           label = "B",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)

fig_ES_scatter <- figS1_2.a / figS1_2.b

ggsave(
  "./figures/figure_S1_2.pdf",
  plot = fig_ES_scatter,
  device = "pdf",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_S1_2.png",
  plot = fig_ES_scatter,
  device = "png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

# Figure S1.3 -------------------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
figS1_3.a <- plot(ES1, type = "true_h0", h0_x_var = "ecological_effect", facet_effort = TRUE)

#Mangroove root epibionts
figS1_3.b <- plot.effect_size_data_nested(ES2, type = "true_h0",
                                          h0_x_var = "ecological_effect",
                                          facet_effort = FALSE)+
  annotate("text",
           label = "B",
           x = -Inf,
           y = Inf,
           vjust = 1.5,
           hjust = -1.0,
           size = 5)

fig_ES_trueH0 <- figS1_3.a

ggsave(
  "./figures/figure_S1_3.pdf",
  plot = fig_ES_trueH0,
  device = "pdf",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_S1_3.png",
  plot = fig_ES_trueH0,
  device = "png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

# Figure S1.4 -------------------------------------------------------------------------------------------------------------------------------------------
#Beach macrofauna
figS1_4.a <- plot(ES1, type = "ordination", reduction_levels = c(0, 0.3, 0.6, 1),
                  label_centroids = FALSE,
                  point_alpha = 0.5, point_size = 1)+
  ggplot2::stat_ellipse()


figS1_4.b <- plot.effect_size_data_nested(ES2, type = "ordination",
                             reduction_levels = c(0, 1/3, 2/3, 1),
                             label_centroids = FALSE,
                             point_alpha = 0.5, point_size = 1)+
  ggplot2::stat_ellipse()

fig_ES_centroid <- figS1_4.a / figS1_4.b

ggsave(
  "./figures/figure_S1_4.pdf",
  plot = fig_ES_centroid,
  device = "pdf",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)

ggsave(
  "./figures/figure_S1_4.png",
  plot = fig_ES_centroid,
  device = "png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 600
)
