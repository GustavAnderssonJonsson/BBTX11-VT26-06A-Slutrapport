library(ggplot2)
library(dplyr)
library(tidyr)

pasi75 <- read.csv("summary_PASI75_limma_by_drug.csv") %>%
  mutate(Endpoint = "PASI75")

pasi50 <- read.csv("summary_PASI50_limma_by_drug.csv") %>%
  mutate(Endpoint = "PASI50")

heat_df <- bind_rows(pasi75, pasi50) %>%
  transmute(
    Vavnad = recode(
      tissue,
      "blood" = "Blod",
      "lesional" = "Lesional hud",
      "nonlesional" = "Icke-lesional hud"
    ),
    Lakemedel = recode(
      drug,
      "adalimumab" = "Adalimumab",
      "ustekinumab" = "Ustekinumab"
    ),
    Endpoint = Endpoint,
    Baseline = glm_mean_auc,
    RandomForest = rf_mean_auc
  ) %>%
  pivot_longer(
    cols = c(Baseline, RandomForest),
    names_to = "Modell",
    values_to = "AUC"
  ) %>%
  mutate(
    Modell = recode(
      Modell,
      "Baseline" = "Baseline",
      "RandomForest" = "RF"
    ),
    Kolumn = paste(Endpoint, Modell, Lakemedel, sep = " - ")
  )

heat_df$Vavnad <- factor(
  heat_df$Vavnad,
  levels = c("Blod", "Lesional hud", "Icke-lesional hud")
)

heat_df$Endpoint <- factor(heat_df$Endpoint, levels = c("PASI75", "PASI50"))
heat_df$Modell <- factor(heat_df$Modell, levels = c("Baseline", "RF"))
heat_df$Lakemedel <- factor(heat_df$Lakemedel, levels = c("Adalimumab", "Ustekinumab"))

heat_df <- heat_df %>%
  arrange(Lakemedel, Endpoint, Modell) %>%
  mutate(Kolumn = factor(Kolumn, levels = unique(Kolumn)))

p <- ggplot(heat_df, aes(x = Kolumn, y = Vavnad, fill = AUC)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", AUC)), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#D73027",
    mid = "white",
    high = "#2C7FB8",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "AUC"
  ) +
  labs(
    title = "AUC för klassificering av behandlingssvar",
    subtitle = "Separata modeller per läkemedel för PASI75 och PASI50",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

print(p)

ggsave(
  filename = "heatmap_auc_pasi75_pasi50_per_drug.png",
  plot = p,
  width = 14,
  height = 5,
  dpi = 300
)