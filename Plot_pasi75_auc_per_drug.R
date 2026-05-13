library(readr)
library(dplyr)
library(tidyr)

summary_table <- read.csv("summary_PASI75_limma_by_drug.csv")

plot_df <- summary_table %>%
  transmute(
    Vavnad = recode(
      tissue,
      "blood" = "Blod",
      "lesional" = "Lesional hud",
      "nonlesional" = "Icke-lesional hud"
    ),
    Drug = recode(
      drug,
      "adalimumab" = "Adalimumab",
      "ustekinumab" = "Ustekinumab"
    ),
    Baseline_AUC = glm_mean_auc,
    Baseline_SD = glm_sd_auc,
    RF_AUC = rf_mean_auc,
    RF_SD = rf_sd_auc
  ) %>%
  pivot_longer(
    cols = c(Baseline_AUC, RF_AUC),
    names_to = "Modell",
    values_to = "AUC"
  ) %>%
  mutate(
    SD = ifelse(Modell == "Baseline_AUC", Baseline_SD, RF_SD),
    Modell = recode(
      Modell,
      "Baseline_AUC" = "Baseline",
      "RF_AUC" = "Random Forest"
    )
  )
plot_df$Vavnad <- factor(
  plot_df$Vavnad,
  levels = c("Blod", "Lesional hud", "Icke-lesional hud")
)

p <- ggplot(
  plot_df,
  aes(
    x = Drug,
    y = AUC,
    fill = Drug,
    alpha = Modell
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = AUC - SD,
      ymax = AUC + SD,
      group = interaction(Drug, Modell)
    ),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.4
  ) +
  facet_wrap(~ Vavnad, nrow = 1) +
  scale_fill_manual(
    name = "Drug",
    values = c(
      "Adalimumab" = "#128EBB",
      "Ustekinumab" = "#EF3B4A"
    )
  ) +
  scale_alpha_manual(
    name = "Model",
    values = c(
      "Baseline" = 0.35,
      "Random Forest" = 1
    ),
    labels = c(
      "Null (baseline-PASI only)",
      "RF (signalvägsgener)"
    )
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = NULL,
    x = NULL,
    y = "AUC (medel ± SD)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 13, face = "bold")
  ) +
  guides(
    alpha = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

print(p)