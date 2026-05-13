library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plot_df <- bind_rows(
  
  # Blood - Adalimumab
  data.frame(
    Vavnad = "Blod",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$blood$adalimumab$lm_rmse_mean,
      topgene_results$top50$blood$adalimumab$rf_rmse_mean,
      topgene_results$top20$blood$adalimumab$rf_rmse_mean,
      topgene_results$top10$blood$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$blood$adalimumab$lm_rmse_sd,
      topgene_results$top50$blood$adalimumab$rf_rmse_sd,
      topgene_results$top20$blood$adalimumab$rf_rmse_sd,
      topgene_results$top10$blood$adalimumab$rf_rmse_sd
    )
  ),
  
  # Blood - Ustekinumab
  data.frame(
    Vavnad = "Blod",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$blood$ustekinumab$lm_rmse_mean,
      topgene_results$top50$blood$ustekinumab$rf_rmse_mean,
      topgene_results$top20$blood$ustekinumab$rf_rmse_mean,
      topgene_results$top10$blood$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$blood$ustekinumab$lm_rmse_sd,
      topgene_results$top50$blood$ustekinumab$rf_rmse_sd,
      topgene_results$top20$blood$ustekinumab$rf_rmse_sd,
      topgene_results$top10$blood$ustekinumab$rf_rmse_sd
    )
  ),
  
  # Lesional - Adalimumab
  data.frame(
    Vavnad = "Lesional hud",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$lesional$adalimumab$lm_rmse_mean,
      topgene_results$top50$lesional$adalimumab$rf_rmse_mean,
      topgene_results$top20$lesional$adalimumab$rf_rmse_mean,
      topgene_results$top10$lesional$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$lesional$adalimumab$lm_rmse_sd,
      topgene_results$top50$lesional$adalimumab$rf_rmse_sd,
      topgene_results$top20$lesional$adalimumab$rf_rmse_sd,
      topgene_results$top10$lesional$adalimumab$rf_rmse_sd
    )
  ),
  
  # Lesional - Ustekinumab
  data.frame(
    Vavnad = "Lesional hud",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$lesional$ustekinumab$lm_rmse_mean,
      topgene_results$top50$lesional$ustekinumab$rf_rmse_mean,
      topgene_results$top20$lesional$ustekinumab$rf_rmse_mean,
      topgene_results$top10$lesional$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$lesional$ustekinumab$lm_rmse_sd,
      topgene_results$top50$lesional$ustekinumab$rf_rmse_sd,
      topgene_results$top20$lesional$ustekinumab$rf_rmse_sd,
      topgene_results$top10$lesional$ustekinumab$rf_rmse_sd
    )
  ),
  
  # Nonlesional - Adalimumab
  data.frame(
    Vavnad = "Icke-lesional hud",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$nonlesional$adalimumab$lm_rmse_mean,
      topgene_results$top50$nonlesional$adalimumab$rf_rmse_mean,
      topgene_results$top20$nonlesional$adalimumab$rf_rmse_mean,
      topgene_results$top10$nonlesional$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$nonlesional$adalimumab$lm_rmse_sd,
      topgene_results$top50$nonlesional$adalimumab$rf_rmse_sd,
      topgene_results$top20$nonlesional$adalimumab$rf_rmse_sd,
      topgene_results$top10$nonlesional$adalimumab$rf_rmse_sd
    )
  ),
  
  # Nonlesional - Ustekinumab
  data.frame(
    Vavnad = "Icke-lesional hud",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener"),
    MedelRMSE = c(
      topgene_results$top50$nonlesional$ustekinumab$lm_rmse_mean,
      topgene_results$top50$nonlesional$ustekinumab$rf_rmse_mean,
      topgene_results$top20$nonlesional$ustekinumab$rf_rmse_mean,
      topgene_results$top10$nonlesional$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      topgene_results$top50$nonlesional$ustekinumab$lm_rmse_sd,
      topgene_results$top50$nonlesional$ustekinumab$rf_rmse_sd,
      topgene_results$top20$nonlesional$ustekinumab$rf_rmse_sd,
      topgene_results$top10$nonlesional$ustekinumab$rf_rmse_sd
    )
  )
)

# ------------------------------------------------------------
# Factors
# ------------------------------------------------------------
  plot_df <- plot_df %>%
    mutate(
      Vavnad = factor(Vavnad, levels = c("Blod", "Lesional hud", "Icke-lesional hud")),
      Lakemedel = factor(Lakemedel, levels = c("Adalimumab", "Ustekinumab")),
      Modell = factor(Modell, levels = c("Baseline", "Topp 50 gener", "Topp 20 gener", "Topp 10 gener")),
      
      Farggrupp = factor(
        if_else(
          Modell == "Baseline",
          "Baseline",
          paste(Lakemedel, Modell, sep = " - ")
        ),
        levels = c(
          "Baseline",
          "Adalimumab - Topp 50 gener",
          "Adalimumab - Topp 20 gener",
          "Adalimumab - Topp 10 gener",
          "Ustekinumab - Topp 50 gener",
          "Ustekinumab - Topp 20 gener",
          "Ustekinumab - Topp 10 gener"
        )
      )
    )

p <- ggplot(plot_df, aes(x = Lakemedel, y = MedelRMSE, fill = Farggrupp)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = MedelRMSE - SDRMSE, ymax = MedelRMSE + SDRMSE),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.4
  ) +
  facet_wrap(~ Vavnad, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Baseline" = "grey60",
      "Adalimumab - Topp 50 gener" = "#128EBB",
      "Adalimumab - Topp 20 gener" = "#5FAFCC",
      "Adalimumab - Topp 10 gener" = "#B8DCE8",
      "Ustekinumab - Topp 50 gener" = "#EF3B4A",
      "Ustekinumab - Topp 20 gener" = "#F47A84",
      "Ustekinumab - Topp 10 gener" = "#F8C3C7"
    ),
    breaks = c(
      "Baseline",
      "Adalimumab - Topp 50 gener",
      "Adalimumab - Topp 20 gener",
      "Adalimumab - Topp 10 gener",
      "Ustekinumab - Topp 50 gener",
      "Ustekinumab - Topp 20 gener",
      "Ustekinumab - Topp 10 gener"
    )
  ) +
  labs(
    title = "Jämförelse mellan baseline-modell och modeller med toppgener",
    x = NULL,
    y = "RMSE (medel ± SD)",
    fill = "Läkemedel och modell"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p)