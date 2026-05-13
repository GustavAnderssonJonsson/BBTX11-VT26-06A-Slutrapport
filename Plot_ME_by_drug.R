library(ggplot2)
library(dplyr)

# ============================================================
# Plot: Jämförelse mellan baseline-modell och Random Forest
# från scriptet:
#   ME_by_drug_with_baseline.R
# ============================================================

plot_df <- bind_rows(
  
  # ----------------------------------------------------------
  # BLOD
  # ----------------------------------------------------------
  data.frame(
    Vavnad = "Blod",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_blood$adalimumab$lm_rmse_mean,
      res_blood$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_blood$adalimumab$lm_rmse_sd,
      res_blood$adalimumab$rf_rmse_sd
    )
  ),
  
  data.frame(
    Vavnad = "Blod",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_blood$ustekinumab$lm_rmse_mean,
      res_blood$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_blood$ustekinumab$lm_rmse_sd,
      res_blood$ustekinumab$rf_rmse_sd
    )
  ),
  
  # ----------------------------------------------------------
  # LESIONAL HUD
  # ----------------------------------------------------------
  data.frame(
    Vavnad = "Lesional hud",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_lesional$adalimumab$lm_rmse_mean,
      res_lesional$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_lesional$adalimumab$lm_rmse_sd,
      res_lesional$adalimumab$rf_rmse_sd
    )
  ),
  
  data.frame(
    Vavnad = "Lesional hud",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_lesional$ustekinumab$lm_rmse_mean,
      res_lesional$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_lesional$ustekinumab$lm_rmse_sd,
      res_lesional$ustekinumab$rf_rmse_sd
    )
  ),
  
  # ----------------------------------------------------------
  # ICKE-LESIONAL HUD
  # ----------------------------------------------------------
  data.frame(
    Vavnad = "Icke-lesional hud",
    Lakemedel = "Adalimumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_nonlesional$adalimumab$lm_rmse_mean,
      res_nonlesional$adalimumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_nonlesional$adalimumab$lm_rmse_sd,
      res_nonlesional$adalimumab$rf_rmse_sd
    )
  ),
  
  data.frame(
    Vavnad = "Icke-lesional hud",
    Lakemedel = "Ustekinumab",
    Modell = c("Baseline", "Random Forest"),
    MedelRMSE = c(
      res_nonlesional$ustekinumab$lm_rmse_mean,
      res_nonlesional$ustekinumab$rf_rmse_mean
    ),
    SDRMSE = c(
      res_nonlesional$ustekinumab$lm_rmse_sd,
      res_nonlesional$ustekinumab$rf_rmse_sd
    )
  )
)

# ============================================================
# Factor order
# ============================================================

plot_df$Vavnad <- factor(
  plot_df$Vavnad,
  levels = c("Blod", "Lesional hud", "Icke-lesional hud")
)

plot_df$Lakemedel <- factor(
  plot_df$Lakemedel,
  levels = c("Adalimumab", "Ustekinumab")
)

plot_df$Modell <- factor(
  plot_df$Modell,
  levels = c("Baseline", "Random Forest")
)

# ============================================================
# Plot
# ============================================================

plot_df$Vavnad <- factor(
  plot_df$Vavnad,
  levels = c("Blod", "Lesional hud", "Icke-lesional hud")
)

plot_df$Lakemedel <- factor(
  plot_df$Lakemedel,
  levels = c("Adalimumab", "Ustekinumab")
)

plot_df$Modell <- factor(
  plot_df$Modell,
  levels = c("Baseline", "Random Forest")
)

p <- ggplot(
  plot_df,
  aes(
    x = Lakemedel,
    y = MedelRMSE,
    fill = Lakemedel,
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
      ymin = MedelRMSE - SDRMSE,
      ymax = MedelRMSE + SDRMSE,
      group = interaction(Lakemedel, Modell)
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
      "Baseline" = "Null (baseline-PASI only)",
      "Random Forest" = "RF (eigengenes)"
    )
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "Mean RMSE ± SD over splits"
  ) +
  
  coord_cartesian(ylim = c(0, 10)) +
  
  theme_bw(base_size = 13) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p)

ggsave(
  filename = "jamforelse_baseline_randomforest_ME_layout_like_example_RMSE.png",
  plot = p,
  width = 11,
  height = 4.5,
  dpi = 300
)

# ============================================================
# Save
# ============================================================

ggsave(
  filename = "jamforelse_baseline_randomforest_ME_per_vavnad_och_lakemedel_RMSE.png",
  plot = p,
  width = 11,
  height = 4.5,
  dpi = 300
)