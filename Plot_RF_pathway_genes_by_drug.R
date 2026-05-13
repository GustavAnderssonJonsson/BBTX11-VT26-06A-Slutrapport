suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

# ============================================================
# RF pathway model per drug arm
# Models:
#   1) Baseline: delta_PASI ~ PASI_baseline
#   2) RF pathway genes: delta_PASI ~ PASI_baseline + pathway genes
#
# Outcome:
#   absolute ΔPASI = PASI_baseline - PASI_followup
#   positive value = improvement

# ------------------------------------------------------------
# 1) Pathway genes
# ------------------------------------------------------------

tnf_genes <- c(
  "TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "TRAF2", "TRAF5", "RIPK1",
  "NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "IKBKG", "CHUK",
  "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK14",
  "JUN", "FOS", "CASP3", "CASP8", "BIRC2", "BIRC3",
  "CXCL8", "CCL2", "CCL5", "ICAM1", "VCAM1"
)

il23_genes <- c(
  "IL12A", "IL12B", "IL23A", "IL12RB1", "IL12RB2", "IL23R",
  "JAK2", "TYK2", "STAT3", "STAT4", "RORC",
  "IL17A", "IL17F", "IL22", "CCR6", "CCL20", "IFNG"
)

psoriasis_genes <- c(
  "S100A7", "S100A8", "S100A9", "DEFB4A",
  "KRT16", "KRT17", "PI3", "SERPINB3", "SERPINB4",
  "LCE3D", "LCE3E"
)

drug_gene_sets <- list(
  adalimumab  = unique(c(tnf_genes, psoriasis_genes)),
  ustekinumab = unique(c(il23_genes, psoriasis_genes))
)

# ------------------------------------------------------------
# 2) Helper functions
# ------------------------------------------------------------

clean_ensembl <- function(x) {
  sub("\\..*$", "", x)
}

map_symbols_to_ensembl <- function(symbols) {
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENSEMBL")
  )
  
  gene_map <- gene_map[!is.na(gene_map$ENSEMBL), , drop = FALSE]
  unique(clean_ensembl(gene_map$ENSEMBL))
}

make_repeated_folds <- function(n, k = 5, repeats = 5, seed = 42) {
  set.seed(seed)
  folds_all <- list()
  
  for (r in seq_len(repeats)) {
    fold_id <- sample(rep(seq_len(k), length.out = n))
    
    for (f in seq_len(k)) {
      folds_all[[paste0("rep", r, "_fold", f)]] <- which(fold_id == f)
    }
  }
  
  folds_all
}

remove_near_zero_variance <- function(X, min_var = 1e-8) {
  vars <- apply(X, 2, var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  keep <- names(vars)[vars > min_var]
  X[, keep, drop = FALSE]
}

calc_regression_metrics <- function(pred, obs) {
  data.frame(
    rmse = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
    cor = suppressWarnings(cor(pred, obs, use = "complete.obs"))
  )
}

# ------------------------------------------------------------
# 3) Prepare dataset
# ------------------------------------------------------------

prepare_pathway_dataset <- function(expr_mat,
                                    sample_info,
                                    tissue_name,
                                    drug_name,
                                    gene_symbols) {
  
  if ("sample_id" %in% colnames(sample_info)) {
    rownames(sample_info) <- sample_info$sample_id
  }
  
  expr_mat <- as.data.frame(expr_mat)
  colnames(expr_mat) <- clean_ensembl(colnames(expr_mat))
  
  selected_ensembl <- map_symbols_to_ensembl(gene_symbols)
  selected_ensembl <- intersect(selected_ensembl, colnames(expr_mat))
  
  if (length(selected_ensembl) < 2) {
    stop("Too few pathway genes matched for ", tissue_name, " / ", drug_name)
  }
  
  expr_selected <- expr_mat[, selected_ensembl, drop = FALSE]
  
  si <- sample_info %>%
    mutate(
      week = as.numeric(as.character(week)),
      PASI = as.numeric(as.character(PASI)),
      drug = tolower(trimws(as.character(drug))),
      sample_id = rownames(sample_info)
    ) %>%
    select(sample_id, patient, drug, week, PASI)
  
  baseline_tbl <- si %>%
    filter(week == 0, drug == drug_name) %>%
    group_by(patient, drug) %>%
    slice(1) %>%
    ungroup() %>%
    rename(
      PASI_baseline = PASI,
      sample_id_baseline = sample_id
    )
  
  followup_tbl <- si %>%
    filter(week == 12, drug == drug_name) %>%
    group_by(patient, drug) %>%
    slice(1) %>%
    ungroup() %>%
    rename(PASI_followup = PASI)
  
  delta_tbl <- baseline_tbl %>%
    inner_join(followup_tbl, by = c("patient", "drug")) %>%
    mutate(delta_PASI = PASI_baseline - PASI_followup) %>%
    filter(
      !is.na(delta_PASI),
      !is.na(PASI_baseline),
      sample_id_baseline %in% rownames(expr_selected)
    )
  
  X <- expr_selected[delta_tbl$sample_id_baseline, , drop = FALSE]
  
  model_df <- delta_tbl %>%
    select(patient, PASI_baseline, PASI_followup, delta_PASI) %>%
    bind_cols(as.data.frame(X))
  
  list(
    tissue = tissue_name,
    drug = drug_name,
    model_df = model_df,
    n_patients = nrow(model_df),
    n_genes = ncol(X)
  )
}

# ------------------------------------------------------------
# 4) Repeated 5-fold CV
# ------------------------------------------------------------

run_repeated_cv_pathway <- function(model_df,
                                    tissue_name,
                                    drug_name,
                                    k = 5,
                                    repeats = 5,
                                    num_trees = 1000,
                                    seed = 42) {
  
  n <- nrow(model_df)
  if (n < k) stop("Too few patients for CV: ", tissue_name, " / ", drug_name)
  
  folds <- make_repeated_folds(n, k = k, repeats = repeats, seed = seed)
  results <- list()
  
  gene_cols <- setdiff(
    colnames(model_df),
    c("patient", "PASI_baseline", "PASI_followup", "delta_PASI")
  )
  
  for (fold_name in names(folds)) {
    
    test_i <- folds[[fold_name]]
    train_i <- setdiff(seq_len(n), test_i)
    
    train_df_raw <- model_df[train_i, , drop = FALSE]
    test_df_raw  <- model_df[test_i, , drop = FALSE]
    
    X_train <- train_df_raw[, gene_cols, drop = FALSE]
    X_test  <- test_df_raw[, gene_cols, drop = FALSE]
    
    # Remove near-zero variance genes inside training fold only
    X_train <- remove_near_zero_variance(X_train)
    kept_genes <- colnames(X_train)
    
    if (length(kept_genes) < 2) next
    
    X_test <- X_test[, kept_genes, drop = FALSE]
    
    # ----------------------------
    # Baseline linear model
    # ----------------------------
    lm_model <- lm(delta_PASI ~ PASI_baseline, data = train_df_raw)
    lm_pred <- predict(lm_model, newdata = test_df_raw)
    
    lm_metrics <- calc_regression_metrics(
      pred = lm_pred,
      obs = test_df_raw$delta_PASI
    ) %>%
      mutate(Model = "Baseline")
    
    # ----------------------------
    # RF pathway gene model
    # ----------------------------
    train_rf <- bind_cols(
      data.frame(
        delta_PASI = train_df_raw$delta_PASI,
        PASI_baseline = train_df_raw$PASI_baseline
      ),
      as.data.frame(X_train)
    )
    
    test_rf <- bind_cols(
      data.frame(
        delta_PASI = test_df_raw$delta_PASI,
        PASI_baseline = test_df_raw$PASI_baseline
      ),
      as.data.frame(X_test)
    )
    
    rf_model <- ranger(
      dependent.variable.name = "delta_PASI",
      data = train_rf,
      num.trees = num_trees,
      importance = "permutation",
      seed = seed
    )
    
    rf_pred <- predict(rf_model, data = test_rf)$predictions
    
    rf_metrics <- calc_regression_metrics(
      pred = rf_pred,
      obs = test_rf$delta_PASI
    ) %>%
      mutate(Model = "RF pathway-gener")
    
    results[[fold_name]] <- bind_rows(lm_metrics, rf_metrics) %>%
      mutate(
        Tissue = tissue_name,
        Drug = drug_name,
        Fold = fold_name,
        n_train = length(train_i),
        n_test = length(test_i),
        n_genes_initial = length(gene_cols),
        n_genes_after_nzv = length(kept_genes)
      )
  }
  
  bind_rows(results)
}

# ------------------------------------------------------------
# 5) Run all tissues and drug arms
# ------------------------------------------------------------

tissue_objects <- list(
  "Blod" = list(expr = expr_blood, sample_info = sample_info_blood),
  "Lesional hud" = list(expr = expr_lesional, sample_info = sample_info_lesional),
  "Icke-lesional hud" = list(expr = expr_nonlesional, sample_info = sample_info_nonlesional)
)

all_results <- list()

for (tissue_name in names(tissue_objects)) {
  
  for (drug_name in names(drug_gene_sets)) {
    
    message("Running ", tissue_name, " / ", drug_name)
    
    dat <- prepare_pathway_dataset(
      expr_mat = tissue_objects[[tissue_name]]$expr,
      sample_info = tissue_objects[[tissue_name]]$sample_info,
      tissue_name = tissue_name,
      drug_name = drug_name,
      gene_symbols = drug_gene_sets[[drug_name]]
    )
    
    cv_res <- run_repeated_cv_pathway(
      model_df = dat$model_df,
      tissue_name = tissue_name,
      drug_name = drug_name,
      k = 5,
      repeats = 5,
      num_trees = 1000,
      seed = 42
    )
    
    all_results[[paste(tissue_name, drug_name, sep = "_")]] <- cv_res
  }
}

cv_results <- bind_rows(all_results)

saveRDS(cv_results, "pathway_all_genes_per_drug_cv_results.rds")
write.csv(cv_results, "pathway_all_genes_per_drug_cv_results.csv", row.names = FALSE)

# ------------------------------------------------------------
# 6) Summary
# ------------------------------------------------------------

summary_df <- cv_results %>%
  group_by(Tissue, Drug, Model) %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse   = sd(rmse, na.rm = TRUE),
    mean_cor  = mean(cor, na.rm = TRUE),
    sd_cor    = sd(cor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Drug_label = recode(
      Drug,
      "adalimumab" = "Adalimumab",
      "ustekinumab" = "Ustekinumab"
    )
  )

print(summary_df)

saveRDS(summary_df, "pathway_all_genes_per_drug_summary.rds")
write.csv(summary_df, "pathway_all_genes_per_drug_summary.csv", row.names = FALSE)

# ------------------------------------------------------------
# 7) Plot RMSE
# ------------------------------------------------------------

plot_df <- summary_df %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c("Blod", "Lesional hud", "Icke-lesional hud")
    ),
    Drug_label = factor(
      Drug_label,
      levels = c("Adalimumab", "Ustekinumab")
    ),
    Model_group = ifelse(
      Model == "Baseline",
      "Null (baseline-PASI only)",
      "RF (signalvägsgener)"
    ),
    Model_group = factor(
      Model_group,
      levels = c("Null (baseline-PASI only)", "RF (signalvägsgener)")
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = Drug_label,
    y = mean_rmse,
    fill = Drug_label,
    alpha = Model_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = mean_rmse - sd_rmse,
      ymax = mean_rmse + sd_rmse,
      group = interaction(Drug_label, Model_group)
    ),
    position = position_dodge(width = 0.7),
    width = 0.2,
    linewidth = 0.4
  ) +
  facet_wrap(~ Tissue, nrow = 1) +
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
      "Null (baseline-PASI only)" = 0.35,
      "RF (signalvägsgener)" = 1
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
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    alpha = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  )

print(p)

ggsave(
  filename = "pathway_all_genes_RMSE_new_layout.png",
  plot = p,
  width = 11,
  height = 4.5,
  dpi = 300
)