library(dplyr)
library(ranger)

# ------------------------------------------------------------
# Select top variable genes
# ------------------------------------------------------------
select_top_variable_features <- function(features, n = 50) {
  features <- as.matrix(features)
  vars <- apply(features, 2, var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  top <- order(vars, decreasing = TRUE)[seq_len(min(n, ncol(features)))]
  features[, top, drop = FALSE]
}

# ------------------------------------------------------------
# Load tissue-specific expression data + metadata
# Expression matrix:
#   rows = sample IDs
#   cols = genes
# ------------------------------------------------------------
load_tissue_expression_data <- function(tissue,
                                        min_matched = 10) {
  expr_obj_name <- paste0("expr_", tissue)
  si_obj_name   <- paste0("sample_info_", tissue)
  
  if (!exists(expr_obj_name, inherits = TRUE)) {
    stop("Object not found in environment: ", expr_obj_name)
  }
  if (!exists(si_obj_name, inherits = TRUE)) {
    stop("Object not found in environment: ", si_obj_name)
  }
  
  expr <- get(expr_obj_name, inherits = TRUE)
  sample_info <- get(si_obj_name, inherits = TRUE)
  
  expr <- as.data.frame(expr)
  
  if (is.null(rownames(sample_info))) {
    if ("sample_id" %in% colnames(sample_info)) {
      rownames(sample_info) <- sample_info$sample_id
    } else {
      stop("sample_info must have rownames or a sample_id column.")
    }
  }
  
  if (is.null(rownames(expr))) {
    stop("Expression matrix must have rownames = sample IDs.")
  }
  
  common <- intersect(rownames(sample_info), rownames(expr))
  if (length(common) < min_matched) {
    stop("Too few matched samples for tissue = ", tissue, " (matched = ", length(common), ").")
  }
  
  list(
    tissue = tissue,
    expr = expr[common, , drop = FALSE],
    sample_info = sample_info[common, , drop = FALSE]
  )
}

# ------------------------------------------------------------
# Build RF models separately for each drug
# Outcome:
#   delta_PASI = PASI_followup - PASI_baseline
#
# Features inside RF:
#   - PASI_baseline
#   - top variable genes
#
# drug is NOT included as feature because data is split by drug
# ------------------------------------------------------------
build_rf_topgenes_by_drug <- function(features,
                                      sample_info,
                                      top_n = 50,
                                      baseline_week = 0,
                                      followup_week = 12,
                                      n_splits = 10,
                                      test_frac = 0.2,
                                      num_trees = 1000,
                                      seed = 42,
                                      importance_mode = "permutation") {
  
  set.seed(seed)
  
  stopifnot(is.matrix(features) || is.data.frame(features))
  features <- as.data.frame(features)
  
  if (is.null(rownames(features))) stop("features must have rownames = sample IDs.")
  if (is.null(rownames(sample_info))) stop("sample_info must have rownames = sample IDs.")
  
  if (!all(rownames(sample_info) %in% rownames(features))) {
    stop("Not all sample_info rownames exist in features rownames. Check alignment.")
  }
  
  # Select top variable genes ONCE at the full dataset level
  # (same pragmatic style as your earlier exploratory top-gene work)
  features_top <- as.data.frame(select_top_variable_features(features, n = top_n))
  
  si <- sample_info %>%
    dplyr::mutate(
      week = as.numeric(as.character(week)),
      PASI = as.numeric(as.character(PASI)),
      drug = tolower(trimws(as.character(drug)))
    ) %>%
    dplyr::mutate(sample_id = rownames(sample_info)) %>%
    dplyr::select(sample_id, patient, drug, week, PASI)
  
  baseline_tbl <- si %>%
    dplyr::filter(week == baseline_week) %>%
    dplyr::group_by(patient, drug) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(
      PASI_baseline = PASI,
      sample_id_baseline = sample_id
    )
  
  followup_tbl <- si %>%
    dplyr::filter(week == followup_week) %>%
    dplyr::group_by(patient, drug) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(PASI_followup = PASI)
  
  outcome_tbl <- baseline_tbl %>%
    dplyr::inner_join(followup_tbl, by = c("patient", "drug")) %>%
    dplyr::mutate(delta_PASI = PASI_followup - PASI_baseline) %>%
    dplyr::filter(!is.na(PASI_baseline), !is.na(PASI_followup), !is.na(delta_PASI))
  
  arm_results <- list()
  
  for (drug_nm in c("adalimumab", "ustekinumab")) {
    arm <- outcome_tbl %>% dplyr::filter(drug == drug_nm)
    
    if (nrow(arm) < 10) {
      warning(sprintf("Skipping %s: only %d patients", drug_nm, nrow(arm)))
      next
    }
    
    X <- features_top[arm$sample_id_baseline, , drop = FALSE]
    rownames(X) <- arm$patient
    
    split_perf <- lapply(seq_len(n_splits), function(s) {
      set.seed(seed + s)
      
      test_n <- max(1, floor(test_frac * nrow(arm)))
      test_i <- sample(seq_len(nrow(arm)), size = test_n)
      train_i <- setdiff(seq_len(nrow(arm)), test_i)
      
      arm_train <- arm[train_i, , drop = FALSE]
      arm_test  <- arm[test_i, , drop = FALSE]
      
      X_train <- X[train_i, , drop = FALSE]
      X_test  <- X[test_i, , drop = FALSE]
      
      train_df <- dplyr::bind_cols(
        data.frame(
          delta_PASI = arm_train$delta_PASI,
          PASI_baseline = arm_train$PASI_baseline
        ),
        as.data.frame(X_train)
      )
      
      test_df <- dplyr::bind_cols(
        data.frame(
          delta_PASI = arm_test$delta_PASI,
          PASI_baseline = arm_test$PASI_baseline
        ),
        as.data.frame(X_test)
      )
      
      rf_fit <- ranger(
        dependent.variable.name = "delta_PASI",
        data = train_df,
        num.trees = num_trees,
        importance = importance_mode,
        seed = seed + s
      )
      
      rf_pred <- predict(rf_fit, data = test_df)$predictions
      rf_rmse <- sqrt(mean((rf_pred - test_df$delta_PASI)^2))
      rf_cor  <- suppressWarnings(cor(rf_pred, test_df$delta_PASI, use = "complete.obs"))
      
      lm_fit <- lm(delta_PASI ~ PASI_baseline, data = arm_train)
      lm_pred <- predict(lm_fit, newdata = arm_test)
      lm_rmse <- sqrt(mean((lm_pred - arm_test$delta_PASI)^2))
      lm_cor  <- suppressWarnings(cor(lm_pred, arm_test$delta_PASI, use = "complete.obs"))
      
      data.frame(
        rf_rmse = rf_rmse,
        rf_cor = rf_cor,
        lm_rmse = lm_rmse,
        lm_cor = lm_cor
      )
    })
    
    sp <- dplyr::bind_rows(split_perf)
    
    train_df_full <- dplyr::bind_cols(
      data.frame(
        patient = arm$patient,
        sample_id_baseline = arm$sample_id_baseline,
        PASI_baseline = arm$PASI_baseline,
        PASI_followup = arm$PASI_followup,
        delta_PASI = arm$delta_PASI
      ),
      as.data.frame(X)
    )
    
    rf_final <- ranger(
      dependent.variable.name = "delta_PASI",
      data = train_df_full %>%
        dplyr::select(-patient, -sample_id_baseline, -PASI_followup),
      num.trees = num_trees,
      importance = importance_mode,
      seed = seed
    )
    
    imp <- sort(rf_final$variable.importance, decreasing = TRUE)
    
    arm_results[[drug_nm]] <- list(
      n_patients = nrow(arm),
      top_n = top_n,
      outcome_tbl = arm,
      train_df = train_df_full,
      model = rf_final,
      split_perf = sp,
      
      rf_rmse_mean = mean(sp$rf_rmse, na.rm = TRUE),
      rf_rmse_sd   = sd(sp$rf_rmse, na.rm = TRUE),
      rf_cor_mean  = mean(sp$rf_cor, na.rm = TRUE),
      rf_cor_sd    = sd(sp$rf_cor, na.rm = TRUE),
      
      lm_rmse_mean = mean(sp$lm_rmse, na.rm = TRUE),
      lm_rmse_sd   = sd(sp$lm_rmse, na.rm = TRUE),
      lm_cor_mean  = mean(sp$lm_cor, na.rm = TRUE),
      lm_cor_sd    = sd(sp$lm_cor, na.rm = TRUE),
      
      importance = data.frame(
        feature = names(imp),
        importance = as.numeric(imp),
        row.names = NULL
      ),
      features = colnames(X)
    )
  }
  
  arm_results
}

# ------------------------------------------------------------
# Run model for one tissue and one top_n
# ------------------------------------------------------------
run_topgenes_for_tissue_by_drug <- function(tissue,
                                            top_n = 50,
                                            seed = 42,
                                            baseline_week = 0,
                                            followup_week = 12,
                                            n_splits = 10,
                                            test_frac = 0.2,
                                            num_trees = 1000,
                                            importance_mode = "permutation") {
  
  dat <- load_tissue_expression_data(tissue = tissue)
  
  res <- build_rf_topgenes_by_drug(
    features = dat$expr,
    sample_info = dat$sample_info,
    top_n = top_n,
    baseline_week = baseline_week,
    followup_week = followup_week,
    n_splits = n_splits,
    test_frac = test_frac,
    num_trees = num_trees,
    seed = seed,
    importance_mode = importance_mode
  )
  
  cat("TISSUE:", tissue, "| Top genes:", top_n, "\n")
  
  for (drug_nm in names(res)) {
    arm <- res[[drug_nm]]
    
    cat(sprintf("\n  %s (n=%d)\n", drug_nm, arm$n_patients))
    
    cat(sprintf("  RF model: delta_PASI ~ PASI_baseline + top %d genes\n", top_n))
    cat(sprintf("    RMSE: %.2f +/- %.2f\n", arm$rf_rmse_mean, arm$rf_rmse_sd))
    cat(sprintf("    Correlation: %.2f +/- %.2f\n", arm$rf_cor_mean, arm$rf_cor_sd))
    
    cat("  Baseline LM: delta_PASI ~ PASI_baseline\n")
    cat(sprintf("    RMSE: %.2f +/- %.2f\n", arm$lm_rmse_mean, arm$lm_rmse_sd))
    cat(sprintf("    Correlation: %.2f +/- %.2f\n", arm$lm_cor_mean, arm$lm_cor_sd))
    
    cat("  Top 10 features:\n")
    print(head(arm$importance, 10))
  }
  
  invisible(res)
}

# ------------------------------------------------------------
# Run for all tissues and top gene sets
# ------------------------------------------------------------
top_gene_sizes <- c(50, 20, 10)

topgene_results <- list()

for (top_n in top_gene_sizes) {
  cat("\n========================================\n")
  cat("RUNNING TOP", top_n, "GENES\n")
  cat("========================================\n")
  
  topgene_results[[paste0("top", top_n)]] <- list(
    blood = run_topgenes_for_tissue_by_drug("blood", top_n = top_n),
    lesional = run_topgenes_for_tissue_by_drug("lesional", top_n = top_n),
    nonlesional = run_topgenes_for_tissue_by_drug("nonlesional", top_n = top_n)
  )
}