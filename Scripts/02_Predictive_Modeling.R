# ==============================================================================
# SCRIPT 2: HIERARCHICAL TreMs MODELING AND FORMULA EXTRACTION
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Functionality:
# 1. Loading LiDAR metrics and field data (TreMs).
# 2. Automated construction of ecological hierarchy (Macro, Forms, Groups, Types).
# 3. Multicollinearity analysis of LiDAR metrics (cut-off threshold = 80%).
# 4. Feature Selection based on Random Forest Importance (Top 5 predictors).
# 5. 10-fold Cross-Validation for 4 models (Random Forest, XGBoost, Linear Model, ZINB).
# 6. Calculation of mean R2, RMSE, and Significance (p-value).
# 7. DIRECT EXTRACTION OF MATHEMATICAL FORMULAS (LM and ZINB) without post-processing.
# 8. Automated export of:
#    - "Complete_Model_Results_with_Pvalues.xlsx"
#    - "Best_Models_CV_Summary.txt" (with Exact Formulas)
#    - Global comparative plot at 600 DPI ("Comparison_4_Models.png")
#    - Individual plots for each target ("Single_Target_Plots/")
# ==============================================================================

rm(list = ls())
gc()

# 1. ENVIRONMENT AND PACKAGES ====
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, MASS, pscl, randomForest, xgboost, caret, dplyr, tidyr, ggplot2, writexl, stringr)

cat("\n===================================================================\n")
cat("    MASTER PIPELINE: HIERARCHICAL TreMs MODELING & FORMULAS    \n")
cat("===================================================================\n\n")

# Output folder selection
dir_output <- "../Outputs"
dir_plot_perf <- file.path(dir_output, "Performance_Plots")
dir_plot_single <- file.path(dir_plot_perf, "Single_Target_Plots")
dir.create(dir_plot_perf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plot_single, showWarnings = FALSE, recursive = TRUE)

time_start_global <- Sys.time()
time_start_prep <- Sys.time()

COL_ID_FIELD <- "treeID"

# 2. DATA LOADING ====
cat("\n>>> PHASE 1: Loading Data Files <<<\n")
file_csv_lidar <- "../Data/Extracted_TLS_Metrics_Full.csv"
file_excel_field <- "../Data/TreMs_Field_Surveys.xlsx"

if(!file.exists(file_csv_lidar)) stop("ERROR: Extracted_TLS_Metrics_Full.csv not found in ../Data/")
if(!file.exists(file_excel_field)) stop("ERROR: TreMs_Field_Surveys.xlsx not found in ../Data/")
df_lidar <- read.csv(file_csv_lidar, stringsAsFactors = FALSE)
df_field <- read_excel(file_excel_field)

# Standardize ID column name
if ("treeID" %in% names(df_lidar)) {
  df_lidar <- df_lidar %>% rename(Tree_ID = treeID)
} else if (names(df_lidar)[1] != "Tree_ID") {
  names(df_lidar)[1] <- "Tree_ID"
}
df_lidar$Tree_ID <- as.character(df_lidar$Tree_ID)

# Search for field ID column
found_id_col <- names(df_field)[grep("^treeID$|^ID$|^ID_Pianta$|^ID_Albero$", names(df_field), ignore.case = TRUE)[1]]
if (is.na(found_id_col)) found_id_col <- names(df_field)[1]
df_field <- df_field %>% rename(!!COL_ID_FIELD := all_of(found_id_col)) %>% mutate(!!COL_ID_FIELD := as.character(get(COL_ID_FIELD)))

# 3. DATA WRANGLING: ECOLOGICAL TreMs HIERARCHY ====
cat("\n>>> PHASE 2: Automated Construction of Ecological Hierarchy <<<\n")

# Identify individual TreMs columns (2-letter codes + numbers, e.g., CV11, DE12)
trems_64 <- grep("^[A-Z]{2}[0-9]+", names(df_field), value = TRUE)
if (length(trems_64) == 0) stop("ERROR: No TreMs columns (e.g., CV11, DE12) found in the Excel file.")

df_field <- df_field %>%
  mutate(across(all_of(trems_64), as.numeric)) %>%
  mutate(across(all_of(trems_64), ~ replace_na(., 0)))

# A. Macro Scale
df_field$Macro_TreMs_Abundance <- rowSums(df_field[, trems_64], na.rm = TRUE)
df_field$Macro_TreMs_Richness  <- rowSums(df_field[, trems_64] > 0, na.rm = TRUE)

# B. Forms (First 2 letters, e.g., CV, IN, DE, EP)
for (fp in unique(substr(trems_64, 1, 2))) {
  col_form <- grep(paste0("^", fp), names(df_field), value = TRUE)
  df_field[[paste0("Form_", fp)]] <- rowSums(df_field[, col_form, drop = FALSE], na.rm = TRUE)
}

# C. Groups (First 3 letters/numbers, e.g., CV1, CV2, DE1)
for (gp in unique(substr(trems_64, 1, 3))) {
  col_group <- grep(paste0("^", gp), names(df_field), value = TRUE)
  df_field[[paste0("Group_", gp)]] <- rowSums(df_field[, col_group, drop = FALSE], na.rm = TRUE)
}

ALL_TARGETS <- c("Macro_TreMs_Abundance", "Macro_TreMs_Richness", 
                 grep("^Form_|^Group_", names(df_field), value = TRUE), 
                 trems_64)

# Merge Dataset
df_all <- inner_join(df_lidar, df_field %>% dplyr::select(all_of(c(COL_ID_FIELD, ALL_TARGETS))), 
                     by = setNames(COL_ID_FIELD, "Tree_ID"))

df_model <- df_all %>% dplyr::select(-any_of(c("Tree_ID", COL_ID_FIELD))) %>% na.omit()
cat("Dataset successfully merged. Total trees analyzed:", nrow(df_model), "\n")

# 4. MULTICOLLINEARITY ANALYSIS ====
cat("\n>>> PHASE 3: Multicollinearity Analysis of LiDAR Metrics <<<\n")
present_lidar_names <- setdiff(names(df_model), ALL_TARGETS)
variances <- sapply(df_model[, present_lidar_names, drop = FALSE], var, na.rm = TRUE)
valid_names <- present_lidar_names[variances > 0 & !is.na(variances)]

cat("Initial LiDAR metrics with variance > 0:", length(valid_names), "\n")

# Correlation matrix
cor_mat <- cor(df_model[, valid_names], use = "complete.obs")
high_corr_idx <- findCorrelation(cor_mat, cutoff = 0.80)

if (length(high_corr_idx) > 0) {
  removed_names <- valid_names[high_corr_idx]
  valid_names <- valid_names[-high_corr_idx]
  cat(sprintf("[!] Discarded %d redundant metrics (|correlation| > 80%%):\n", length(removed_names)))
  cat(paste(removed_names, collapse = ", "), "\n")
} else {
  cat("No metrics exceeded the 80% correlation threshold.\n")
}
cat("=> Final non-redundant LiDAR metrics:", length(valid_names), "\n")

time_end_prep <- Sys.time()

# 5. AUXILIARY FUNCTIONS (CROSS-VALIDATION & FORMULA EXTRACTION) ====

# Function to extract the exact mathematical formula
get_math_formula <- function(model, target_name = "Y", model_type = "LM") {
  if (model_type == "LM" && !is.null(model)) {
    co <- tryCatch(coef(model), error = function(e) NULL)
    if (is.null(co)) return("Error extracting linear coefficients")
    co <- co[!is.na(co)]
    clean_names <- gsub("`", "", names(co))
    terms <- sapply(2:length(co), function(j) sprintf("%+.4f*%s", co[j], clean_names[j]))
    return(sprintf("%s = %.4f %s", target_name, co[1], paste(terms, collapse = " ")))
    
  } else if (model_type == "ZINB" && !is.null(model)) {
    co_count <- tryCatch(model$coefficients$count, error = function(e) NULL)
    co_zero  <- tryCatch(model$coefficients$zero, error = function(e) NULL)
    if (is.null(co_count) || is.null(co_zero)) return("Error extracting ZINB coefficients")
    
    co_count <- co_count[!is.na(co_count)]
    co_zero  <- co_zero[!is.na(co_zero)]
    
    names_c <- gsub("`", "", names(co_count))
    term_count <- sapply(2:length(co_count), function(j) sprintf("%+.4f*%s", co_count[j], names_c[j]))
    str_count <- sprintf("Count ln(%s) = %.4f %s", target_name, co_count[1], paste(term_count, collapse = " "))
    
    names_z <- gsub("`", "", names(co_zero))
    term_zero <- sapply(2:length(co_zero), function(j) sprintf("%+.4f*%s", co_zero[j], names_z[j]))
    str_zero <- sprintf(" | Zeroes logit(p) = %.4f %s", co_zero[1], paste(term_zero, collapse = " "))
    
    return(paste0(str_count, str_zero))
  } else {
    return("Non-parametric (Machine Learning)")
  }
}

# Function for comparative 10-Fold Cross-Validation
cv_total_comparison <- function(data_cv, target_name) {
  set.seed(123)
  folds <- createFolds(data_cv[[target_name]], k = 10)
  model_names <- c("RF", "XGB", "LM", "ZINB")
  
  results_r2   <- matrix(NA, nrow = 10, ncol = length(model_names), dimnames = list(NULL, model_names))
  results_rmse <- matrix(NA, nrow = 10, ncol = length(model_names), dimnames = list(NULL, model_names))
  
  for (i in 1:10) {
    train_data <- data_cv[-folds[[i]], ]; test_data <- data_cv[folds[[i]], ]
    
    # 1. Random Forest
    m_rf <- tryCatch(randomForest(as.formula(paste(target_name, "~ .")), data = train_data, ntree = 150), error = function(e) NULL)
    
    # 2. XGBoost
    m_xgb <- tryCatch({
      X_tr <- as.matrix(train_data %>% dplyr::select(-all_of(target_name)))
      xgb.train(params = list(objective = "count:poisson", eta = 0.1), 
                data = xgb.DMatrix(X_tr, label = train_data[[target_name]]), 
                nrounds = 50, verbose = 0)
    }, error = function(e) NULL)
    
    # 3. Linear Model
    m_lm <- tryCatch(lm(as.formula(paste(target_name, "~ .")), data = train_data), error = function(e) NULL)
    
    # 4. ZINB
    top_var_zero <- names(train_data)[2]
    m_zinb <- tryCatch(zeroinfl(as.formula(paste(target_name, "~ . |", top_var_zero)), data = train_data, dist = "negbin"), error = function(e) NULL)
    
    # Internal metric calculation function
    calc_metrics <- function(mod, test, m_type = "default") {
      if (is.null(mod)) return(c(NA, NA))
      p <- tryCatch({
        if (m_type == "xgb") predict(mod, xgb.DMatrix(as.matrix(test %>% dplyr::select(-all_of(target_name)))))
        else pmax(0, predict(mod, test, type = "response"))
      }, error = function(e) return(c(NA, NA)))
      
      if (any(is.na(p)) || length(unique(p)) <= 1 || sd(p, na.rm = TRUE) < 0.0001) return(c(NA, NA))
      
      actual <- test[[target_name]]
      r2_val <- cor(actual, p)^2
      rmse_val <- sqrt(mean((actual - p)^2))
      return(c(r2_val, rmse_val))
    }
    
    m_rf_res   <- calc_metrics(m_rf, test_data)
    results_r2[i, "RF"] <- m_rf_res[1]; results_rmse[i, "RF"] <- m_rf_res[2]
    
    m_xgb_res  <- calc_metrics(m_xgb, test_data, "xgb")
    results_r2[i, "XGB"] <- m_xgb_res[1]; results_rmse[i, "XGB"] <- m_xgb_res[2]
    
    m_lm_res   <- calc_metrics(m_lm, test_data)
    results_r2[i, "LM"] <- m_lm_res[1]; results_rmse[i, "LM"] <- m_lm_res[2]
    
    m_zinb_res <- calc_metrics(m_zinb, test_data, "stats")
    results_r2[i, "ZINB"] <- m_zinb_res[1]; results_rmse[i, "ZINB"] <- m_zinb_res[2]
  }
  
  return(list(
    R2   = colMeans(results_r2, na.rm = TRUE),
    RMSE = colMeans(results_rmse, na.rm = TRUE)
  ))
}

# 6. MASTER LOOP FOR MODELING & FORMULA EXTRACTION ====

time_start_training <- Sys.time()

cat("\n>>> PHASE 4: 10-Fold CV Execution & Mathematical Formula Extraction <<<\n")
final_comparison_table <- data.frame()

for (target in ALL_TARGETS) {
  if (!target %in% names(df_model) || sum(df_model[[target]] > 0) < 15) next
  cat(sprintf("\nAnalyzing 4 models for: %s (Trees with presence: %d)...\n", target, sum(df_model[[target]] > 0)))
  
  # A. Feature Selection via Random Forest
  rf_sel <- tryCatch(randomForest(as.formula(paste(target, "~ .")), data = df_model[, c(target, valid_names)], ntree = 100), error = function(e) NULL)
  if (is.null(rf_sel)) next
  
  imp_mat <- importance(rf_sel)
  top_5 <- rownames(imp_mat)[order(imp_mat[, 1], decreasing = TRUE)][1:min(5, nrow(imp_mat))]
  
  # Compute Spearman Correlations with Top 5
  correlations_list <- sapply(top_5, function(var) {
    valore_cor <- cor(df_model[[target]], df_model[[var]], method = "spearman", use = "complete.obs")
    sprintf("%s (%+.3f)", var, valore_cor)
  })
  string_correlations <- paste(correlations_list, collapse = ", ")
  
  # B. 10-Fold Cross-Validation
  res_cv_target <- cv_total_comparison(df_model[, c(target, top_5)], target)
  
  # C. Fit models on full dataset for p-values and exact formulas
  full_data <- df_model[, c(target, top_5)]
  full_formula <- as.formula(paste(target, "~ ."))
  top_v <- top_5[1]
  formula_zinb <- as.formula(paste(target, "~ . |", top_v))
  
  m_lm_full   <- tryCatch(lm(full_formula, data = full_data), error = function(e) NULL)
  m_zinb_full <- tryCatch(zeroinfl(formula_zinb, data = full_data, dist = "negbin"), error = function(e) NULL)
  m_zinb_null <- tryCatch(zeroinfl(as.formula(paste(target, "~ 1 | 1")), data = full_data, dist = "negbin"), error = function(e) NULL)
  
  # Calculate P-Values
  p_val_lm <- NA
  if (!is.null(m_lm_full)) {
    f_stat <- summary(m_lm_full)$fstatistic
    if (!is.null(f_stat)) p_val_lm <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  }
  
  p_val_zinb <- NA
  if (!is.null(m_zinb_full) && !is.null(m_zinb_null)) {
    ll_full <- logLik(m_zinb_full); ll_null <- logLik(m_zinb_null)
    df_diff <- attr(ll_full, "df") - attr(ll_null, "df")
    if (df_diff > 0) p_val_zinb <- pchisq(2 * as.numeric(ll_full - ll_null), df = df_diff, lower.tail = FALSE)
  }
  
  # D. Construction of Results DataFrame with Exact Formulas
  df_res <- data.frame(
    Target = target,
    Model = names(res_cv_target$R2),
    R2 = as.numeric(res_cv_target$R2),
    RMSE = as.numeric(res_cv_target$RMSE),
    P_Value = NA_real_,
    Used_Variables = paste(top_5, collapse = ", "),
    Spearman_Correlations = string_correlations,
    Exact_Formula = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Assign P-Values
  df_res$P_Value[df_res$Model == "LM"]   <- p_val_lm
  df_res$P_Value[df_res$Model == "ZINB"] <- p_val_zinb
  
  # Assign Exact Mathematical Formulas
  df_res$Exact_Formula[df_res$Model == "LM"]   <- get_math_formula(m_lm_full, target, "LM")
  df_res$Exact_Formula[df_res$Model == "ZINB"] <- get_math_formula(m_zinb_full, target, "ZINB")
  df_res$Exact_Formula[df_res$Model %in% c("RF", "XGB")] <- "Non-parametric (Machine Learning)"
  
  # Print summary to console
  res_print <- df_res %>% arrange(desc(R2)) %>% filter(!is.na(R2))
  res_print$R2 <- round(res_print$R2, 3)
  res_print$RMSE <- round(res_print$RMSE, 3)
  res_print$Significance <- ifelse(is.na(res_print$P_Value), "NA (ML)", 
                                   ifelse(res_print$P_Value < 0.001, "< 0.001 ***", 
                                          ifelse(res_print$P_Value < 0.05, paste0(round(res_print$P_Value, 3), " *"), 
                                                 as.character(round(res_print$P_Value, 3)))))
  print(res_print[, c("Target", "Model", "R2", "RMSE", "Significance")], row.names = FALSE)
  
  final_comparison_table <- rbind(final_comparison_table, df_res)
}

time_end_training <- Sys.time()
time_start_prediction <- Sys.time()
if (exists("rf_sel")) {
  # Prediction simulation on the entire dataset to calculate pure computational cost
  dummy_pred <- predict(rf_sel, df_model) 
}
time_end_prediction <- Sys.time()

# 7. SAVING RESULTS (EXCEL & TXT REPORT) ====
cat("\n>>> PHASE 5: Saving Complete Results & Formulas <<<\n")

time_start_reporting <- Sys.time()

# A. Export Excel Table
file_excel_out <- file.path(dir_output, "Complete_Model_Results_with_Pvalues.xlsx")
write_xlsx(final_comparison_table, file_excel_out)
cat(sprintf("[✔] Excel file saved: %s\n", basename(file_excel_out)))

# B. Generate TXT Report with Exact Formulas
best_models <- final_comparison_table %>%
  group_by(Target) %>%
  slice_max(R2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Category_Order = case_when(
      grepl("^Macro_", Target) ~ 1,
      grepl("^Form_", Target)  ~ 2,
      grepl("^Group_", Target) ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(Category_Order, Target)

txt_lines <- unlist(lapply(1:nrow(best_models), function(i) {
  p_val_text <- ifelse(is.na(best_models$P_Value[i]), "NA (ML)", 
                       ifelse(best_models$P_Value[i] < 0.001, "< 0.001", sprintf("%.4f", best_models$P_Value[i])))
  
  sprintf("TARGET: %-25s | BEST MODEL: %-15s | R2 CV: %.3f | RMSE CV: %.3f | P-VALUE: %-7s | VARIABLES: %s | CORRELATIONS: %s | EXACT FORMULA: %s", 
          best_models$Target[i], 
          best_models$Model[i], 
          best_models$R2[i], 
          best_models$RMSE[i],
          p_val_text,
          best_models$Used_Variables[i],
          best_models$Spearman_Correlations[i],
          best_models$Exact_Formula[i])
}))

file_txt_out <- file.path(dir_output, "Best_Models_CV_Summary.txt")
writeLines(txt_lines, file_txt_out)
cat(sprintf("[✔] TXT file with exact formulas saved: %s\n", basename(file_txt_out)))

# 8. SCIENTIFIC PLOTS GENERATION ====
cat("\n>>> PHASE 6: Generation of Comparative Plots (600 DPI) <<<\n")

# A. Global Comparative Plot for 4 Models
df_plot <- final_comparison_table %>% 
  group_by(Target) %>% 
  filter(max(R2, na.rm = TRUE) >= 0.09) %>%
  ungroup() %>%
  mutate(Sig_Label = case_when(
    is.na(P_Value) ~ "", P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ ""
  ))

if (nrow(df_plot) > 0) {
  plot_4_models <- ggplot(df_plot, aes(x = reorder(Target, R2), y = R2, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
    geom_text(aes(label = Sig_Label), position = position_dodge(width = 0.9), hjust = -0.2, size = 4, color = "black") +
    coord_flip() + 
    theme_minimal() + 
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) + 
    labs(
      title = "Comparison of 4 Predictive Models (10-fold CV)", 
      subtitle = "Asterisks denote statistical significance (* p<0.05, ** p<0.01, *** p<0.001; NA for non-parametric ML)",
      x = "Tree-Related Microhabitats (Target Variables)", y = "Mean Cross-Validated R-Squared"
    ) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold", size = 11),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
  
  file_plot_global <- file.path(dir_plot_perf, "Comparison_4_Models.png")
  ggsave(file_plot_global, plot = plot_4_models, width = 12, height = 10, dpi = 600, bg = "white")
  cat(sprintf("[✔] Global comparative plot saved: %s\n", basename(file_plot_global)))
}

# B. Individual Plots for each Target
cat("Generating individual plots for each target...\n")
for (tgt in unique(final_comparison_table$Target)) {
  df_single <- final_comparison_table %>% filter(Target == tgt) %>%
    mutate(Sig_Label = case_when(
      is.na(P_Value) ~ "", P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ ""
    ))
  
  max_r2 <- max(df_single$R2, na.rm = TRUE)
  best_m <- df_single$Model[which.max(df_single$R2)]
  
  p_single <- ggplot(df_single, aes(x = reorder(Model, R2), y = R2, fill = Model)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_text(aes(label = Sig_Label), hjust = -0.2, size = 6, color = "black") +
    coord_flip() + theme_minimal() + scale_fill_viridis_d(option = "turbo") +
    labs(
      title = paste("Model Performance:", tgt),
      subtitle = paste("Best Performing Model:", best_m),
      x = "Predictive Algorithm", y = expression("10-fold CV R"^2)
    ) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, max(max_r2 * 1.25, 0.05)))
  
  ggsave(file.path(dir_plot_single, paste0("Plot_", tgt, ".png")), plot = p_single, width = 8, height = 5, dpi = 300, bg = "white")
}

time_end_reporting <- Sys.time()
time_end_global <- Sys.time()

cat("\n===================================================================\n")
cat("🎉 MODELING AND FORMULA EXTRACTION COMPLETED SUCCESSFULLY!\n")
cat("All results and plots have been saved in:", dir_output, "\n")
cat("===================================================================\n")

# --- TIME CALCULATION ---
duration_prep       <- difftime(time_end_prep, time_start_prep, units = "secs")
duration_training   <- difftime(time_end_training, time_start_training, units = "mins")
duration_prediction <- difftime(time_end_prediction, time_start_prediction, units = "secs")
duration_reporting  <- difftime(time_end_reporting, time_start_reporting, units = "secs")
duration_total      <- difftime(time_end_global, time_start_global, units = "mins")

# --- PRINT CONSOLE TABLE ---
cat("\n=================================================================================\n")
cat("📊 COMPUTATIONAL EFFICIENCY REPORT\n")
cat("=================================================================================\n")
cat(sprintf(" %-65s %6.2f seconds\n", "Data prep, Multicollinearity analysis, top 5 Predictors:", as.numeric(duration_prep)))
cat(sprintf(" %-65s %6.2f minutes\n",  "Model Training:", as.numeric(duration_training)))
cat(sprintf(" %-65s %6.4f seconds\n", "Prediction time (per dataset):", as.numeric(duration_prediction)))
cat(sprintf(" %-65s %6.2f seconds\n", "Reporting, Excel & Plotting:", as.numeric(duration_reporting)))
cat(" ---------------------------------------------------------------------------------\n")
cat(sprintf(" %-65s %6.2f minutes\n",  "TOTAL R PIPELINE TIME:", as.numeric(duration_total)))
cat("=================================================================================\n\n")
