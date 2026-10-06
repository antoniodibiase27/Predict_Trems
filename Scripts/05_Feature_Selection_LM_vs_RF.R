# ==============================================================================
# FEATURE SELECTION TEST: Linear Model vs Random Forest
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Description: 
# This script evaluates the impact of the LiDAR metrics selection method 
# on the predictive performance of 4 models (RF, XGB, LM, ZINB) for TreMs estimation. 
# Specifically:
# 1. Selects the top 5 variables using the p-values of a Linear Model (LM).
# 2. Calculates R2 and RMSE via 10-fold Cross-Validation.
# 3. Automatically compares these results with those of a previous 
#    pipeline (where selection was based on Random Forest importance).
# 
# OUTPUT:
# Excel result tables, textual reports, and a dumbbell plot 
# showing which feature selection method maximizes the R2.
# ==============================================================================
rm(list = ls())
gc()

# 1. ENVIRONMENT AND PACKAGE LOADING
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, MASS, pscl, randomForest, xgboost, caret, dplyr, tidyr, ggplot2, writexl, ggalt)

# 2. FOLDER DEFINITION AND FILE LOADING
cat("\n--- PHASE 1: SETUP AND DATA LOADING ---\n")
dir_output <- "../Outputs"
setwd(dir_output)
dir.create("Plot_Performance_LM_Test", showWarnings = FALSE)
dir.create("Plot_Performance", showWarnings = FALSE)

COL_ID_TERRA <- "treeID"

df_lidar <- read.csv("../Data/Extracted_TLS_Metrics_Full.csv") %>% rename(ID_Albero = treeID) %>% mutate(ID_Albero = as.character(ID_Albero))
df_terra <- read_excel("../Data/TreMs_Field_Surveys.xlsx") %>% mutate(!!COL_ID_TERRA := as.character(get(COL_ID_TERRA)))
df_rf <- read_excel("../Data/Results_Model_Pvalue.xlsx")

# ==============================================================================
# 3. DATA WRANGLING AND DATA PREPARATION
# ==============================================================================
cat("\n--- PHASE 2: DATA PREPARATION ---\n")
tipi_64 <- grep("^[A-Z]{2}[0-9]+", names(df_terra), value = TRUE)
df_terra <- df_terra %>% mutate(across(all_of(tipi_64), as.numeric)) %>% mutate(across(all_of(tipi_64), ~ replace_na(., 0)))
df_terra$Macro_TreMs_Abundance <- rowSums(df_terra[, tipi_64], na.rm = TRUE)
df_terra$Macro_TreMs_Richness <- rowSums(df_terra[, tipi_64] > 0, na.rm = TRUE)

for(fp in unique(substr(tipi_64, 1, 2))) { df_terra[[paste0("Form_", fp)]] <- rowSums(df_terra[, grep(paste0("^", fp), names(df_terra)), drop=FALSE], na.rm=TRUE) }
for(gp in unique(substr(tipi_64, 1, 3))) { df_terra[[paste0("Group_", gp)]] <- rowSums(df_terra[, grep(paste0("^", gp), names(df_terra)), drop=FALSE], na.rm=TRUE) }

TUTTI_I_TARGETS <- c("Macro_TreMs_Abundance", "Macro_TreMs_Richness", grep("Form_|Group_", names(df_terra), value=T), tipi_64)

df_all <- inner_join(df_lidar, df_terra %>% dplyr::select(all_of(c(COL_ID_TERRA, TUTTI_I_TARGETS))), by = setNames(COL_ID_TERRA, "ID_Albero"))
df_model <- df_all %>% dplyr::select(-any_of(c("ID_Albero", COL_ID_TERRA))) %>% na.omit()

nomi_lidar_presenti <- setdiff(names(df_model), TUTTI_I_TARGETS)
varianze <- sapply(df_model[, nomi_lidar_presenti], var)
nomi_validi <- nomi_lidar_presenti[varianze > 0 & !is.na(varianze)]

# ==============================================================================
# 4. MULTICOLLINEARITY ANALYSIS
# ==============================================================================
cat("\n--- PHASE 3: MULTICOLLINEARITY ANALYSIS ---\n")
num_iniziali <- length(nomi_validi)
cat("Number of initial LiDAR metrics analyzed:", num_iniziali, "\n")

cor_mat <- cor(df_model[, nomi_validi], use = "complete.obs")
high_corr_idx <- findCorrelation(cor_mat, cutoff = 0.80)

if(length(high_corr_idx) > 0) {
  nomi_rimossi <- nomi_validi[high_corr_idx]
  nomi_validi <- nomi_validi[-high_corr_idx]
  cat("\n[!] MULTICOLLINEARITY FOUND:\n")
  cat("- Number of discarded metrics (correlation > 80%):", length(nomi_rimossi), "\n")
  cat("=> Number of final usable metrics:", length(nomi_validi), "\n")
} else {
  cat("\nNo variables exceeded the 80% correlation threshold.\n")
}
cat("--------------------------------------------\n")

# ==============================================================================
# 5. CROSS-VALIDATION FUNCTION
# ==============================================================================
cv_total_comparison <- function(data_cv, target_name) {
  set.seed(123)
  folds <- createFolds(data_cv[[target_name]], k = 10)
  nomi_modelli <- c("RF", "XGB", "LM", "ZINB")
  
  risultati_r2 <- matrix(NA, nrow = 10, ncol = length(nomi_modelli), dimnames = list(NULL, nomi_modelli))
  risultati_rmse <- matrix(NA, nrow = 10, ncol = length(nomi_modelli), dimnames = list(NULL, nomi_modelli))
  
  for(i in 1:10) {
    train_data <- data_cv[-folds[[i]], ]; test_data <- data_cv[folds[[i]], ]
    
    m_rf <- tryCatch(randomForest(as.formula(paste(target_name, "~ .")), data=train_data, ntree=150), error=function(e) NULL)
    m_xgb <- tryCatch({
      X_tr <- as.matrix(train_data %>% dplyr::select(-all_of(target_name)))
      xgb.train(params=list(objective="count:poisson", eta=0.1), data=xgb.DMatrix(X_tr, label=train_data[[target_name]]), nrounds=50, verbose=0)
    }, error=function(e) NULL)
    m_lm <- tryCatch(lm(as.formula(paste(target_name, "~ .")), data=train_data), error=function(e) NULL)
    m_zinb <- tryCatch(zeroinfl(as.formula(paste(target_name, "~ . |", names(train_data)[2])), data=train_data, dist="negbin"), error=function(e) NULL)
    
    calc_metrics <- function(mod, test, m_type="default") {
      if(is.null(mod)) return(c(NA, NA))
      p <- tryCatch({
        if(m_type=="xgb") predict(mod, xgb.DMatrix(as.matrix(test %>% dplyr::select(-all_of(target_name)))))
        else pmax(0, predict(mod, test, type="response"))
      }, error=function(e) return(c(NA, NA)))
      
      if(any(is.na(p)) || length(unique(p)) <= 1 || sd(p, na.rm=T) < 0.0001) return(c(NA, NA))
      
      actual <- test[[target_name]]
      r2_val <- cor(actual, p)^2
      rmse_val <- sqrt(mean((actual - p)^2))
      return(c(r2_val, rmse_val))
    }
    
    metrics_rf <- calc_metrics(m_rf, test_data)
    risultati_r2[i, "RF"] <- metrics_rf[1]; risultati_rmse[i, "RF"] <- metrics_rf[2]
    metrics_xgb <- calc_metrics(m_xgb, test_data, "xgb")
    risultati_r2[i, "XGB"] <- metrics_xgb[1]; risultati_rmse[i, "XGB"] <- metrics_xgb[2]
    metrics_lm <- calc_metrics(m_lm, test_data)
    risultati_r2[i, "LM"] <- metrics_lm[1]; risultati_rmse[i, "LM"] <- metrics_lm[2]
    metrics_zinb <- calc_metrics(m_zinb, test_data, "stats")
    risultati_r2[i, "ZINB"] <- metrics_zinb[1]; risultati_rmse[i, "ZINB"] <- metrics_zinb[2]
  }
  
  return(list(R2 = colMeans(risultati_r2, na.rm = TRUE), RMSE = colMeans(risultati_rmse, na.rm = TRUE)))
}

# ==============================================================================
# 6. MASTER LOOP (SELECTION WITH LINEAR MODEL)
# ==============================================================================
cat("\n--- PHASE 4: EXECUTING MODELS WITH LM SELECTION ---\n")
tabella_confronto_finale <- data.frame()

for (target in TUTTI_I_TARGETS) {
  if(!target %in% names(df_model) || sum(df_model[[target]] > 0) < 15) next
  cat(sprintf("\nAnalyzing 4 models: %s... \n", target))
  
  lm_sel <- tryCatch(lm(as.formula(paste(target, "~ .")), data = df_model[, c(target, nomi_validi)]), error=function(e) NULL)
  if(is.null(lm_sel)) next
  
  summary_lm <- summary(lm_sel)$coefficients
  if(nrow(summary_lm) < 2) next 
  
  p_values <- summary_lm[-1, 4] 
  top_5 <- names(sort(p_values, decreasing = FALSE))[1:min(5, length(p_values))]
  
  correlazioni_lista <- sapply(top_5, function(var) {
    valore_cor <- cor(df_model[[target]], df_model[[var]], method = "spearman", use = "complete.obs")
    sprintf("(%+.3f)", valore_cor) %>% paste(var, .)
  })
  stringa_correlazioni <- paste(correlazioni_lista, collapse = ", ")
  
  res_cv_target <- cv_total_comparison(df_model[, c(target, top_5)], target)
  stringa_variabili <- paste(top_5, collapse = " + ")
  
  df_res <- data.frame(
    Target = target, 
    Model = names(res_cv_target$R2), 
    R2 = as.numeric(res_cv_target$R2), 
    RMSE = as.numeric(res_cv_target$RMSE),
    P_Value = NA,
    Variabili_Usate = paste(top_5, collapse = ", "),
    Correlazioni_Spearman = stringa_correlazioni,
    Formula = NA
  )
  
  df_res$Formula[df_res$Model == "LM"] <- paste(target, "=", stringa_variabili)
  df_res$Formula[df_res$Model == "ZINB"] <- paste("Count: ln(", target, ") =", stringa_variabili, "| Zeros (Logit):", top_5[1])
  df_res$Formula[df_res$Model %in% c("RF", "XGB")] <- "Non-parametric"
  
  dati_completi <- df_model[, c(target, top_5)]
  formula_completa <- as.formula(paste(target, "~ ."))
  top_v <- top_5[1] 
  formula_zinb <- as.formula(paste(target, "~ . |", top_v))
  
  m_lm_full <- tryCatch(lm(formula_completa, data=dati_completi), error=function(e) NULL)
  if(!is.null(m_lm_full)) {
    f_stat <- summary(m_lm_full)$fstatistic
    if(!is.null(f_stat)) df_res$P_Value[df_res$Model == "LM"] <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE)
  }
  
  m_zinb_full <- tryCatch(zeroinfl(formula_zinb, data=dati_completi, dist="negbin"), error=function(e) NULL)
  m_zinb_null <- tryCatch(zeroinfl(as.formula(paste(target, "~ 1 | 1")), data=dati_completi, dist="negbin"), error=function(e) NULL)
  if(!is.null(m_zinb_full) && !is.null(m_zinb_null)) {
    ll_full <- logLik(m_zinb_full); ll_null <- logLik(m_zinb_null)
    df_diff <- attr(ll_full, "df") - attr(ll_null, "df")
    if(df_diff > 0) df_res$P_Value[df_res$Model == "ZINB"] <- pchisq(2 * as.numeric(ll_full - ll_null), df = df_diff, lower.tail = FALSE)
  }
  
  res_print <- df_res %>% arrange(desc(R2)) %>% filter(!is.na(R2))
  res_print$R2 <- round(res_print$R2, 3)
  res_print$RMSE <- round(res_print$RMSE, 3)
  res_print$Significativita <- ifelse(is.na(res_print$P_Value), "NA (ML)", 
                                      ifelse(res_print$P_Value < 0.001, "< 0.001 ***", 
                                             ifelse(res_print$P_Value < 0.05, paste0(round(res_print$P_Value, 3), " *"), 
                                                    as.character(round(res_print$P_Value, 3)))))
  
  print(res_print[, c("Target", "Model", "R2", "RMSE", "Significativita")], row.names = FALSE)
  tabella_confronto_finale <- rbind(tabella_confronto_finale, df_res)
}

# ==============================================================================
# 7. LM RESULTS EXPORT AND PLOT
# ==============================================================================
cat("\n--- PHASE 5: DATA EXPORT AND LM PLOT ---\n")
write_xlsx(tabella_confronto_finale, "Full_Results_LM_Selection.xlsx")

df_plot <- tabella_confronto_finale %>% 
  group_by(Target) %>% 
  filter(max(R2, na.rm=TRUE) >= 0.0900) %>% 
  ungroup() %>%
  mutate(Sig_Label = case_when(
    is.na(P_Value) ~ "", P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", TRUE ~ ""
  ))

plot_4_modelli <- ggplot(df_plot, aes(x = reorder(Target, R2), y = R2, fill = Model)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.2) +
  geom_text(aes(label = Sig_Label), position = position_dodge(width = 0.9), hjust = -0.2, size = 4, color = "black") +
  coord_flip() + 
  theme_minimal() + 
  scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) + 
  labs(
    title = "Comparison 4 models (LM Selection Test)", 
    subtitle = "Top 5 features selected via Linear Model p-values instead of Random Forest",
    x = "TreMs (Target variables)", y = "Mean CV R-Squared"
  ) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 11),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave("Plot_Performance_LM_Test/Comparison_4_Models_LM_Selection.png", plot=plot_4_modelli, width=12, height=10, dpi=600)

vincitori <- tabella_confronto_finale %>% group_by(Target) %>% slice_max(R2, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(Categoria_Ordine = case_when(grepl("^Macro_", Target) ~ 1, grepl("^Form_", Target) ~ 2, grepl("^Group_", Target) ~ 3, TRUE ~ 4)) %>%
  arrange(Categoria_Ordine, Target)

writeLines(unlist(lapply(1:nrow(vincitori), function(i) {
  p_val_text <- ifelse(is.na(vincitori$P_Value[i]), "NA", sprintf("%.4f", vincitori$P_Value[i]))
  sprintf("TARGET: %-25s | BEST MODEL: %-15s | R2 CV: %.3f | RMSE CV: %.3f | P-VALUE: %-6s | VARIABLES: %s | CORRELATIONS: %s", 
          vincitori$Target[i], vincitori$Model[i], vincitori$R2[i], vincitori$RMSE[i], p_val_text, vincitori$Variabili_Usate[i], vincitori$Correlazioni_Spearman[i])
})), "Summary_Best_CV_Models_LM_Selection.txt")

# ==============================================================================
# 8. FINAL METHODOLOGICAL COMPARISON: RF VS LM
# ==============================================================================
cat("\n--- PHASE 6: FINAL EXTRACTION COMPARISON (RF vs LM) ---\n")

best_rf <- df_rf %>%
  group_by(Target) %>%
  slice_max(R2, n = 1, with_ties = FALSE) %>%
  dplyr::select(Target, R2_RF = R2, Model_RF = Model) %>%
  ungroup()

best_lm <- tabella_confronto_finale %>%
  group_by(Target) %>%
  slice_max(R2, n = 1, with_ties = FALSE) %>%
  dplyr::select(Target, R2_LM = R2, Model_LM = Model) %>%
  ungroup()

df_compare <- inner_join(best_rf, best_lm, by = "Target") %>%
  mutate(
    Delta_R2 = R2_RF - R2_LM,
    Vincitore = case_when(
      Delta_R2 > 0.01 ~ "Random Forest",
      Delta_R2 < -0.01 ~ "Linear Model",
      TRUE ~ "Tie"
    )
  ) %>%
  filter(R2_RF > 0.01 | R2_LM > 0.01) %>%
  arrange(R2_RF)

vittorie_rf <- sum(df_compare$Vincitore == "Random Forest")
vittorie_lm <- sum(df_compare$Vincitore == "Linear Model")
pareggi <- sum(df_compare$Vincitore == "Tie")
totale <- nrow(df_compare)

media_r2_rf <- mean(df_compare$R2_RF, na.rm = TRUE)
media_r2_lm <- mean(df_compare$R2_LM, na.rm = TRUE)
differenza_media_perc <- ((media_r2_rf - media_r2_lm) / media_r2_lm) * 100

cat("\n===========================================================================\n")
cat("                       FINAL METHODOLOGICAL VERDICT                        \n")
cat("===========================================================================\n\n")
cat(sprintf("Total targets analyzed (with R2 > 0.01): %d\n", totale))
cat(sprintf("- Times Random Forest wins (higher R2): %d (%.1f%%)\n", vittorie_rf, (vittorie_rf/totale)*100))
cat(sprintf("- Times Linear Model wins  (higher R2): %d (%.1f%%)\n", vittorie_lm, (vittorie_lm/totale)*100))
cat(sprintf("- Ties (difference < 0.01): %d (%.1f%%)\n\n", pareggi, (pareggi/totale)*100))

cat(sprintf("Absolute Mean R-Squared (RF Selection): %.3f\n", media_r2_rf))
cat(sprintf("Absolute Mean R-Squared (LM Selection): %.3f\n\n", media_r2_lm))

if (media_r2_rf > media_r2_lm) {
  cat(">>> CONCLUSION FOR THE PAPER:\n")
  cat(sprintf("Generally, feature selection via Random Forest is significantly SUPERIOR.\n"))
  cat(sprintf("By replacing it with a Linear Model (LM), the overall performance of the framework \ndrops by an average of %.1f%%.\n\n", differenza_media_perc))
  cat("Ecological reasoning: LM discards non-linear LiDAR metrics (e.g., roughness, complexity),\ndepriving subsequent algorithms of key variables, whereas RF retains them, maximizing \nthe R-squared. The use of RF in the feature selection phase is fully justified.\n")
} else {
  cat(">>> CONCLUSION:\n")
  cat("Surprisingly, Linear Model selection proved equal or superior.\nThis indicates that the relationships between forest structure and TreMs in your dataset are \npredominantly linear. You might consider using LM for computational simplicity.\n")
}
cat("===========================================================================\n")

# ==============================================================================
# 9. DUMBBELL PLOT (RF VS LM COMPARISON)
# ==============================================================================
df_plot_dumb <- df_compare %>%
  mutate(Target = factor(Target, levels = Target))

plot_dumb <- ggplot(df_plot_dumb) +
  geom_segment(aes(y = Target, yend = Target, x = R2_LM, xend = R2_RF), color = "grey60", linewidth = 1) +
  geom_point(aes(y = Target, x = R2_LM, color = "Linear Model Selection"), size = 3.5) +
  geom_point(aes(y = Target, x = R2_RF, color = "Random Forest Selection"), size = 3.5) +
  scale_color_manual(values = c("Linear Model Selection" = "#E74C3C", "Random Forest Selection" = "#1ABC9C")) +
  theme_minimal() +
  labs(
    title = "Impact of Feature Selection Method on Predictive Accuracy",
    subtitle = "Comparing maximum R² obtained via Random Forest vs Linear Model variable selection",
    x = "Maximum R-Squared (10-fold CV)", y = "TreMs (Target variables)", color = "Feature Selection Method"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank()
  )

ggsave("Plot_Performance/Comparison_RF_vs_LM_Selection.png", plot = plot_dumb, width = 10, height = 8, dpi = 500, bg = "white")

cat("\n[✔] Processing complete!\n")
cat("- All files and plots have been saved in the selected directory.\n")