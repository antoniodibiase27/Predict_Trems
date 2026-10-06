# ==============================================================================
# SCRIPT 04: BASELINE COMPARISON (TRADITIONAL ALLOMETRIC VS ADVANCED TLS METRICS)
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Description:
# This script compares the predictive performance of traditional allometric baseline 
# models against the advanced TLS-derived structural metrics framework.
# 1. Defines traditional baseline variables (DBH, Height, Slenderness, etc.).
# 2. Automatically merges TLS metrics, field surveys, and inventory tree heights.
# 3. Executes mirrored 10-fold cross-validation across winning algorithms.
# 4. Exports comparative tables (Excel) and structural gain plots (500 DPI).
# ==============================================================================
rm(list = setdiff(ls(), c("df_model", "df_full_max", "ALL_TARGETS")))
gc()

if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, caret, randomForest, xgboost, pscl, ggplot2, writexl, tidyr, stringr)

# ==============================================================================
# 1. DEFINITION OF TRADITIONAL METRICS (BASELINE)
# ==============================================================================
METRICHE_TRADIZIONALI <- c(
  "DBH_cm",           
  "Height",           
  "Slenderness",      
  "Biomass_vox",      
  "crown_diam_x",     
  "crown_diam_y",     
  "projected_area"    
)

# ==============================================================================
# 2. DATA LOADING AND INTEGRATION OF HEIGHT/SLENDERNESS FROM SPECIFIC EXCEL
# ==============================================================================
if(!exists("df_model")) {
  cat("\n[!] Dataset 'df_model' not found in memory. Recreating it.\n")
  
  # 1. Uniformiamo subito le colonne identificative in "ID_Tree" per tutti i dataframe
  df_lidar <- read.csv("../Data/Extracted_TLS_Metrics_Full.csv") %>% 
    rename(ID_Tree = treeID) %>% mutate(ID_Tree = as.character(ID_Tree))
  
  df_terra <- read_excel("../Data/TreMs_Field_Surveys.xlsx") %>% 
    rename(ID_Tree = treeID) %>% mutate(ID_Tree = as.character(ID_Tree))
  
  path_piante_lom <- "../Data/TreMs_Field_Surveys.xlsx" # (Verifica che sia il file giusto per le altezze)
  
  df_lom <- read_excel(path_piante_lom) %>% 
    rename(ID_Tree = treeID) %>% 
    mutate(ID_Tree = as.character(ID_Tree)) %>% 
    rename(Height_LOM = any_of(c("Height_m", "Altezza", "H"))) %>% 
    dplyr::select(ID_Tree, Height_LOM)
  
  # 2. Processiamo i target ecologici di df_terra
  tipi_64 <- grep("^[A-Z]{2}[0-9]+", names(df_terra), value = TRUE)
  df_terra <- df_terra %>% mutate(across(all_of(tipi_64), as.numeric)) %>% mutate(across(all_of(tipi_64), ~ tidyr::replace_na(., 0)))
  df_terra$Macro_TreMs_Abundance <- rowSums(df_terra[, tipi_64], na.rm = TRUE)
  df_terra$Macro_TreMs_Richness <- rowSums(df_terra[, tipi_64] > 0, na.rm = TRUE)
  
  for(fp in unique(substr(tipi_64, 1, 2))) { df_terra[[paste0("Form_", fp)]] <- rowSums(df_terra[, grep(paste0("^", fp), names(df_terra)), drop=FALSE], na.rm=TRUE) }
  for(gp in unique(substr(tipi_64, 1, 3))) { df_terra[[paste0("Group_", gp)]] <- rowSums(df_terra[, grep(paste0("^", gp), names(df_terra)), drop=FALSE], na.rm=TRUE) }
  
  ALL_TARGETS <- c("Macro_TreMs_Abundance", "Macro_TreMs_Richness", grep("Form_|Group_", names(df_terra), value=T), tipi_64)
  
  # 3. Ora i join avvengono in modo pulito e automatico usando l'unica chiave "ID_Tree"
  df_all <- inner_join(df_lidar, df_terra %>% dplyr::select(all_of(c("ID_Tree", ALL_TARGETS))), by = "ID_Tree") %>%
    inner_join(df_lom, by = "ID_Tree")
  
  df_model <- df_all %>% 
    mutate(
      Height = Height_LOM,                  
      Slenderness = Height / DBH_cm         
    ) %>%
    dplyr::select(-any_of(c("ID_Tree", "Height_LOM"))) %>% 
    na.omit()
  
  cat("[✔] Dataset successfully merged and recalculated!\n")
}

metriche_scelte <- intersect(METRICHE_TRADIZIONALI, names(df_model))
cat("\n[*] Traditional metrics verified for Baseline:\n")
print(metriche_scelte)

if(!exists("df_full_max")) {
  linee <- readLines("../Data/Model_Performance_Report.txt")
  
  df_full_max <- data.frame(Target = character(), Best_Model = character(), R2_Full_TLS = numeric(), stringsAsFactors = FALSE)
  
  for (linea in linee) {
    if (grepl("TARGET:", linea)) {
      target <- trimws(sub(".*TARGET:\\s*(.*?)\\s*\\|.*", "\\1", linea))
      best_model <- trimws(sub(".*BEST MODEL:\\s*(.*?)\\s*\\|.*", "\\1", linea))
      r2_val <- as.numeric(trimws(sub(".*R2 CV:\\s*(.*?)\\s*\\|.*", "\\1", linea)))
      
      if (!is.na(r2_val) && r2_val > 0.09) {
        df_full_max <- rbind(df_full_max, data.frame(Target = target, Best_Model = best_model, R2_Full_TLS = r2_val))
      }
    }
  }
}

targets_da_analizzare <- df_full_max$Target

# ==============================================================================
# 3. DYNAMIC CROSS-VALIDATION BASED ON WINNING ALGORITHM
# ==============================================================================
cat("\nExecuting mirrored Baseline Cross-Validation on models (RF, ZINB, XGB, LM)... Please wait...\n")
risultati_baseline <- data.frame()

set.seed(123)
for (target in targets_da_analizzare) {
  if (!target %in% names(df_model)) next
  
  modello_vincente <- df_full_max$Best_Model[df_full_max$Target == target]
  folds <- createFolds(df_model[[target]], k = 10)
  r2_folds <- c()
  
  for (i in 1:10) {
    train_data <- df_model[-folds[[i]], c(target, metriche_scelte)]
    test_data  <- df_model[folds[[i]], c(target, metriche_scelte)]
    pred <- NULL
    
    if (grepl("RF", modello_vincente, ignore.case = TRUE) || grepl("Random", modello_vincente, ignore.case = TRUE)) {
      mod <- tryCatch(randomForest(as.formula(paste(target, "~ .")), data = train_data, ntree = 100), error = function(e) NULL)
      if (!is.null(mod)) pred <- predict(mod, test_data)
      
    } else if (grepl("XGB", modello_vincente, ignore.case = TRUE)) {
      X_train <- as.matrix(train_data[, metriche_scelte]); y_train <- train_data[[target]]
      X_test  <- as.matrix(test_data[, metriche_scelte])
      mod <- tryCatch(xgboost(data = X_train, label = y_train, nrounds = 50, objective = "reg:squarederror", verbose = 0), error = function(e) NULL)
      if (!is.null(mod)) pred <- predict(mod, X_test)
      
    } else if (grepl("ZINB", modello_vincente, ignore.case = TRUE)) {
      train_data[[target]] <- round(train_data[[target]])
      mod <- tryCatch(zeroinfl(as.formula(paste(target, "~ .")), data = train_data, dist = "negbin"), error = function(e) NULL)
      if (!is.null(mod)) pred <- predict(mod, test_data)
      
    } else {
      mod <- tryCatch(lm(as.formula(paste(target, "~ .")), data = train_data), error = function(e) NULL)
      if (!is.null(mod)) pred <- predict(mod, test_data)
    }
    
    if (!is.null(pred)) {
      pred <- pmax(0, pred)
      actual <- test_data[[target]]
      if (sd(pred) > 0 && length(unique(pred)) > 1) {
        r2_folds <- c(r2_folds, cor(actual, pred)^2)
      } else {
        r2_folds <- c(r2_folds, 0)
      }
    } else {
      r2_folds <- c(r2_folds, 0)
    }
  }
  
  risultati_baseline <- rbind(risultati_baseline, data.frame(Target = target, R2_Baseline = mean(r2_folds, na.rm = TRUE)))
}

# ==============================================================================
# 4. DATA MERGING AND STRUCTURAL GAIN (DELTA) CALCULATION
# ==============================================================================
df_confronto <- df_full_max %>%
  left_join(risultati_baseline, by = "Target") %>%
  mutate(
    R2_Baseline = tidyr::replace_na(R2_Baseline, 0),
    Delta_R2 = R2_Full_TLS - R2_Baseline,
    Percent_Improvement = (Delta_R2 / (R2_Baseline + 0.001)) * 100
  ) %>%
  arrange(desc(R2_Full_TLS))

cat("\n--- UPDATED COMPARATIVE TABLE (Baseline LOM vs TLS Model) ---\n")
print(df_confronto %>% select(Target, Best_Model, R2_Baseline, R2_Full_TLS, Delta_R2), row.names = FALSE)

dir.create("Plot_Performance", showWarnings = FALSE)
write_xlsx(df_confronto, "Baseline_vs_Best_TLS_Model_Comparison.xlsx")
cat("\n[✔] Excel report successfully saved.\n")

# ==============================================================================
# 5. PUBLICATION-READY COMPARATIVE PLOT (500 DPI - UNTITLED)
# ==============================================================================
df_confronto <- df_confronto %>%
  mutate(Target_Label = sprintf("%s (%s)", Target, gsub("_", " ", Best_Model)))

df_plot <- df_confronto %>%
  select(Target_Label, `Baseline (Traditional metrics)` = R2_Baseline, `Best Model (Advanced TLS metrics)` = R2_Full_TLS) %>%
  tidyr::pivot_longer(cols = -Target_Label, names_to = "Model_Type", values_to = "R2")

plot_confronto <- ggplot(df_plot, aes(x = reorder(Target_Label, R2, FUN = max), y = R2, fill = Model_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", linewidth = 0.3) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("#7F8C8D", "#1ABC9C")) + 
  labs(x = "TreMs Target (with its best performing algorithm)", y = "Mean CV R-Squared", fill = "Predictor Set") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 9, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank()
  )

print(plot_confronto)

ggsave("Plot_Performance/Baseline_vs_Best_TLS_Model_Comparison.png", plot = plot_confronto, width = 11, height = 8, dpi = 500, bg = "white")
cat("[✔] Comparative plot saved at 500 DPI in 'Plot_Performance/Baseline_vs_Best_TLS_Model_Comparison.png'\n")
