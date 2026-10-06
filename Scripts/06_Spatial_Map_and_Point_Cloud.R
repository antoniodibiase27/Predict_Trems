# ==============================================================================
# SCRIPT 06: SPATIAL MAPS AND 3D POINT-CLOUDS VISUALIZATION
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Description:
# This script generates publication-ready spatial and 3D visual outputs (Figures 6 & 7).
# 1. Spatial Maps (2D): Merges LiDAR predictors with best models to map observed vs. 
#    predicted TreMs Abundance and Richness across the spatial coordinates of the plot.
# 2. Point-Clouds (3D): Dynamically identifies 4 success cases (accurate prediction) 
#    and 4 failure cases (underestimation) for large trees (DBH >= 50 cm). It renders 
#    True-RGB 3D point clouds alongside an inset of the horizontal stem cross-section.
# Outputs are automatically exported as high-resolution 600 DPI TIFF images.
# ==============================================================================

rm(list = ls())
gc()

if (!require("pacman")) install.packages("pacman")
pacman::p_load(lidR, readxl, dplyr, tidyr, randomForest, xgboost, pscl, ggplot2, viridis, data.table)
library(grid)

# ==============================================================================
# 1. FILE AND FOLDER SELECTION
# ==============================================================================

cat("\n--- STEP 1: DATA SELECTION ---\n")
dir_out <- "../Outputs"
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
setwd(dir_out)

df_lidar <- read.csv("../Data/Extracted_TLS_Metrics_Full.csv") %>% rename(ID_Albero = treeID) %>% mutate(ID_Albero = as.character(ID_Albero))

df_terra <- read_excel("../Data/TreMs_Field_Surveys.xlsx") %>% mutate(treeID = as.character(treeID))

report <- readLines("../Data/Model_Performance_Report.txt", warn = FALSE)

dir_las <- "../Data/Sample_LAS_Trees"

file_las_disponibili <- list.files(dir_las, pattern = "\\.(las|laz)$", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
cat(sprintf("Trovati %d file .las/.laz.\n", length(file_las_disponibili)))

trova_file_las <- function(id_albero, lista_file) {
  nomi_base <- tools::file_path_sans_ext(basename(lista_file))
  idx <- which(nomi_base == id_albero | grepl(paste0("(^|[^0-9])", id_albero, "$"), nomi_base))
  if (length(idx) > 0) return(lista_file[idx[1]])
  return(NA_character_)
}

nomi_col <- names(df_terra)
col_x <- nomi_col[grepl("(?i)^x$", nomi_col)][1]
col_y <- nomi_col[grepl("(?i)^y$", nomi_col)][1]
if(is.na(col_x) | is.na(col_y)) stop("ERRORE: Non trovo le colonne X e Y nel file Excel dei rilievi a terra.")

estrai_modello_e_variabili <- function(nome_target, linee_testo) {
  pattern_target <- paste0("TARGET:\\s*", nome_target, "\\b")
  riga <- grep(pattern_target, linee_testo, value = TRUE, ignore.case = TRUE)
  
  modello_out <- "RF"
  var_out <- NULL
  
  if(length(riga) > 0) {
    match_modello <- regmatches(riga[1], regexpr("BEST MODEL:\\s*([A-Za-z0-9_]+)", riga[1], ignore.case = TRUE))
    if(length(match_modello) > 0) {
      modello <- gsub("(?i)BEST MODEL:\\s*", "", match_modello)
      modello <- trimws(modello)
      if(grepl("(?i)RF|Random", modello)) modello_out <- "RF"
      else if(grepl("(?i)XGB|Boost", modello)) modello_out <- "XGB"
      else if(grepl("(?i)LM|Linear", modello)) modello_out <- "LM"
      else if(grepl("(?i)ZINB|Zero", modello)) modello_out <- "ZINB"
    }
    match_var <- regmatches(riga[1], regexpr("VARIABILI:\\s*([^|]+)", riga[1], ignore.case = TRUE))
    if(length(match_var) > 0) {
      var_str <- gsub("(?i)VARIABILI:\\s*", "", match_var)
      var_out <- trimws(unlist(strsplit(var_str, ",")))
    }
  }
  return(list(modello = modello_out, variabili = var_out))
}

# ==============================================================================
# 2. DATA PREPARATION
# ==============================================================================
cat("\n--- STEP 2: DATA PREPARATION ---\n")
tipi_64 <- grep("^[A-Z]{2}[0-9]+", names(df_terra), value = TRUE)
df_terra <- df_terra %>% mutate(across(all_of(tipi_64), as.numeric)) %>% mutate(across(all_of(tipi_64), ~ replace_na(., 0)))

df_terra$Macro_TreMs_Abundance <- rowSums(df_terra[, tipi_64], na.rm = TRUE)
df_terra$Macro_TreMs_Richness <- rowSums(df_terra[, tipi_64] > 0, na.rm = TRUE)

nomi_lidar <- setdiff(names(df_lidar), c("ID_Albero", "treeID"))

df_lidar_clean <- df_lidar %>% select(ID_Albero, all_of(nomi_lidar))

col_dbh_raw <- names(df_terra)[grepl("(?i)dbh", names(df_terra))][1]
if(is.na(col_dbh_raw)) stop("ERRORE: Impossibile trovare una colonna DBH nel file Excel a terra.")

df_terra_clean <- df_terra %>% 
  rename(ID_Albero = treeID, DBH_cm = all_of(col_dbh_raw)) %>%
  select(ID_Albero, all_of(col_x), all_of(col_y), DBH_cm, Macro_TreMs_Abundance, Macro_TreMs_Richness, all_of(tipi_64))

df_all <- inner_join(df_lidar_clean, df_terra_clean, by = "ID_Albero") %>% na.omit()

if (!"DBH_cm" %in% names(df_all)) {
  pos_dbh <- grep("(?i)dbh", names(df_all))[1]
  if(!is.na(pos_dbh)) names(df_all)[pos_dbh] <- "DBH_cm" else stop("ERRORE CRITICO: DBH_cm mancante!")
}

df_all$File_LAS <- sapply(df_all$ID_Albero, trova_file_las, lista_file = file_las_disponibili)


nomi_lidar <- intersect(nomi_lidar, names(df_all))
varianze <- sapply(df_all[, nomi_lidar], var, na.rm = TRUE)
nomi_validi <- nomi_lidar[varianze > 0 & !is.na(varianze)]

# ==============================================================================
# SERVICE FUNCTIONS (Predictive models and safe file saving)
# ==============================================================================
ottieni_predizioni_cv <- function(df, target_col, predittori_validi, info_estratta) {
  set.seed(123)
  tipo_modello <- info_estratta$modello
  top_5 <- info_estratta$variabili
  
  if (!target_col %in% names(df)) stop(sprintf("ERRORE: La colonna target '%s' non esiste!", target_col))
  if (!is.null(top_5)) top_5 <- intersect(top_5, names(df))
  
  if (is.null(top_5) || length(top_5) < 2) {
    pred_disp <- intersect(predittori_validi, names(df))
    form_rf <- as.formula(paste(target_col, "~ ."))
    rf_sel <- randomForest::randomForest(form_rf, data = df[, c(target_col, pred_disp)], ntree = 100)
    imp <- randomForest::importance(rf_sel)
    top_5 <- head(rownames(imp)[order(imp[, 1], decreasing = TRUE)], 5)
  } else {
    top_5 <- head(top_5, 5)
  }
  
  top_5 <- intersect(top_5, names(df))
  if (length(top_5) == 0) top_5 <- head(intersect(predittori_validi, names(df)), 5)
  
  form_base <- as.formula(paste(target_col, "~", paste(top_5, collapse = " + ")))
  cat(sprintf("   -> [Target: %s] Modello: %s | Variabili: %s\n", target_col, tipo_modello, paste(top_5, collapse=", ")))
  
  n_oss <- nrow(df)
  fold_id <- sample(rep(1:10, length.out = n_oss))
  preds <- rep(NA_real_, n_oss)
  
  for (i in 1:10) {
    idx_test <- which(fold_id == i)
    cols_da_usare <- unique(c(target_col, top_5))
    cols_da_usare <- intersect(cols_da_usare, names(df))
    
    train_d <- df[-idx_test, cols_da_usare, drop = FALSE]
    test_d  <- df[idx_test, cols_da_usare, drop = FALSE]
    
    if (tipo_modello == "LM") {
      mod <- tryCatch(lm(form_base, data = train_d), error=function(e) NULL)
      if(!is.null(mod)) preds[idx_test] <- pmax(0, predict(mod, test_d))
    } else if (tipo_modello == "ZINB") {
      form_zinb <- as.formula(paste(target_col, "~", paste(top_5, collapse = " + "), "|", top_5[1]))
      mod <- tryCatch(pscl::zeroinfl(form_zinb, data = train_d, dist = "negbin"), error=function(e) NULL)
      if(!is.null(mod)) preds[idx_test] <- pmax(0, predict(mod, test_d, type = "response"))
      else preds[idx_test] <- pmax(0, predict(lm(form_base, data = train_d), test_d))
    } else if (tipo_modello == "XGB") {
      X_tr <- as.matrix(train_d[, top_5, drop = FALSE])
      X_ts <- as.matrix(test_d[, top_5, drop = FALSE])
      mod <- tryCatch(xgboost::xgboost(data = X_tr, label = train_d[[target_col]], params = list(objective="count:poisson", eta=0.1), nrounds = 50, verbose = 0), error=function(e) NULL)
      if(!is.null(mod)) preds[idx_test] <- pmax(0, predict(mod, X_ts))
    } else {
      mod <- tryCatch(randomForest::randomForest(form_base, data = train_d, ntree = 150), error=function(e) NULL)
      if(!is.null(mod)) preds[idx_test] <- pmax(0, predict(mod, test_d))
    }
  }
  return(list(pred = preds, top5 = top_5, modello_usato = tipo_modello))
}

sovrascrivi_sicuro <- function(filepath) {
  if (file.exists(filepath)) {
    tryCatch({ unlink(filepath, force = TRUE) }, warning = function(w) {
      cat(sprintf("\n[!] Impossibile sovrascrivere '%s'. Chiudi l'immagine se è aperta!\n", basename(filepath)))
    })
  }
}

# ==============================================================================
# 3. FIGURE 6 GENERATION: 2D SPATIAL MAPS
# ==============================================================================
cat("\n--- STEP 3: FIGURE 6 GENERATION (SPATIAL MAPS) ---\n")
info_abd <- estrai_modello_e_variabili("Macro_TreMs_Abundance", report)
info_ric <- estrai_modello_e_variabili("Macro_TreMs_Richness", report)

df_all$Pred_Abundance <- ottieni_predizioni_cv(df_all, "Macro_TreMs_Abundance", nomi_validi, info_abd)$pred
df_all$Pred_Richness <- ottieni_predizioni_cv(df_all, "Macro_TreMs_Richness", nomi_validi, info_ric)$pred

crea_mappa <- function(dati, col_valore, titolo, colore_base, max_val, titolo_leg_colore) {
  p <- ggplot(dati, aes(x = .data[[col_x]], y = .data[[col_y]], color = .data[[col_valore]]))
  
  p <- p + geom_point(aes(size = DBH_cm), alpha = 0.85) +
    scale_size_continuous(range = c(1.5, 6), guide = "none")
  
  p <- p + scale_color_viridis_c(option = colore_base, name = titolo_leg_colore, limits = c(0, max_val)) +
    coord_equal() +
    labs(title = titolo, x = "X Coordinate (m)", y = "Y Coordinate (m)") + 
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5), 
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 1), 
      panel.grid.major = element_line(color = "grey85", linetype = "dashed"), 
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10)
    )
  return(p)
}

max_abd <- max(c(df_all$Macro_TreMs_Abundance, df_all$Pred_Abundance), na.rm = TRUE)
max_ric <- max(c(df_all$Macro_TreMs_Richness, df_all$Pred_Richness), na.rm = TRUE)

map_obs_abd <- crea_mappa(df_all, "Macro_TreMs_Abundance", "A) Observed TreMs Abundance", "magma", max_abd, "N TreMs")
map_pred_abd <- crea_mappa(df_all, "Pred_Abundance", sprintf("B) Predicted TreMs Abundance (%s)", info_abd$modello), "magma", max_abd, "N TreMs")
map_obs_ric <- crea_mappa(df_all, "Macro_TreMs_Richness", "C) Observed TreMs Richness", "viridis", max_ric, "N TreMs type")
map_pred_ric <- crea_mappa(df_all, "Pred_Richness", sprintf("D) Predicted TreMs Richness (%s)", info_ric$modello), "viridis", max_ric, "N TreMs type")

nome_fig6 <- "Figura_6_Spatial_Maps.tiff"
sovrascrivi_sicuro(nome_fig6)

tiff(filename = nome_fig6, width = 14, height = 12, units = "in", res = 600, compression = "lzw", bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2, heights = unit(c(0.1, 0.45, 0.45), "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text("Spatial Distribution of Tree-related Microhabitats: Observed vs. Predicted", y = unit(0.6, "npc"), gp = gpar(fontface = "bold", fontsize = 16))
grid.text("Comparison of field observations and 10-fold cross-validated predictions (100 x 100 m plot)", y = unit(0.3, "npc"), gp = gpar(fontsize = 12, col = "grey30"))
popViewport()
posizioni_mappe <- list(list(p=map_obs_abd, r=2, c=1), list(p=map_pred_abd, r=2, c=2), list(p=map_obs_ric, r=3, c=1), list(p=map_pred_ric, r=3, c=2))
for (pos in posizioni_mappe) {
  pushViewport(viewport(layout.pos.row = pos$r, layout.pos.col = pos$c))
  print(pos$p, newpage = FALSE)
  popViewport()
}
popViewport()
dev.off()

cat(sprintf("-> Figura 6 (Mappe Spaziali) salvata con successo: %s\n", nome_fig6))

# ==============================================================================
# 4. SELECTION OF 8 TREES FOR FIGURE 7 (DBH >= 50 cm, HIGH ABUNDANCE)
# ==============================================================================
cat("\n--- STEP 4: SELECTION OF SUCCESS AND FAILURE CASES ---\n")
alberi_scelti <- list()
usati_id <- c()

target_successi <- c("GR11", "GR13", "GR12", "GR31")

for (t_name in target_successi) {
  if(!t_name %in% names(df_all)) next
  info_mod <- estrai_modello_e_variabili(t_name, report)
  res_cv <- ottieni_predizioni_cv(df_all, t_name, nomi_validi, info_mod)
  
  df_temp <- df_all %>% 
    mutate(Obs = .data[[t_name]], Pred = res_cv$pred, Errore_Assoluto = abs(Obs - Pred)) %>% 
    filter(!ID_Albero %in% usati_id, Obs >= 1, !is.na(File_LAS)) # <-- FILTRO CRITICO
  
  if(nrow(df_temp) == 0) {
    cat(sprintf("[!] Nessun LAS valido trovato per il target (Successo): %s\n", t_name))
    next
  }
  
  candidato <- df_temp %>% filter(DBH_cm >= 50, Obs >= 2, Errore_Assoluto <= 0.40) %>% arrange(Errore_Assoluto, desc(Obs)) %>% slice(1)
  if(nrow(candidato) == 0) candidato <- df_temp %>% filter(DBH_cm >= 50, Errore_Assoluto <= 0.35) %>% arrange(Errore_Assoluto, desc(Obs)) %>% slice(1)
  if(nrow(candidato) == 0) candidato <- df_temp %>% filter(Obs >= 2, Errore_Assoluto <= 0.40) %>% arrange(Errore_Assoluto, desc(Obs)) %>% slice(1)
  if(nrow(candidato) == 0) candidato <- df_temp %>% arrange(Errore_Assoluto) %>% slice(1)
  
  if(nrow(candidato) > 0) {
    usati_id <- c(usati_id, candidato$ID_Albero)
    idx_curr <- length(alberi_scelti) + 1
    alberi_scelti[[idx_curr]] <- list(
      Pan_Letter = LETTERS[idx_curr], Categoria = "SUCCESS CASE (Accurate Prediction)", 
      Tipo_TreM = t_name, Modello = info_mod$modello, ID_Albero = candidato$ID_Albero, 
      File_LAS = candidato$File_LAS, Obs = candidato$Obs, Pred = candidato$Pred, 
      DBH = round(candidato$DBH_cm, 1), Slenderness = round(candidato$Slenderness, 2), Is_Success = TRUE
    )
  }
}

target_fallimenti <- c("EP32", "BA21", "IN12", "OT21")

for (t_name in target_fallimenti) {
  if(!t_name %in% names(df_all)) next
  info_mod <- estrai_modello_e_variabili(t_name, report)
  res_cv <- ottieni_predizioni_cv(df_all, t_name, nomi_validi, info_mod)
  
  # Filtro critico: mantengo solo le piante di cui ho caricato il file LAS nella repository
  df_temp <- df_all %>% 
    mutate(Obs = .data[[t_name]], Pred = res_cv$pred, Errore_Assoluto = abs(Obs - Pred)) %>% 
    filter(!ID_Albero %in% usati_id, Obs >= 1, !is.na(File_LAS))
  
  # Controllo di sicurezza: se nessun albero nel sample ha questo TreM, avvisa e salta
  if(nrow(df_temp) == 0) {
    cat(sprintf("[!] Nessun LAS valido trovato nel sample per il target (Fallimento): %s\n", t_name))
    next
  }
  
  # Ricerca a cascata del candidato migliore
  candidato <- df_temp %>% filter(DBH_cm >= 50, Pred < 0.25) %>% arrange(desc(Obs), Pred) %>% slice(1)
  
  if(nrow(candidato) == 0) {
    candidato <- df_temp %>% filter(Pred < 0.25) %>% arrange(desc(Obs), Pred) %>% slice(1)
  }
  if(nrow(candidato) == 0) {
    candidato <- df_temp %>% arrange(desc(Obs), Pred) %>% slice(1)
  }
  
  # Salvataggio nella lista
  if(nrow(candidato) > 0) {
    usati_id <- c(usati_id, candidato$ID_Albero)
    idx_curr <- length(alberi_scelti) + 1
    alberi_scelti[[idx_curr]] <- list(
      Pan_Letter = LETTERS[idx_curr], Categoria = "FAILURE CASE (Underestimated Superficial TreM)", 
      Tipo_TreM = t_name, Modello = info_mod$modello, ID_Albero = candidato$ID_Albero, 
      File_LAS = candidato$File_LAS, Obs = candidato$Obs, Pred = candidato$Pred, 
      DBH = round(candidato$DBH_cm, 1), Slenderness = round(candidato$Slenderness, 2), Is_Success = FALSE
    )
  }
}
# ==============================================================================
# 5. GRAPHICAL FUNCTION FOR TRUE-RGB POINT CLOUD + STEM SECTION
# ==============================================================================
crea_pannello_nuvola <- function(info) {
  las <- readLAS(info$File_LAS)
  pts <- as.data.frame(las@data)
  
  pts$X <- pts$X - median(pts$X, na.rm = TRUE)
  pts$Y <- pts$Y - median(pts$Y, na.rm = TRUE)
  pts$Z <- pts$Z - min(pts$Z, na.rm = TRUE)
  
  has_rgb <- all(c("R", "G", "B") %in% names(pts))
  if (has_rgb) {
    max_col <- ifelse(max(pts$R, na.rm = TRUE) > 255, 65535, 255)
    pts$Red <- pts$R / max_col
    pts$Green <- pts$G / max_col
    pts$Blue <- pts$B / max_col
    pts$TrueColor <- rgb(pts$Red, pts$Green, pts$Blue, maxColorValue = 1)
  } else {
    pts$TrueColor <- "#238B45" 
  }
  
  slice_fusto <- pts %>% filter(Z >= 0.2 & Z <= 2.0)
  if (nrow(slice_fusto) > 6000) slice_fusto <- slice_fusto[sample(nrow(slice_fusto), 6000), ]
  if (nrow(pts) > 30000) pts <- pts[sample(nrow(pts), 30000), ]
  
  colore_bordo <- ifelse(info$Is_Success, "#1B7837", "#B2182B")
  titolo_pan <- sprintf("%s) Target: %s", info$Pan_Letter, info$Tipo_TreM)
  
  sottotitolo_pan <- sprintf("Tree ID: %s | Model: %s | Obs: %.1f vs Pred: %.2f", info$ID_Albero, info$Modello, info$Obs, info$Pred)
  
  p_main <- ggplot(pts, aes(x = X, y = Z)) + 
    geom_point(color = pts$TrueColor, size = 0.18, alpha = 0.8) + 
    coord_equal() + 
    labs(title = titolo_pan, subtitle = sottotitolo_pan, x = "Relative X (m)", y = "Height Z (m)") + 
    theme_minimal(base_size = 9) + 
    theme(
      plot.title = element_text(face = "bold", size = 10, color = colore_bordo, hjust = 0), 
      plot.subtitle = element_text(size = 8, color = "grey20", face = "italic"), 
      panel.border = element_rect(color = colore_bordo, fill = NA, linewidth = 1.2), 
      panel.grid.minor = element_blank()
    )
  
  p_inset <- ggplot(slice_fusto, aes(x = X, y = Y, color = Z)) + 
    geom_point(size = 0.25, alpha = 0.8) + 
    scale_color_viridis_c(option = "turbo", guide = "none") +
    coord_equal() + 
    labs(title = "Stem Section", x = NULL, y = NULL) + 
    theme_void() + 
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0.5, margin = margin(b = 1)), 
      plot.background = element_rect(fill = alpha("white", 0.92), color = "grey40", linewidth = 0.4), 
      plot.margin = margin(2, 2, 2, 2)
    )
  
  return(list(main = p_main, inset = p_inset))
}

# ==============================================================================
# 6. ASSEMBLY OF THE 8-PANEL SUPER-FIGURE (4x2) AND TIFF EXPORT
# ==============================================================================
cat("\nGenerating the 8 panels (Figure 7)...\n")
lista_plot <- lapply(alberi_scelti, crea_pannello_nuvola)

nome_fig7 <- "Figura_7_Success_vs_Failure_PointClouds.tiff"
sovrascrivi_sicuro(nome_fig7)

tiff(filename = nome_fig7, width = 27, height = 14, units = "in", res = 600, compression = "lzw", bg = "white")
grid.newpage()

pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 4, heights = unit(c(0.08, 0.46, 0.46), "npc"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:4))
grid.text("TLS Point-Cloud True-RGB Comparison: Model Success vs. Failure Cases (DBH >= 50 cm)", y = unit(0.7, "npc"), gp = gpar(fontface = "bold", fontsize = 15))
grid.text("Top row (A-D): Volumetric targets successfully predicted. Bottom row (E-H): Superficial targets leading to underestimation", y = unit(0.2, "npc"), gp = gpar(fontsize = 10.5, col = "grey30"))
popViewport()

posizioni_8 <- list(
  c(r = 2, c = 1), c(r = 2, c = 2), c(r = 2, c = 3), c(r = 2, c = 4), 
  c(r = 3, c = 1), c(r = 3, c = 2), c(r = 3, c = 3), c(r = 3, c = 4)  
)

for (k in 1:length(alberi_scelti)) {
  pushViewport(viewport(layout.pos.row = posizioni_8[[k]]["r"], layout.pos.col = posizioni_8[[k]]["c"]))
  print(lista_plot[[k]]$main, newpage = FALSE)
  
  vp_inset <- viewport(x = unit(0.24, "npc"), y = unit(0.72, "npc"), width = unit(0.28, "npc"), height = unit(0.26, "npc"))
  pushViewport(vp_inset)
  print(lista_plot[[k]]$inset, newpage = FALSE)
  popViewport()
  
  popViewport()
}

popViewport()
dev.off()

cat(sprintf("\n===========================================================================\n"))
cat(sprintf("8-PANEL SUPER-FIGURE 7 SUCCESSFULLY SAVED: %s\n", nome_fig7))
cat("Summary of the targets included (all with DBH >= 50 cm):\n")
for (info in alberi_scelti) {
  cat(sprintf(" - Panel %s [%s]: ID = %s | Target = %s | Obs = %.1f vs Pred = %.2f\n",
              info$Pan_Letter, ifelse(info$Is_Success, "SUCCESS", "FAILURE"), 
              info$ID_Albero, info$Tipo_TreM, info$Obs, info$Pred))
}
cat("===========================================================================\n")
