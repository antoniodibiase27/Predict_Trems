# ==============================================================================
# SCRIPT 00: MASTER LAUNCHER - REPOSITORY PIPELINE
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)

# Description:
# This script allows to execute the entire analytical pipeline of the paper:

# "Predicting Tree-related Microhabitats abundance and richness using 
#     Terrestrial Laser Scanner data in an old-growth forest: a comparison 
#     between parametric linear models and machine learning"

# Authors: di Biase et al.
# Journal: Ecological Informatics
# ==============================================================================

cat("\n===========================================================================\n")
cat("   WELCOME TO THE PROJECT REPOSITORY (Di Biase et al.)                     \n")
cat("===========================================================================\n")
cat("Please choose which part of the analysis you want to reproduce:\n\n")

master_choice <- menu(
  c(
    "[Script 01] Main Predictive Modeling (RF, XGB, LM, ZINB & Formula Extraction)",
    "[Script 02] Generate Dendrogram and Semantic Network (Requires 'Risultati_Completi.xlsx')",
    "[Script 03] Baseline vs TLS (Comparison with traditional allometric metrics)",
    "[Script 04] Feature Selection Test (Random Forest vs Linear Model approaches)",
    "[Script 05] Utility: Export Word tables and PNG plots for the manuscript",
    "[Script 06] Spatial Maps 2D and 3D Point-Clouds (Success/Failure cases)",
    "[ALL SCRIPTS] Run the entire pipeline sequentially (Scripts 01 to 06)"
  ),
  title = "Type the number of the operation to execute (0 to exit):"
)

# Dynamic execution based on user choice
script_folder <- getwd() 

if (master_choice == 1) {
  cat("\n>>> Running Script 01: Main Predictive Modeling...\n")
  source(file.path(script_folder, "01_Main_Predictive_Modeling.R"))
  
} else if (master_choice == 2) {
  cat("\n>>> Running Script 02: Clustering and Semantic Network...\n")
  source(file.path(script_folder, "02_Clustering_and_Semantic_Network.R"))
  
} else if (master_choice == 3) {
  cat("\n>>> Running Script 03: Baseline vs TLS...\n")
  source(file.path(script_folder, "03_Baseline_vs_TLS_Comparison.R"))
  
} else if (master_choice == 4) {
  cat("\n>>> Running Script 04: Feature Selection Test...\n")
  source(file.path(script_folder, "04_Feature_Selection_RF_vs_LM.R"))
  
} else if (master_choice == 5) {
  cat("\n>>> Running Script 05: Export Tables and Figures...\n")
  source(file.path(script_folder, "05_Export_Tables_and_Figures.R"))
  
} else if (master_choice == 6) {
  cat("\n>>> Running Script 06: Spatial Maps and 3D Point-Clouds...\n")
  source(file.path(script_folder, "06_Spatial_Maps_and_3D_PointClouds.R"))
  
} else if (master_choice == 7) {
  cat("\n===========================================================================\n")
  cat(">>> RUNNING FULL PIPELINE (Scripts 01 to 06)\n")
  cat("NOTE: The scripts are interactive. You will be prompted to select the \n")
  cat("      required data files and output folders during the execution.\n")
  cat("===========================================================================\n")
  
  cat("\n---> [1/6] Running Script 01...\n")
  source(file.path(script_folder, "01_Main_Predictive_Modeling.R"))
  
  cat("\n---> [2/6] Running Script 02...\n")
  source(file.path(script_folder, "02_Clustering_and_Semantic_Network.R"))
  
  cat("\n---> [3/6] Running Script 03...\n")
  source(file.path(script_folder, "03_Baseline_vs_TLS_Comparison.R"))
  
  cat("\n---> [4/6] Running Script 04...\n")
  source(file.path(script_folder, "04_Feature_Selection_RF_vs_LM.R"))
  
  cat("\n---> [5/6] Running Script 05...\n")
  source(file.path(script_folder, "05_Export_Tables_and_Figures.R"))
  
  cat("\n---> [6/6] Running Script 06...\n")
  source(file.path(script_folder, "06_Spatial_Maps_and_3D_PointClouds.R"))
  
  cat("\n===========================================================================\n")
  cat(">>> FULL PIPELINE COMPLETED SUCCESSFULLY!\n")
  cat("===========================================================================\n")
  
} else {
  cat("\nOperation cancelled. Exiting the launcher.\n")
}