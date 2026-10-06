# Predicting Tree-related Microhabitats abundance and richness using Terrestrial Laser Scanner data in an old-growth forest: a comparison between parametric linear models and machine learning
Dataset and R Scripts associated with the manuscript submitted to Ecological Informatics

**Authors:** Antonio di Biase*, César Alvites, Pierdomenico Spina, Renzo Motta, Giovanni Santopuoli

        * corresponding author: a.dibiase5@studenti.unimol.it
## 1. Description
This repository contains the data and R scripts necessary to fully replicate the statistical and machine learning modeling of Tree-related Microhabitats (TreMs) extracted from Terrestrial Laser Scanner (TLS) point clouds.

In accordance with the transparency and reproducibility standards of Ecological Informatics, all analytical workflows are provided. Due to storage limitations, the raw multi-gigabyte point clouds (totaling 889 trees) are not uploaded in their entirety. Instead, we provide the complete final tabular datasets required for modeling, alongside a subset of raw compressed .laz files to demonstrate the visual assessment and mapping methodology.

## 2. Repository Structure
The repository is organized into two main folders:

### Data/

* #### TreMs_Field_Surveys.xlsx: 
The ground-truth dataset containing tree IDs and the abundance/richness of TreMs observed during field surveys.

* #### Extracted_TLS_Metrics_Full.csv: 
The complete dataset containing the morpho-structural metrics already extracted for all 889 trees.

* #### Model_Performance_Report.txt:
The summary report containing the best cross-validated predictive models, exact mathematical formulas, and selected LiDAR features. Providing this pre-computed file allows users to run downstream graphical scripts (e.g., spatial mapping) independently.

* **FOLDER: Sample_LAZ_Trees/** 
A folder containing 40 representative individual tree point clouds (.laz format). 
These files are used to test the feature extraction algorithm and the generation of 3D True-RGB success/failure prediction cases.


### Scripts/
A suite of sequential R scripts handling the entire pipeline. To facilitate the review process, an interactive Master Launcher is provided.

* #### 00_MASTER_LAUNCHER.R:
An interactive R script acting as a control panel to easily run any part of the analysis, or the entire pipeline sequentially.

* #### 01_Metrics_Extraction_lidR.R:
Extracts geometric and volumetric LiDAR metrics from .laz point clouds. Features an automated cross-section "rescue" logic to ensure robust DBH calculation even in presence of understory noise.

* #### 02_Predictive_Modeling.R:
Performs automated ecological hierarchy construction, multicollinearity filtering, RF feature selection, 10-fold CV (RF, XGB, LM, ZINB), and exact formula extraction for linear and zero-inflated models.

* #### 03_Clustering_and_Network.R: 
Generates hierarchical clustering and the ecological network graph.

* #### 04_Baseline_vs_TLS_Comparison.R:
Benchmarks advanced TLS metrics against traditional allometric variables.

* #### 05_Feature_Selection_LM_vs_RF.R:
Compares non-linear vs. linear variable selection impact on predictive accuracy.

* #### 06_Spatial_Maps_and_Point_Cloud.R:
Generates 2D spatial distribution maps and renders 8-panel True-RGB 3D point-cloud comparisons to visually assess model success vs. failure cases based on TreM volumetric properties.

## 3. Workflow and Reproducibility Instructions
All analyses were performed in R (version 4.5.1). Ensure that all necessary packages listed in the scripts (e.g., randomForest, xgboost, pscl, caret, lidR, ggplot2, viridis) are installed.

To reproduce the analyses:

1. Create a dedicated folder on your local machine (e.g., TreMs_Modeling) containing the Data/, Scripts/, and Outputs/ subfolders.
2. Open the Scripts/ folder in your R environment and run the 00_MASTER_LAUNCHER.R file.
3. An interactive menu will appear in the R console. Type the number corresponding to the analysis you wish to reproduce, or type 7 to run the entire pipeline sequentially.
4. **Note on Interactivity and Working Directory**: The scripts are designed to be fully interactive. At the beginning of each execution, you will be prompted via pop-up windows to select the necessary input files from the Data/ folder. Crucially, you will also be prompted to choose the output directory: ensure you select the Outputs/ folder so that all generated TIFF plots and summary tables are correctly stored in one place.
5. Alternatively, each script can be used independently by running it directly in R, without using the Master Launcher.

---
**License**: CC-BY 4.0  
**Contact**: For any questions regarding the dataset or code, please contact the corresponding author (Antonio di Biase: a.dibiase5@studenti.unimol.it).
