# ==============================================================================
# SCRIPT 1: TLS POINT CLOUD METRICS EXTRACTION
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Description:
# This script processes Terrestrial Laser Scanning (TLS) point clouds of individual
# segmented trees (.las/.laz format). It extracts advanced geometric, structural, 
# and volumetric metrics required for ecological modeling (e.g., TreMs prediction).
# 
# Key Features:
# 1. Extracts stem diameter by slicing the point cloud 
#    at standard breast height (1.2m - 1.4m). If noise or gaps prevent computation, 
#    a 'rescue' logic automatically expands the slice to (1.1m - 1.5m).
# 2. Hull Convex Area & DBSCAN Clustering applied to stem cross-sections to isolate 
#    the main trunk from noise.
# 3. Automatic 2D graphical export of stem cross-sections for visual validation.
# 4. Calculation of 28 advanced structural metrics:
#    - Height percentiles (p05, p25, p50, p75, p95), Roughness, Skewness, Kurtosis
#    - Voxel-based Volume (Biomass index), 3D Envelope Volume, 2D Projected Area
#    - Crown metrics (Depth, Diameter X/Y), Stratification density
#    - PCA-based Shape Indices (Linearity, Planarity, Sphericity, Asymmetry)
#    - Slenderness index (Height/DBH)
# 5. Export to "COMPLETE_FINAL_Metrics.csv"
# ==============================================================================

# 1. LOAD REQUIRED LIBRARIES
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lidR, dplyr, geometry, pracma, tools, ggplot2, tidyverse, cowplot, dbscan, moments)

start_time <- Sys.time()

# =========================================================
# 2. FOLDER CONFIGURATION
# =========================================================
cat("\n>>> PHASE 1: Folder Preparation <<<\n")

folder_path <- "../Data/Sample_LAS_Trees"
plot_folder <- "../Outputs/DBH_Graph" # or output folder
if (!dir.exists(plot_folder)) {dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)}

if(!dir.exists(folder_path)) stop("ERROR: '../Data/Sample_LAS_Trees' folder not found! Check directory structure.")

las_files <- list.files(folder_path, pattern = "\\.las$|\\.laz$", full.names = TRUE, ignore.case = TRUE)
if(length(las_files) == 0) stop("NO LAS/LAZ FILES FOUND IN THE DIRECTORY!")
cat("Found", length(las_files), "files to analyze.\n")

# =========================================================
# 3. METRIC EXTRACTION ALGORITHM
# =========================================================
cat("\n>>> PHASE 2: Comprehensive Point Cloud Analysis <<<\n")

extract_metrics_and_plot <- function(file_path, plot_folder) {
  
  las <- tryCatch({ readLAS(file_path) }, error = function(e) return(NULL))
  if (is.empty(las)) return(NULL)
  
  tree_id <- file_path_sans_ext(basename(file_path))
  
  # -------------------------------------------------------
  # A. CLEAN DBH CALCULATION (Fallback Logic: 1.10m - 1.50m)
  # -------------------------------------------------------
  df <- las@data
  dbh_cm <- NA
  
  # Define stem slices in priority order
  slice_attempts <- list(
    c(1.20, 1.40), # 1st Attempt: Standard DBH range
    c(1.10, 1.50)  # 2nd Attempt: Rescue/Emergency slice
  )
  
  for (slice in slice_attempts) {
    z_min <- slice[1]
    z_max <- slice[2]
    
    slice_raw <- df[df$Z >= z_min & df$Z <= z_max, ]
    
    if (nrow(slice_raw) >= 10) {
      # 10cm Distance Filter (DBSCAN) to remove noise
      db <- dbscan::dbscan(slice_raw[, c("X", "Y")], eps = 0.10, minPts = 5)
      valid_mask <- db$cluster > 0
      
      if (sum(valid_mask) > 0) {
        counts <- sort(table(db$cluster[valid_mask]), decreasing = TRUE)
        main_cluster_id <- as.numeric(names(counts)[1])
        slice_clean <- slice_raw[db$cluster == main_cluster_id, ]
        
        if (nrow(slice_clean) >= 3) {
          # Convex Hull computation for stem geometry
          hull_idx <- chull(slice_clean$X, slice_clean$Y)
          hull_coords <- slice_clean[hull_idx, c("X", "Y")]
          hull_coords <- rbind(hull_coords, hull_coords[1, ]) 
          dbh_cm <- max(as.matrix(dist(hull_coords)), na.rm = TRUE) * 100
          
          # Validation Plot Generation
          min_x <- min(slice_clean$X); min_y <- min(slice_clean$Y)
          p_data <- slice_clean; p_data$Xr <- p_data$X - min_x; p_data$Yr <- p_data$Y - min_y
          h_data <- hull_coords; h_data$Xr <- h_data$X - min_x; h_data$Yr <- h_data$Y - min_y
          raw_data <- slice_raw; raw_data$Xr <- raw_data$X - min_x; raw_data$Yr <- raw_data$Y - min_y
          
          p <- ggplot() +
            geom_point(data = raw_data, aes(Xr, Yr), color = "grey85", shape = 3, size = 1.5) +
            geom_point(data = p_data, aes(Xr, Yr), color = "blue", size = 2, alpha = 0.9) +
            geom_path(data = h_data, aes(Xr, Yr), color = "red", linewidth = 1.2) +
            labs(title = paste0("Tree ID: ", tree_id), 
                 subtitle = paste0("Calculated DBH: ", round(dbh_cm, 1), " cm (Slice: ", z_min, "m - ", z_max, "m)"), 
                 x = "Relative X (m)", y = "Relative Y (m)") +
            coord_fixed() + theme_minimal()
          
          ggsave(file.path(plot_folder, paste0(tree_id, ".png")), p, width = 5, height = 5)
          
          # Exit loop successfully once DBH is calculated
          break
        }
      }
    }
  }
  
  # -------------------------------------------------------
  # B. COMPREHENSIVE METRICS CALCULATION
  # -------------------------------------------------------
  
  # 1. Base Statistics
  H_max <- max(las$Z, na.rm = TRUE)
  metrics <- list(
    Tree_ID = tree_id,
    DBH_cm = dbh_cm,
    Height = H_max,
    Mean_height = mean(las$Z, na.rm = TRUE),
    Roughness = sd(las$Z, na.rm = TRUE),
    p05 = quantile(las$Z, 0.05, na.rm = TRUE),
    p25 = quantile(las$Z, 0.25, na.rm = TRUE),
    p50 = quantile(las$Z, 0.50, na.rm = TRUE),
    p75 = quantile(las$Z, 0.75, na.rm = TRUE),
    p95 = quantile(las$Z, 0.95, na.rm = TRUE)
  )
  
  # 2. Point Density and Stratification
  n_tot <- npoints(las)
  metrics$Dens_point <- n_tot / area(las)
  metrics$Dens_above_p50 <- sum(las$Z > metrics$p50) / n_tot
  metrics$Dens_0_2m <- sum(las$Z >= 0 & las$Z < 2) / n_tot
  metrics$Dens_2_5m <- sum(las$Z >= 2 & las$Z < 5) / n_tot
  metrics$Dens_5m_up <- sum(las$Z >= 5) / n_tot
  
  # 3. Crown Architecture
  metrics$Crown_diam_x <- diff(range(las$X, na.rm = TRUE))
  metrics$Crown_diam_y <- diff(range(las$Y, na.rm = TRUE))
  metrics$Crown_diameter <- mean(c(metrics$Crown_diam_x, metrics$Crown_diam_y))
  metrics$Crown_depth <- diff(range(las$Z, na.rm = TRUE))
  
  # 4. Geometry (3D Convex Hull Volume, 2D Area)
  tryCatch({
    xyz <- as.matrix(las@data[, c("X", "Y", "Z")])
    xy  <- as.matrix(las@data[, c("X", "Y")])
    metrics$Envelope_vol <- convhulln(xyz, options = "FA")$vol
    metrics$Projected_area <- convhulln(xy, options = "FA")$area
  }, error = function(e) {
    metrics$Envelope_vol <- NA; metrics$Projected_area <- NA
  })
  
  # 5. Asymmetry and Inclination
  tryCatch({
    pca <- prcomp(las@data[, c("X", "Y")], center = TRUE, scale. = FALSE)
    metrics$Asymmetry <- ifelse(length(pca$sdev) >= 2, pca$sdev[1] / pca$sdev[2], NA)
    fit <- lm(Z ~ X + Y, data = las@data)
    horiz_slope <- sqrt(coef(fit)["X"]^2 + coef(fit)["Y"]^2)
    metrics$Inclination_deg <- atan(horiz_slope) * 180 / pi
  }, error = function(e) { metrics$Asymmetry <- NA; metrics$Inclination_deg <- NA })
  
  # 6. Shape Indices (Eigenvalues)
  tryCatch({
    coords <- scale(as.matrix(las@data[, c("X", "Y", "Z")]), center = T, scale = F)
    eig <- sort(eigen(cov(coords))$values, decreasing = T)
    metrics$Linearity <- (eig[1] - eig[2]) / eig[1]
    metrics$Planarity <- (eig[2] - eig[3]) / eig[1]
    metrics$Sphericity <- eig[3] / eig[1]
  }, error = function(e) {
    metrics$Linearity <- NA; metrics$Planarity <- NA; metrics$Sphericity <- NA
  })
  
  # 7. SLENDERNESS INDEX (H/D Ratio)
  if(!is.na(dbh_cm) && dbh_cm > 0) {
    metrics$Slenderness <- H_max / (dbh_cm / 100)
  } else {
    metrics$Slenderness <- NA
  }
  
  # 8. VOXEL-BASED BIOMASS VOLUME (10cm Voxel Size)
  tryCatch({
    res <- 0.10 # 10 cm voxel
    voxels <- las@data %>%
      mutate(Xv = floor(X/res), Yv = floor(Y/res), Zv = floor(Z/res)) %>%
      distinct(Xv, Yv, Zv)
    metrics$Biomass_vox <- nrow(voxels) * (res^3)
  }, error = function(e) metrics$Biomass_vox <- NA)
  
  # 9. SKEWNESS & KURTOSIS (Z Distribution)
  metrics$Z_Skewness <- tryCatch(skewness(las$Z), error = function(e) NA)
  metrics$Z_Kurtosis <- tryCatch(kurtosis(las$Z), error = function(e) NA)
  
  return(as.data.frame(metrics))
}

# =========================================================
# 4. EXECUTION AND LOOP
# =========================================================
cat("\n>>> PHASE 3: Data Processing <<<\n")

output_dir <- "../Outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

pb <- txtProgressBar(min = 0, max = length(las_files), style = 3)
results_list <- list()

for (i in 1:length(las_files)) {
  res <- extract_metrics_and_plot(las_files[i], plot_folder)
  if (!is.null(res)) results_list[[i]] <- res
  setTxtProgressBar(pb, i)
}
close(pb)

if (length(results_list) > 0) {
  Extracted_Data <- bind_rows(results_list)
  
  # Export to CSV format
  outfile <- file.path(output_dir, "COMPLETE_FINAL_Metrics.csv")
  write.csv(Extracted_Data, outfile, row.names = FALSE)
  
  cat("\n\n------------------------------------------------")
  cat("\nANALYSIS SUCCESSFULLY COMPLETED!")
  cat("\nResults saved to:", outfile)
  cat("\nCheck the 'Plots_Complete' folder to review cross-sections and DBH calculations.")
  cat("\n------------------------------------------------\n")
} else {
  stop("Error: No data could be extracted.")
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat(sprintf("\n=======================================================\n"))
cat(sprintf(" SCRIPT EXECUTION COMPLETED!\n"))
cat(sprintf(" Total execution time: %.2f minutes.\n", total_time))
cat(sprintf("=======================================================\n"))
