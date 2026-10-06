# ==============================================================================
# SCRIPT 03: CLUSTERING AND SEMANTIC NETWORK (TLS vs TreMs)
# ==============================================================================
# Author: Antonio di Biase (a.dibiase5@studenti.unimol.it)
# Description:
# This script performs a meta-analysis of the predictive modeling results.
# 1. Hierarchical Clustering: Analyzes the co-occurrence of LiDAR metrics across 
#    all predictive models to identify structural functional groups (Ward.D2).
# 2. Interactive UI (Shiny): Allows the user to semantically validate and name 
#    the identified LiDAR clusters through a pop-up interface.
# 3. Ecological Network Graph: Generates a bipartite network (TreMs and TLS metrics) 
#    to visualize the complex interrelations between forest microhabitats and 
#    3D structural traits.
# ==============================================================================

rm(list = ls())
gc()

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, stringr, ggplot2, ggdendro, dendextend, pheatmap, reshape2, igraph, ggraph, ggforce, shiny, miniUI, tidygraph, cluster, ggnewscale)

# Create output folders explicitly
dir_clust <- "../Outputs/Plot_Funzionali_LiDAR"
dir_net <- "../Outputs/Plot_Network"
dir.create(dir_clust, showWarnings = FALSE, recursive = TRUE) 
dir.create(dir_net, showWarnings = FALSE, recursive = TRUE)          

# ==============================================================================
# PHASE 1: HIERARCHICAL CLUSTERING (From Excel file)
# ==============================================================================
cat("\n--- PHASE 1: HIERARCHICAL CLUSTERING ---\n")

# Dynamic File Loading (checks Outputs first, then Data)
excel_out <- "../Outputs/Complete_Model_Results_with_Pvalues.xlsx"
excel_data <- "../Data/Results_Model_Pvalue.xlsx"

if (file.exists(excel_out)) {
  file_excel_res <- excel_out
} else if (file.exists(excel_data)) {
  file_excel_res <- excel_data
} else {
  stop("ERROR: The model results Excel file could not be found.")
}

df_results <- read_excel(file_excel_res)

cat("Processing extracted variables...\n")
if ("Variabili_Usate" %in% names(df_results)) {
  df_results <- df_results %>% rename(Used_Variables = Variabili_Usate)
}

df_target_vars <- df_results %>% 
  group_by(Target) %>% 
  summarise(
    Variables = first(Used_Variables), 
    Max_R2 = max(R2, na.rm = TRUE), 
    .groups = 'drop'
  )

all_metrics <- unique(unlist(str_split(df_target_vars$Variables, ",\\s*")))
all_metrics <- all_metrics[!is.na(all_metrics) & all_metrics != "NA"]

presence_matrix <- matrix(0, nrow = length(all_metrics), ncol = nrow(df_target_vars), dimnames = list(all_metrics, df_target_vars$Target))

for(i in 1:nrow(df_target_vars)) {
  tgt <- df_target_vars$Target[i]
  vars <- unlist(str_split(df_target_vars$Variables[i], ",\\s*"))
  if(length(vars) > 0 && !all(is.na(vars)) && vars[1] != "NA") { presence_matrix[vars, tgt] <- 1 }
}

clean_matrix <- presence_matrix[rowSums(presence_matrix) > 0, colSums(presence_matrix) > 0, drop = FALSE]

if(nrow(clean_matrix) >= 2) {
  dist_mat <- dist(clean_matrix, method = "binary")
  hc <- hclust(dist_mat, method = "ward.D2")
  
  heights <- hc$height
  heights_diff <- diff(heights)
  
  cut_index <- which.max(tail(heights_diff, n = min(10, length(heights_diff))))
  DYNAMIC_CLUSTER_NUM <- length(heights) - (length(heights) - length(tail(heights_diff, n = min(10, length(heights_diff)))) + cut_index) + 1
  DYNAMIC_CLUSTER_NUM <- max(2, min(8, DYNAMIC_CLUSTER_NUM))
  
  cat(sprintf("\n[!] The algorithm automatically identified %d optimal Clusters.\n", DYNAMIC_CLUSTER_NUM))
  
  dendro <- as.dendrogram(hc) %>% 
    color_branches(k = DYNAMIC_CLUSTER_NUM) %>% 
    set("branches_lwd", 2.3) %>% 
    set("labels_cex", 0.9)
  
  png(file.path(dir_clust, "1_Dendrogramma_Cluster.png"), width = 4666, height = 3666, res = 500)
  par(mar = c(5, 1, 3, 15))
  plot(dendro, horiz = TRUE, main = "")
  title(xlab = "Dissimilarity (Ward.D2)", line = 2.5, cex.lab = 1.1)
  
  x_max_plot <- grconvertX(1, "ndc", "user")
  length_prop <- 0.65 
  
  x_left_box <- 0
  x_right_box <- x_max_plot * length_prop
  
  branch_colors <- get_leaves_branches_col(dendro)
  cluster_blocks <- rle(branch_colors)
  current_y <- 0.5 
  
  for (i in 1:length(cluster_blocks$lengths)) {
    metrics_in_cluster <- cluster_blocks$lengths[i]
    cluster_color <- cluster_blocks$values[i]
    y_bottom <- current_y; y_top <- current_y + metrics_in_cluster
    bg_color <- adjustcolor(cluster_color, alpha.f = 0.15)
    
    rect(xleft = x_left_box, ybottom = y_bottom, xright = x_right_box, ytop = y_top, col = bg_color, border = cluster_color, lwd = 1.8, xpd = NA) 
    current_y <- y_top
  }
  dev.off()
  
  assignments <- cutree(hc, k = DYNAMIC_CLUSTER_NUM)
  
  df_sil <- tryCatch({
    sil <- cluster::silhouette(assignments, dist_mat)
    data.frame(
      Metric = labels(dist_mat),
      Alternative_Cluster_ID = as.character(sil[, "neighbor"]),
      Silhouette_Score = round(sil[, "sil_width"], 3)
    )
  }, error = function(e) {
    data.frame(Metric = labels(dist_mat), Alternative_Cluster_ID = NA, Silhouette_Score = NA)
  })
  
  df_final <- data.frame(
    Metric = names(assignments), 
    Cluster = as.factor(assignments), 
    Frequency = rowSums(clean_matrix)
  ) %>%
    left_join(df_sil, by = "Metric")
  
  bar_plot <- ggplot(df_final, aes(x = reorder(Metric, Frequency), y = Frequency, fill = Cluster)) + 
    geom_bar(stat = "identity", color = "black", alpha = 0.8) + 
    coord_flip() + theme_minimal() + scale_fill_brewer(palette = "Set1") + 
    labs(x = "TLS Metrics", y = "Number of selections") + 
    theme(legend.position = "bottom")
  ggsave(file.path(dir_clust, "2_Frequency_BarPlot.png"), plot = bar_plot, width = 10, height = 8, dpi = 500, bg = "white")
  
  # ==============================================================================
  # PHASE 2: SHINY UI FOR SEMANTIC VALIDATION
  # ==============================================================================
  cat("\nStarting the semantic validation interface...\n")
  
  ui_cluster_data <- df_final %>% arrange(Cluster, desc(Frequency)) %>% group_by(Cluster) %>% summarise(Included_Metrics = paste(Metric, collapse = ", "), Total_Freq = sum(Frequency), .groups = 'drop')
  
  ui_cluster_data$Suggested_Name <- paste("Cluster", ui_cluster_data$Cluster)
  for(i in 1:nrow(ui_cluster_data)) {
    if(grepl("(?i)DBH", ui_cluster_data$Included_Metrics[i])) ui_cluster_data$Suggested_Name[i] <- "Allometry and Dimension"
    if(grepl("(?i)planarity", ui_cluster_data$Included_Metrics[i])) ui_cluster_data$Suggested_Name[i] <- "Complexity and Anomalies"
    if(grepl("(?i)crown", ui_cluster_data$Included_Metrics[i])) ui_cluster_data$Suggested_Name[i] <- "Irregularity and Spatial Envelope"
    if(grepl("(?i)dens_0_2m", ui_cluster_data$Included_Metrics[i])) ui_cluster_data$Suggested_Name[i] <- "Basal Density"
  }
  
  ui <- miniPage(
    gadgetTitleBar("Semantic Validation of TLS Clusters", right = miniTitleBarButton("done", "Confirm and Execute", primary = TRUE)),
    miniContentPanel(
      h4("Define Structural Families"),
      p("Check the included metrics and confirm the Cluster category. Clicking 'Confirm' will generate the plots and save the tables."), hr(),
      uiOutput("input_panel")
    )
  )
  
  server <- function(input, output, session) {
    output$input_panel <- renderUI({
      lapply(1:nrow(ui_cluster_data), function(i) {
        wellPanel(
          h5(strong(paste("CLUSTER ID:", ui_cluster_data$Cluster[i], "| Selections:", ui_cluster_data$Total_Freq[i]))),
          p(em("Metrics:", ui_cluster_data$Included_Metrics[i])),
          textInput(inputId = paste0("cluster_name_", ui_cluster_data$Cluster[i]), label = "Validated Name:", value = ui_cluster_data$Suggested_Name[i], width = "100%")
        )
      })
    })
    
    observeEvent(input$done, {
      validated_names <- sapply(as.character(ui_cluster_data$Cluster), function(c_id) { input[[paste0("cluster_name_", c_id)]] })
      name_map <- setNames(validated_names, ui_cluster_data$Cluster)
      
      df_final_agg <- df_final %>% mutate(
        Validated_Name = name_map[as.character(Cluster)],
        Alternative_Family = ifelse(is.na(Alternative_Cluster_ID), NA, name_map[as.character(Alternative_Cluster_ID)])
      )
      
      cat("\n--- SAVING TABLES AND PLOTS ---\n")
      
      assignment_table <- df_final_agg %>% select(Metric, Cluster_ID = Cluster, Functional_Family = Validated_Name, Frequency, Most_Similar_Alternative = Alternative_Family, Silhouette_Score) %>% arrange(Functional_Family, desc(Frequency))
      write.csv(assignment_table, file.path(dir_clust, "3_Clustering_table.csv"), row.names = FALSE)
      
      # --- HEATMAP ---
      cat("Generating Heatmap...\n")
      row_annotations <- df_final_agg %>% select(Metric, Validated_Name) %>% rename(Cluster = Validated_Name)
      rownames(row_annotations) <- row_annotations$Metric; row_annotations$Metric <- NULL
      num_categories <- length(unique(row_annotations$Cluster))
      color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
      heatmap_colors <- setNames(color_palette[1:num_categories], unique(row_annotations$Cluster))
      
      png(file.path(dir_clust, "4_Heatmap_Co_Occurrence.png"), width = 12, height = 8, units = "in", res = 500)
      pheatmap(clean_matrix, 
               color = c("white", "#2C3E50"), 
               breaks = c(-0.5, 0.5, 1.5), 
               annotation_row = row_annotations, 
               annotation_colors = list(Cluster = heatmap_colors), 
               cluster_rows = TRUE, 
               cluster_cols = TRUE, 
               treeheight_row = 100, 
               fontsize_row = 10, 
               fontsize_col = 8, 
               border_color = "white",
               cellwidth = 12,        
               cellheight = 12,
               legend = FALSE)
      dev.off()
      
      # --- BOXPLOT ---
      cat("Generating Dynamic Boxplot...\n")
      boxplot_data <- df_final_agg %>% select(Metric, Cluster_ID = Cluster, Frequency, Cluster_Name = Validated_Name)
      level_order <- boxplot_data %>% group_by(Cluster_Name) %>% summarise(Median_Val = median(Frequency)) %>% arrange(desc(Median_Val)) %>% pull(Cluster_Name)
      boxplot_data$Cluster_Name <- factor(boxplot_data$Cluster_Name, levels = rev(level_order))
      
      importance_plot <- ggplot(boxplot_data, aes(x = Cluster_Name, y = Frequency, fill = Cluster_Name)) + 
        geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black", linewidth = 0.6) + 
        geom_jitter(width = 0.15, size = 2.5, color = "black", alpha = 0.8) + 
        coord_flip() + theme_classic() + scale_fill_viridis_d(option = "turbo", direction = -1) + 
        labs(x = "", y = "Selection Frequency in the Models") + 
        theme(legend.position = "none")
      ggsave(file.path(dir_clust, "5_Hierarchy_Importance_LiDAR.png"), plot = importance_plot, width = 10, height = 6, dpi = 500, bg = "white") 
      
      # ==============================================================================
      # --- NETWORK ---
      # ==============================================================================
      cat("\n--- File Request for NETWORK ---\n")
      
      # Dynamic file loading for text report
      txt_out <- "../Outputs/Best_Models_CV_Summary.txt"
      txt_data <- "../Data/Model_Performance_Report.txt"
      
      if (file.exists(txt_out)) {
        file_txt <- txt_out
      } else if (file.exists(txt_data)) {
        file_txt <- txt_data
      } else {
        stop("ERROR: Model Summary report not found in Outputs or Data folders.")
      }
      
      file_field <- "../Data/TreMs_Field_Surveys.xlsx"
      
      lines_txt <- readLines(file_txt)
      df_field <- read_excel(file_field)
      
      trems_64 <- grep("^[A-Z]{2}[0-9]+", names(df_field), value = TRUE)
      if(length(trems_64) > 0) {
        df_field <- df_field %>% mutate(across(all_of(trems_64), as.numeric)) %>% mutate(across(all_of(trems_64), ~ tidyr::replace_na(., 0)))
        df_field$Macro_TreMs_Abundance <- rowSums(df_field[, trems_64], na.rm = TRUE)
        df_field$Macro_TreMs_Richness <- rowSums(df_field[, trems_64] > 0, na.rm = TRUE)
        
        for(fp in unique(substr(trems_64, 1, 2))) { df_field[[paste0("Form_", fp)]] <- rowSums(df_field[, grep(paste0("^", fp), names(df_field)), drop=FALSE], na.rm=TRUE) }
        for(gp in unique(substr(trems_64, 1, 3))) { df_field[[paste0("Group_", gp)]] <- rowSums(df_field[, grep(paste0("^", gp), names(df_field)), drop=FALSE], na.rm=TRUE) }
      }
      
      edges_net <- data.frame()
      parsed_targets <- c()
      
      for (line in lines_txt) {
        if (grepl("TARGET:", line)) {
          target <- trimws(sub(".*TARGET: (.*?) \\|.*", "\\1", line))
          parsed_targets <- c(parsed_targets, target)
          
          # Supports both English and Italian output formats
          if (grepl("CORRELATIONS:|CORRELAZIONI:", line, ignore.case = TRUE)) {
            corr_str <- sub(".*(?:CORRELATIONS|CORRELAZIONI): (.*?)( \\|.*|$)", "\\1", line, ignore.case = TRUE)
            elements <- unlist(strsplit(corr_str, ","))
            
            for (elem in elements) {
              elem <- trimws(elem)
              metric_name <- trimws(sub("\\(.*", "", elem))
              val_str <- sub(".*\\((.*?)\\).*", "\\1", elem)
              
              num_val <- as.numeric(val_str)
              if (is.na(num_val)) {
                corr_strength <- 0.05; corr_sign <- "Positive"
              } else {
                corr_strength <- abs(num_val); corr_sign <- ifelse(num_val < 0, "Negative", "Positive") 
              }
              
              edges_net <- rbind(edges_net, data.frame(from = metric_name, to = target, Correlation = corr_strength, Sign = corr_sign))
            }
          }
        }
      }
      
      trems_frequencies <- data.frame(Clean_Target = unique(parsed_targets), Abundance = NA)
      
      for(j in 1:nrow(trems_frequencies)) {
        t_clean <- trems_frequencies$Clean_Target[j]
        
        if (t_clean %in% names(df_field)) {
          trems_frequencies$Abundance[j] <- sum(df_field[[t_clean]] > 0, na.rm = TRUE)
        } else {
          real_col <- names(df_field)[grepl(paste0(t_clean, "$"), names(df_field), ignore.case = TRUE)]
          if(length(real_col) > 0) {
            trems_frequencies$Abundance[j] <- sum(df_field[[real_col[1]]] > 0, na.rm = TRUE)
          } else {
            trems_frequencies$Abundance[j] <- 1 
          }
        }
      }
      
      # 3. Graph Creation
      graph_pro <- as_tbl_graph(edges_net, directed = FALSE) %>%
        mutate(Type = ifelse(name %in% boxplot_data$Metric, "TLS metrics", "TreMs")) %>%
        left_join(boxplot_data, by = c("name" = "Metric")) %>%
        left_join(trems_frequencies, by = c("name" = "Clean_Target")) %>%
        mutate(
          Cluster_Name = ifelse(Type == "TreMs", "Microhabitat (TreMs)", 
                                ifelse(is.na(Cluster_ID), "Unclassified", name_map[as.character(Cluster_ID)])),
          node_label = ifelse(Type == "TreMs", 
                              toupper(gsub("^Group_|^Form_|^Macro_Trems_", "", name, ignore.case = TRUE)), 
                              name)
        )
      
      # --- HARMONIC AND CONTRASTING COLOR PALETTE ---
      base_colors <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999")
      network_colors <- setNames(base_colors[1:length(name_map)], name_map)
      network_colors["Microhabitat (TreMs)"] <- "grey35"
      network_colors["Unclassified"] <- "grey80" 
      
      category_order <- c("Microhabitat (TreMs)", setdiff(names(network_colors), c("Microhabitat (TreMs)", "Unclassified")), "Unclassified")
      
      # --- DYNAMIC LEGEND CONTROLS ---
      present_categories <- intersect(category_order, unique(igraph::V(graph_pro)$Cluster_Name))
      
      legend_shapes <- ifelse(present_categories == "Microhabitat (TreMs)", 15, 19)
      legend_sizes <- ifelse(present_categories == "Microhabitat (TreMs)", 5.5, 5.5)
      
      create_range_labels <- function(breaks) {
        b <- na.omit(breaks)
        if(length(b) < 2) return(as.character(breaks))
        lbls <- c()
        for (i in 1:(length(b)-1)) {
          lbls <- c(lbls, sprintf("%d - %d", as.integer(b[i]), as.integer(b[i+1]-1)))
        }
        lbls <- c(lbls, sprintf("\u2265 %d", as.integer(b[length(b)])))
        res <- rep(NA, length(breaks))
        res[!is.na(breaks)] <- lbls
        return(res)
      }
      
      # 4. Generate Network Plot
      plot_final <- ggraph(graph_pro, layout = "stress") +
        
        geom_mark_hull(aes(x, y, group = Cluster_Name, fill = Cluster_Name, label = Cluster_Name), 
                       data = function(df) df %>% filter(Type == "TLS metrics" & Cluster_Name != "Unclassified"), 
                       alpha = 0.15, color = NA, show.legend = FALSE, concavity = 5, expand = unit(3, "mm"), label.fontsize = 12) +
        
        geom_edge_link(aes(edge_width = Correlation, edge_color = Sign), alpha = 0.55) + 
        
        geom_node_point(aes(filter = Type == "TreMs", color = Cluster_Name, shape = Type, size = Abundance), stroke = 0.5, alpha = 0.70) +
        scale_size_continuous(
          range = c(4, 14), 
          name = "TreMs Abundance", 
          breaks = c(15, 30, 50, 100), 
          labels = create_range_labels, 
          guide = guide_legend(
            order = 3, 
            override.aes = list(shape = 15) 
          )
        ) +
        
        ggnewscale::new_scale("size") +
        
        geom_node_point(aes(filter = Type == "TLS metrics", color = Cluster_Name, shape = Type, size = Frequency), stroke = 0.5, alpha = 0.80) +
        scale_size_continuous(
          range = c(1, 10), 
          name = "TLS Metric Frequency", 
          labels = create_range_labels, 
          guide = guide_legend(order = 4)
        ) +
        
        geom_node_text(aes(label = node_label), repel = TRUE, size = 3.2, fontface = "bold", point.padding = unit(3, "mm"), min.segment.length = 0, max.overlaps = Inf, bg.color = "white", bg.r = 0.15) + 
        
        scale_shape_manual(values = c("TLS metrics" = 19, "TreMs" = 15), name = "Type") +
        scale_color_manual(values = network_colors, breaks = category_order, name = "Categories") + 
        scale_fill_manual(values = network_colors, breaks = category_order, name = "Categories") + 
        scale_edge_width_continuous(range = c(0.2, 2.5), name = "Strength of correlation") +
        scale_edge_color_manual(values = c("Positive" = "grey20", "Negative" = "red"), name = "Direction of correlation") +
        
        theme_void() + 
        theme( 
          legend.position = "right",
          legend.box = "vertical", 
          legend.box.just = "left",
          legend.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
          legend.text = element_text(size = 12),                
          legend.title = element_text(size = 14, face = "bold"),
          legend.key.height = unit(0.8, "cm"),
          legend.background = element_blank(), 
          legend.box.background = element_rect(fill = alpha("grey60", 0.3), color = NA), 
          legend.box.margin = ggplot2::margin(15, 15, 15, 15) 
        ) +
        
        guides(
          shape = guide_legend(order = 1, override.aes = list(size = 5)), 
          
          color = guide_legend(
            order = 2, 
            override.aes = list(shape = legend_shapes, size = legend_sizes)
          ), 
          
          fill = guide_legend(
            order = 2, 
            override.aes = list(shape = legend_shapes, size = legend_sizes)
          ), 
          
          edge_width = guide_legend(order = 5), 
          edge_color = guide_legend(order = 6)
        )
      
      ggsave(file.path(dir_net, "6_Ecological_Network.png"), plot = plot_final, width = 14, height = 11, dpi = 500, bg = "white") 
      
      cat("\nFULL PIPELINE SUCCESSFULLY COMPLETED!\n")
      stopApp(invisible()) 
    }) 
  } 
  
  runGadget(ui, server, viewer = dialogViewer("Semantic Validation", width = 800, height = 700))
}

