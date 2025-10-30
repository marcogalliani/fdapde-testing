# = ========================================================================== =
# - Script: plot_results.R
# - Desc: Provides visualization utilities for simulation results. Includes
#         quantitative summaries (e.g., RMSE, time) and qualitative comparisons
#         of reconstructed loadings, both at locations and on high-resolution
#         grids, for different functional PCA approaches.
# = ========================================================================== =


## Function: plot_quantitative_results
# - Args:
#   * loaded_results: list containing aggregated quantitative results and model info.
#       Expected fields include:
#         - $model_names, $model_labels, $model_colors
#         - $execution_time, $rmse (with nested lists for various metrics)
# - Desc:
#   Generates boxplots summarizing execution times and RMSE-based performance metrics
#   for all models. Separate figures are produced for reconstruction accuracy,
#   score orthogonality deviation, and component-wise RMSE for loadings and scores.
plot_quantitative_results <- function(loaded_results) {
  ## Get models details
  model_names <- loaded_results$model_names
  model_labels <- loaded_results$model_labels
  model_colors <- loaded_results$model_colors
  
  library(scales)
  library(patchwork)
  
  gcv_scores <- loaded_results$gcv_scores
  mse <- loaded_results$mse
  
  n_comp <- ncol(gcv_scores[[1]])
  
  global_gcv_max <- 3 * do.call(max, lapply(gcv_scores, max))
  global_gcv_min <- do.call(min, lapply(gcv_scores, min))
  global_mse_max <- do.call(max, lapply(mse, max))
  global_mse_min <- do.call(min, lapply(mse, min))
  global_y_max <- max(global_gcv_max, global_mse_max)
  global_y_min <- min(global_gcv_min, global_mse_min)
  
  plot_list <- list()
  
  metric_labels <- c(setNames(paste0("fPC", 1:n_comp), paste0("fPC", 1:n_comp)),
                     mse = expression(paste("MSE(",lambda,")",sep="")),
                     sum_gcv = expression(paste("GCV(",lambda,")",sep="")))
  
  for (model_idx in seq_along(model_names)) {
    model_name <- model_names[model_idx]
    base_color <- test_options$model_colors[model_idx]
    
    # Define color shades for fPCs (light → dark)
    fpc_colors <- colorRampPalette(c("grey90", base_color))(n_comp + 2)[-1]
    
    # Define color mapping
    curve_colors <- c(
      setNames(fpc_colors[1:n_comp], paste0("fPC", 1:n_comp)),
      mse = base_color,
      sum_gcv = base_color
    )
    
    # Define line types
    line_types <- c(
      setNames(rep(c("dashed", "dotdash", "longdash", "dotted"), 
                   length.out = n_comp), paste0("fPC", 1:n_comp)),
      sum_gcv = "solid",
      mse = "twodash"
    )
    
    # Define point shapes (different marker for each line)
    point_shapes <- c(
      setNames(rep(c(21, 22, 23, 24, 25, 4, 8), length.out = n_comp), paste0("fPC", 1:n_comp)),
      mse = 19,         # solid circle
      sum_gcv = 17      # solid triangle
    )
    
    current_gcv_matrix <- gcv_scores[[model_name]]
    
    df_plot <- as.data.frame(current_gcv_matrix)
    colnames(df_plot) <- paste0("fPC", 1:n_comp)
    df_plot$lambda <- test_options$regularization$lambda_grid
    df_plot$mse <- mse[[model_name]]
    df_plot$sum_gcv <- rowSums(current_gcv_matrix)
    
    df_long <- df_plot %>%
      pivot_longer(cols = c(starts_with("fPC"), "mse", "sum_gcv"),
                   names_to = "Metric",
                   values_to = "GCV_Score")
    
    # Find minima for each line
    min_gcv_points <- df_long %>%
      group_by(Metric) %>%
      filter(GCV_Score == min(GCV_Score))
    
    plot_list[[model_name]] <- ggplot(df_long, aes(x = lambda, y = GCV_Score,
                                                   color = Metric, linetype = Metric)) +
      geom_line(linewidth = 1) +
      geom_point(data = min_gcv_points,
                 aes(x = lambda, y = GCV_Score, shape = Metric),
                 size = 3, stroke = 1.5, fill = "white") +
      scale_x_log10(labels = label_scientific()) +
      scale_y_log10(limits = c(global_y_min, global_y_max),
                    labels = label_scientific()) +
      scale_color_manual(values = curve_colors) +
      scale_linetype_manual(values = line_types, labels = metric_labels) +
      scale_shape_manual(values = point_shapes) +
      labs(
        title = test_options$model_labels[model_idx],
        x = expression(lambda),
        y = NULL,
        color = NULL,
        linetype = NULL,
        shape = NULL
      ) +
      guides(color = "none", shape="none") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.key.width = unit(2.5, "cm"),
        legend.spacing.x = unit(0.8, "cm"),
        legend.position = "top"
      )
  }  
  # Combine plots
  final_plot <- Reduce(`+`, plot_list) +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  )
  print(final_plot)
}

## Function: plot_qualitative_results
# - Args:
#   * quantitative_results: list containing quantitative summaries (for selecting representative samples)
#   * qualitative_results: list containing qualitative outputs for each model, including:
#       - $loadings, $loadings_locs, $loadings_HR
#       - $loadings_true, $loadings_true_locs, $loadings_true_HR
#       - $domain, $grid, $locations, $nodes, and $boundary
# - Desc:
#   Produces qualitative visual comparisons of true and reconstructed loadings across
#   models and principal components. For each method, three representative replicates
#   are selected based on RMSE quantiles (min, median, max). Each component’s fields
#   are plotted at nodes, evaluation locations, and high-resolution grids, both with
#   and without isolines. Results are arranged in labeled grids for clarity.
plot_qualitative_results <- function(quantitative_results, qualitative_results) {
  ## Extract components
  scores <- qualitative_results$scores
  loadings <- qualitative_results$loadings
  loadings_locs <- qualitative_results$loadings_locs
  loadings_HR <- qualitative_results$loadings_HR
  
  loadings_true <- qualitative_results$loadings_true
  loadings_true_locs <- qualitative_results$loadings_true_locs
  loadings_true_HR <- qualitative_results$loadings_true_HR
  
  ## Metadata
  domain <- qualitative_results$domain
  boundary <- domain$boundary
  n_comp <- ncol(loadings_true)
  
  model_names <- quantitative_results$model_names
  model_labels <- qualitative_results$model_labels
  
  ## Labels for grid
  labels_cols <- c("True", model_labels)
  labels_rows <- paste0("f", 1:n_comp)
  
  ## Plot storage
  plot_list <- list()
  plot_list_locs <- list()
  plot_list_HR_clean <- list()
  plot_list_HR <- list()
    
  
  ## Build grid row by row (one per component)
  for (i in seq_len(n_comp)) {
    ## Value limits for this component
    limits <- range(loadings_true_HR[, i])
    breaks <- seq(limits[1], limits[2], length = 10)
    
    ## Start each row with TRUE component plots
    row_plots <- list()
    row_plots_locs <- list()
    row_plots_HR_clean <- list()
    row_plots_HR <- list()
    
    ## True component
    row_plots[[1]] <- plot.field_points(
      qualitative_results$nodes, loadings_true[, i],
      boundary = boundary, size = 1.5, LEGEND = FALSE
    ) + std_plot_settings_fields()
    
    row_plots_locs[[1]] <- plot.field_points(
      qualitative_results$locations, loadings_true_locs[, i],
      boundary = boundary, size = 1.5, LEGEND = FALSE
    ) + std_plot_settings_fields()
    
    row_plots_HR_clean[[1]] <- plot.field_tile(
      qualitative_results$grid, loadings_true_HR[, i],
      boundary = boundary, LEGEND = FALSE
    ) + std_plot_settings_fields()
    
    row_plots_HR[[1]] <- plot.field_tile(
      qualitative_results$grid, loadings_true_HR[, i],
      boundary = boundary, limits = limits, breaks = breaks, LEGEND = FALSE
    ) + std_plot_settings_fields()
    
    ## Now append estimated components from each model
    for (m in seq_along(model_names)) {
      name_model <- model_names[m]
      
      row_plots[[m + 1]] <- plot.field_points(
        qualitative_results$nodes, loadings[[name_model]][[1]][, i],
        boundary = boundary, size = 1.5
      ) + std_plot_settings_fields()
      
      row_plots_locs[[m + 1]] <- plot.field_points(
        qualitative_results$locations, loadings_locs[[name_model]][[1]][, i],
        boundary = boundary, size = 1.5
      ) + std_plot_settings_fields()
      
      row_plots_HR_clean[[m + 1]] <- plot.field_tile(
        qualitative_results$grid, loadings_HR[[name_model]][[1]][, i],
        boundary = boundary
      ) + std_plot_settings_fields()
      
      row_plots_HR[[m + 1]] <- plot.field_tile(
        qualitative_results$grid, loadings_HR[[name_model]][[1]][, i],
        boundary = boundary, limits = limits, breaks = breaks
      ) + std_plot_settings_fields()
    }
    
    ## Combine the row
    plot_list[[i]] <- arrangeGrob(grobs = row_plots, ncol = length(labels_cols))
    plot_list_locs[[i]] <- arrangeGrob(grobs = row_plots_locs, ncol = length(labels_cols))
    plot_list_HR_clean[[i]] <- arrangeGrob(grobs = row_plots_HR_clean, ncol = length(labels_cols))
    plot_list_HR[[i]] <- arrangeGrob(grobs = row_plots_HR, ncol = length(labels_cols))
  }
  
  ## Combine rows vertically
  final_grid <- arrangeGrob(grobs = plot_list, nrow = n_comp)
  final_grid_locs <- arrangeGrob(grobs = plot_list_locs, nrow = n_comp)
  final_grid_HR_clean <- arrangeGrob(grobs = plot_list_HR_clean, nrow = n_comp)
  final_grid_HR <- arrangeGrob(grobs = plot_list_HR, nrow = n_comp)
  
  ## Add labels
  grid_final <- labled_plots_grid(final_grid, "All Models", labels_cols, labels_rows)
  grid_final_locs <- labled_plots_grid(final_grid_locs, "All Models", labels_cols, labels_rows)
  grid_final_HR_clean <- labled_plots_grid(final_grid_HR_clean, "All Models", labels_cols, labels_rows)
  grid_final_HR <- labled_plots_grid(final_grid_HR, "All Models", labels_cols, labels_rows)
  
  ## Display
  grid.arrange(grid_final_locs)
  grid.arrange(grid_final)
  grid.arrange(grid_final_HR_clean)
  grid.arrange(grid_final_HR)
}
