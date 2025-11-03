# = ========================================================================== =
# - Script: plot_results.R
# - Desc: Provides visualization utilities for simulation results. Includes
#         quantitative summaries (e.g., RMSE, time) and qualitative comparisons
#         of reconstructed loadings, both at locations and on high-resolution
#         grids, for different functional PCA approaches.
# = ========================================================================== =

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
plot_qualitative_results <- function(qualitative_results) {
  ## Extract components
  scores <- qualitative_results$scores
  loadings <- qualitative_results$loadings
  loadings_locs <- qualitative_results$loadings_locs
  loadings_HR <- qualitative_results$loadings_HR
  
  ## Metadata
  domain <- qualitative_results$domain
  boundary <- domain$boundary
  n_comp <- ncol(scores[[1]])
  
  model_names <- qualitative_results$model_names
  model_labels <- qualitative_results$model_labels
  
  ## Labels for grid
  labels_cols <- model_labels
  labels_rows <- paste0("f", 1:n_comp)
  
  ## Plot storage
  plot_list <- list()
  plot_list_locs <- list()
  plot_list_HR_clean <- list()
  plot_list_HR <- list()
    
  
  ## Build grid row by row (one per component)
  for (i in seq_len(n_comp)) {
    ## Value limits for this component
    limits <- range(loadings_HR[[1]][, i],na.rm =T)
    breaks <- seq(limits[1], limits[2], length = 10)
    
    ## Start each row with TRUE component plots
    row_plots <- list()
    row_plots_locs <- list()
    row_plots_HR_clean <- list()
    row_plots_HR <- list()
    
    
    ## Now append estimated components from each model
    for (m in seq_along(model_names)) {
      name_model <- model_names[m]
      
      row_plots[[m]] <- plot.field_points(
        qualitative_results$nodes, loadings[[name_model]][, i],
        boundary = boundary, size = 1.5
      ) + std_plot_settings_fields()
      
      row_plots_locs[[m]] <- plot.field_points(
        qualitative_results$locations, loadings_locs[[name_model]][, i],
        boundary = boundary, size = 1.5
      ) + std_plot_settings_fields()
      
      row_plots_HR_clean[[m]] <- plot.field_tile(
        qualitative_results$grid, loadings_HR[[name_model]][, i],
        boundary = boundary
      ) + std_plot_settings_fields()
      
      row_plots_HR[[m]] <- plot.field_tile(
        qualitative_results$grid, loadings_HR[[name_model]][, i],
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