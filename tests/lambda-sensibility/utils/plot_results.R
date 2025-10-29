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
  ## Time ----
  times <- loaded_results$execution_time
  indexes <- which(!is.nan(colSums(times[, model_names])))
  plot <- plot.grouped_boxplots(
    times[, c("Group", names(indexes))],
    values_name = "Time [seconds]",
    group_name = "",
    group_labels = "",
    subgroup_name = "Approaches",
    subgroup_labels = model_labels[indexes],
    subgroup_colors = model_colors[indexes]
  ) + std_plot_settings() + ggtitle("Time")
  print(plot)
  ## RMSE ----
  ## Overall measures
  rmses <- loaded_results$rmse
  names <- c("reconstruction_locs", "scores_orth")
  titles <- c("Reconstruction at locations", "Deviation from scores orthogonality")

  for (i in seq_along(names)) {
    name <- names[i]
    title <- titles[i]
    indexes <- which(!is.nan(colSums(rmses[[name]][, model_names])))
    plot <- plot.grouped_boxplots(
      rmses[[name]][, c("Group", names(indexes))],
      values_name = "RMSE",
      group_name = "",
      group_labels = "",
      subgroup_name = "Approaches",
      subgroup_labels = model_names[indexes],
      subgroup_colors = model_colors[indexes]
    ) + std_plot_settings() + ggtitle(title)
    print(plot)
  }

  ## Component-by-component measures
  names <- c("loadings_locs", "scores")
  titles <- c("Loadings at locations", "Scores")

  for (i in seq_along(names)) {
    name <- names[i]
    title <- titles[i]
    indexes <- which(!is.nan(colSums(rmses[[name]][, model_names])))
    plot <- plot.grouped_boxplots(
      rmses[[name]][, c("Group", names(indexes))],
      values_name = "RMSE",
      subgroup_name = "Approaches",
      subgroup_labels = model_names[indexes],
      subgroup_colors = model_colors[indexes]
    ) + std_plot_settings() + ggtitle(title)
    print(plot)
  }

  ## Angles ----
  angles <- loaded_results$angles
  ## Component by component measures
  names <- c("components_m")
  titles <- c("Angle between true and estimated loading")
  for (i in seq_along(names)) {
    name <- names[i]
    title <- titles[i]
    indexes <- which(!is.nan(colSums(angles[[name]][, model_names])))
    plot <- plot.grouped_boxplots(
      angles[[name]][, c("Group", names(indexes))],
      values_name = "Angle",
      subgroup_name = "Approaches",
      subgroup_labels = model_names[indexes],
      subgroup_colors = model_colors[indexes]
    ) + std_plot_settings() + ggtitle(title)
    print(plot)
  } 
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
  ## Debugging helpers
  # quantitative_results <- loaded_qnt_results
  # qualitative_results  <- loaded_qlt_results

  ## Get fitted quantities
  scores <- qualitative_results$scores
  loadings <- qualitative_results$loadings
  loadings_locs <- qualitative_results$loadings_locs
  loadings_HR <- qualitative_results$loadings_HR

  ## Get true loadings
  loadings_true <- qualitative_results$loadings_true
  loadings_true_locs <- qualitative_results$loadings_true_locs
  loadings_true_HR <- qualitative_results$loadings_true_HR
  
  ## Get infos
  domain <- qualitative_results$domain
  n_comp <- ncol(loadings_true)
  boundary <- domain$boundary

  ## Models info
  model_names <- quantitative_results$model_names
  model_labels <- qualitative_results$model_labels

  ## Labels for plot grids
  labels_cols <- c("True", "Quantile 0", "Quantile 0.5", "Quantile 1")
  labels_rows <- paste0("f", 1:n_comp)

  ## Room for plots
  plots <- list()
  plots_locs <- list()
  plots_HR_clean <- list() # high-resolution without isolines
  plots_HR <- list() # high-resolution with isolines

  ## Generate figures for each model
  for (m in seq_along(model_names)) { # m <- 1

    ## Model details
    name_model <- model_names[m]
    label_model <- model_labels[m]

    ## Room for plots
    plot_list <- list()
    plot_list_locs <- list()
    plot_list_HR_clean <- list()
    plot_list_HR <- list()

    ## Select representative replicates (min, median, max RMSE)
    indexes <- tapply(
      quantitative_results$rmse$loadings_locs[[name_model]],
      quantitative_results$rmse$loadings_locs$Group,
      function(x) {
        sapply(quantile(x, c(0, 0.5, 1), na.rm = TRUE), function(q) which.min(abs(x - q)))
      }
    )

    ## Compute limits
    limits_HR <- apply(do.call(rbind, unlist(loadings_HR, recursive = FALSE)), 2, range)

    for (i in seq_len(n_comp)) {
      limits <- range(c(loadings_true_HR[, i], limits_HR[, i]))
      breaks <- seq(limits[1], limits[2], length = 10)

      ## True components
      plot_list[[4 * (i - 1) + 1]] <- plot.field_points(
        qualitative_results$nodes, loadings_true[, i],
        boundary = boundary, size = 1.5, LEGEND = FALSE
      ) + std_plot_settings_fields()

      plot_list_locs[[4 * (i - 1) + 1]] <- plot.field_points(
        qualitative_results$locations, loadings_true_locs[, i],
        boundary = boundary, size = 1.5, LEGEND = FALSE
      ) + std_plot_settings_fields()

      plot_list_HR_clean[[4 * (i - 1) + 1]] <- plot.field_tile(
        qualitative_results$grid, loadings_true_HR[, i],
        boundary = boundary, LEGEND = FALSE
      ) + std_plot_settings_fields()

      plot_list_HR[[4 * (i - 1) + 1]] <- plot.field_tile(
        qualitative_results$grid, loadings_true_HR[, i],
        boundary = boundary, limits = limits, breaks = breaks, LEGEND = FALSE
      ) + std_plot_settings_fields()

      ## Reconstructed components for quantile-based replicates
      for (j in 1:3) {
        plot_list[[4 * (i - 1) + j + 1]] <- plot.field_points(
          qualitative_results$nodes,
          loadings[[name_model]][[indexes[[i]][j]]][, i],
          boundary = boundary, size = 1.5
        ) + std_plot_settings_fields()

        plot_list_locs[[4 * (i - 1) + j + 1]] <- plot.field_points(
          qualitative_results$locations,
          loadings_locs[[name_model]][[indexes[[i]][j]]][, i],
          boundary = boundary, size = 1.5
        ) + std_plot_settings_fields()

        plot_list_HR_clean[[4 * (i - 1) + j + 1]] <- plot.field_tile(
          qualitative_results$grid,
          loadings_HR[[name_model]][[indexes[[i]][j]]][, i],
          boundary = boundary
        ) + std_plot_settings_fields()

        plot_list_HR[[4 * (i - 1) + j + 1]] <- plot.field_tile(
          qualitative_results$grid,
          loadings_HR[[name_model]][[indexes[[i]][j]]][, i],
          boundary = boundary, limits = limits, breaks = breaks
        ) + std_plot_settings_fields()
      }
    }

    ## Arrange labeled grids for each display mode
    plots[[m]] <- labled_plots_grid(arrangeGrob(grobs = plot_list, nrow = n_comp), label_model, labels_cols, labels_rows)
    plots_locs[[m]] <- labled_plots_grid(arrangeGrob(grobs = plot_list_locs, nrow = n_comp), label_model, labels_cols, labels_rows)
    plots_HR_clean[[m]] <- labled_plots_grid(arrangeGrob(grobs = plot_list_HR_clean, nrow = n_comp), label_model, labels_cols, labels_rows)
    plots_HR[[m]] <- labled_plots_grid(arrangeGrob(grobs = plot_list_HR, nrow = n_comp), label_model, labels_cols, labels_rows)
  }

  ## Display plots sequentially
  for (m in seq_along(model_names)) grid.arrange(plots_locs[[m]])
  for (m in seq_along(model_names)) grid.arrange(plots[[m]])
  for (m in seq_along(model_names)) grid.arrange(plots_HR_clean[[m]])
  for (m in seq_along(model_names)) grid.arrange(plots_HR[[m]])
}
