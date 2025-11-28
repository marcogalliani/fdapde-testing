# = ========================================================================== =
# - Script: plotting_utils.R
# - Desc: Utility functions for plotting, including standard plot settings
#         and point visualization over 2D domains or meshes.
# = ========================================================================== =

TEXT_SZ <- 20

## Function: std_plot_settings
# - Args:
#   * NONE
# - Desc:
#   Returns a ggplot2 theme with standard visual settings for general-purpose
#   plots, using a clean black-and-white background and centered bold titles.
std_plot_settings <- function() {
  
  ## Create standard theme
  standard_plot_settings <- theme_bw() +
    theme(
      text = element_text(size = TEXT_SZ),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = TEXT_SZ,
        hjust = 0.5,
        vjust = 1
      ),
      legend.position = "top",
      legend.text = element_text(size = TEXT_SZ),
      axis.title.x = element_text(size = TEXT_SZ),
      axis.title.y = element_text(size = TEXT_SZ),
      axis.text.x = element_text(size = TEXT_SZ),
      axis.text.y = element_text(size = TEXT_SZ),
      legend.text.position = "top"
    )
}


## Function: std_plot_settings_fields
# - Args:
#   * NONE
# - Desc:
#   Returns a ggplot2 theme optimized for visualizing spatial fields.
#   Removes axes, ticks, and grid lines, keeping only essential elements
#   such as color legend and plot title.
std_plot_settings_fields <- function() {
  
  ## Create theme for field visualization
  standard_plot_settings_fields <- theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(
        color = "black",
        face = "bold",
        size = 14,
        hjust = 0.5,
        vjust = 1
      ),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(angle = 45, hjust = 1)
    )
}


## Function: plot.points
# - Args:
#   * locations: matrix or data.frame with 2 columns (x, y) for point coordinates
#   * boundary: optional SpatialPolygons or data.frame defining a boundary
#   * group: optional vector assigning each point to a group
#   * group_name: string, name of the legend for groups (default: "Groups")
#   * group_colors: optional named vector of colors for each group
#   * group_labels: optional vector of labels for each group
#   * size: numeric, point size
#   * LEGEND: logical, whether to display the legend (default: FALSE)
# - Desc:
#   Plots 2D points (locations) with optional grouping and boundary outline.
#   Automatically assigns colors and labels when not provided, and ensures
#   fixed aspect ratio for spatial accuracy.
plot.points <- function(locations, boundary = NULL, group = NULL, 
                        group_name = "Groups", group_colors = NULL, group_labels = NULL,
                        size = 1, LEGEND = FALSE) {
  
  ## Assemble data
  if (is.null(group)) {
    group <- rep(0, nrow(locations))
  }
  data <- data.frame(locations, group)
  colnames(data) <- c("x", "y", "Group")
  
  ## Grouping
  group_levels <- unique(data$Group)
  if (is.null(group_labels)) {
    group_labels <- group_levels
  }
  
  ## Define association between group labels and colors
  if (is.null(group_colors)) {
    if (length(group_labels) > 1) {
      group_colors <- rainbow(length(group_labels))
    } else {
      group_colors <- c("#000000")
    }
  }
  names(group_colors) <- group_labels
  
  ## Refactor categorical variables
  data <- data %>%
    mutate(Group = factor(Group, levels = group_levels, labels = group_labels))
  
  ## Build plot
  plot <- ggplot() +
    geom_point(data = data, aes(x = x, y = y, color = Group), size = size) +
    coord_fixed() +
    scale_color_manual(
      name = group_name,
      values = group_colors,
      labels = group_labels
    )
  
  ## Add boundary if provided
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(
        data = fortify(boundary),
        aes(x = long, y = lat),
        fill = "transparent",
        color = "black",
        linewidth = 1
      )
  }
  
  ## Add or remove legend
  if (!LEGEND) {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}

## Function: plot.field_points
# - Args:
#   * locations: matrix or data.frame with 2 columns (x, y) for coordinates
#   * f: numeric vector of values to be plotted
#   * boundary: optional SpatialPolygons or data.frame defining a boundary
#   * size: numeric, point size
#   * limits: optional numeric range for color scaling
#   * colormap: string, Viridis palette option (default "D")
#   * discrete: logical, whether to use discrete colormap
#   * LEGEND: logical, whether to display the legend
# - Desc:
#   Plots a spatial field defined at scattered locations using colored points.
#   Supports both continuous and discrete color scales, with optional boundary overlay.
plot.field_points <- function(locations, f, boundary = NULL,
                              size = 1, limits = NULL, colormap = "D",
                              discrete = FALSE, LEGEND = FALSE) {
  
  ## Handle null input
  if (is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  ## Assemble data
  data <- data.frame(locations, value = f)
  colnames(data) <- c("x", "y", "value")
  
  ## Build base plot
  plot <- ggplot() +
    geom_point(data = data, aes(x = x, y = y, color = value), size = size) +
    coord_fixed()
  
  ## Apply color scale
  if (!discrete) {
    if (is.null(limits)) {
      plot <- plot + scale_color_viridis(option = colormap)
    } else {
      plot <- plot + scale_color_viridis(option = colormap, limits = limits)
    }
  } else {
    plot <- plot + scale_color_viridis_d(option = colormap)
  }
  
  ## Add boundary if provided
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat),
                   fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## Add or remove legend
  if (!LEGEND) {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}


## Function: plot.field_tile
# - Args:
#   * nodes: matrix or data.frame with 2 columns (x, y) for coordinates
#   * f: numeric vector of field values
#   * boundary: optional boundary polygon
#   * limits: optional numeric range for color scaling
#   * breaks: optional numeric vector for contour levels
#   * colormap: string, Viridis palette option
#   * discrete: logical, whether to use discrete colormap
#   * ISOLINES: logical, whether to draw contour isolines
#   * LEGEND: logical, whether to display the legend
# - Desc:
#   Plots a spatial field on a regular grid using colored tiles.
#   Supports contour overlays, color limits, and boundary visualization.
plot.field_tile <- function(nodes, f, boundary = NULL,
                            limits = NULL, breaks = NULL, colormap = "D",
                            discrete = FALSE, ISOLINES = FALSE, LEGEND = FALSE) {
  
  ## Handle null input
  if (is.null(f)) {
    return(ggplot() + theme_void())
  }
  
  ## Assemble data
  data <- data.frame(nodes, value = f)
  colnames(data) <- c("x", "y", "value")
  data <- na.omit(data)
  
  ## Build base plot
  plot <- ggplot() +
    geom_tile(data = data, aes(x = x, y = y, fill = value)) +
    coord_fixed()
  
  ## Add contour isolines if requested
  if (!is.null(breaks) || ISOLINES) {
    color <- "black"
    limits_real <- range(data$value)
    if (is.null(breaks)) {
      breaks <- seq(limits_real[1], limits_real[2], length = 10)
    }
    breaks_initial <- breaks
    h <- breaks[2] - breaks[1]
    if (limits_real[1] < min(breaks)) {
      breaks <- c(sort(seq(min(breaks), limits_real[1] - h, by = -h)[-1]), breaks)
    }
    if (limits_real[2] > max(breaks)) {
      breaks <- c(breaks, seq(max(breaks), limits_real[2] + h, by = h)[-1])
    }
    if (length(breaks) > 2 * length(breaks_initial)) {
      breaks <- breaks_initial
      color <- "red"
    }
    plot <- plot +
      geom_contour(data = data, aes(x = x, y = y, z = value),
                   color = color, breaks = breaks)
  }
  
  ## Apply color scale
  if (!discrete) {
    if (!is.null(breaks)) {
      h <- breaks[2] - breaks[1]
      limits <- limits + c(-h, h)
    }
    if (is.null(limits)) {
      plot <- plot + scale_fill_viridis(option = colormap)
    } else {
      plot <- plot + scale_fill_viridis(option = colormap, limits = limits)
    }
  } else {
    plot <- plot + scale_fill_viridis_d(option = colormap)
  }
  
  ## Add boundary if provided
  if (!is.null(boundary)) {
    plot <- plot +
      geom_polygon(data = fortify(boundary), aes(x = long, y = lat),
                   fill = "transparent", color = "black", linewidth = 1)
  }
  
  ## Add or remove legend
  if (!LEGEND) {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}


## Function: plot.grouped_boxplots
# - Args:
#   * data: data.frame with columns Group | Model1 | ... | ModelN
#   * group_name: string, label for the x-axis
#   * subgroup_name: string, label for legend entries
#   * subgroup_colors: optional named vector of colors for subgroups
#   * values_name: string, label for y-axis
#   * limits: optional y-axis range
#   * DIVIDERS: logical, whether to draw vertical separators between groups
#   * LEGEND: logical, whether to display the legend
# - Desc:
#   Creates grouped boxplots comparing models (columns) within multiple groups.
plot.grouped_boxplots <- function(data,
                                  group_name = "Components", group_labels = NULL,
                                  subgroup_name = "Models", subgroup_labels = NULL, subgroup_colors = NULL,
                                  values_name = "Score", limits = NULL,
                                  DIVIDERS = TRUE, LEGEND = TRUE) {
  
  ## Data integrity check
  if (!("Group" %in% names(data))) stop("The dataframe must contain a column named 'Group'")
  
  ## Reshape data
  data <- data %>%
    pivot_longer(cols = -Group, names_to = "SubGroup", values_to = "Score")
  
  ## Extract labels
  groups_levels <- unique(data$Group)
  subgroup_levels <- unique(data$SubGroup)
  if (is.null(group_labels)) group_labels <- groups_levels
  if (is.null(subgroup_labels)) subgroup_labels <- subgroup_levels
  
  ## Assign colors
  if (is.null(subgroup_colors)) subgroup_colors <- rainbow(length(subgroup_labels))
  names(subgroup_colors) <- subgroup_labels
  
  ## Refactor categorical variables
  data <- data %>%
    mutate(Group = factor(Group, levels = groups_levels, labels = group_labels)) %>%
    mutate(SubGroup = factor(SubGroup, levels = subgroup_levels, labels = subgroup_labels))
  
  ## Frame color (Darker)
  outline_colors <- adjustcolor(subgroup_colors, red.f = 0.6, green.f = 0.6, blue.f = 0.6)
  names(outline_colors) <- subgroup_labels
  
  ## Medians for the bars 
  data_medians <- data %>%
    group_by(Group, SubGroup) %>%
    summarise(Median = median(Score, na.rm = TRUE), .groups = "drop")
  
  box_width <- 0.75
  pos_dodge <- position_dodge(width = box_width)

  ## Build plot
  plot <- ggplot(data, aes(x = Group, y = Score, fill = SubGroup)) +
    ##Underlying Bars 
    geom_col(data = data_medians, 
             aes(y = Median, group = SubGroup),
             position = pos_dodge, width = box_width, 
             alpha = 0.3, color = NA) +
    # Update Boxplots to use the explicit position and width
    stat_boxplot(geom = "errorbar", position = pos_dodge, width = box_width, aes(color = SubGroup)) + 
    geom_boxplot(na.rm = TRUE, position = pos_dodge, width = box_width, aes(color = SubGroup)) +
    labs(x = group_name, y = values_name) +
    scale_fill_manual(name = NULL, values = subgroup_colors) +
    scale_color_manual(name = NULL, values = outline_colors)
  
  ## Add limits
  if (!is.null(limits)) {
    plot <- plot + scale_y_continuous(limits = limits, expand = expansion(mult = c(0, 0.05)))
  }
  
  ## Add dividers
  if (DIVIDERS && length(group_labels) > 1) {
    plot <- plot +
      geom_vline(xintercept = seq(1.5, length(unique(group_labels)) - 0.5, 1),
                 lwd = 0.2, colour = "grey")
  }
  
  ## Legend
  if (!LEGEND) {
    plot <- plot + guides(fill = "none", color = "none")
  } else {
    plot <- plot + 
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(nrow = 1, label.position = "right"), 
             color = guide_legend(nrow = 1, label.position = "right"))
  }
  
  return(plot)
}


## Function: plot.multiple_lines
# - Args:
#   * data: data.frame with columns x | y1 | ... | yN
#   * x_name: string, label for x-axis
#   * subgroup_name: string, label for legend
#   * subgroup_labels, subgroup_colors: optional aesthetics
#   * limits: optional y-axis range
#   * LOGX, LOGY, LOGLOG: logical flags for logarithmic scales
#   * NORMALIZED: logical, normalize curves by their first value
#   * LEGEND: logical, whether to display the legend
# - Desc:
#   Plots multiple lines (e.g., performance curves) on the same axes with
#   optional logarithmic and normalized scaling.
plot.multiple_lines <- function(data,
                                x_name = "Components",
                                x_breaks = NULL,
                                subgroup_name = "Models", subgroup_labels = NULL, subgroup_colors = NULL,
                                values_name = "Score", limits = NULL,
                                LOGX = FALSE, LOGY = FALSE, LOGLOG = FALSE,
                                NORMALIZED = FALSE, LEGEND = TRUE) {
  
  ## Rename first column
  columns_names <- colnames(data)
  columns_names[1] <- "x"
  colnames(data) <- columns_names
  
  ## Normalize if requested
  if (NORMALIZED) {
    for (name in columns_names[-1])
      data[, name] <- data[, name] / min(data[, name])
  }
  
  ## Log-log consistency
  if (LOGLOG) {
    LOGX <- TRUE
    LOGY <- TRUE
  }
  
  ## Auto-generate breaks if requested
  if (is.logical(x_breaks) && x_breaks) x_breaks <- unique(data$x)
  
  ## Reshape data
  data <- data %>%
    pivot_longer(cols = -x, names_to = "SubGroup", values_to = "Score")
  
  ## Extract labels
  subgroup_levels <- unique(data$SubGroup)
  if (is.null(subgroup_labels)) subgroup_labels <- subgroup_levels
  
  ## Assign colors
  if (is.null(subgroup_colors)) subgroup_colors <- rainbow(length(subgroup_labels))
  names(subgroup_colors) <- subgroup_labels
  
  ## Refactor categories
  data <- data %>%
    mutate(SubGroup = factor(SubGroup, levels = subgroup_levels, labels = subgroup_labels))
  
  ## Build plot
  plot <- ggplot(data, aes(x = x, y = Score, color = SubGroup)) +
    geom_line(linewidth = 1) +
    labs(x = x_name, y = values_name) +
    scale_color_manual(name = subgroup_name, values = subgroup_colors)
  
  ## Logarithmic reference lines for normalized log-log plots
  if (LOGLOG && NORMALIZED) {
    x <- seq(min(data$x), max(data$x), length = 10)
    plot <- plot +
      geom_line(data = data.frame(x = x, y = x / x[1]), aes(x = x, y = y),
                linetype = "dashed", color = "grey", linewidth = 0.3) +
      geom_line(data = data.frame(x = x, y = (x / x[1])^2), aes(x = x, y = y),
                linetype = "dashed", color = "grey", linewidth = 0.3) +
      geom_line(data = data.frame(x = x, y = (x / x[1])^3), aes(x = x, y = y),
                linetype = "dashed", color = "grey", linewidth = 0.3)
  }
  
  ## Axis scaling
  if (LOGX) {
    plot <- plot + scale_x_log10(breaks = x_breaks)
  } else {
    plot <- plot + scale_x_continuous(breaks = x_breaks)
  }
  if (LOGY) {
    plot <- plot + scale_y_log10(limits = limits)
  } else if (!is.null(limits)) {
    plot <- plot + scale_y_continuous(limits = limits)
  }
  
  ## Add or remove legend
  if (!LEGEND) {
    plot <- plot + guides(color = "none")
  }
  
  return(plot)
}


## Function: labled_plots_grid
# - Args:
#   * plot: ggplot or grid object
#   * title: optional string title
#   * labels_cols, labels_rows: optional lists of strings labeling rows/cols
#   * height, width: numeric dimensions
# - Desc:
#   Adds row and column labels (and optional title) to a grid of ggplots.
labled_plots_grid <- function(plot, title = NULL, labels_cols = NULL,
                              labels_rows = NULL, height = 15, width = 18) {
  
  ## Compute grid dimensions
  n_row <- max(plot$layout$t)
  n_col <- max(plot$layout$l)
  
  ## Add column labels
  if (!is.null(labels_cols)) {
    labels_grobs_cols <- lapply(labels_cols, function(lab)
      textGrob(lab, gp = gpar(fontsize = TEXT_SZ, fontface = "bold")))
    labels_grobs_cols <- arrangeGrob(grobs = labels_grobs_cols, nrow = 1)
  }
  
  ## Add row labels
  add <- 0
  if (!is.null(labels_rows)) {
    labels_grobs_rows <- list()
    if (!is.null(labels_cols)) {
      add <- 1
      labels_grobs_rows[[1]] <- textGrob(" ", gp = gpar(fontsize = TEXT_SZ, fontface = "bold"))
    }
    for (row in 1:length(labels_rows) + add) {
      label_row <- labels_rows[[row - add]]
      labels_grobs_rows[[row]] <- textGrob(label_row, gp = gpar(fontsize = TEXT_SZ, fontface = "bold"), rot = 90)
    }
    labels_grobs_rows <- arrangeGrob(grobs = labels_grobs_rows, ncol = 1, heights = c(1, rep(height, n_row)))
  }
  
  ## Add title
  if (!is.null(title)) {
    title_grob <- textGrob(title, gp = gpar(fontsize = TEXT_SZ, fontface = "bold"))
  }
  
  ## Combine all components
  if (!is.null(labels_cols)) plot <- arrangeGrob(labels_grobs_cols, plot, heights = c(1, height * n_row))
  if (!is.null(labels_rows)) plot <- arrangeGrob(labels_grobs_rows, plot, widths = c(1, width * n_col))
  if (!is.null(title)) plot <- arrangeGrob(title_grob, plot, heights = c(1, add + height * n_row))
  
  return(plot)
}


## Function: plot.grouped_violins
# - Args:
#   * data: data.frame with columns Group | Model1 | ... | ModelN
#   * group_name, subgroup_name: strings, axis and legend labels
#   * subgroup_colors: optional colors for subgroups
#   * values_name: string, y-axis label
#   * limits: optional numeric y-axis range
#   * DIVIDERS, LEGEND: logical flags
#   * show_boxplot: logical, overlay boxplots inside violins
# - Desc:
#   Draws grouped violin plots comparing model distributions across groups,
#   with optional embedded boxplots and group separators.
plot.grouped_violins <- function(data,
                                 group_name = "Components", group_labels = NULL,
                                 subgroup_name = "Models", subgroup_labels = NULL, subgroup_colors = NULL,
                                 values_name = "Score", limits = NULL,
                                 DIVIDERS = TRUE, LEGEND = TRUE, show_boxplot = TRUE) {
  
  ## Data integrity check
  if (!("Group" %in% names(data))) stop("The dataframe must contain a column named 'Group'")
  
  ## Reshape data
  data <- data %>%
    pivot_longer(cols = -Group, names_to = "SubGroup", values_to = "Score")
  
  ## Extract labels
  groups_levels <- unique(data$Group)
  subgroup_levels <- unique(data$SubGroup)
  if (is.null(group_labels)) group_labels <- groups_levels
  if (is.null(subgroup_labels)) subgroup_labels <- subgroup_levels
  
  ## Assign colors
  if (is.null(subgroup_colors)) subgroup_colors <- rainbow(length(subgroup_labels))
  names(subgroup_colors) <- subgroup_labels
  
  ## Refactor variables
  data <- data %>%
    mutate(Group = factor(Group, levels = groups_levels, labels = group_labels)) %>%
    mutate(SubGroup = factor(SubGroup, levels = subgroup_levels, labels = subgroup_labels))
  
  ## Build violin plot
  plot <- ggplot(data, aes(x = Group, y = Score, fill = SubGroup)) +
    geom_violin(trim = FALSE, width = 1.5, position = position_dodge(width = 0.8),
                na.rm = TRUE, alpha = 0.8) +
    labs(x = group_name, y = values_name) +
    scale_fill_manual(name = subgroup_name, values = subgroup_colors)
  
  ## Optionally overlay boxplots
  if (show_boxplot) {
    plot <- plot +
      geom_boxplot(width = 0.15, position = position_dodge(width = 0.8),
                   outlier.shape = NA, alpha = 0.6)
  }
  
  ## Apply y-limits
  if (!is.null(limits)) {
    plot <- plot + scale_y_continuous(
      limits = limits
    )
  }
  
  ## Add groups divider if required
  if (DIVIDERS) {
    if (length(group_labels) > 1) {
      plot <- plot +
        geom_vline(
          xintercept = seq(1.5, length(unique(group_labels)) - 0.5, 1),
          lwd = 0.2, colour = "grey"
        )
    }
  }
  
  ## Legend control
  if (!LEGEND) {
    plot <- plot + guides(fill = "none")
  }
  
  return(plot)
}



# = ========================================================================== =
# - Script: aggregated_plots.R
# - Author: Pietro Donelli
# - Date: 2025-10-24
# - Desc: Aggregates quantitative results across varying options and produces
#         comparative plots (boxplots, lines, log-x, log-log, normalized).
#         Supports multi-dimensional option grids with labeled faceting.
# = ========================================================================== =


## Function: plot.aggregated_data
# - Args:
#   * loaded_results: list carrying model metadata and (optionally) varying_options
#       - $model_names, $model_labels, $model_colors
#       - $varying_options: character vector of option names used for grouping
#   * data_plot_orig: data.frame with columns:
#       - "Group" (aggregation variable = first varying option)
#       - one column per model in loaded_results$model_names
#       - one column per varying option (added upstream)
#   * title_prefix: string prefix for figure titles
#   * values_names: unused here (kept for API parity; pass NULL)
#   * order: optional integer vector to reorder varying_options (default: 1:k)
#   * limits: optional y-axis limits c(ymin, ymax) for value scales
#   * plots_catalog: list of toggles {boxplots, lines, logx, loglog, normalized}
# - Desc:
#   Builds per-group panels of aggregated results over the first varying option,
#   conditioning on all remaining varying options. For each combination, it can
#   render boxplots and/or line plots (linear, log-x, log-log, normalized).
plot.aggregated_data <- function(loaded_results, data_plot_orig, title_prefix, values_names, order = NULL, limits = NULL, plots_catalog = NULL) {
  
  ## Get model-related properties
  model_names  <- loaded_results$model_names
  model_labels <- loaded_results$model_labels
  model_colors <- loaded_results$model_colors
  
  ## Get varying options and apply ordering
  varying_options <- loaded_results$varying_options
  if (is.null(order)) order <- seq_along(varying_options)
  varying_options <- varying_options[order]
  name_varying_options <- varying_options
  
  ## Build options grid (unique sorted values for each varying option)
  options_grid <- list()
  for (name_ao in varying_options) {
    options_grid[[name_ao]] <- unique(data_plot_orig[, name_ao])
    options_grid[[name_ao]] <- sort(options_grid[[name_ao]])
  }
  
  ## Plots catalog defaults
  if (is.null(plots_catalog)) {
    plots_catalog <- list(
      boxplots = TRUE,
      lines = FALSE,
      logx = FALSE,
      loglog = FALSE,
      normalized = FALSE
    )
  }
  
  ## Detect groups (values of the first varying option will become x-axis)
  groups <- sort(unique(data_plot_orig$Group))
  
  ## If there is more than one group, plot them sequentially
  for (group in groups) {  # group <- groups[1]
    
    ## Select the specific group
    data_plot <- data_plot_orig[data_plot_orig$Group == group, ]
    
    ## Generate title
    title <- paste(
      title_prefix,
      name_varying_options[1],
      ifelse(length(groups) > 1, paste0("- ", group), "")
    )
    
    ## Rooms for plots
    boxplot_list <- list()
    plot_list <- list()
    plot_logx_list <- list()
    plot_loglog_list <- list()
    plot_loglog_normalized_list <- list()
    
    name_aggregation_option <- varying_options[1]
    group_name <- name_varying_options[1]
    
    options_grid_selected <- options_grid
    options_grid_selected[[name_aggregation_option]] <- NULL
    names_options_selected  <- names(options_grid_selected)
    labels_options_selected <- name_varying_options[-1]
    
    ## Handle 1D case (only one varying option)
    if (length(options_grid_selected) == 0) {
      combinations_options <- data.frame(dummy = 1)
      labels_rows <- ""
      labels_cols <- ""
    } else {
      mg <- do.call(expand.grid, options_grid_selected)
      combinations_options <- do.call(data.frame, lapply(mg, as.vector))
      colnames(combinations_options) <- names_options_selected
      
      labels_rows <- if (length(labels_options_selected) >= 1) paste(labels_options_selected[1], "=", options_grid_selected[[1]]) else ""
      labels_cols <- if (length(labels_options_selected) >= 2) paste(labels_options_selected[2], "=", options_grid_selected[[2]]) else ""
    }
    
    for (j in 1:nrow(combinations_options)) {  # j <- 1
      
      ## Data preparation
      if (length(names_options_selected) == 0) {
        ## 1D case: use all data
        data_plot_trimmed <- data_plot[, c(name_aggregation_option, model_names)]
      } else {
        ## ND case: filter by combination
        condition <- rep(TRUE, nrow(data_plot))
        for (k in seq_along(names_options_selected)) {
          condition <- condition & (data_plot[[names_options_selected[k]]] == combinations_options[j, k])
        }
        data_plot_trimmed <- data_plot[condition, c(name_aggregation_option, model_names)]
      }
      colnames(data_plot_trimmed)[1] <- "Group"
      
      ## Skip if data is empty
      if (nrow(data_plot_trimmed) == 0) next
      
      ## Remove models with all-NaN results
      valid_models <- model_names[!apply(data_plot_trimmed[, model_names], 2, function(x) all(is.nan(x)))]
      if (length(valid_models) == 0) next
      
      ## Aggregate with median per x (Group)
      data_plot_aggregated <- aggregate(. ~ Group, data = data_plot_trimmed[, c("Group", valid_models)], FUN = median)
      
      ## Boxplots
      if (isTRUE(plots_catalog$boxplots)) {
        boxplot_list[[j]] <- plot.grouped_boxplots(
          data_plot_trimmed[, c("Group", valid_models)],
          values_name = NULL,
          group_name = group_name,
          subgroup_name = "Approaches",
          subgroup_labels = model_labels[match(valid_models, model_names)],
          subgroup_colors = model_colors[match(valid_models, model_names)],
          limits = limits,
          LEGEND = TRUE
        ) + std_plot_settings()
      }
      
      ## Lines (linear x)
      if (isTRUE(plots_catalog$lines)) {
        plot_list[[j]] <- plot.multiple_lines(
          data_plot_aggregated[, c("Group", valid_models)],
          values_name = NULL,
          x_name = group_name,
          x_breaks = TRUE,
          subgroup_name = "Approaches",
          subgroup_labels = model_labels[match(valid_models, model_names)],
          subgroup_colors = model_colors[match(valid_models, model_names)],
          LEGEND = FALSE,
          limits  = limits,
          NORMALIZED = FALSE,
          LOGX = FALSE
        ) + std_plot_settings()
      }
      
      ## Lines (log-x)
      if (isTRUE(plots_catalog$logx)) {
        plot_logx_list[[j]] <- plot.multiple_lines(
          data_plot_aggregated[, c("Group", valid_models)],
          values_name = NULL,
          x_name = group_name,
          x_breaks = TRUE,
          subgroup_name = "Approaches",
          subgroup_labels = model_labels[match(valid_models, model_names)],
          subgroup_colors = model_colors[match(valid_models, model_names)],
          LEGEND = FALSE,
          limits = limits,
          NORMALIZED = FALSE,
          LOGX = TRUE
        ) + std_plot_settings()
      }
      
      ## Lines (log-log)
      if (isTRUE(plots_catalog$loglog)) {
        plot_loglog_list[[j]] <- plot.multiple_lines(
          data_plot_aggregated[, c("Group", valid_models)],
          values_name = NULL,
          x_name = group_name,
          x_breaks = TRUE,
          subgroup_name = "Approaches",
          subgroup_labels = model_labels[match(valid_models, model_names)],
          subgroup_colors = model_colors[match(valid_models, model_names)],
          LEGEND = FALSE,
          limits = NULL,
          NORMALIZED = FALSE,
          LOGLOG = TRUE
        ) + std_plot_settings()
      }
      
      ## Lines (log-log, normalized)
      if (isTRUE(plots_catalog$normalized)) {
        plot_loglog_normalized_list[[j]] <- plot.multiple_lines(
          data_plot_aggregated[, c("Group", valid_models)],
          values_name = NULL,
          x_name = group_name,
          x_breaks = TRUE,
          subgroup_name = "Approaches",
          subgroup_labels = model_labels[match(valid_models, model_names)],
          subgroup_colors = model_colors[match(valid_models, model_names)],
          LEGEND = FALSE,
          limits = NULL,
          NORMALIZED = TRUE,
          LOGLOG = TRUE
        ) + std_plot_settings()
      }
    }
    
    ## Handle layout safely if 1D: use 1 column
    ncols <- if (length(labels_cols) == 0) 1 else length(labels_cols)
    
    if (isTRUE(plots_catalog$boxplots)) {
      boxplot <- arrangeGrob(grobs = boxplot_list, ncol = ncols)
      boxplot <- labled_plots_grid(boxplot, title, labels_cols, labels_rows)
      grid.arrange(boxplot)
    }
    
    if (isTRUE(plots_catalog$lines)) {
      plot <- arrangeGrob(grobs = plot_list, ncol = ncols)
      plot <- labled_plots_grid(plot, title, labels_cols, labels_rows)
      grid.arrange(plot)
    }
    
    if (isTRUE(plots_catalog$logx)) {
      plot_logx <- arrangeGrob(grobs = plot_logx_list, ncol = ncols)
      plot_logx <- labled_plots_grid(plot_logx, title, labels_cols, labels_rows)
      grid.arrange(plot_logx)
    }
    
    if (isTRUE(plots_catalog$loglog)) {
      plot_loglog <- arrangeGrob(grobs = plot_loglog_list, ncol = ncols)
      plot_loglog <- labled_plots_grid(plot_loglog, title, labels_cols, labels_rows)
      grid.arrange(plot_loglog)
    }
    
    if (isTRUE(plots_catalog$normalized)) {
      plot_loglog_normalized <- arrangeGrob(grobs = plot_loglog_normalized_list, ncol = ncols)
      plot_loglog_normalized <- labled_plots_grid(plot_loglog_normalized, title, labels_cols, labels_rows)
      grid.arrange(plot_loglog_normalized)
    }
  }
}