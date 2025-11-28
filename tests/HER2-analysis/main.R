rm(list = ls())

# [CONFIG] ----
LATEX_WIDTH_IN <- 345 / 72.27

## Load libraries ----
invisible(suppressMessages(sapply(c(
  # discretization
  "fdaPDE", "femR",
  # smoothing comparison
  "mgcv",
  # algebraic utils
  "pracma",
  # data manipulation
  "MASS", "tidyr", "dplyr",
  # visualization
  "ggplot2", "viridis", "stringr", "RColorBrewer", "grid", "gridExtra","ggrastr","cowplot", "tikzDevice",
  # meshing
  "RTriangle",
  # json
  "jsonlite",
  # sampling
  "sf", "sp", "raster"
), require, character.only = TRUE)))

## Load general utility functions
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/options.R")
source("src/utils/mesh_utils.R")
source("src/utils/domain_utils.R")
source("src/utils/plotting_utils.R")
source("src/utils/error_metrics.R")
source("src/utils/load_results_utils.R")
suppressMessages(sapply(list.files("src/data-generation", pattern = "\\.R$", full.names = TRUE), source))


## Load configuration file
path_this <- get_script_path()
source(paste0(path_this, "config.R"))

## Load test-specific functions
source("src/wrappers/fPCA.R")
source(paste0("tests/", test_suite, "/utils/generate_options.R"))

## Define and create work directories for the test suite
path_list <- create_paths(test_suite)

## Read arguments passed from the terminal
args <- commandArgs(trailingOnly = TRUE)

## Parse the arguments, if any
if (length(args) == 0) {
  ## Switch to interactive mode
  INTERACTIVE <- TRUE

  ## Load the option-generation function
  source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))

  ## Select the test you're interested in
  name_main_test <- name_main_test_default

  ## Update directories according to the new test
  path_list$queue <- paste0(path_list$queue, name_main_test, "/")
  path_list$logs <- paste0(path_list$logs, name_main_test, "/")
  mkdir(c(path_list$queue, path_list$logs))

  ## Generate all the options for that test
  generate_options(test_suite, name_main_test, path_list$queue)

  ## Read all the available options
  file_options_list <- sort(list.files(path_list$queue), decreasing = FALSE)
  file_options <- NULL
} else {
  ## Switch to non-interactive mode
  INTERACTIVE <- FALSE

  ## Set the requested configuration
  name_main_test <- args[1]
  file_options <- args[2]

  ## Update directories according to the new test
  path_list$queue <- paste0(path_list$queue, name_main_test, "/")
  path_list$logs <- paste0(path_list$logs, name_main_test, "/")
  mkdir(c(path_list$queue, path_list$logs))
}
## Select the test option
if (is.null(file_options)) {
  file_options_list
  file_options <- file_options_list[1] ## <====== INPUT HERE
}

## Load selected options
test_options <- fromJSON(paste0(path_list$queue, file_options))

## Update and create work directories for the test selected
path_list <- update_paths(path_list, name_main_test, test_options)

# [DATA] ----
# Data content
# - locations:      (data.frame)    n_locations x 2
# - counts:         (dgCMatrix)     n_genes x n_locations
load("tests/HER2-analysis/data/preprocessed_data.RData")
gene_names <- names(counts[,1])

## Visualise expert labeling
exp_labeling <- c(
  "3"= "In-situ cancer",
  "6"= "Invasive cancer",
  "4"= "Connective tissue",
  "1"= "Adipose tissue",
  "5"= "Immune infiltrates",
  "2"= "Breast glands",
  "7"= "Unlabeled"
)

exp_lab_cols <- c(
  "In-situ cancer"= "#FE7F29",
  "Invasive cancer"="#ED2023",
  "Connective tissue"="#3F48CC",
  "Adipose tissue"="#71E3DF",
  "Immune infiltrates"="#FFF204",
  "Breast glands"="#0ED145",
  "Unlabeled"= "black"
)

plot_data <- data.frame(
  x = locations$x,
  y = locations$y,
  numeric_label = as.factor(true_labels$true_label)
)

# Use the 'exp_labeling' vector to create a new column with the descriptive names
# We use 'as.character' to make sure we index the vector correctly
plot_data$tissue_type <- exp_labeling[as.character(plot_data$numeric_label)]

# Ensure the new column is a factor for a discrete color scale
# Ordering by 'exp_labeling' ensures the legend matches your intended order
plot_data$tissue_type <- factor(plot_data$tissue_type, levels = unname(exp_labeling))

exp_lab_plot <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(aes(color = tissue_type), size=5) +
  scale_color_manual(values=exp_lab_cols) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.width = unit(2, "lines"),
    text = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1)) +
  labs(color = NULL)



pdf(paste0(path_list$images,"exp_labels.pdf"),width = 1.8*LATEX_WIDTH_IN, height = 1*LATEX_WIDTH_IN)
print(exp_lab_plot)
dev.off()


## Function: plot_stacked_fields
# - Args:
#   * locations: data.frame with 2 columns (x, y) for all layers
#   * f_list: A NAMED list of numeric vectors. Names become labels (e.g., "Gene 1")
#   * boundary: optional data.frame (x, y) for a single boundary at the bottom
#   * vertical_shift_prop: Proportional y-shift (e.g., 0.5 = 50% of data's y-range)
#   * point_size: numeric, diameter of the points (now width of tile)
#   * limits: optional numeric range for color scaling (applied to all layers)
#   * colormap: string, Viridis palette option (default "D")
#   * y_squash_factor: numeric, controls the "perspective" (foreshortening)
#   * label_size: numeric, size of text labels
#   * label_nudge_x: numeric, horizontal nudge for labels
# - Desc:
#   Creates a stacked-layer plot by applying a vertical shift and y-axis
#   foreshortening (perspective) to each layer's coordinates.
plot_stacked_fields <- function(locations, f_list, boundary = NULL,
                                vertical_shift_prop = 0.5,
                                point_size = 1, # Renamed from 'size' to avoid confusion with tile parameters
                                limits = NULL, colormap = "D",
                                y_squash_factor = 0.4,
                                label_size = 5, label_nudge_x = -2) {
  
  master_data <- list()
  annotation_data <- list()
  
  # Calculate y-range of the *original* data for squashing and shift calculation
  y_range_orig <- max(locations$y) - min(locations$y)
  min_y_orig <- min(locations$y)
  
  # Calculate the actual vertical shift *between* the squashed layers
  # A prop of 0.5 means layers will overlap by 50% of the squashed height
  vertical_shift <- (y_range_orig * y_squash_factor) * vertical_shift_prop
  
  layer_names <- names(f_list)
  min_x <- min(locations$x)
  
  for (i in 1:length(layer_names)) {
    current_name <- layer_names[i]
    current_f <- f_list[[current_name]]
    
    # Calculate the total vertical offset for this layer
    current_y_offset <- (i - 1) * vertical_shift
    
    # Apply y-squashing and vertical offset
    y_plot_data <- (locations$y - min_y_orig) * y_squash_factor + current_y_offset
    x_plot_data <- locations$x
    
    temp_data <- data.frame(
      x = x_plot_data,
      y = y_plot_data,
      value = current_f,
      layer = current_name
    )
    master_data[[current_name]] <- temp_data
    
    # Calculate label position (squash median Y, then apply offset)
    y_label_plot <- (median(locations$y) - min_y_orig) * y_squash_factor + current_y_offset
    x_label_plot <- min_x
    
    annotation_data[[current_name]] <- data.frame(
      x = x_label_plot,
      y = y_label_plot,
      label = current_name
    )
  }
  
  # Combine all data into single data frames
  plot_data <- do.call(rbind, master_data)
  label_data <- do.call(rbind, annotation_data)
  
  # Set layer factor order (important for ggplot drawing order for correct overlaps)
  plot_data$layer <- factor(plot_data$layer, levels = layer_names)
  
  plot <- ggplot()
  
  # 1. Add the boundary (bottom-most layer)
  if (!is.null(boundary)) {
    squashed_boundary <- boundary
    squashed_boundary$y <- (squashed_boundary$y - min_y_orig) * y_squash_factor
    
    plot <- plot +
      geom_polygon(data = squashed_boundary, aes(x = x, y = y),
                   fill = "grey80", color = "black", linewidth = 1)
  }
  
  # 2. Add all the "points" (now tiles) from all layers
  #    The height of the tile is also squashed by y_squash_factor
  #    We want the height to be proportional to the original point size,
  #    but then squashed. width remains proportional to original point size.
  plot <- plot +
    geom_tile(data = plot_data, 
              aes(x = x, y = y, fill = value), 
              width = point_size * 0.75, # Tile width
              height = point_size * 0.75 * y_squash_factor, # Tile height is squashed
              color = "black", linewidth = 0.1) # Outline for each tile
  
  # 3. Add the text labels
  plot <- plot +
    geom_text(data = label_data, 
              aes(x = x, y = y, label = label),
              hjust = 1, # Right-align text
              nudge_x = label_nudge_x, # Nudge to the left
              size = label_size)
  
  # 4. Apply styling
  plot <- plot +
    coord_fixed() + # Essential for spatial data
    theme_void() +  # Remove all axes, gridlines, etc.
    guides(fill = "none") # Hide the legend (fill now used instead of color)
  
  # 5. Apply color scale (using scale_fill_viridis now)
  if (is.null(limits)) {
    plot <- plot + scale_fill_viridis(option = colormap, na.value = "transparent")
  } else {
    plot <- plot + scale_fill_viridis(option = colormap, limits = limits, na.value = "transparent")
  }
  
  return(plot)
}

gene_list <- lapply(1:5, function(i) counts[i,])
names(gene_list) <- names(counts[1:5,1])

plot_stacked_fields(locations, 
              gene_list, 
              boundary = ,
              vertical_shift_prop = 1.1, limits = NULL, 
              colormap = "magma",
              label_size = 10)

# [MESH] ----
# Create a planar straight line graph object
p <- pslg(P = locations)

# Create a regular mesh of the spatial domain
triangulation <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(triangulation$H)) triangulation$H <- matrix(numeric(0), ncol = 2)

# Create a regular mesh of the spatial domain
triangulation <- triangulate(triangulation, a = 0.1, q = 30, D = TRUE)
plot(triangulation)
mesh <- Mesh(triangulation)

domain <- list(femr_mesh = mesh, 
                boundary = NULL,
                fdapde_mesh = create.mesh.2D(nodes = mesh$nodes()))

# [MODELS FIT] ----
model_list <- list()
for (model_name in c("mv","sequential", "subspace")) {
    file_model <- paste0(path_list$results, "fitted_model_", model_name, ".RData")
    if (file.exists(file_model) && !FORCE_FIT) {
        cat("- Loading fitted model:", model_name, "... \n")
        load(file_model)
        model_list[[model_name]] <- get(paste0("model_", model_name))
    } else{
        ## Initialize empty model
        cat("- Fitting model:", model_name, "... ")    
        ## Fit the model
        model_list[[model_name]] <- 
            fit_model(model_name, 
                        domain,
                        data = list(X=counts, locations=locations), path_list, test_options)
        ## Save fitted model
      assign(paste("model_", model_name, sep = ""), model_list[[model_name]])
      save(
        list = paste("model_", model_name, sep = ""),
        file = file_model
      )
    }
}

# [VISUAL] ----
## Construct HR grid ----
## 1) Generate a rectangular grid
grid_step <- 1/6
seed_point <- SpatialPoints(data.frame(x = 11, y = -11))

bbox <- SpatialPoints(locations)@bbox

xmin <- bbox[1,1]
xmax <- bbox[1,2]
ymin <- bbox[2,1]
ymax <- bbox[2,2]

x <- seq(xmin, xmax, by = grid_step)
y <- seq(ymin, ymax, by = grid_step)
grid <- expand.grid(x = x, y = y)

## 2) Adapt the rectangular grid to the convex hull of the locations
locations_sf <- st_as_sf(as.data.frame(locations), coords = c("x", "y"), crs = 4326)
grid_sf <- st_as_sf(grid, coords = c("x", "y"), crs = 4326)

# Compute convex hull
convex_hull <- st_convex_hull(st_union(st_as_sf(as.data.frame(locations), coords = c("x", "y"), crs = 4326)))
# Keep only grid points within the convex hull
grid <- as.data.frame(st_coordinates(grid_sf[st_within(grid_sf, convex_hull, sparse = FALSE), ]))


##
gene_list <- lapply(1:3, function(i) evaluate_field(grid, model_list[["subspace"]]$results$loadings[, i],domain$fdapde_mesh))

names(gene_list) <- paste0("f",1:3)

names(grid) <- c("x","y")

plot_stacked_fields(grid, 
              rev(gene_list), 
              boundary = NULL,
              vertical_shift_prop = 1.1, limits = NULL, 
              colormap = "magma",
              label_size = 10)

## Principal functions ----
### At locations ----
n_comp <- 3

plot_list <- list()
for (i in seq_len(n_comp)) {
  row_plots <- list()
  for (name_model in c("mv","sequential", "subspace")) {
      row_plots[[name_model]] <- plot.field_points(
              locations, model_list[[name_model]]$results$loadings_locs[, i],
              boundary = NULL, 
              size = 1,
              colormap = "magma"
          ) + std_plot_settings_fields()
  }
  plot_list <- c(plot_list,row_plots)
}
final_grid <- labled_plots_grid(arrangeGrob(grobs = unlist(plot_list), nrow = n_comp), NULL,c("mv","seq", "sub"),paste0("fPC",1:n_comp))

pdf(paste0(path_list$images,"fPCs.pdf"),height = 5, width = 5)
grid.arrange(final_grid)
dev.off()


tikz(paste0(path_list$images,"fPCs.tex"), width = LATEX_WIDTH_IN, height=LATEX_WIDTH_IN)
grid.arrange(final_grid)
dev.off()



### At HR grid ----
plot_list <- list()
for (i in seq_len(n_comp)) {
    row_plots <- list()
    for (name_model in c("subspace")) {
        sign <- 1
        if(i %in% c(3)){ sign <- -1}
        row_plots[[name_model]] <- plot.field_points(
                grid, sign*evaluate_field(grid, model_list[[name_model]]$results$loadings[, i],domain$fdapde_mesh),
                boundary = as(convex_hull, "Spatial"), 
                colormap = "magma"
            ) + std_plot_settings_fields()
    }
    plot_list <- c(plot_list,row_plots)
}
final_grid <- arrangeGrob(grobs = plot_list, nrow = n_comp)
grid.arrange(final_grid)

png(paste0(path_list$images,"sub_fpcs.png"),height = 10, width = 5, units="in", res=150)
grid.arrange(final_grid)
dev.off()


## first PC
pdf(paste0(path_list$images,"fPC1.pdf"),height = 5, width = 5)
plot.field_points(
                grid, evaluate_field(grid, model_list[["subspace"]]$results$loadings[, 1],domain$fdapde_mesh),
                boundary = as(convex_hull, "Spatial"), 
                colormap = "magma"
            ) + std_plot_settings_fields()
dev.off()


## Mean ----
plot.field_points(
                locations, model_list$subspace$results$center_locs,
                boundary = NULL, size = 13
            ) + std_plot_settings_fields()



# [VARIANCE EXPLAINED] ----
test_options$model_options$n_comp <- 10

model_sub <- fit_model(
      "subspace",
      domain,
      data = list(X = counts, locations = locations),
      path_list = path_list,
      test_options = test_options
    )

var_explained_sub <- apply(model_sub$results$scores, 2, var)

pdf(paste0(path_list$images,"var_explained.pdf"),width = 10,height = 15)
barplot(var_explained_sub, main="Variance Explained",cex.main=2.3, cex.axis=1.5, cex.names=1.5, names.arg=paste0("fPC",1:10))
dev.off()



# [GCV CURVES] ----
lambda_grid <- test_options$regularization$lambda_grid

# Combine the scores into tidy data frames
gcv_sub_df <- as.data.frame(model_list$subspace$results$gcv_scores)
gcv_seq_df <- as.data.frame(model_list$sequential$results$gcv_scores)

# Add lambda (not log)
gcv_sub_df$lambda <- lambda_grid
gcv_seq_df$lambda <- lambda_grid

# Reshape to long format
gcv_sub_long <- gcv_sub_df %>%
  pivot_longer(
    cols = -lambda,
    names_to = "Component",
    values_to = "Subscore"
  )

gcv_seq_long <- gcv_seq_df %>%
  pivot_longer(
    cols = -lambda,
    names_to = "Component",
    values_to = "Seqscore"
  )

# Merge all three datasets
gcv_long <- gcv_sub_long %>%
  left_join(gcv_seq_long, by = c("lambda", "Component")) 

gcv_long <- gcv_long %>%
  mutate(
    Component = factor(
      Component,
      levels = unique(Component),
      labels = paste0("fPC", seq_along(unique(Component)))
    )
  )

# Find minima for red points
min_points_sub <- gcv_long %>%
  group_by(Component) %>%
  slice_min(Subscore, n = 1)

min_points_seq <- gcv_long %>%
  group_by(Component) %>%
  slice_min(Seqscore, n = 1)

# --- Compute helper data frames ---

# Subspace mean (for single "sub" line)
gcv_sub_mean <- gcv_long %>%
  group_by(lambda) %>%
  summarise(Subscore = sum(Subscore, na.rm = TRUE))

# Labels for sequential component lines (rightmost point)
seq_labels <- gcv_long %>%
  group_by(Component) %>%
  filter(lambda == max(lambda)) %>%
  slice_tail(n = 1)

library(ggrepel)

gcv_curves_plot <- ggplot() +
  # Multiple sequential lines (one per component)
  geom_line(
    data = gcv_long,
    aes(x = lambda, y = Seqscore, group = Component, color = "seq"),
    linewidth = 2,
    alpha = 0.5
  ) +
  # Label each sequential component
  geom_text_repel(
    data = seq_labels,
    aes(x = lambda, y = Seqscore, label = Component),
    size = 1,
    color = "coral3",
    direction = "y",
    hjust = 0,
    nudge_x = 0.15,
    segment.color = "grey70",
    segment.size = 0.3,
    box.padding = 0.2,
    min.segment.length = 0
  ) +
  # Add the averaged subspace line (bold aquamarine)
  geom_line(
    data = gcv_sub_mean,
    aes(x = lambda, y = Subscore, color = "sub"),
    linewidth = 1
  ) +
  geom_point(
    data = gcv_sub_mean %>%
      slice_min(Subscore, n = 1),
    aes(x = lambda, y = Subscore, color = "sub"),
    shape = 21, size = 3, fill = "#4DAF4A", stroke = 1
  ) +
  # --- ADDED THIS BLOCK ---
  # Add markers for the minimum of each sequential line
  geom_point(
    data = min_points_seq,
    aes(x = lambda, y = Seqscore, color = "seq"),
    shape = 21, size = 2, fill = "#E41A1C", stroke = 1
  ) +
  # -------------------------
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_log10(limits = c(min(gcv_long[,c("Subscore","Seqscore")]),max(gcv_sub_mean$Subscore))) +
  scale_color_manual(
    values = c(
      "sub" = "#4DAF4A",
      "seq" = "#E41A1C"),
    labels = c(
      "sub" = "subspace",
      "seq" = "sequential"
    )
  ) +
  labs(
    x = paste0("$\\","lambda$"),
    y = "GCV score",
    colour=NULL
  ) +
  std_plot_settings() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

pdf(paste0(path_list$images,"gcv_curves.pdf"), width = 12, height = 15)
gcv_curves_plot
dev.off()



# Combine GCV curves plot and barplot
library(gridExtra)
library(gridGraphics)


par(mar = c(3.5, 3.5, 2, 1),  # Shrink margins significantly
    mgp = c(2.0, 0.7, 0),     # Pull labels closer to the axis
    las = 1)
# Convert the base R barplot to a grid object
barplot_grob <- as_grob(
  ~barplot(
    var_explained_sub, 
    main = "Variance Explained", 
    names.arg = paste0("f", 1:10), 
    ylim = c(0, max(var_explained_sub) * 1.2)
  )
)

# Combine the GCV curves plot and the barplot
combined_plot <- grid.arrange(
  gcv_curves_plot, 
  barplot_grob, 
  ncol = 2
)


# Save the combined plot
tikz(paste0(path_list$images, "combined_gcv_variance.tex"), width = LATEX_WIDTH_IN, height = 0.5*LATEX_WIDTH_IN)
grid.draw(combined_plot)
dev.off()


# [SUBSPACE FIT] ----
model_fPCA <- model_list[["subspace"]]

sm_mean.grid <- evaluate_field(grid, model_fPCA$results$center,domain$fdapde_mesh)

## Smooth mean
pdf(paste0(path_list$images,"smooth_mean.pdf"),height = 5, width = 5)
plot.field_tile(
                grid, sm_mean.grid,
                boundary = as(convex_hull, "Spatial"), 
                colormap = "magma"
    ) + 
    std_plot_settings_fields() +
    ggtitle("Smooth mean")
dev.off()
 
## fPCs
fPCs.grid <- sapply(1:3, function(i) evaluate_field(grid, model_fPCA$results$loadings[, i],domain$fdapde_mesh))

tikz(paste0(path_list$images,"smooth_mean.tex"),height = LATEX_WIDTH_IN, width = LATEX_WIDTH_IN)
plot.field_tile(
                grid, sm_mean.grid,
                boundary = as(convex_hull, "Spatial"), 
                colormap = "magma"
    ) + 
    std_plot_settings_fields() +
    ggtitle("Smooth mean")
dev.off()

# [ERB22 analysis] ----
idx.HER2 <- 22

plot_HER2 <- plot.field_points(
                locations, counts[idx.HER2,],
                colormap = "magma", size=2
    ) + 
    std_plot_settings_fields() +
    ggtitle("ERBB2")


pdf(paste0(path_list$images,"erbb2.pdf"), width = 5, height = 5)
plot_HER2 
dev.off()

## Reconstruction
erb22_reconstruction <- 
  sweep(model_fPCA$results$scores[,1:3] %*% t(fPCs.grid), 2, 
        sm_mean.grid, FUN = "+")[idx.HER2,]
  
plot_HER2_reconstruction <- 
  plot.field_tile(
                grid, erb22_reconstruction,
                boundary = as(convex_hull, "Spatial"), 
                colormap = "magma"
    ) + 
    std_plot_settings_fields() +
    ggtitle("ERBB2 reconstruction")

pdf(paste0(path_list$images,"erbb2_rec.pdf"),height = 5, width = 5)
plot_HER2_reconstruction
dev.off()


# Combine ERBB2 plots using patchwork
library(patchwork)

combined_plot <- plot_HER2 + plot_HER2_reconstruction + plot_layout(ncol = 2)

tikz(paste0(path_list$images, "erbb2_raw_and_rec.tex"), width = LATEX_WIDTH_IN, height = 0.5*LATEX_WIDTH_IN)
print(combined_plot)
dev.off()

# scores1
gene_names[which.max(model_fPCA$results$scores[,1])]

#boxplot of scores1
pdf(paste0(path_list$images,"scores_1.pdf"),height = 5, width = 7)
boxplot(model_fPCA$results$scores[,1], 
  horizontal = T,
  frame.plot=F,
  main = "Score1", cex.main=2)
points(max(model_fPCA$results$scores[,1]), 1, 
  col = "red", pch = 19, cex = 2)
text(max(model_fPCA$results$scores[,1]), 1.1,
  labels = "ERB22", col = "red", pos = 2, cex=1.5)
dev.off()

# loadings
pdf("ERBB2_scores.pdf", width = 7, height = 5)
barplot(model_fPCA$results$scores[idx.HER2,1:3],
        names.arg = paste0("fPC", 1:3), 
        main = "ERB22 scores", cex.main=3,
        ylim=c(-10,20))
abline(h=0)
dev.off()



# scores1
gene_names[which.max(model_fPCA$results$scores[,1])]

#boxplot of scores1
tikz(paste0(path_list$images,"scores_1.tex"),height = 0.5*LATEX_WIDTH_IN, width = LATEX_WIDTH_IN)
par(mfrow=c(1,2))

boxplot(model_fPCA$results$scores[,1], 
  horizontal = T,
  frame.plot=F,
  main = "Score1")
points(max(model_fPCA$results$scores[,1]), 1, 
  col = "red", pch = 19)
text(max(model_fPCA$results$scores[,1]), 1.3,
  labels = "ERB22", col = "red", pos = 2)

barplot(model_fPCA$results$scores[idx.HER2,1:3],
        names.arg = paste0("fPC", 1:3), 
        main = "ERBB2",
        ylim=c(-10,20))
abline(h=0)
dev.off()

par(mfrow=c(1,1))

# [CROSS-VALIDATION] ----
cross_validate_model <- function(model_name, domain, counts, locations, path_list, test_options, K = 5, seed = 123, verbose = TRUE) {
  set.seed(seed)
  
  n <- nrow(counts)
  folds <- sample(rep(1:K, length.out = n))
  
  results_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    if (verbose) message(sprintf("Fold %d / %d ...", k, K))
    
    # --- Split data ---
    test_idx <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    X_train <- counts[train_idx, , drop = FALSE]    
    X_test <- counts[test_idx, , drop = FALSE]
    
    # --- Fit model on training data ---
    model_k <- fit_model(
      model_name,
      domain,
      data = list(X = X_train, locations = locations),
      path_list = path_list,
      test_options = test_options
    )
    # normalise loadings
    loadings_locs <- model_k$results$loadings_locs
    

    # --- Get denoised matrix ---
    X_hat_test <- sweep(X_test,2,model_k$results$center_locs)%*%loadings_locs%*%t(loadings_locs)
    X_hat_test <- sweep(X_hat_test,2,model_k$results$center_locs,FUN="+")
    
    # --- Compute error metrics ---
    mse <- mean((X_test - X_hat_test)^2, na.rm = TRUE)
    mae <- mean(abs(X_test - X_hat_test), na.rm = TRUE)
    
    results_list[[k]] <- list(
      fold = k,
      test_idx = test_idx,
      mse = mse,
      mae = mae
    )
  }
  
  # --- Aggregate metrics ---
  cv_summary <- data.frame(
    Fold = seq_len(K),
    MSE = sapply(results_list, `[[`, "mse"),
    MAE = sapply(results_list, `[[`, "mae")
  )
  
  cv_summary$MSE_mean <- mean(cv_summary$MSE)
  cv_summary$MAE_mean <- mean(cv_summary$MAE)
  
  if (verbose) {
    message("Cross-validation complete.")
    print(cv_summary)
  }
  
  return(list(
    folds = folds,
    results = results_list,
    summary = cv_summary
  ))
}

# cv_results <- list()
# for (model_name in c("smv","tpsPCA","sequential","subspace")) {
#     cv_results[[model_name]] <- cross_validate_model(model_name, domain, counts, locations, path_list, test_options, K = 5, seed = 1412, verbose = TRUE)
# }