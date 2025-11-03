rm(list = ls())

# [CONFIG] ----

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
  "ggplot2", "viridis", "stringr", "RColorBrewer", "grid", "gridExtra",
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

# [MESH] ----
# Create a planar straight line graph object
p <- pslg(P = locations)

# Create a regular mesh of the spatial domain
triangulation <- triangulate(p, Y = FALSE, D = TRUE)
if (is.null(triangulation$H)) triangulation$H <- matrix(numeric(0), ncol = 2)

# Create a regular mesh of the spatial domain
triangulation <- triangulate(triangulation, a = 0.1, q = 20, D = TRUE)
plot(triangulation)
mesh <- Mesh(triangulation)

domain <- list(femr_mesh = mesh, 
                boundary = NULL,
                fdapde_mesh = create.mesh.2D(nodes = mesh$nodes()))

# [MODELS FIT] ----
model_list <- list()
for (model_name in test_options$model_names) {
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


## Principal functions ----
### At HR grid ----
n_comp <- test_options$model_options$n_comp

plot_list <- list()
for (i in seq_len(n_comp)) {
  row_plots <- list()
  for (name_model in c("sequential", "subspace")) {
      row_plots[[name_model]] <- plot.field_points(
              locations, model_list[[name_model]]$results$loadings_locs[, i],
              boundary = NULL, size = 2.5
          ) + std_plot_settings_fields()
  }
  plot_list[[i]] <- arrangeGrob(grobs = row_plots, ncol = length(test_options$model_names))
}
final_grid <- arrangeGrob(grobs = plot_list, nrow = n_comp)
grid.arrange(final_grid)

### At locations ----
plot_list <- list()
for (i in seq_len(n_comp)) {
    row_plots <- list()
    for (name_model in test_options$model_names) {
        row_plots[[name_model]] <- plot.field_points(
                locations, model_list[[name_model]]$results$loadings_locs[, i],
                boundary = NULL, size = 2.5
            ) + std_plot_settings_fields()
    }
    plot_list[[i]] <- arrangeGrob(grobs = row_plots, ncol = length(test_options$model_names))
}
final_grid <- arrangeGrob(grobs = plot_list, nrow = n_comp)
grid.arrange(final_grid)


## Mean ----
plot.field_points(
                locations, model_list$subspace$results$center_locs,
                boundary = NULL, size = 13
            ) + std_plot_settings_fields()



# [VARIANCE EXPLAINED] ----
var_explained_sub <- apply(model_list$subspace$results$scores, 2, var)
barplot(var_explained_sub)

var_explained_seq <- apply(model_list$sequential$results$scores, 2, var)
barplot(var_explained_seq)

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

# Sequential sum (for total "seq_sum" line)
gcv_seq_sum <- gcv_long %>%
  group_by(lambda) %>%
  summarise(Seqscore_sum = sum(Seqscore, na.rm = TRUE))

# Labels for sequential component lines (rightmost point)
seq_labels <- gcv_long %>%
  group_by(Component) %>%
  filter(lambda == max(lambda)) %>%
  slice_tail(n = 1)

library(ggrepel)

ggplot() +
  # Multiple sequential lines (one per component)
  geom_line(
    data = gcv_long,
    aes(x = lambda, y = Seqscore, group = Component, color = "seq"),
    linewidth = 0.8,
    alpha = 0.5
  ) +
  # Label each sequential component
  geom_text_repel(
    data = seq_labels,
    aes(x = lambda, y = Seqscore, label = Component),
    size = 3,
    color = "coral3",
    direction = "y",
    hjust = 0,
    nudge_x = 0.15,
    segment.color = "grey70",
    segment.size = 0.3,
    box.padding = 0.2,
    min.segment.length = 0
  ) +
  # Add the summed sequential line (bold coral)
  geom_line(
    data = gcv_seq_sum,
    aes(x = lambda, y = Seqscore_sum, color = "seq_sum"),
    linewidth = 1.5,
    linetype = "solid"
  ) +
  # Add the averaged subspace line (bold aquamarine)
  geom_line(
    data = gcv_sub_mean,
    aes(x = lambda, y = Subscore, color = "sub"),
    linewidth = 1.3
  ) +
  # Highlight minima
  geom_point(
    data = gcv_seq_sum %>%
      slice_min(Seqscore_sum, n = 1),
    aes(x = lambda, y = Seqscore_sum, color = "seq_sum"),
    shape = 21, size = 3, fill = "white", stroke = 1
  ) +
  geom_point(
    data = gcv_sub_mean %>%
      slice_min(Subscore, n = 1),
    aes(x = lambda, y = Subscore, color = "sub"),
    shape = 21, size = 3, fill = "white", stroke = 1
  ) +
  scale_x_log10(labels = scales::label_scientific(), limits = c(1e-3, 1e3)) +
  scale_y_log10(limits = c(1e-0, 1e3)) +
  scale_color_manual(
    values = c(
      "sub" = "aquamarine4",
      "seq" = "coral3",
      "seq_sum" = "firebrick3"
    ),
    labels = c(
      "sub" = "Subspace (avg)",
      "seq" = "Sequential (per fPC)",
      "seq_sum" = "Sequential (sum)"
    )
  ) +
  labs(
    x = expression(lambda),
    y = "GCV score",
    color = "Approach"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

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

cv_results <- list()
for (model_name in c("smv","tpsPCA","sequential","subspace")) {
    cv_results[[model_name]] <- cross_validate_model(model_name, domain, counts, locations, path_list, test_options, K = 5, seed = 1412, verbose = TRUE)
}