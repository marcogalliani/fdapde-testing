# = ========================================================================== =
# - Test: Example data decomposition — Aggregate results
# - Desc: Loads quantitative results for all options of a selected test and
#         returns them as a structured list, then aggregates the data in plots
# - Args:
#     [1] name_main_test : name of the main test to scan (e.g., "test1")
# = ========================================================================== =

rm(list = ls())
graphics.off()
options(warn = -1)

# README ----
# Assumes all batches for each option have already been run and saved.
# For each option JSON in the queue/<test>/, this script loads the corresponding
# quantitative results (via load_quantitative_results) and collects them.

# Configuration ----

## Load libraries ----
invisible(suppressMessages(sapply(c(
  # discretization
  "fdaPDE", "femR",
  # algebraic utils
  "pracma",
  # data manipulation
  "MASS", "tidyr", "dplyr",
  # keep plotting deps (not used here) to avoid missing imports in utils
  "ggplot2", "viridis", "stringr", "RColorBrewer", "grid", "gridExtra",
  # json
  "jsonlite",
  # sampling
  "sf", "sp", "raster"
), require, character.only = TRUE)))

## Load functions ----
source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/options.R")
source("src/utils/load_results_utils.R")
source("src/utils/plotting_utils.R")

## Load configuration file
path_this <- get_script_path()
source(paste0(path_this, "config.R"))

## Load test-specific helpers
source(paste0("tests/", test_suite, "/utils/generate_options.R"))

## Create suite directories ----
path_list <- create_paths(test_suite)

# Loader ----

## Generate options ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  INTERACTIVE <- TRUE
  name_main_test <- "test1"
} else {
  INTERACTIVE <- FALSE
  name_main_test <- args[1]
}

## Prepare queue/log dirs for this test
path_list$queue <- paste0(path_list$queue, name_main_test, "/")
path_list$logs <- paste0(path_list$logs, name_main_test, "/")
mkdir(c(path_list$queue, path_list$logs))

## Generate all option files for this test
generate_options(test_suite, name_main_test, path_list$queue)


## Load results ----

## Load results (all option combinations, all batches)
loaded_results <- load_all_quantitiative_results(path_list, name_main_test)


# Plot aggregated results ----

# loaded_results$varying_options
if (is.null(order) || length(order) != length(loaded_results$varying_options)) {
  order <- 1:length(loaded_results$varying_options) # default
}


### Time complexity ----

## Set plots parameters
data_plot <- loaded_results$execution_time
title_prefix <- "Execution times w.r.t the"
values_name <- "Time [seconds]"
limits <- c(0, max(data_plot[loaded_results$model_names]))

WIDTH <- 15
HEIGHT <- 7.5

plots_catalog <- list(
  boxplots = TRUE,
  lines = TRUE,
  logx = TRUE,
  loglog = TRUE,
  normalized = TRUE
)


## Open a pdf where to save the plots
pdf(paste(path_list$images, name_main_test, "/time_complexity.pdf", sep = ""), width = WIDTH, height = HEIGHT)
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits, plots_catalog = plots_catalog
)
dev.off()



### RMSE ----

## Open a pdf where to save the plots
pdf(paste(path_list$images, name_main_test, "/rmse.pdf", sep = ""), width = WIDTH, height = HEIGHT)


#### Reconstruction error at locations ----

## Set plots parameters
data_plot <- loaded_results$rmse$reconstruction_locs
title_prefix <- "RMSE[Reconstruction] at locations w.r.t the"
values_name <- "RMSE"
limits <- c(0, max(data_plot[loaded_results$model_names]))

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits
)

#### Scores orthogonality ----

## Set plots parameters
data_plot <- loaded_results$rmse$scores_orth
title_prefix <- "RMSE[Scores orthogonality] at locations w.r.t the"
values_name <- "RMSE"
limits <- c(0, max(data_plot[loaded_results$model_names]))

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits
)

#### Scores ----

## Set plots parameters
data_plot <- loaded_results$rmse$scores
title_prefix <- "RMSE[Scores] w.r.t the"
values_name <- "RMSE"
limits <- c(0, max(data_plot[loaded_results$model_names]))

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits
)

#### Loadings at locations ----

## Set plots parameters
data_plot <- loaded_results$rmse$loadings_locs
title_prefix <- "RMSE[Loadings] at locations w.r.t the"
values_name <- "RMSE"
limits <- c(0, max(data_plot[loaded_results$model_names]))

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits
)

## Close pdf
dev.off()

### Angles ----

## Open a pdf where to save the plots
pdf(paste(path_list$images, name_main_test, "/angles.pdf", sep = ""), width = WIDTH, height = HEIGHT)

## Set plots parameters
data_plot <- loaded_results$angles$components_m
data_plot[loaded_results$model_names] <- log10(data_plot[loaded_results$model_names])
title_prefix <- "Angles between true and estimated functions at locations"
values_name <- "angle"
limits <- c(0, 10)

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names, 
  order = order, limits = limits
)
## Close pdf
dev.off()

### Regularization ----

## Open a pdf where to save the plots
pdf(paste(path_list$images, name_main_test, "/regularization.pdf", sep = ""), width = WIDTH, height = HEIGHT)

## Set plots parameters
data_plot <- loaded_results$lambdas
data_plot[loaded_results$model_names] <- log10(data_plot[loaded_results$model_names])
title_prefix <- "lambda selected w.r.t the"
values_name <- "lambda"
limits <- c(-12, 1)

## Plot aggregated results
plot.aggregated_data(
  loaded_results, data_plot, title_prefix, values_names,
  order = order, limits = limits
)

## Close pdf
dev.off()
