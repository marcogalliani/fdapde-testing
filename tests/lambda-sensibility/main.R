# = ========================================================================== =
# - Test: Example data decomposition
# - Desc: Automates a full test cycle: picks an option, runs all methods,
#         and saves results, logs, and plots for later analysis.
# - Args (when calling it from terminal):
#   [1] name_main_test : name of the main test to run (e.g., "test1")
#   [2] file_options : JSON file containing the test options to use
# = ========================================================================== =

rm(list = ls())
graphics.off()
options(warn = -1)


# README ----

## Welcome!
# Hi there, this script aims to make your life easier when it comes
# to testing your brand-new method in different scenarios.

## IMPORTANT!
# Before starting:
# - Remember to set your working directory to the root directory.
# - Updated the config.R file in the directory of the test
# - Type "make help" in your terminal to discover all the shortcuts!


# Configuration ----

## Load libraries ----

invisible(suppressMessages(sapply(c(
  # discretization
  "fdaPDE", "femR",
  # algebraic utils
  "pracma",
  # data manipulation
  "MASS", "tidyr", "dplyr",
  # visualization
  "ggplot2", "viridis", "stringr", "RColorBrewer", "grid", "gridExtra",
  # json
  "jsonlite",
  # sampling
  "sf", "sp", "raster"
), require, character.only = TRUE)))


## Load functions ----

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
source(paste0("tests/", test_suite, "/utils/fit_and_evaluate.R"))
source(paste0("tests/", test_suite, "/utils/adjust_results.R"))
source(paste0("tests/", test_suite, "/utils/models_evaluation.R"))
source(paste0("tests/", test_suite, "/utils/plot_results.R"))
source(paste0("tests/", test_suite, "/utils/load_qualitative_results.R"))


## Create suite directories ----

## Define and create work directories for the test suite
path_list <- create_paths(test_suite)


# Test ----

## Options ----

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

## Log file
file_log_global <- paste0(path_list$logs, "log.txt")
file_log_specific <- paste0(path_list$logs, "log_", test_options$name_test, ".txt")
if (!INTERACTIVE) {
  ## Console
  cat("- Running:", test_suite, "/", test_options$name_test, "\n")
  ## Global log
  sink(file_log_global, append = TRUE)
  cat("- Running:", test_suite, "/", test_options$name_test, "\n")
  sink()
  ## Test specific log
  sink(file_log_specific, append = TRUE)
}

## Options visualization
cat.script_title(paste("Test:", TEST_SUITE))
cat.section_title("Options")
cat.json(test_options)


## Load domain and locations ----

## Generate the domain and its mesh
domain <- generate_domain(
  test_options$domain_and_locations$name_mesh,
  test_options$dimensions$n_nodes
)

## Sample the locations
locations <- generate_locations(
  domain,
  test_options$domain_and_locations$locs_eq_nodes,
  test_options$dimensions$n_locs
)

## Plot domain and locations
plot.points(
  locations = locations,
  boundary = domain$boundary,
  group_colors = "darkblue", size = 1
) + std_plot_settings_fields() + ggtitle("Domain and locations")


## Select generators ----

## Data generation utilities
source("src/data-generation/fpca/generate_2D_fpca_data.R")

## Select the desired generators
generate_data <- generate_2D_fpca_data
generate_loadings_true <- function(locs, i) {
  translated_laplacian_eigenfunction(locs, i, x_t = 0.2, y_t = 0)
}
generate_mean_true <- function(locs) {
  if (test_options$data$mean) {
    return(log_mean_generator(locs))
  } else {
    return(locs[, 1] * 0)
  }
}


## Recursive fit ----

## Fit the models n_reps times
if (RUN$tests) {
  for (batch_idx in 1:test_options$test_options$n_reps) {
    cat(paste0("\nBatch ", batch_idx, ":\n"))

    ## Create batch directory
    path_list$batch <- paste0(path_list$results, "batch_", batch_idx, "/")
    mkdir(path_list$batch)

    ### Generate data ----
    cat("- Generate data\n")

    ## File names where the results should be found
    file_model_vect <- paste0(path_list$batch, "batch_", batch_idx, "_fitted_model_", test_options$model_names, ".RData")

    ## Generate data only if necessary (no fit found of fit is forced)
    if (any(!file.exists(file_model_vect)) || FORCE_FIT || FORCE_EVALUATE) {
      data <- generate_data(
        domain = domain, locations = locations,
        test_options = test_options,
        loadings_generator = generate_loadings_true,
        mean_generator = generate_mean_true,
        seed = 4 * batch_idx + test_options$noise$seed
      )
    } else {
      cat("Skipped, data are not necessary!\n")
    }

    ## Fit and evaluation ----
    cat("- Fit models\n")

    fit_and_evaluate_models(
      path_list = path_list,
      data = data,
      domain = domain,
      batch_index = batch_idx,
      test_options = test_options
    )

    ## Save data for qualitative results analysis ----
    if (batch_idx == 1) {
      save(data, file = paste0(path_list$data, test_options$name_test, ".RData"))
    }
  }
} else {
  cat("Skipped, relying on the saved results!\n")
}

# Test finalization ----

## Remove options file from the queue ----
file.remove(paste0(path_list$queue, file_options))

## Close the log file
if (!INTERACTIVE) {
  sink()
}
