# = ========================================================================== =
# - Test: Example data decomposition - Inspect results
# - Desc: orchestrates quantitative and qualitative analyses for a selected test
#         configuration: loads options and data, prepares directories, and renders
#         summary plots to PDF.
# - Args (when calling it from terminal):
#   [1] name_main_test : name of the main test to run (e.g., "test1")
#   [2] file_options : JSON file containing the test options to use
# = ========================================================================== =

rm(list = ls())
graphics.off()
options(warn = -1)


# README ----
# Assumes all batches for each option have already been run and saved.
# This script loads the quantitative and qualitative results for the
# desired test and option.



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
source("src/utils/error_metrics.R")
source("src/utils/plotting_utils.R")
source("src/utils/load_results_utils.R")
sapply(list.files("src/data-generation", pattern = "\\.R$", full.names = TRUE), source)

## Load configuration file
path_this <- get_script_path()
source(paste0(path_this, "config.R"))

## Load test-specific functions
source("src/wrappers/fPCA.R")
source(paste0("tests/", test_suite, "/utils/plot_results.R"))
source(paste0("tests/", test_suite, "/utils/load_qualitative_results.R"))

## Create suite directories ----
## Define and create work directories for the test suite
path_list <- create_paths(test_suite)

# Results analysis ----
cat.section_title("Results analysis")

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
  file_options <- file_options_list[7] ## <====== input here
}

## Load selected options
test_options <- fromJSON(paste0(path_list$queue, file_options))

## Update and create work directories for the test selected
path_list <- update_paths(path_list, name_main_test, test_options)


## Load data ----
load(paste0(path_list$data, test_options$name_test, ".RData"))


## Create sub-directory ----
path_list$images <- paste0(path_list$images, test_options$name_test, "/")
mkdir(path_list$images)


## Quantitative analysis ----
cat.subsection_title("Quantitative analysis")

## Load data
loaded_qnt_results <- load_quantitative_results(test_options, path_list)

## Plot
pdf(file = paste(path_list$images, test_options$name_test, "_quantitative.pdf", sep = ""))
plot_quantitative_results(loaded_qnt_results)
dev.off()

## Quantitative analysis ----
cat.subsection_title("Qualitative analysis")

## Load data
loaded_qlt_results <- load_qualitative_results(test_options, data, path_list)

## Plot
pdf(file = paste(path_list$images, test_options$name_test, "_qualitative.pdf", sep = ""), width = 15,height=10)
plot_qualitative_results(loaded_qnt_results, loaded_qlt_results)
dev.off()

## Optionally open the results directory
if (INTERACTIVE) {
  open(path_list$images)
}
