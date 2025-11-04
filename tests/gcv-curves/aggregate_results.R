rm(list = ls())
graphics.off()
options(warn = -1)


# README ----
# Assumes all batches for each option have already been run and saved.
# This script loads the quantitative results for all the options and plot aggregated results

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
sapply(list.files("src/data-generation", pattern = "\\.R$", full.names = TRUE), source)

## Load configuration file
path_this <- get_script_path()
source(paste0(path_this, "config.R"))

## Load test-specific functions
source("src/wrappers/fPCA.R")
source(paste0("tests/", test_suite, "/utils/load_results.R"))
source(paste0("tests/", test_suite, "/utils/plot_results.R"))

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

## Quantitative analysis ----
cat.subsection_title("Quantitative analysis")

## Concatenate quantitive plots (GCV curves)
source(paste0("tests/", test_suite, "/utils/plot_results.R"))




varying_option <- "NSR"
NSR_vect <- numeric(length(file_options_list))

plot_list <- list()
for(i in 1:length(file_options_list)){
    ## Update directories according to the new test
    path_list <- create_paths(test_suite)
    path_list$queue <- paste0(path_list$queue, name_main_test, "/")
    path_list$logs <- paste0(path_list$logs, name_main_test, "/")
    mkdir(c(path_list$queue, path_list$logs))
    ## Load selected options
    test_options <- fromJSON(paste0(path_list$queue, file_options_list[[i]]))
    path_list <- update_paths(path_list, name_main_test, test_options)
    ## Load data
    NSR_vect[i] <- test_options$noise$NSR
    loaded_qnt_results <- load_quantitative_results(test_options, path_list)
    plot_list[[file_options_list[[i]]]] <- plot_quantitative_results(loaded_qnt_results)
}


## gridExtra
# Custom function to extract the legend grob
get_legend <- function(my_ggplot) {
  # Convert the ggplot to a gtable object
  tmp <- ggplotGrob(my_ggplot)
  
  # Find the position of the legend grob ("guide-box")
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  # Extract the legend grob
  legend <- tmp$grobs[[leg]]
  return(legend)
}

all_single_plots <- unlist(lapply(plot_list, as.list), 
                        recursive = TRUE)
shared_legend <- get_legend(all_single_plots[[1]])
all_single_plots <- lapply(all_single_plots, function(p) {
  p + theme(plot.title = element_blank(),
  axis.text.x = element_text(size = 8),
  legend.position='none')
})



grid <- arrangeGrob(grobs=all_single_plots, ncol=4)
labled_grid <- labled_plots_grid(grid, NULL, test_options$model_labels, paste0("NSR=",NSR_vect))


pdf(paste(path_list$images, "/aggregated_gcv_curves.pdf", sep = ""), width = 10, height = 15)

grid.arrange(labled_grid,shared_legend,heights = c(20, 1))
dev.off()

## Optionally open the results directory
if (INTERACTIVE) {
  open(path_list$images)
}
