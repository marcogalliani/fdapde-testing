# = ========================================================================== =
# - Script: generate_options.R
# - Desc: Generates JSON option files for test configurations.
# = ========================================================================== =

## Function: generate_options(test_suite, name_main_test, path_queue)
# - Args:
#   * test_suite: name of the calling test suite (used for directory structure)
#   * name_main_test: identifier of the specific test to generate options for
#   * path_queue: directory where JSON files will be written
# - Desc: 
#   Defines model parameters, expands selected grid options, and writes the 
#   resulting combinations to JSON files ready for execution.
generate_options <- function(test_suite, name_main_test, path_queue) {
  
  ## Create the directory (if it does not exist yet)
  mkdir(c(path_queue))
  
  ## Names of the models you want to compare
  # - model_names: used for indexing (no spaces, please)
  # - model_labels: used for plotting
  model_names <- c(
    "mv", "sequential", "subspace"
  )
  model_labels <- c(
    "mv", "sequential", "subspace"
  )
  
  ## Define the color palette
  model_colors <- brewer.pal(length(model_labels), "Set1")
  
  ## Options that you want to be common across tests 
  lambda_grid <- pracma::logspace(-2, 3, n = 20)
  seed <- 1412
  
  switch(
    name_main_test,
    test1 = {
      
      ## Set the desired options
      options <- list(
        model_names = model_names,
        model_labels = model_labels,
        model_colors = model_colors,
        cpp_script = "fPCA-2D",
        test_options = list(
          varying_options = c("n_nodes")
        ),
        domain_and_locations = list(
          name_mesh = "HER2"
        ),
        dimensions = list(
          n_nodes = c(1000),       # vectors to be combined
          n_nodes_HR_grid = 1000
        ),
        model_options = list(
          n_comp = 5,
          mean = TRUE
        ),
        regularization = list(
          lambda_grid = lambda_grid
        )
      )
      
      ## File naming policy
      name_fun <- function(opts_i, comb_row) {
        paste(
          name_main_test,
          ## Include all the varying options!
          "nn", sprintf("%04d", comb_row$n_nodes),
          sep = "_"
        )
      }
      
      ## Expand ONLY the varying options
      options_list <- explode_options(
        options,
        by = options$test_options$varying_options,
        name_fun = name_fun
      )
      
      ## Write JSON files
      write_options_json(
        options_list,
        dir = path_queue,
        name_field = "name_test"
      )
    },
    
    {
      stop(paste("The test", name_main_test, "does not exist"))
    }
  )
}