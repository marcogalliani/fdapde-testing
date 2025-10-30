## Function: extract_new_results
# - Args:
#   * results_evaluation: nested list containing results per model
#   * names_models: character vector with model identifiers
#   * name_result: character or vector specifying which result(s) to extract
# - Desc:
#   Extracts specific evaluation metrics from each model’s result structure.
#   Handles nested lists, converts time units to seconds, and replaces NULLs
#   with NaN values for consistency.
extract_new_results <- function(results_evaluation, names_models, name_result) {
  
  ## Room for new results
  new_results <- list()
  for (name_model in names_models) {
    ## Router for reading the data
    if (length(name_result) == 1) {
      new_results[[name_model]] <- results_evaluation[[name_model]][[name_result]]
    } else if (length(name_result) == 2) {
      new_results[[name_model]] <- results_evaluation[[name_model]][[name_result[1]]][[name_result[2]]]
    } else {
      stop()
    }
    ## Convert time units if necessary
    if (length(new_results[[name_model]]) == 1 &&
        "units" %in% names(attributes(new_results[[name_model]]))) {
      new_results[[name_model]] <- format_time(new_results[[name_model]])
    }
    ## Replace NULL with NaN
    if (is.null(new_results[[name_model]])) {
      new_results[[name_model]] <- c(NaN)
    }
  }
  return(new_results)
}

## Quantitative analysis ----
load_quantitative_results <- function(test_options, path_list){
  cat(paste0("\nLoading quantitative results for ", test_options$name_test, " ...\n"))
  batch_index <- 1 # test on GCV curves use a single batch
  ## Get model names, labels, and colors
  model_names  <- test_options$model_names
  model_labels <- test_options$model_labels
  model_colors <- test_options$model_colors
  ## Load the first batch
  batch_index <- 1
  ok <- tryCatch({
    path_batch <- file.path(path_list$results, paste0("batch_", batch_index))
    load(file.path(path_batch, paste0("batch_", batch_index, "_results_evaluation.RData")))
    TRUE
  }, error = function(e) {
    cat(sprintf("Error in test %s - batch %d: %s\n", test_options$name_test, batch_index, conditionMessage(e)))
    FALSE
  })
  if (!ok) next
  
  ## Load gcv scores and mse
  gcv_scores <- extract_new_results(results_evaluation, model_names, "gcv_scores")
  mse <- extract_new_results(results_evaluation, model_names, "mse")
  cat(sprintf("- Batch %d loaded\n", batch_index))

  return(list(
    gcv_scores = gcv_scores,
    mse = mse,
    model_names = model_names,
    model_labels = model_labels,
    model_colors = model_colors,
    varying_options = test_options$test_options$varying_options
  ))
}


## Function: load_qualitative_results
# - Args:
#   * test_options: list containing model and simulation parameters, including:
#       - $model_names: vector of model identifiers
#       - $model_labels: vector of model display names
#       - $model_options$n_comp: number of fPCs
#       - $dimensions$n_nodes_HR_grid: number of high-resolution grid points
#       - $test_options$n_reps: number of simulation repetitions (batches)
#   * data: list containing reference quantities and domain info:
#       - $domain: includes $fdapde_mesh and $boundary
#       - $locations: observation points
#       - $loadings_generator: function generating true fPC loadings
#       - $loadings_true, $loadings_true_locs
#   * path_list: list of paths where batch results and meshes are stored
#       - $results: directory containing batch subfolders
# - Desc:
#   Loads model outputs from each batch (scores, loadings, loadings at locations
#   and high-resolution grid) for all models. If available, it also evaluates
#   functional loadings on a regular grid using the domain mesh. Returns all
#   loaded data organized by model and batch, along with the true loadings and
#   evaluation domains.
load_qualitative_results <- function(test_options, data, path_list) {
  cat("\nLoading results for qualitative analysis ...\n")

  ## Number of components
  n_comp <- test_options$model_options$n_comp

  ## Generators
  loadings_generator <- data$loadings_generator

  ## Locations and grid
  nodes <- data$domain$fdapde_mesh$nodes
  locations <- data$locations
  grid <- spsample(
    data$domain$boundary,
    test_options$dimensions$n_nodes_HR_grid,
    "regular"
  )@coords

  ## Room for solutions
  scores <- list()
  loadings <- list()
  loadings_locs <- list()
  loadings_HR <- list()

  ## True fPCs at HR grid
  loadings_true_HR <- matrix(0, nrow = nrow(grid), ncol = n_comp)
  for (m in seq_len(n_comp)) {
    loadings_true_locs <- loadings_generator(locations, m)
    norm <- norm_l2(loadings_true_locs)
    loadings_true_HR[, m] <- loadings_generator(grid, m) / norm
  }

  ## Load batches
  n_reps <- test_options$test_options$n_reps
  for (batch_index in seq_len(n_reps)) {
    tryCatch(
      {
        path_batch <- paste0(path_list$results, "/", "batch_", batch_index, "/")

        for (name_model in test_options$model_names) {
          ## Load model file
          model_file <- paste0(
            path_batch,
            "batch_", batch_index, "_fitted_model_", name_model, ".RData"
          )

          load(model_file)
          model <- get(paste0("model_", name_model))

          ## Store results
          scores[[name_model]][[batch_index]] <- model$results$scores
          loadings_locs[[name_model]][[batch_index]] <- model$results$loadings_locs

          if ("loadings" %in% names(model$results)) {
            loadings[[name_model]][[batch_index]] <- model$results$loadings
            loadings_HR[[name_model]][[batch_index]] <- evaluate_field(
              grid, model$results$loadings, data$domain$fdapde_mesh
            )
          }
        }
      },
      error = function(e) {
        cat(
          paste(
            "Error in test ", test_options$name_test,
            " - batch ", batch_index, ": ",
            conditionMessage(e), "\n",
            sep = ""
          )
        )
      }
    )
    cat(paste("- Batch", batch_index, "loaded\n"))
  }

  ## Return results
  return(list(
    model_names = test_options$model_names,
    model_labels = test_options$model_labels,
    scores = scores,
    loadings_true = data$loadings_true,
    loadings_true_locs = data$loadings_true_locs,
    loadings_true_HR = loadings_true_HR,
    loadings = loadings,
    loadings_locs = loadings_locs,
    loadings_HR = loadings_HR,
    domain = data$domain,
    nodes = nodes,
    locations = locations,
    grid = grid
  ))
}