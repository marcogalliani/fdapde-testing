# = ========================================================================== =
# - Script: load_qualitative_results.R
# - Desc: Loads model outputs from all simulation batches for qualitative analysis.
#         Reconstructs functional principal components (fPCs) at multiple spatial
#         resolutions — nodes, observed locations, and a high-resolution grid —
#         for each model and repetition. Also computes true fPCs at the same grid
#         for visual comparison.
# = ========================================================================== =


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
