# = ========================================================================== =
# - Script: mnodels_evaluation.R
# - Desc: Utilities for adjusting model outputs and evaluating performance
#         via RMSE/IRMSE and angular measures on loadings and subspaces.
# = ========================================================================== =


# - Function: evaluate_results
# - Args:
#   * model: fitted model object with $results and $model_traits fields
#   * data: list from data generator with true quantities
# - Desc:
#   Computes RMSE/IRMSE metrics for centering, loadings, scores, and
#   reconstructions; and angular measures for components and subspaces.
evaluate_results <- function(model, data) {
  
  ## Number of computed components ----
  n_comp <- data$dimensions$n_comp
  
  ## Room for results ----
  rmse <- list()
  irmse <- list()
  angles <- list()
  lambdas <- numeric(n_comp)
  
  ## Execution time ----
  execution_time <- model$results$execution_time
  
  ## Lambdas ----
  lambdas <- as.vector(model$results$lambda)
  
  ## RMSE at locations ----
  
  ## Centering
  if (!is.null(model$results$X_mean_locs)) {
    norm <- ifelse(RMSE(data$X_mean_true_locs) == 0, 1, RMSE(data$X_mean_true_locs))
    rmse$centering_locs <- RMSE(model$results$X_mean_locs - data$X_mean_true_locs) / norm
  }
  
  ### Loadings & scores ----
  for (h in 1:n_comp) {
    norm <- RMSE(data$loadings_true_locs[, h])
    rmse$loadings_locs[h] <- RMSE(model$results$loadings_locs[, h] - data$loadings_true_locs[, h]) / norm
    norm <- RMSE(data$scores_true[, h])
    rmse$scores[h] <- RMSE(model$results$scores[, h] - data$scores_true[, h]) / norm
  }
  scores_normalized <- model$results$scores
  scores_norms <- apply(scores_normalized, MARGIN = 2, function(x) { sqrt(sum(x^2)) })
  scores_normalized <- sweep(scores_normalized, MARGIN = 2, scores_norms, FUN = "/")
  rmse$scores_orth <- RMSE(diag(n_comp) - t(scores_normalized) %*% scores_normalized)
  
  ### Data reconstruction ----
  norm <- RMSE(data$X_true_locs)
  rmse$reconstruction_locs <- RMSE(model$results$X_hat_locs - data$X_true_locs) / norm
  
  ### RMSE at nodes (if possible) ----
  if (model$model_traits$has_interpolator) {
    norm <- ifelse(RMSE(data$X_mean_true) == 0, 1, RMSE(data$X_mean_true))
    rmse$centering <- RMSE(model$results$X_mean - data$X_mean_true) / norm
    for (h in 1:n_comp) {
      norm <- RMSE(data$loadings_true[, h])
      rmse$loadings[h] <- RMSE(model$results$loadings[, h] - data$loadings_true[, h]) / norm
    }
    norm <- RMSE(data$X_true)
    rmse$reconstruction <- RMSE(model$results$X_hat - data$X_true) / norm
  }
  
  ## IRMSE (if possible) ----
  if (model$model_traits$is_functional) {
    if (!is.null(model$results$X_mean)) {
      norm <- IRMSE(data$X_mean_true, model$R0())
      norm <- ifelse(norm == 0, 1, norm)
      irmse$centering <- IRMSE(model$results$X_mean - data$X_mean_true, model$R0()) / norm
    }
    for (h in 1:n_comp) {
      norm <- IRMSE(data$loadings_true[, h], model$R0())
      irmse$loadings[h] <- IRMSE(model$results$loadings[, h] - data$loadings_true[, h], model$R0()) / norm
    }
    norm <- IRMSE(t(data$X_true), model$R0())
    irmse$reconstruction <- IRMSE(t(model$results$X_hat - data$X_true), model$R0()) / norm
  }
  
  ## Angles ----
  for (h in 1:n_comp) {
    angles$subspaces_m[h] <- 180 * subspace(model$results$loadings_locs[, 1:h], data$loadings_true_locs[, 1:h]) / pi
    angles$components_m[h] <- 180 * subspace(model$results$loadings_locs[, h], data$loadings_true_locs[, h]) / pi
    if (model$model_traits$is_functional) {
      angles$components_f[h] <- 180 * angle_between_functions(
        model$results$loadings[, h],
        data$loadings_true[, h],
        model$R0()
      ) / pi
    }
    if (h < n_comp) {
      for (j in (h + 1):n_comp) {
        angles$orthogonality_m[h + j - 2] <- 180 * subspace(model$results$loadings_locs[, h], model$results$loadings_locs[, j]) / pi
        if (model$model_traits$is_functional) {
          angles$orthogonality_f[h + j - 2] <- 180 * angle_between_functions(
            model$results$loadings[, h],
            model$results$loadings[, j],
            model$R0()
          ) / pi
        }
      }
    }
  }
  
  return(list(
    execution_time = execution_time,
    lambdas = lambdas,
    rmse = rmse,
    irmse = irmse,
    angles = angles
  ))
}