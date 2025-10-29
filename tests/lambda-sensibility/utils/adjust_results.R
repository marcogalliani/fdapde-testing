# = ========================================================================== =
# - Script: adjust_results.R
# - Desc: Provides normalization and alignment utilities for model outputs.
#         Includes functions to adjust loadings and scores by their l2 norms,
#         and to align estimated components with true ones by correcting signs
#         and scaling factors while preserving data reconstruction consistency.
# = ========================================================================== =

# - Function: adjust_norms
# - Args:
#   * FF: matrix of loadings evaluated at locations (n_locs x n_comp)
#   * SS: matrix of scores (n_stat_units x n_comp)
# - Desc:
#   Adjusts loadings and scores by component-wise norms to keep reconstruction
#   invariant. Returns adjusted loadings at locations, scores, and the norms.
adjust_norms <- function(FF, SS) {
  ## Number of components
  n_comp <- ncol(SS)
  ## Component norms (placeholder as in original code)
  f_norms <- ones(ncol(FF))
  for (h in 1:n_comp) {
    f_norms[h] <- norm_l2(FF[, h])
    FF[, h] <- FF[, h] / f_norms[h]
    SS[, h] <- SS[, h] * f_norms[h]
  }
  return(list(loadings_locs = FF, scores = SS, norms = f_norms))
}

# - Function: adjust_results
# - Args:
#   * F_hat_locs: estimated loadings at locations (n_locs x n_comp)
#   * S_hat: estimated scores (n_stat_units x n_comp)
#   * F_true_locs: true loadings at locations (n_locs x n_comp)
#   * Fs_hat_evaluated: optional list of loadings evaluated on nodes (for functional metrics)
# - Desc:
#   Aligns the sign of estimated components to match the true ones, normalizes
#   results at locations, and adjusts node-evaluated loadings accordingly.
adjust_results <- function(model, data) {
  ## Get quantities to be adjusted
  F_hat_locs <- model$results$loadings_locs
  S_hat <- model$results$scores
  if ("loadings" %in% names(model$results)) {
    Fs_hat_evaluated <- list(loadings = model$results$loadings)
  } else {
    Fs_hat_evaluated <- NULL
  }

  ## Get reference
  F_true_locs <- data$loadings_true_locs

  ## Number of components
  n_comp <- ncol(S_hat)

  ## Change signs to match the true ones
  for (h in 1:n_comp) {
    if (RMSE(F_hat_locs[, h] + F_true_locs[, h]) < RMSE(F_hat_locs[, h] - F_true_locs[, h])) {
      F_hat_locs[, h] <- -F_hat_locs[, h]
      S_hat[, h] <- -S_hat[, h]
      if (!is.null(Fs_hat_evaluated)) {
        for (i in 1:length(Fs_hat_evaluated)) {
          Fs_hat_evaluated[[i]][, h] <- -Fs_hat_evaluated[[i]][, h]
        }
      }
    }
  }

  ## Normalize results at locations
  adjusted_results <- adjust_norms(F_hat_locs, S_hat)

  ## Adjust data at nodes accordingly
  for (h in 1:n_comp) {
    if (!is.null(Fs_hat_evaluated)) {
      for (i in 1:length(Fs_hat_evaluated)) {
        Fs_hat_evaluated[[i]][, h] <- Fs_hat_evaluated[[i]][, h] / adjusted_results$norms[h]
      }
    }
  }

  ## Update the model results
  model$results$loadings_locs <- adjusted_results$loadings_locs
  model$results$scores <- adjusted_results$scores
  if ("loadings" %in% names(model$results)) {
    model$results$loadings <- Fs_hat_evaluated$loadings
  }


  return(model)
}
