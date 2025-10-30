# = ========================================================================== =
# - Script: fit_and_evaluate_models.R
# - Desc: Fits one or more models via external C++ tools and evaluates results.
#         Handles data/mesh export, parameter JSON creation, caching of fits,
#         and storage of evaluation outputs for each batch.
# = ========================================================================== =


# - Function: fit_and_evaluate_models
# - Args:
#   * path_list: list of directories (expects $cpp_script, $batch, $tmp_data, $tmp_results)
#   * data: list with fields $X (matrix) and $locations (matrix/data.frame)
#   * domain: list with $fdapde_mesh (fdaPDE mesh)
#   * batch_index: integer, identifier for the current batch
#   * test_options: nested list with $model_names, $regularization, etc.
# - Desc:
#   Exports data/mesh, prepares parameters for C++ executables, runs model
#   fitting if required, optionally evaluates fitted models, and saves results.
fit_and_evaluate_models <- function(path_list, 
                                    data,
                                    domain,
                                    batch_index,
                                    test_options){
  
  # Room for results ----
  results_evaluation <- list()
  
  # Paths ----
  path_batch <- path_list$batch
  
  # Load results if available ----
  ## Reload previously saved evaluation results (if present)
  if (file.exists(paste0(path_batch, "batch_", batch_index, "_results_evaluation.RData"))) {
    load(paste0(path_batch, "batch_", batch_index, "_results_evaluation.RData"))
  }
  
  
  # Fit and Evaluate ----
  for (model_name in test_options$model_names) { # model_name <- test_options$model_names[1]
    
    ## Initialize empty model
    model <- NULL
    
    ## File name where the results should be found
    file_model <- paste(path_batch, "batch_", batch_index, "_fitted_model_", model_name, ".RData", sep = "")
    
    ## Fit the model only if necessary (no fit found or fit is forced)
    if (file.exists(file_model) && !FORCE_FIT) {
      if (FORCE_EVALUATE) {
        cat("- Loading fitted model:", model_name, "... \n")
        load(file_model)
      }
    } else {
      cat("- Fitting model:", model_name, "... ")
      
      ## Fit the calibrated model to get gcv scores
      model <- fit_model(model_name, domain, data, path_list, test_options)

      ## Adjust results
      model <- adjust_results(model, data)
        
      ## Fit #lambda_grid uncalibrated models
      ## Write JSON arguments for the C++ solver ----
      cpp_script_arguments <- list()
      cpp_script_arguments$path_list <- list(
        mesh = paste0(path_list$tmp_data, "mesh/"),
        data = path_list$tmp_data,
        results = path_list$tmp_results
      )
      cpp_script_arguments$options$solver <- model_name
      cpp_script_arguments$options$n_comp <- test_options$model_options$n_comp
      ## write cpp params
      file_name_params <- paste0(
        test_options$name_test, "_", model_name, "_batch",
        test_options$batch_index, "_params.json"
      )
      ## compute mse
      mse <- numeric(length(test_options$regularization$lambda_grid))
      for(j in 1:length(test_options$regularization$lambda_grid)){
        ## update the json for each value of lambda
        cpp_script_arguments$options$lambda_grid <- I(test_options$regularization$lambda_grid[j])
        write_json(
          path = paste0(path_list$cpp_script, file_name_params),
          cpp_script_arguments,
          auto_unbox = TRUE,
          pretty = TRUE,
          digits = 10
        )
        ## run the cpp script
        system(paste0("cd ", path_list$cpp_script, " && ", "./fit_model_uncal ", file_name_params),
                ignore.stdout = IGNORE_CPP_OUTPUT)
        ## evaluate mse        
        X_hat <- as.matrix(read.csv(paste(path_list$tmp_results, "reconstruction_at_locs.csv", sep = "")))
        mse[j] <- mean((X_hat-data$X_true_locs)^2)
      }
      model$results$mse <- mse
      ## Save fitted model
      assign(paste("model_", model_name, sep = ""), model)
      save(
        index_batch = batch_index,
        list = paste("model_", model_name, sep = ""),
        file = file_model
      )
      rm(list = paste("model_", model_name, sep = ""))
    }
    
    if (!is.null(model)) {
      ## Model evaluation ----
      results_evaluation[[model_name]] <- evaluate_results(model, data)
    }
  }
  
  # Save results of the evaluation ----
  save(
    index_batch = batch_index,
    results_evaluation,
    file = paste(path_batch, "batch_", batch_index, "_results_evaluation.RData", sep = "")
  )
  cat(paste("- Batch", batch_index, "completed.\n"))
}