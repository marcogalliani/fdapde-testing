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
      
      ## Fit the model
      model <- list()
      # Paths ----
      path_cpp_script  <- path_list$cpp_script
      path_batch       <- path_list$batch
      path_tmp_data    <- path_list$tmp_data
      path_mesh        <- paste0(path_list$tmp_data, "mesh/")
      mkdir(path_mesh)
      path_tmp_results <- path_list$tmp_results
      
      # Write data for C++ scripts ----
      
      ## Data matrix and locations ----
      write.csv(format(data$X, digits = 16), file = paste0(path_tmp_data, "X.csv"))
      write.csv(format(data$locations, digits = 16), file = paste0(path_tmp_data, "locs.csv"))
      
      ## Mesh ----
      mesh <- domain$fdapde_mesh
      write.csv(format(mesh$nodes, digits = 16), paste0(path_mesh, "points.csv"))
      write.csv(format(mesh$triangles, digits = 16), paste0(path_mesh, "elements.csv"))
      write.csv(format(1 * mesh$nodesmarkers, digits = 16), paste0(path_mesh, "boundary.csv"))
      write.csv(format(mesh$neighbors, digits = 16), paste0(path_mesh, "neigh.csv"))
      write.csv(format(mesh$edges, digits = 16), paste0(path_mesh, "edges.csv"))
      ## Write JSON arguments for the C++ solver ----
      cpp_script_arguments <- list()
      cpp_script_arguments$path_list <- list(
        mesh = paste0(path_list$tmp_data, "mesh/"),
        data = path_list$tmp_data,
        results = path_list$tmp_results
      )
      cpp_script_arguments$options$solver <- model_name
      cpp_script_arguments$options$n_comp <- test_options$model_options$n_comp
      cpp_script_arguments$options$lambda_grid <- I(test_options$regularization$lambda)
      ## write cpp params
      file_name_params <- paste0(
        test_options$name_test, "_", model_name, "_batch",
        test_options$batch_index, "_params.json"
      )
      write_json(
        path = paste0(path_list$cpp_script, file_name_params),
        cpp_script_arguments,
        auto_unbox = TRUE,
        pretty = TRUE,
        digits = 10
      ) 
      ## run the cpp script
      start.time <- Sys.time()
      system(paste0("cd ", path_list$cpp_script, " && ", "./fit_model_uncal ", file_name_params),
              ignore.stdout = IGNORE_CPP_OUTPUT)
      end.time <- Sys.time()
      cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
      ## read results
      path_tmp_results <- path_list$tmp_results
      model$results$loadings <- as.matrix(read.csv(paste0(path_tmp_results, "loadings.csv")))
      model$results$loadings_locs <- as.matrix(read.csv(paste0(path_tmp_results, "loadings_locs.csv")))
      model$results$scores <- as.matrix(read.csv(paste0(path_tmp_results, "scores.csv")))
      model$results$X_hat <- as.matrix(read.csv(paste0(path_tmp_results, "reconstruction.csv")))
      model$results$X_hat_locs <- as.matrix(read.csv(paste0(path_tmp_results, "reconstruction_at_locs.csv")))
      model$results$lambda <- as.matrix(read.csv(paste0(path_tmp_results, "lambda.csv")))
      # Add flags ----
      model$model_traits$is_functional <- FALSE
      model$model_traits$has_interpolator <- FALSE
      ## Execution time
      model$results$execution_time <- end.time - start.time
      
      ## Adjust results
      model <- adjust_results(model, data)
      
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