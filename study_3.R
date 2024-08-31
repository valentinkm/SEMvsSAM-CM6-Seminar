# study_3.R

# Load necessary scripts and libraries
source("gen_pop_data.R")
source("calc_metrics.R")
source("run_analysis.R")
source("gen_mat.R")
source("gen_pop_lavsyntax.R")
library(furrr)
library(parallel)
library(dplyr)
library(purrr)
library(tidyr)
library(progress)
library(lavaan)

# Generate seeds
parallel_seeds <- function(n, seed = NULL) {
  if (is.null(seed))
    stop("seed must be provided.")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  purrr::accumulate(seq_len(n - 1), function(s, x) parallel::nextRNGStream(s), 
                    .init = .Random.seed)
}

# Generate parameters grid with seeds
n_reps <- 100
params <- expand.grid(
  seed = parallel_seeds(n_reps, seed = 42),
  model_type = c("3.1", "3.2", "3.1_negative", "3.2_negative"),
  N = c(50, 100, 250, 400),
  reliability = c(0.3, 0.5, 0.7),
  method = c("SEM", "gSAM", "lSAM_ML", "lSAM_ULS")
)


# Set population values
B_true <- c(
  'f3~f1' = 0.1, 'f3~f2' = 0.1, 'f3~f4' = 0.1,
  'f4~f1' = 0.1, 'f4~f2' = 0.1,
  'f5~f3' = 0.1, 'f5~f4' = 0.1, "f5~f1" = 0.1
)
true_values <- list(
  B = B_true
)

sam_model_syntax <- "
    # Measurement 
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f3 =~ y7 + y8 + y9
    f4 =~ y10 + y11 + y12
    f5 =~ y13 + y14 + y15
    
    # Structural 
    f3 ~ f1 + f2 + f4
    f4 ~ f1 + f2
    f5 ~ f3 + f4 + f1
"

sem_model_syntax <- "
  # Structural 
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9
  f4 =~ y10 + y11 + y12
  f5 =~ y13 + y14 + y15
  
  # lower bound factor variances:
  f1 ~~ v_f1*f1
  f2 ~~ v_f2*f2
  f3 ~~ v_f3*f3
  f4 ~~ v_f4*f4
  f5 ~~ v_f5*f5
  
  v_f1 > 0.01
  v_f2 > 0.01
  v_f3 > 0.01
  v_f4 > 0.01
  v_f5 > 0.01
  
  # Measurement 
  f3 ~ f1 + f2 + f4
  f4 ~ f1 + f2
  f5 ~ f3 + f4 + f1
  
  # lower bound residual variances:
  y1 ~~ v1*y1
  y2 ~~ v2*y2
  y3 ~~ v3*y3
  y4 ~~ v4*y4
  y5 ~~ v5*y5
  y6 ~~ v6*y6
  y7 ~~ v7*y7
  y8 ~~ v8*y8
  y9 ~~ v9*y9
  y10 ~~ v10*y10
  y11 ~~ v11*y11
  y12 ~~ v12*y12
  y13 ~~ v13*y13
  y14 ~~ v14*y14
  y15 ~~ v15*y15
  
  v1 > 0.01
  v2 > 0.01
  v3 > 0.01
  v4 > 0.01
  v5 > 0.01
  v6 > 0.01
  v7 > 0.01
  v8 > 0.01
  v9 > 0.01
  v10 > 0.01
  v11 > 0.01
  v12 > 0.01
  v13 > 0.01
  v14 > 0.01
  v15 > 0.01
"


# Study function
simulate_inner <- function(model_type, N, reliability, method) {
  safe_quiet_run_analysis <- safely(quietly(run_analysis))
  data <- gen_pop_model_data(model_type, N, reliability)$data
  
  # Choose the appropriate model syntax based on the method
  model_syntax <- if(method == "SEM") sem_model_syntax else sam_model_syntax
  
  fit_result <- safe_quiet_run_analysis(data, model_syntax, method)
  sanity_check_estimates <- run_sanity_check(model_type, model_syntax)
  
  warnings_detected <- fit_result$result$warnings
  improper_solution <- any(grepl("some estimated ov variances are negative", warnings_detected))
  
  if (!is.null(fit_result$result$result) && lavInspect(fit_result$result$result, "converged")) {
    PT <- parTable(fit_result$result$result)
    estimated_paths <- PT[PT$op == "~", "est"]
    names(estimated_paths) <- paste0(PT[PT$op == "~", "lhs"], "~", PT[PT$op == "~", "rhs"])
    
    # Calculate performance metrics
    coverage <- calculate_coverage(fit_result$result$result, true_values)
    relative_bias <- calculate_relative_bias(estimated_paths, true_values)
    relative_rmse <- calculate_relative_rmse(estimated_paths, true_values)
    relative_bias_se <- calculate_relative_bias_se(fit_result$result$result, model_type, model_syntax)  # Updated here
    
    return(list(
      Converged = 1, NonConverged = 0,
      EstimatedPaths = list(estimated_paths),
      SanityCheckEstimates = list(sanity_check_estimates),
      Coverage = coverage,
      RelativeBias = relative_bias,
      RelativeRMSE = relative_rmse,
      RelativeBiasSE = relative_bias_se,  # Added here
      RelativeBiasList = list(relative_bias),
      RelativeRMSEList = list(relative_rmse),
      RelativeBiasSEList = list(relative_bias_se),  # Added here
      ImproperSolution = improper_solution,
      Warnings = toString(fit_result$result$warnings),
      Messages = toString(fit_result$result$messages),
      Errors = if (is.null(fit_result$error)) NA_character_ else toString(fit_result$error$message)
    ))
  } else {
    return(list(
      Converged = 0, NonConverged = 1,
      EstimatedPaths = list(setNames(rep(NA, length(true_values$B)), names(true_values$B))),
      SanityCheckEstimates = list(sanity_check_estimates),
      Coverage = NA,
      RelativeBias = NA,
      RelativeRMSE = NA,
      RelativeBiasSE = NA,  # Added here
      RelativeBiasList = list(NA),
      RelativeRMSEList = list(NA),
      RelativeBiasSEList = list(NA),  # Added here
      ImproperSolution = improper_solution,
      Warnings = toString(fit_result$result$warnings),
      Messages = toString(fit_result$result$messages),
      Errors = if (is.null(fit_result$error)) NA_character_ else toString(fit_result$error$message)
    ))
  }
}

simulate_mid <- function(seed, design) {
  future_pmap(design, simulate_inner,
              .options = furrr_options(seed = rep.int(list(seed), nrow(design))))
}

simulate_outer <- function(chunk, chunk_params) {
  future_pmap(chunk_params, simulate_mid,
              .options = furrr_options(seed = NULL, scheduling = 1))
}

run_study_3 <- function(params, true_values) {
  # Run the simulations and analysis in parallel
  results_df <- params %>%
    mutate(chunk = row_number() %% 100) %>%
    nest(.by = c(chunk, seed), .key = 'design') %>%
    nest(.by = chunk, .key = 'chunk_params') %>%
    mutate(results = future_map2(chunk, chunk_params, simulate_outer, .options = furrr_options(seed=NULL))) %>%
    unnest(c(chunk_params, results)) %>%
    unnest(c(design, results)) %>%
    mutate(
      Converged = map_dbl(results, ~ .x$Converged),
      NonConverged = map_dbl(results, ~ .x$NonConverged),
      Warnings = map_chr(results, ~ .x$Warnings),
      Messages = map_chr(results, ~ .x$Messages),
      Errors = map_chr(results, ~ .x$Errors),
      EstimatedPaths = map(results, ~ .x$EstimatedPaths[[1]]),
      SanityCheck = map(results, ~ .x$SanityCheckEstimates[[1]]),
      Coverage = map_dbl(results, ~ .x$Coverage),
      RelativeBias = map_dbl(results, ~ .x$RelativeBias),
      RelativeRMSE = map_dbl(results, ~ .x$RelativeRMSE),
      RelativeBiasSE = map_dbl(results, ~ .x$RelativeBiasSE),  # Added here
      RelativeBiasList = map(results, ~ .x$RelativeBiasList[[1]]),
      RelativeRMSEList = map(results, ~ .x$RelativeRMSEList[[1]]),
      RelativeBiasSEList = map(results, ~ .x$RelativeBiasSEList[[1]]),  # Added here
      ImproperSolution = map_lgl(results, ~ .x$ImproperSolution)
    )
  
  print('Parallel computation done')
  
  # Remove NAs from the lists before calculating MCSE
  relative_bias_list <- unlist(results_df$RelativeBiasList)
  relative_rmse_list <- unlist(results_df$RelativeRMSEList)
  relative_bias_se_list <- unlist(results_df$RelativeBiasSEList)  # Added here
  
  relative_bias_list <- na.omit(relative_bias_list)
  relative_rmse_list <- na.omit(relative_rmse_list)
  relative_bias_se_list <- na.omit(relative_bias_se_list)  # Added here
  
  # Calculate summary statistics
  summary_stats <- results_df %>%
    group_by(model_type, N, reliability, method) %>%
    summarise(
      ConvergenceRate = mean(Converged),
      NonConvergenceCount = sum(NonConverged),
      n_converged = sum(Converged),
      MeanCoverage = mean(Coverage, na.rm = TRUE),
      MeanRelativeBias = mean(RelativeBias, na.rm = TRUE),
      MeanRelativeRMSE = mean(RelativeRMSE, na.rm = TRUE),
      MeanRelativeBiasSE = mean(RelativeBiasSE, na.rm = TRUE),  # Added here
      MCSE_RelativeBias = calculate_mcse_bias(relative_bias_list),
      MCSE_RelativeRMSE = calculate_mcse_rmse(relative_rmse_list),
      MCSE_RelativeBiasSE = calculate_mcse_bias(relative_bias_se_list),  # Added here
      ImproperSolutionsCount = sum(ImproperSolution, na.rm = TRUE),  # Update to count
      .groups = 'drop'
    ) %>%
    arrange(model_type, N, reliability, method)
  
  list(Summary = summary_stats, DetailedResults = results_df)
}


cat("Starting simulation study 3...\n")
execution_time <- system.time({
  simulation_results <- run_study_3(params, true_values)
})
cat("Simulation study completed. Saving results...\n")

# Print execution time
print(execution_time)

# Save with timestamp
timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
filename <- paste0("simulation_results_study3r", timestamp, ".rda")
save_results(simulation_results, filename)

cat("Results saved to:", file.path(results_dir, filename), "\n")
print(simulation_results$Summary)
