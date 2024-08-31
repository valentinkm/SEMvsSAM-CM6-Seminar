library(dplyr)
library(stringr)
library(purrr)
source("gen_mat.R")


standardize_lavaan_warning <- function(warning) {
  if (is.na(warning) || warning == "") return(NA_character_)
  
  warnings <- str_split(warning, ", ")[[1]]
  
  standardized_warnings <- sapply(warnings, function(w) {
    if (str_detect(w, "The variance-covariance matrix of the estimated parameters \\(vcov\\)")) {
      w <- str_replace(w, "\\(= .+\\)", "(= [VALUE])")
      w <- str_replace(w, "is (close to|smaller than) zero", "is [RELATION] zero")
      return("lavaan WARNING: The variance-covariance matrix of the estimated parameters (vcov) does not appear to be positive definite! The smallest eigenvalue is smaller than or close to zero. This may be a symptom that the model is not identified.")
    } else if (str_detect(w, "some estimated lv variances are negative")) {
      return("lavaan WARNING: some estimated lv variances are negative")
    } else {
      return(w)
    }
  })
  
  paste(unique(standardized_warnings), collapse = ", ")
}

process_study_warnings <- function(detailed_results, study_number) {
  warnings_summary <- detailed_results %>%
    mutate(
      Warnings = map_chr(Warnings, standardize_lavaan_warning),
      Errors = map_chr(Errors, standardize_lavaan_warning),
      MessageType = case_when(
        !is.na(Warnings) ~ "Warning",
        !is.na(Errors) ~ "Error",
        TRUE ~ NA_character_
      ),
      Message = coalesce(Warnings, Errors)
    ) %>%
    filter(!is.na(MessageType)) %>%
    mutate(
      method = case_when(
        str_detect(method, "lSAM_ULS") ~ "lSAM-ULS",
        TRUE ~ method
      ),
      Study = paste("Study", study_number)
    ) %>%
    group_by(Study, model_type, N, reliability, method, seed, MessageType, Message) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(count = pmin(count, 1)) %>%
    group_by(Study, model_type, N, reliability, method, MessageType, Message) %>%
    summarise(count = sum(count), .groups = 'drop') %>%
    arrange(desc(count))
  
  return(warnings_summary)
}

filter_and_summarize_study3 <- function(detailed_results) {
  # Adjust improper solutions and track estimation general issues
  detailed_results <- detailed_results %>%
    mutate(ImproperSolution = case_when(
      grepl("variances are negative", Warnings, fixed = TRUE) ~ 1,
      TRUE ~ as.numeric(ImproperSolution)
    ),
    OtherEstimationIssue = case_when(
      is.na(Warnings) | Warnings == "" ~ 0,
      grepl("variances are negative", Warnings, fixed = TRUE) ~ 0,
      TRUE ~ 1
    ))
  
  # Define the grouping variables specific to Study 3
  grouping_vars <- c("model_type", "N", "reliability", "method")
  
  # Create a summary before filtering
  summary <- detailed_results %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      Total = n(),
      Converged = sum(Converged == 1),
      NotConverged = sum(Converged == 0),
      ProperSolutions = sum(ImproperSolution == 0),
      ImproperSolutions = sum(ImproperSolution == 1),
      OtherEstimationIssues = sum(OtherEstimationIssue == 1),
      .groups = "drop"
    )
  
  # Filter out improper and not converged rows
  filtered_results <- detailed_results %>%
    filter(
      ImproperSolution != 1,
      Converged != 0,
      OtherEstimationIssue != 1
    )
  
  # Save the summary
  saveRDS(summary, file = "results/convergence_rate_study3.rds")
  
  return(filtered_results)
}



get_true_values <- function(study, model_type, for_aggregation = FALSE) {
  MLIST <- gen_mat(model_type)
  BETA <- MLIST$beta
  
  required_params <- c("f3~f1", "f3~f2", "f3~f4", "f4~f1", "f4~f2", "f5~f3", "f5~f4", "f5~f1")
  
  true_values <- numeric(length(required_params))
  names(true_values) <- required_params
  
  for (param in required_params) {
    parts <- strsplit(param, "~")[[1]]
    i <- as.numeric(substr(parts[1], 2, 2))
    j <- as.numeric(substr(parts[2], 2, 2))
    true_values[param] <- BETA[i, j]
  }
  
  return(true_values)
}

process_row <- function(row) {
  true_values <- get_true_values(3, row$model_type)
  estimated_paths <- row$EstimatedPaths
  
  tibble(
    model_type = row$model_type,
    N = row$N,
    reliability = row$reliability,
    method = row$method,
    parameter = names(estimated_paths),
    true_value = true_values[names(estimated_paths)],
    estimated_value = unlist(estimated_paths),
    is_misspecified = FALSE,
    Coverage = row$Coverage
    # RelativeBiasSE = row$RelativeBiasSE  # Include the relative SE bias in the processed row
  )
}


calculate_metrics <- function(group, study_number) {
  K <- nrow(group)
  is_misspecified <- first(group$is_misspecified)
  true_value <- first(group$true_value)
  
  results <- group %>%
    summarise(
      true_value = first(true_value),
      mean_estimate = if(n() > 0) mean(estimated_value, na.rm = TRUE) else NA_real_,
      var_estimate = if(n() > 0) var(estimated_value, na.rm = TRUE) else NA_real_,
      n = n()
    )
  
  # Absolute metrics (for all studies and parameters)
  results <- results %>%
    mutate(
      bias = mean_estimate - true_value,
      mse = mean((group$estimated_value - true_value)^2),
      rmse = sqrt(mse),
      mcse_bias = sqrt(var_estimate / K),
      mcse_rmse = sqrt(var(sqrt((group$estimated_value - true_value)^2)) * (K - 1) / K)
    )
  
  # Conditional calculation of relative metrics
  if (abs(true_value) > 1e-10) {  # Avoid division by zero or very small numbers
    results <- results %>%
      mutate(
        relative_bias = bias / true_value,
        relative_rmse = sqrt((bias^2 + var_estimate) / true_value^2),
        mcse_relative_bias = sqrt(var_estimate / (K * true_value^2)),
        mcse_relative_rmse = sqrt(var(sqrt(((group$estimated_value - true_value) / true_value)^2)) * (K - 1) / K)
      )
  } else {
    results <- results %>%
      mutate(
        relative_bias = NA_real_,
        relative_rmse = NA_real_,
        mcse_relative_bias = NA_real_,
        mcse_relative_rmse = NA_real_
      )
  }
  
  results
}

process_study_parameterwise <- function(study_data) {
  processed_data <- study_data %>%
    rowwise() %>%
    do(process_row(.)) %>%
    ungroup()
  
  grouped_data <- processed_data %>%
    filter(!is.na(parameter)) %>% 
    group_by(model_type, N, reliability, method, parameter)
  
  results <- grouped_data %>%
    summarise(
      n = n(),
      TrueValue = first(true_value),
      SumEstimate = sum(estimated_value, na.rm = TRUE),
      SumSquaredEstimate = sum(estimated_value^2, na.rm = TRUE),
      # SumRelativeBiasSE = sum(RelativeBiasSE, na.rm = TRUE),  # Sum of RelativeBiasSE
      MeanCoverage = mean(Coverage, na.rm = TRUE),
      .groups = 'drop'
    )
  
  final_results <- results %>%
    mutate(
      MeanEstimate = SumEstimate / n,
      VarEstimate = (SumSquaredEstimate - (SumEstimate^2 / n)) / (n - 1),
      Bias = MeanEstimate - TrueValue,
      RMSE = sqrt(VarEstimate + Bias^2),
      MCSE_Bias = sqrt(VarEstimate / n),
      MCSE_RMSE = sqrt(VarEstimate / (2 * n)),
      RelativeBias = if_else(abs(TrueValue) > 1e-10, Bias / TrueValue, NA_real_),
      RelativeRMSE = if_else(abs(TrueValue) > 1e-10, RMSE / abs(TrueValue), NA_real_),
      MCSE_RelativeBias = if_else(abs(TrueValue) > 1e-10, MCSE_Bias / abs(TrueValue), NA_real_),
      MCSE_RelativeRMSE = if_else(abs(TrueValue) > 1e-10, MCSE_RMSE / abs(TrueValue), NA_real_),
      Coverage = MeanCoverage
      # MeanRelativeBiasSE = SumRelativeBiasSE / n  # Calculate mean relative SE bias
    ) %>%
    select(-SumEstimate, -SumSquaredEstimate, -MeanCoverage)
  
  return(final_results)
}


aggregate_results <- function(paramwise_results) {
  grouping_vars <- c("model_type", "N", "reliability", "method")
  exclude_from_mean <- c("n", "total_cases", "included_cases", "TrueValue")
  numeric_cols <- names(paramwise_results)[sapply(paramwise_results, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c(grouping_vars, exclude_from_mean))
  
  summary_exprs <- list()
  for (col in exclude_from_mean) {
    if (col %in% names(paramwise_results)) {
      summary_exprs[[col]] <- expr(first(!!sym(col)))
    }
  }
  for (col in numeric_cols) {
    if (col %in% c("Bias", "RelativeBias")) {
      # Use mean of absolute values for Bias and RelativeBias
      summary_exprs[[paste0("mean_", col)]] <- expr(mean(abs(!!sym(col)), na.rm = TRUE))
    } else {
      summary_exprs[[paste0("mean_", col)]] <- expr(mean(!!sym(col), na.rm = TRUE))
    }
  }
  summary_exprs$n_parameters <- expr(n())
  
  aggregated_results <- paramwise_results %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(!!!summary_exprs, .groups = "drop") %>%
    rename(average_coverage = mean_Coverage)
  
  return(aggregated_results)
}


process_study3 <- function() {
  print("Reading results file")
  study_3 <- readRDS("results/simulation_results_study3.rda")
  detailed_results_3 <- study_3$DetailedResults
  rm(study_3)
  gc()
  
  print("Processing warnings and errors")
  summary_warnings_3 <- process_study_warnings(detailed_results_3, 3)
  
  print("Filter and summarize the results for Study 3")
  detailed_results_3 <- filter_and_summarize_study3(detailed_results_3)
  
  print("Calculate parameter-wise metrics")
  study3_paramwise <- process_study_parameterwise(detailed_results_3)
  
  print("Calculate metrics aggregated across parameters")
  study3_aggregated <- aggregate_results(study3_paramwise)
  
  saveRDS(study3_paramwise, file = "results/parameter_wise_summary_study3.rds")
  saveRDS(study3_aggregated, file = "results/aggregated_summary_study3.rds")
  cat("Study 3 results saved to: \n
      results/parameter_wise_summary_study3.rds \n
      results/aggregated_summary_study3.rds")
  
  rm(detailed_results_3, study3_paramwise, study3_aggregated)
  gc()
  
  print("Study 3 processing completed.")
  
  return(summary_warnings_3)
}

# Run the processing for Study 3
process_study3()
