# Load necessary libraries
library(dplyr)
library(lavaan)
library(purrr)


# Function to generate the true standard errors using the population covariance matrix
get_true_se_for_model <- function(model_type, model_syntax) {
  # Generate the population model syntax and implied covariance matrix
  popmodel <- gen_pop_model_syntax(gen_mat(model_type = model_type))
  fit0 <- sem(model = popmodel, do.fit = FALSE)
  COV <- lavInspect(fit0, what = "implied")$cov[, ]
  
  # Fit the model to this implied covariance matrix
  sanity_check_fit <- sem(sample.cov = COV, model = model_syntax, sample.nobs = 10^6)
  
  # Extract the standard errors
  param_estimates <- parameterEstimates(sanity_check_fit) %>%
    filter(op %in% c("~", "=~")) %>%
    select(lhs, rhs, op, se)
  
  # Create a named vector of SEs, where names are in the format 'lhs ~ rhs' or 'lhs =~ rhs'
  true_se <- param_estimates %>%
    mutate(param = paste(lhs, op, rhs)) %>%
    select(param, se) %>%
    pull(se, param)
  
  return(true_se)
}

# Function to calculate relative bias in SE
calculate_relative_bias_se <- function(fit, model_type, model_syntax) {
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(NA)
  }
  
  # Extract parameter estimates and their standard errors from the fitted model
  param_estimates <- parameterEstimates(fit) %>%
    filter(op == "~")
  
  # Generate the true SEs using the renamed function
  true_se <- get_true_se_for_model(model_type, model_syntax)
  
  # Ensure we only calculate relative bias for regressors in the model
  common_params <- intersect(paste(param_estimates$lhs, "~", param_estimates$rhs), names(true_se))
  aligned_true_se <- true_se[common_params]
  
  relative_bias_se <- map_dbl(common_params, function(param) {
    row <- param_estimates %>%
      filter(lhs == unlist(strsplit(param, " ~ "))[1], rhs == unlist(strsplit(param, " ~ "))[2])
    
    if (nrow(row) > 0 && !is.na(aligned_true_se[param])) {
      (row$se - aligned_true_se[param]) / aligned_true_se[param]
    } else {
      NA
    }
  })
  
  # Return mean relative bias in SE across all parameters
  mean(relative_bias_se, na.rm = TRUE)
}


# Function to calculate coverage of confidence intervals
calculate_coverage <- function(fit, true_values) {
  if (is.null(fit) || !lavInspect(fit, "converged")) {
    return(NA)
  }
  
  # Extract parameter estimates with confidence intervals
  param_estimates <- parameterEstimates(fit)
  
  # Filter for regression paths (op == "~")
  param_estimates <- param_estimates %>%
    filter(op == "~")
  
  # Calculate coverage for each parameter in true_values$B
  coverage <- map_dbl(names(true_values$B), function(param) {
    parts <- unlist(strsplit(param, "~"))
    row <- param_estimates %>%
      filter(lhs == parts[1], rhs == parts[2])
    if (nrow(row) > 0) {
      row$ci.lower <= true_values$B[param] && row$ci.upper >= true_values$B[param]
    } else {
      NA
    }
  })
  
  # Return mean coverage across all parameters
  mean(coverage, na.rm = TRUE)
}

# Function to calculate empirical relative bias
calculate_relative_bias <- function(estimated_paths, true_values) {
  if (all(is.na(estimated_paths))) {
    return(NA)
  }
  
  # Align true values with estimated paths
  common_params <- intersect(names(estimated_paths), names(true_values$B))
  aligned_true_values <- true_values$B[common_params]
  aligned_estimated_paths <- estimated_paths[common_params]
  
  bias <- (aligned_estimated_paths - aligned_true_values) / aligned_true_values
  mean(bias, na.rm = TRUE)
}

# Function to calculate empirical relative RMSE
calculate_relative_rmse <- function(estimated_paths, true_values) {
  if (all(is.na(estimated_paths))) {
    return(NA)
  }
  
  # Align true values with estimated paths
  common_params <- intersect(names(estimated_paths), names(true_values$B))
  aligned_true_values <- true_values$B[common_params]
  aligned_estimated_paths <- estimated_paths[common_params]
  
  rmse <- sqrt(mean((aligned_estimated_paths - aligned_true_values)^2, na.rm = TRUE)) / mean(aligned_true_values)
  rmse
}

# Function to calculate Monte Carlo Standard Error for Bias
calculate_mcse_bias <- function(bias_list) {
  bias_list <- na.omit(bias_list)  # Remove NAs
  if (length(bias_list) == 0) {
    return(NA)
  }
  K <- length(bias_list)
  S_T2 <- var(bias_list)  # This is already S_T^2
  sqrt(S_T2 / K)
}

# Function to calculate Monte Carlo Standard Error for RMSE
calculate_mcse_rmse <- function(rmse_list) {
  rmse_list <- na.omit(rmse_list)  # Remove NAs
  if (length(rmse_list) == 0) {
    return(NA)
  }
  K <- length(rmse_list)
  rmse_mean <- mean(rmse_list)
  sqrt(sum((rmse_list - rmse_mean)^2) / (K * (K - 1)))
}