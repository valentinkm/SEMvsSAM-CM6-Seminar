#run_analysis.R

# Load necessary libraries
library(lavaan)
library(dplyr)
library(purrr)

# Function to run analysis using SEM or SAM
run_analysis <- function(data, model_syntax, method = "SEM", b = 5) {
  if (method == "SEM") {
    fit <- sem(model_syntax, data = as.data.frame(data))
  } else if (method == "gSAM") {
    fit <- sam(model_syntax, data = as.data.frame(data), sam.method = "global")
  } else if (method %in% c("lSAM_ML", "lSAM_ULS")) {
    mm.list <- NULL
    if (b == 5) {
      mm.list <- list(
        "f1" = "f1",
        "f2" = "f2",
        "f3" = "f3",
        "f4" = "f4",
        "f5" = "f5"
      )
    } else if (b == 3) {
      mm.list <- list(
        "exo" = c("f1", "f2"),
        "endo" = c("f3", "f4", "f5")
      )
    } else {
      stop("Invalid number of measurement blocks specified")
    }
    
    struc_args <- list(estimator = ifelse(method == "lSAM_ML", "ML", "ULS"))
    
    fit <- sam(
      model_syntax, 
      data = as.data.frame(data), 
      sam.method = "local", 
      mm.list = mm.list, 
      struc.args = struc_args
    )
  } else {
    stop("Unknown method specified")
  }
  return(fit)
}



# Function to fit model to population matrix
run_sanity_check <- function(model_type, model_syntax) {
  popmodel <- gen_pop_model_syntax(gen_mat(model_type = model_type))
  fit0 <- sem(model = popmodel, do.fit = FALSE)
  COV <- lavInspect(fit0, what = "implied")$cov[,]
  sanity_check_fit <- sem(sample.cov = COV, model = model_syntax, sample.nobs = 10^6)
  sanity_check_estimates <- coef(sanity_check_fit)
  names(sanity_check_estimates) <- paste0(names(sanity_check_estimates), "_pop")
  return(sanity_check_estimates)
}

# Function to check the sanity check results
check_sanity <- function(sanity_check_estimates, true_values) {
  comparison <- compare_sanity_check(sanity_check_estimates, true_values)
  max_difference <- max(comparison$Differences$Difference, na.rm = TRUE)
  return(list(
    MaxDifference = max_difference,
    Alarm = comparison$Alarm
  ))
}

# Function to compare sanity check estimates with true values
compare_sanity_check <- function(sanity_check_estimates, true_values, threshold = 0.1) {
  true_values_flat <- unlist(true_values)
  
  # Align true values with the sanity check estimates
  common_params <- intersect(names(sanity_check_estimates), names(true_values_flat))
  aligned_true_values <- true_values_flat[common_params]
  aligned_sanity_check_estimates <- sanity_check_estimates[common_params]
  
  differences <- abs(aligned_sanity_check_estimates - aligned_true_values)
  differences_df <- data.frame(
    Parameter = names(differences),
    TrueValue = aligned_true_values,
    SanityCheckEstimate = aligned_sanity_check_estimates,
    Difference = differences
  )
  
  # Check if any difference exceeds the threshold
  alarm <- any(differences > threshold)
  
  return(list(
    Differences = differences_df,
    Alarm = alarm
  ))
}
