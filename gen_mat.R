# gen_mat.R

# Load necessary libraries
library(lavaan)
library(Matrix)

# Generate a set of matrices that serve as basis for the population model
gen_mat <- function(model_type, nfactors = 5, nvar.factor = 3, lambda = 0.70, 
                    beta_value = 0.1, psi.cor = 0.3, reliability = 0.80, 
                    rho = 0.80, R_squared = 0.1) {
  
  if (model_type %in% c("2.1", "2.2_exo", "2.2_endo", "2.2_both")) {
    # For study 2 models
    # Lambda matrix
    LAMBDA <- matrix(0, nrow = 15, ncol = 5)
    LAMBDA[1:3, 1] <- c(1, 0.9, 0.8)
    LAMBDA[4:6, 2] <- c(1, 0.9, 0.8)
    LAMBDA[7:9, 3] <- c(1, 0.9, 0.8)
    LAMBDA[10:12, 4] <- c(1, 0.9, 0.8)
    LAMBDA[13:15, 5] <- c(1, 0.9, 0.8)
    
    # Adjust Lambda for misspecifications
    if (model_type %in% c("2.2_exo", "2.2_both")) {
      LAMBDA[6, 3] <- lambda
    }
    if (model_type %in% c("2.2_endo", "2.2_both")) {
      LAMBDA[12, 5] <- lambda
    }
    
    # Beta matrix
    BETA <- matrix(0, nrow = nfactors, ncol = nfactors)
    
    # Number of predictors for each endogenous variable
    predictors_endogenous_4 <- 2  # for eta_4
    predictors_endogenous_5 <- 3  # for eta_5
    
    if (R_squared == 0.1) {
      R_squared_4 <- 0.1
      R_squared_5 <- 0.1
    } else if (R_squared == 0.4) {
      R_squared_4 <- 0.4
      R_squared_5 <- 0.4
    }
    
    beta_value_4 <- sqrt(R_squared_4 / predictors_endogenous_4)
    beta_value_5 <- sqrt(R_squared_5 / predictors_endogenous_5)
    
    BETA[4, 1] <- BETA[4, 2] <- beta_value_4
    BETA[5, 2] <- BETA[5, 3] <- BETA[5, 4] <- beta_value_5
    
    # Psi matrix
    PSI <- diag(nfactors)
    
    # Theta matrix
    THETA <- diag(15)
    diag(THETA) <- apply(LAMBDA, 1, function(row) {
      sum(row^2) / reliability - sum(row^2)
    })
    
    if (model_type %in% c("2.2_exo", "2.2_both")) {
      THETA[8, 9] <- THETA[9, 8] <- 0.6 * min(diag(THETA)[c(8, 9)])
    }
    if (model_type %in% c("2.2_endo", "2.2_both")) {
      THETA[14, 15] <- THETA[15, 14] <- 0.6 * min(diag(THETA)[c(14, 15)])
    }
  } else {
    # For study 1 and 3 models
    # 1. LAMBDA
    fac <- matrix(c(1, rep(lambda, times = (nvar.factor - 1L))), nvar.factor, 1L)
    LAMBDA <- bdiag(rep(list(fac), nfactors))
    
    if (model_type %in% c("1.2", "3.1")) {
      # misspecification in the measurement part: cross-loadings
      i.cross <- (0:(nfactors-1)) * nvar.factor + ceiling(nvar.factor / 2)
      for (j in 1:ncol(LAMBDA)) {
        LAMBDA[i.cross[j], c(2:nfactors, 1)[j]] <- rho * lambda
      }
    } else if (model_type == "3.1_negative") {
      # model with negative cross-loadings
      i.cross <- (0:(nfactors-1)) * nvar.factor + ceiling(nvar.factor / 2)
      for (j in 1:ncol(LAMBDA)) {
        LAMBDA[i.cross[j], c(2:nfactors, 1)[j]] <- -rho * lambda
      }
    }
    
    LAMBDA <- as.matrix(LAMBDA)
    
    # 2. BETA
    BETA <- matrix(0, nrow = nfactors, ncol = nfactors)
    BETA[3:5, 1] <- BETA[3:4, 2] <- BETA[5, 3] <- BETA[3, 4] <- beta_value
    
    # Add the new paths as per your requirements
    if (model_type %in% c("1.1", "1.2", "1.3", "1.4", "3.1", "3.2", "3.1_negative", "3.2_negative")) {
      BETA[5, 4] <- beta_value  # Additional path from f4 to f5
    }
    if (model_type == "1.3") {
      BETA[4, 3] <- beta_value  # Path from f3 to f4 instead of f4 to f3
      BETA[3, 4] <- 0           # Ensure the path f4 to f3 is removed
    }
    if (model_type == "1.4") {
      BETA[4, 3] <- beta_value  # Path from f3 to f4
      BETA[3, 4] <- beta_value  # Path from f4 to f3
    }
    
    BETA.model <- BETA  # Structural part to be fitted
    VAL <- BETA[BETA != 0]  # true values
    
    # 3. PSI
    PSI <- matrix(0, nrow = nfactors, ncol = nfactors)
    PSI[1, 1] <- PSI[2, 2] <- 1  # the exogenous latent variables
    RES <- (1 - beta_value^2)
    PSI[lav_matrix_diag_idx(nfactors)[-c(1:2)]] <- RES
    
    # Calculate Sigma_eta
    IB_inv <- solve(diag(nfactors) - BETA)
    Sigma_eta <- IB_inv %*% PSI %*% t(IB_inv)
    
    # 4. THETA
    tmp <- diag(LAMBDA %*% Sigma_eta %*% t(LAMBDA))
    theta_diag <- tmp / reliability - tmp
    stopifnot(all(theta_diag > 0))
    THETA <- matrix(0, nrow(LAMBDA), nrow(LAMBDA))
    diag(THETA) <- theta_diag
    
    if (model_type %in% c("1.3", "3.2", "3.2_negative")) {
      pairs <- cbind(c(2, 5, 8, 11, 14), c(3, 6, 9, 12, 15))
      for (pair in 1:nrow(pairs)) {
        i <- pairs[pair, 1]
        j <- pairs[pair, 2]
        THETA[i, j] <- THETA[j, i] <- -0.6 * min(theta_diag[c(i, j)])
      }
    }
  }
  
  MLIST <- list(lambda = LAMBDA, theta = THETA, psi = PSI, beta = BETA)
  
  return(MLIST)
}

# Utility function to get diagonal indices
lav_matrix_diag_idx <- function(n) {
  return(seq(1, n^2, by = n + 1))
}

# # Test the function with models 1.1 to 1.4, 3.1, 3.2, and return the matrices
# test_models_study1 <- function() {
#   models <- c("1.1", "1.2", "1.3", "1.4", "3.1", "3.2", "3.1_negative", "3.2_negative")
#   for (model in models) {
#     cat("Testing model:", model, "\n")
#     MLIST <- gen_mat(model, nfactors = 5, nvar.factor = 3, lambda = 0.70,
#                      beta_value = 0.1, psi.cor = 0.3, reliability = 0.80,
#                      rho = 0.80, R_squared = 0.1)
#     cat("Lambda matrix:\n")
#     print(MLIST$lambda)
#     cat("Theta matrix:\n")
#     print(MLIST$theta)
#     cat("Psi matrix:\n")
#     print(MLIST$psi)
#     cat("Beta matrix:\n")
#     print(MLIST$beta)
#     cat("\n")
#   }
# }
# 
# test_models_study1()

# Uncomment the following to test study 2 models
# test_all_models_2 <- function() {
#   # models <- c("2.1", "2.2_exo", "2.2_endo", "2.2_both")
#   for (model in models) {
#     cat("Testing model:", model, "\n")
#     MLIST <- gen_mat(model, nfactors = 5, nvar.factor = 3, lambda = 0.70,
#                      beta_value = 0.1, psi.cor = 0.3, reliability = 0.70,
#                      rho = 0.80, R_squared = 0.4)
#     cat("Beta matrix:\n")
#     print(MLIST$beta)
#     cat("\n")
#   }
# }
# # 
# test_all_models_2()
