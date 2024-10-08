# gen_pop_lavsyntax.R
library(lavaan)
library(Matrix)

# Function to round numbers to 3 digits
round_three <- function(x) {
  return(format(round(x, 3), nsmall = 3))
}

# Function to generate lavaan model syntax from model matrices
gen_pop_model_syntax <- function(MLIST, ov.prefix = "y", lv.prefix = "f", include.values = TRUE) {
  
  LAMBDA <- MLIST$lambda
  THETA  <- MLIST$theta
  PSI    <- MLIST$psi
  BETA   <- MLIST$beta
  
  if (ov.prefix == lv.prefix) {
    stop("lavaan ERROR: ov.prefix can not be the same as lv.prefix")
  }
  
  header <- "# syntax generated by gen_pop_model_syntax()"
  
  # LAMBDA
  if (!is.null(LAMBDA)) {
    IDXV <- row(LAMBDA)[(LAMBDA != 0)]
    IDXF <- col(LAMBDA)[(LAMBDA != 0)]
    
    unique_factors <- unique(IDXF)
    IDXV <- as.integer(unlist(sapply(unique_factors, function(j) {
      ji <- IDXV[which(IDXF == j)]
      j1 <- which(abs(LAMBDA[ji, j] - 1) < .Machine$double.eps)
      if (length(j1) > 0) {
        ji[c(1, j1)] <- ji[c(j1, 1)]
      }
      return(ji)
    })))
    
    IDXF <- rep(unique_factors, times = sapply(unique_factors, function(j) length(which(IDXF == j))))
    
    nel <- length(IDXF)
    lambda.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- LAMBDA[IDXV[i], IDXF[i]]
      if (include.values) {
        lambda.txt[i] <- sprintf("%s%d =~ %s*%s%d", lv.prefix, IDXF[i], round_three(value), ov.prefix, IDXV[i])
      } else {
        lambda.txt[i] <- sprintf("%s%d =~ %s%d", lv.prefix, IDXF[i], ov.prefix, IDXV[i])
      }
    }
  } else {
    lambda.txt <- character(0L)
  }
  
  # THETA
  if (!is.null(THETA)) {
    IDX1 <- row(THETA)[(THETA != 0) & upper.tri(THETA, diag = TRUE)]
    IDX2 <- col(THETA)[(THETA != 0) & upper.tri(THETA, diag = TRUE)]
    nel <- length(IDX1)
    theta.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- THETA[IDX1[i], IDX2[i]]
      if (include.values) {
        theta.txt[i] <- sprintf("%s%d ~~ %s*%s%d", ov.prefix, IDX1[i], round_three(value), ov.prefix, IDX2[i])
      } else {
        theta.txt[i] <- sprintf("%s%d ~~ %s%d", ov.prefix, IDX1[i], ov.prefix, IDX2[i])
      }
    }
  } else {
    theta.txt <- character(0L)
  }
  
  # PSI
  if (!is.null(PSI)) {
    IDX1 <- row(PSI)[(PSI != 0) & upper.tri(PSI, diag = TRUE)]
    IDX2 <- col(PSI)[(PSI != 0) & upper.tri(PSI, diag = TRUE)]
    nel <- length(IDX1)
    psi.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- PSI[IDX1[i], IDX2[i]]
      if (include.values) {
        psi.txt[i] <- sprintf("%s%d ~~ %s*%s%d", lv.prefix, IDX1[i], round_three(value), lv.prefix, IDX2[i])
      } else {
        psi.txt[i] <- sprintf("%s%d ~~ %s%d", lv.prefix, IDX1[i], lv.prefix, IDX2[i])
      }
    }
  } else {
    psi.txt <- character(0L)
  }
  
  # BETA
  if (!is.null(BETA)) {
    IDX1 <- row(BETA)[(BETA != 0)]
    IDX2 <- col(BETA)[(BETA != 0)]
    nel <- length(IDX1)
    beta.txt <- character(nel)
    for (i in seq_len(nel)) {
      value <- BETA[IDX1[i], IDX2[i]]
      if (include.values) {
        beta.txt[i] <- sprintf("%s%d ~ %s*%s%d", lv.prefix, IDX1[i], round_three(value), lv.prefix, IDX2[i])
      } else {
        beta.txt[i] <- sprintf("%s%d ~ %s%d", lv.prefix, IDX1[i], lv.prefix, IDX2[i])
      }
    }
  } else {
    beta.txt <- character(0L)
  }
  
  # Assemble
  syntax <- paste(c(header, lambda.txt, theta.txt, psi.txt, beta.txt, ""),
                  collapse = "\n")
  
  return(syntax)
}

# Test the function with different models
# test_models <- function() {
#     models <- c("1.1", "1.2", "1.3", "1.4")
#   # models <- c(2.1", "2.2_exo", "2.2_endo", "2.2_both")
#   # models <- c("3.1","3.2", "3.1_negative", "3.2_negative")
#   for (model in models) {
#     cat("Testing model:", model, "\n")
#     MLIST <- gen_mat(model, nfactors = 5, nvar.factor = 3, lambda = 0.70,
#                      beta_value = 0.1, psi.cor = 0.3, reliability = 0.7,
#                      rho = 0.80, R_squared = 0.1)
#     syntax <- gen_pop_model_syntax(MLIST)
#     cat(syntax, "\n\n")
#   }
# }
# test_models()

# # # get model 2.1 syntax for low and medium R-squared:
# MLIST <- gen_mat("2.1", nfactors = 5, nvar.factor = 3, lambda = 0.70,
#                  beta_value = 0.1, psi.cor = 0.3, reliability = 0.80,
#                  rho = 0.80, R_squared = 0.1)
# syntax <- gen_pop_model_syntax(MLIST)
# cat(syntax, "\n\n")
# 
# MLIST <- gen_mat("2.1", nfactors = 5, nvar.factor = 3, lambda = 0.70,
#                  beta_value = 0.1, psi.cor = 0.3, reliability = 0.80,
#                  rho = 0.80, R_squared = 0.4)
# syntax <- gen_pop_model_syntax(MLIST)
# cat(syntax, "\n\n")