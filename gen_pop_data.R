# gen_pop_data.R

# Load necessary libraries
library(lavaan)
library(Matrix)
source("gen_mat.R")
source("gen_pop_lavsyntax.R")

# Data generation function
gen_pop_model_data <- function(model_type, N, reliability, R_squared = 0.1) {
  # Generate matrices and syntax
  matrices <- gen_mat(model_type, reliability = reliability, R_squared = R_squared)
  syntax <- gen_pop_model_syntax(matrices)
  
  # Generate data using lavaan's simulateData function
  data <- simulateData(model = syntax, sample.nobs = N)
  
  return(list(data = data, syntax = syntax, matrices = matrices))
}
