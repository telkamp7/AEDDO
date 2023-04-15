# This script contains the negative log-likelihood (nll) functions for the Compound Poisson
# Gamma model considered in this master thesis project.

# The nll using the case totals
nll <- function(theta, x){
  r <- 1/theta[2]
  p <- 1/(theta[1]*theta[2]+1)
  -sum(dnbinom(x = x, size = r, prob = p, log = TRUE))
}

# The nll using age specific phi parameters with 11 agegroups
nll.age11 <- function(theta, data, formula){
  # Extract counts
  y <- data$y
  # Construct the design matrix
  designMatrix <- model.matrix(object = formula, data = data)
  # Construct regression parameters
  beta <- theta[1:ncol(designMatrix)]
  #... and model parameters
  phi <- theta[-(1:ncol(designMatrix))]
  # # Extract agegroups
  ageGroup <- data$ageGroup
  
  # Define parameters
  lambda <- exp(designMatrix %*% beta - log(data$n))
  
  # Construct the size and probability for the negative binomial distribution
  r <- 1/phi
  p <- 1/(lambda*phi+1)
  
  # Initilize the log-likelihood
  ll <- 0
  for(i in 1:nrow(data)){
    ll = ll + dnbinom(x = y[i], size = r[ageGroup[i]], prob = p[i], log = TRUE) 
  }
  
  # Return the negative log-likelihood
  -ll
}

# The nll using the same phi for all agegroups
nll.age6 <- function(theta, data, formula){
  # Extract counts
  y <- data$y
  # Construct the design matrix
  designMatrix <- model.matrix(object = formula, data = data)
  # Construct regression parameters
  beta <- theta[1:ncol(designMatrix)]
  #... and model parameters
  phi <- theta[-(1:ncol(designMatrix))]
  # # Extract agegroups
  ageGroup <- data$ageGroup
  
  # Define parameters
  lambda <- exp(designMatrix %*% beta - log(data$n))
  
  # Construct the size and probability for the negative binomial distribution
  r <- 1/phi
  p <- 1/(lambda*phi+1)
  
  # Initilize the log-likelihood
  ll <- 0
  for(i in 1:nrow(data)){
    ll = ll + dnbinom(x = y[i], size = r, prob = p[i], log = TRUE) 
  }
  
  # Return the negative log-likelihood
  -ll
}

nll.model <- function(theta, data, designMatrix){
  # Extract counts
  y <- data$y
  # Construct regression parameters
  beta <- theta[1:ncol(designMatrix)]
  #... and model parameters
  phi <- theta[-(1:ncol(designMatrix))]
  # Define parameters
  lambda <- exp(designMatrix %*% beta - log(data$n))
  # Construct the size and probability for the negative binomial distribution
  r <- 1/phi
  p <- 1/(lambda*phi+1)
  # Return the negative log-likelihood
  -sum( dnbinom(x = y, size = r, prob = p, log = TRUE) )
}

nll.age6.window <- function(theta, data, designMatrix){
  # Extract counts
  y <- data$y
  # # Construct the design matrix
  # designMatrix <- model.matrix(object = formula, data = data)
  # Construct regression parameters
  beta <- theta[1:ncol(designMatrix)]
  #... and model parameters
  phi <- theta[-(1:ncol(designMatrix))]
  # # Extract agegroups
  ageGroup <- data$ageGroup
  
  # Define parameters
  lambda <- exp(designMatrix %*% beta - log(data$n))
  
  # Construct the size and probability for the negative binomial distribution
  r <- 1/phi
  p <- 1/(lambda*phi+1)
  
  # Initilize the log-likelihood
  # ll <- 0
  # for(i in 1:nrow(data)){
  #   ll = ll + dnbinom(x = y[i], size = r, prob = p[i], log = TRUE) 
  # }
  # -ll
  
  # Return the negative log-likelihood
  -sum( dnbinom(x = y, size = r, prob = p, log = TRUE) )
}
