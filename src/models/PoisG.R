
# Import libraries
library(readr)
library(dplyr)

# Load in the processed data
dat <- read_rds(file = "../../data/processed/dat.rds")

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date) %>%
  reframe(y = sum(cases), n = sum(n))

?dnbinom

# Construct the negative log-likelihood function
nll <- function(theta, x){
  r <- 1/theta[2]
  p <- 1/(theta[1]*theta[2]+1)
  -sum(dnbinom(x = x, size = r, prob = p, log = TRUE))
}

# Optimize the parameters
opt <- nlminb(start = c(1,1), objective = nll, x = y$y, lower = c(0,0))

# Extract the optimized parameters
lambda <- opt$par[1]
beta <- opt$par[2]

# Reconstruct r and p
r <- 1/beta
p <- 1/(lambda*beta+1)

# Calculate the expectation from the geometric distribution
PoisG.mu <- r*(1-p)/p

u <- (y$y * beta + 1)/(lambda*beta+1)

plot(u)

# Agegroups ---------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

nll.age <- function(theta, data, formula){
  # Extract counts
  y <- data$y
  # Construct the design matrix
  designMatrix <- model.matrix(object = formula, data = data)
  # Construct regression parameters
  beta <- theta[1:ncol(designMatrix)]
  # # Extract agegroups
  ageGroup <- data$ageGroup
  
  # Define parameters
  lambda <- exp(designMatrix %*% beta - log(data$n))
  phi <- theta[-(1:ncol(designMatrix))]
  
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

formula <- y ~ -1 + ageGroup

# Optimize the parameters
opt.age <- nlminb(start = rep(1,22), objective = nll.age, data = y.age, formula = formula, lower = rep(1e-6,22))

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Define parameters
par.PoisG <- tibble(ageGroup = levels(y.age$ageGroup),
                    beta = opt.age$par[1:n.ageGroup],
                    phi = opt.age$par[(n.ageGroup+1):(n.ageGroup*2)])

# Construct parameters and calculate the expectation for the negative binomial distribution
par.PoisG %>%
  # Construct the size and probability
  mutate(r = 1/phi, p = 1/( beta*phi + 1 )) %>%
  # Calculate the expectation
  reframe(r, p, mu = r*(1-p)/p)

# Make inference on individual random effects
y.age <- y.age %>%
  mutate(u = ( y*beta.age[ageGroup]+1 )/( lambda.age[ageGroup] * beta.age[ageGroup] + 1 ))

# Save the results in a list
PoisG <- list(par = par.PoisG, fn = nll.age, results = y.age)

# Write the results
write_rds(x = PoisG, file = "PoissonGamma.rds")
