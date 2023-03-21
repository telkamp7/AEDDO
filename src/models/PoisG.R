
# Import libraries
library(readr)
library(dplyr)
library(numDeriv)

# Load in the processed data
# dat <- read_rds(file = "../../data/processed/dat.rds") # 11-agegroups
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6-agegroups

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date) %>%
  reframe(y = sum(cases), n = sum(n))

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

# Agegroups 11 ---------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

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

formula <- y ~ -1 + ageGroup

# Optimize the parameters
opt.age <- nlminb(start = c(rep(14,nlevels(y.age$ageGroup)),rep(0.5,nlevels(y.age$ageGroup))),
                  objective = nll.age11,
                  data = y.age,
                  formula = formula,
                  lower = rep(1e-6,nlevels(y.age$ageGroup)*2))

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Define parameters
par.PoisG <- tibble(ageGroup = levels(y.age$ageGroup),
                    beta = opt.age$par[1:n.ageGroup],
                    phi = opt.age$par[(n.ageGroup+1):(n.ageGroup*2)])

# Construct parameters and calculate the expectation for the negative binomial distribution
# par.PoisG %>%
#   # Construct the size and probability
#   mutate(r = 1/phi, p = 1/( beta*phi + 1 )) %>%
#   # Calculate the expectation
#   reframe(r, p, mu = r*(1-p)/p)

# library(ggplot2)
# y.age %>%
#   mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 ) / ( lambda*phi+1 ) ) %>%
#   ggplot(mapping = aes(x = Date, y = u)) +
#   geom_point() +
#   facet_wrap(facets = vars(ageGroup))

# Make inference on individual random effects
y.age <- y.age %>%
  full_join(y = par.PoisG, join_by(ageGroup)) %>%
  mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 )/( lambda * phi + 1 )) %>%
  select(Date:n, u)

# y.age %>%
#   group_by(ageGroup) %>%
#   summarize(mean(y), mean(u))

# Save the results in a list
PoisG <- list(par = par.PoisG, fn = nll.age, results = y.age)

# Write the results
write_rds(x = PoisG, file = "PoissonGamma.rds")


# 6 agegroups -----------------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

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

formula <- y ~ -1 + ageGroup

# Optimize the parameters
opt.age <- nlminb(start = c(rep(14,nlevels(y.age$ageGroup)),rep(0.5,1)),
                  objective = nll.age6,
                  data = y.age,
                  formula = formula,
                  lower = rep(1e-6,nlevels(y.age$ageGroup)+1))

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Calculate uncertainty directly in R
H <- hessian(nll.age6, opt.age$par, data = y.age, formula = formula)
se.theta <- sqrt(diag(solve(H)))

# Construct parameter table
par.tbl.PoisG <- tibble(Parameter = c(paste0("$\\beta_{", levels(y.age$ageGroup), "}$"), "$\\phi$"),
                        Estimate = opt.age$par,
                        `Std. Error` = se.theta)

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Define parameters
par.inf.PoisG <- tibble(ageGroup = levels(y.age$ageGroup),
                    beta = opt.age$par[1:n.ageGroup],
                    phi = opt.age$par[-c(1:n.ageGroup)])

# Construct parameters and calculate the expectation for the negative binomial distribution
# par.PoisG %>%
#   # Construct the size and probability
#   mutate(r = 1/phi, p = 1/( beta*phi + 1 )) %>%
#   # Calculate the expectation
#   reframe(r, p, mu = r*(1-p)/p)

# library(ggplot2)
# y.age %>%
#   mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 ) / ( lambda*phi+1 ) ) %>%
#   ggplot(mapping = aes(x = Date, y = u)) +
#   geom_point() +
#   facet_wrap(facets = vars(ageGroup))

# Make inference on individual random effects
results <- y.age %>%
  full_join(y = par.inf.PoisG, join_by(ageGroup)) %>%
  mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 )/( lambda * phi + 1 )) %>%
  select(Date:n, u)

# y.age %>%
#   group_by(ageGroup) %>%
#   summarize(mean(y), mean(u))

# Save the results in a list
PoisG <- list(par.tbl = par.tbl.PoisG,
              par.inf = par.inf.PoisG,
              fn = nll.age6,
              results = results)

# Write the results
write_rds(x = PoisG, file = "PoissonGamma.rds")

