
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

# Construct the negative log-likelihood function
nll <- function(theta, x){
  p <- 1/(theta[1]*theta[2]+1)
  -sum(dgeom(x = x, prob = p, log = TRUE))
}

# Optimize the parameters
opt <- nlminb(start = c(1,1), objective = nll, x = y$y)

# Extract the optimized parameters
lambda <- opt$par[1]
phi <- opt$par[2]

# Reconstruct p
p <- 1/(lambda*phi+1)

# Calculate the expectation from the geometric distribution
(1-p)/p

u <- (y$y + 1)/(phi/(lambda*phi+1))

plot(log(u))

mean(y$y)


# Agegroups ---------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

nll.age <- function(theta, data){
  # Extract counts
  y <- data$y
  # Extract agegroups
  ageGroup <- data$ageGroup
  # Extract counts agegroups
  n.ageGroup <- n_distinct(data$ageGroup)
  
  # Define parameters
  lambda <- theta[1:n.ageGroup]
  phi <- theta[(n.ageGroup+1):(n.ageGroup*2)]
  
  # Construct the probability for the geometric distribution
  p <- 1/(lambda*phi+1)
  
  # Initilize the log-likelihood
  ll <- 0
  for(i in 1:nrow(data)){
    ll = ll + dgeom(x = y[i], prob = p[ageGroup[i]], log = TRUE) 
  }
  
  # Return the negative log-likelihood
  -ll
}

# Optimize the parameters
opt.age <- nlminb(start = rep(1,22), objective = nll.age, data = y.age)




as.numeric(as.character(ageGroup))
