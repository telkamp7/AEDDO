# This script tries to replicate the outbreak simulation study, as seen in https://doi.org/10.1002/sim.5595

# Import the libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Construct the function for the mean mu_t and linear predictor including trend and 
# seasonality determined by Fourier terms.
mu_t <- function(t, theta, beta, gamma1, gamma2, m, trend){
  # Start construction of linear predictor
  intExponent <- theta
  # Add a trend, if the scenario request it
  if(trend == 1){
    intExponent <- intExponent + beta * t
  }
  # Add a seasonal component, it the scenario request it
  if(m != 0){
    j <- 1:m
    
    intExponent <- intExponent + sum( gamma1 * cos( (2*pi*j*t)/52 ) + gamma2 * sin( (2*pi*j*t)/52 ) )
    
  }
  # Return the exponentation
  return( exp(intExponent) )
}

# Make the function to simulate from the negative binomial distribution with
# mean mu and variance phi*mu
simFromNB <- function(t, theta, beta, gamma1, gamma2, m, trend, phi){
  # Exctract the mu for this timepoint
  mu <- mu_t(t, theta, beta, gamma1, gamma2, m, trend)
  # Compute the variance of the distribution in this timepoint
  variance = phi * mu
  # Infer the size in the NB
  size = (mu + mu^2)/variance
  # Simulate the data
  rnbinom(n = 1, mu = mu, size = size)
}

simOutbreak <- function(k, baselineData, weeks, period){
  
  # Randomly generate a number from a Poisson distribution with mean k
  randN <- rpois(n = 1, lambda = k)
  # Calculate the standard deviation of the baseline data
  sdBaseline <- sd(baselineData)
  # Generate the outbreak size
  outbreakSize <- randN * sdBaseline
  # Randomly select a start for the oubreak
  outbreakStart <- sample(x = weeks, size = 1 + (3 * period == "baseline")  , replace = FALSE)
  
  
  finalData <- baselineData
  
  dlnorm(x = 0:10, meanlog = 0, sdlog = 0.5) * 10
  
  sum(rlnorm(10, meanlog = 0, sdlog = 0.5))
  
  
}


dlnorm(x = 0:100, meanlog = 0, sdlog = 0.5)


# Determine paramters
theta <- c(0.1, -2, 1.5, 0.5, 2.5, 3.75, 5)
beta <- c(0.0025, 0.005, 0.003, 0.002, 0.001, 0.001, 0.0001)
gamma1 <- c(0.6, 0.1, 0.2, 0.5, 1, 0.1, 0.05)
gamma2 <- c(0.6, 0.3, -0.4, 0.5, 0.1, -0.1, 0.01)
phi <- c(1.5, 2, 1, 5, 3, 1.1, 1.2)
m <- 0:2
trend <- 0:1

# Allocate space for scenarios
scenarios <- tibble()

# Generate the 42 scenarios
for(i in 1:7){
  clusterScenario <- expand_grid(theta = theta[i], beta = c(0,beta[i]), m = m, phi = phi[i])
  
  clusterScenario$gamma1 <- gamma1[i]
  clusterScenario[clusterScenario$m == 0,]$gamma1 <- 0
  
  clusterScenario$gamma2 <- gamma2[i]
  clusterScenario[clusterScenario$m == 0,]$gamma2 <- 0
  
  clusterScenario$trend <- rep(trend, each = 3)
  
  scenarios <- bind_rows(scenarios, clusterScenario)
  
}

# Inspect the scenarios and see if they match the paper
scenarios %>%
  select(theta, beta, gamma1, gamma2, phi, m, trend) %>%
  print(n=42)

# Assign the number of replicates
nRep <- 100
# ... and the size of each scenario (in weeks)
n <- 624

# Now we assign the weeks used for training the adaptive re-weighting schemes
trainWeeks <- 1:312
# ... and the weeks for baseline
baselineWeeks <- 313:575
# ... and the curWeeks
curWeeks <- 576:624


simulationsPar <- expand_grid(scenario = 1:42, t = 1:n) %>%
  mutate(par = map(scenario, function(x) scenarios[x,]))

# realizations

# test <- simulations %>%
#   mutate(par = map(scenario, function(x) scenarios[x,]))

realizations <- simulationsPar %>%
  mutate(realization = map2_dbl(.x = t, .y = par, .f = function(x, y) simFromNB(t = x, theta = y$theta, beta = y$beta, gamma1 = y$gamma1, gamma2 = y$gamma2, m = y$m, trend = y$trend, phi = y$phi)))


# test <- test %>%
#   filter(scenario %in% c(8,10,12,17)) %>%
#   mutate(realization = map2_dbl(.x = t, .y = par, .f = function(x, y) simFromNB(t = x, theta = y$theta, beta = y$beta, gamma1 = y$gamma1, gamma2 = y$gamma2, m = y$m, trend = y$trend, phi = y$phi)))

realizations %>%
  filter(scenario %in% c(8, 10, 12, 17)) %>%
  ggplot(mapping = aes(x = t, y = realization)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(scenario), scales = "free_y")


subsetScenarios %>%
  ggplot(mapping = aes(x = t, y = realization)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(scenario), scales = "free_y")

simulations <- simulations %>%
  mutate(par = map(scenario, function(x) scenarios[x,])) %>%
  unnest(par) %>%
  group_by(scenario) %>%
  nest(data = t:trend)

subsetSimulations <- simulations %>% 
  filter(scenario %in% c(8, 10, 12, 17)) %>%
  mutate(realization = map(data, function(x) do.call(what = "simFromNB", args = x)))




prob <- diff(plnorm(q = 0:9, meanlog = 0, sdlog = 0.5))

startOutbreak <- 510
PoissonOutbreakSize <- 20


floor(rlnorm(PoissonOutbreakSize, meanlog = 0, sdlog = 0.5)) + startOutbreak



tmp <- sample(x = startOutbreak+0:(length(prob)-1), size = PoissonOutbreakSize, replace = TRUE, prob = prob)
table(tmp)


round(dlnorm(x= 0:10, meanlog = 0, sdlog = 0.5),2)

simulations <- tibble()
for(j in 1:nrow(scenarios)){
  
  for(t in 1:n){
    scenarioPar <- scenarios[j, ] %>%
      mutate(t = t) %>%
      as.list()
    
    
    do.call("mu_t", args = scenarioPar)
    
  }
}



mu_t(t = 1, theta = 0.1, beta = 0, gamma1 = 0, gamma2 = 0, m = 0, trend = 0)




