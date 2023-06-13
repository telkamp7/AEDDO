# This script tries to replicate the outbreak simulation study, as seen in https://doi.org/10.1002/sim.5595

# Import the libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(surveillance)
library(TMB)

source("../models/aeddo.R")

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

simOutbreak <- function(baselineData){
  
  realization <- baselineData$realization
  
  # Randomly draw k
  kBaseline <- sample(x = c(2, 3, 5, 10), size = 4, replace = TRUE)
  # Count the number of outbreaks
  nOutbreaks <- 4
  # Calculate the standard deviation of the baseline data
  sdBaseline <- sd(realization)
  # Randomly select a start for the oubreak
  outbreakStarts <- sample(x = baselineWeeks, size = nOutbreaks, replace = FALSE)
  
  logs <- tibble()
  # Loop over outbreaks
  for(i in 1:nOutbreaks){
    # Randomly generate a number from a Poisson distribution with mean k
    outbreakSizes <- rpois(n = 1, lambda = kBaseline[i] * sdBaseline)
    # Distribute the cases in the individual outbreaks
    cases <- floor(rlnorm(outbreakSizes, meanlog = 0, sdlog = 0.5)) + outbreakStarts[i]
    # Count up the cases
    tmp <- plyr::count(cases)
    # Add them to the baseline data
    realization[tmp$x] <- realization[tmp$x] + tmp$freq
    logs <- bind_rows(
      logs,
      tibble(Start = outbreakStarts[i],
             End = max(c(Start,tmp$x)),
             k = kBaseline[i],
             Size = outbreakSizes,
             Distributed = list(tmp))
      )
  }
  
  # Randomly draw k
  kCurrent <- sample(x = 1:10, size = 1)
  # Count the number of outbreaks
  nOutbreak <- 1
  # Randomly generate a number from a Poisson distribution with mean k
  outbreakSize <- rpois(n = nOutbreak, lambda = kCurrent * sdBaseline)
  # Randomly select a start for the oubreak
  outbreakStart <- sample(x = curWeeks, size = nOutbreak, replace = FALSE)
  
  # Distribute the cases in the individual outbreaks
  cases <- floor(rlnorm(outbreakSize, meanlog = 0, sdlog = 0.5)) + outbreakStart
  cases[which(cases > max(curWeeks))] <- max(curWeeks)
  # Count up the cases
  tmp <- plyr::count(cases)
  # Add them to the baseline data
  realization[tmp$x] <- realization[tmp$x] + tmp$freq
  logs <- bind_rows(
    logs,
    tibble(Start = outbreakStart,
           End = max(c(Start,tmp$x)),
           k = kCurrent,
           Size = outbreakSize,
           Distributed = list(tmp))
  )
  
  return(list(realization = tibble(t = baselineData$t, realization = realization), logs = logs))
  
}

# Determine paramters
theta <- c(0.1, -2, 1.5, 0.5, 2.5, 3.75, 5)
beta <- c(0.0025, 0.005, 0.003, 0.002, 0.001, 0.001, 0.0001)
gamma1 <- c(0.6, 0.1, 0.2, 0.5, 1, 0.1, 0.05)
gamma2 <- c(0.6, 0.3, -0.4, 0.5, 0.1, -0.1, 0.01)
phi <- c(1.5, 2, 1, 5, 3, 1.1, 1.2)
m <- 0:1
trend <- 0:1

# Allocate space for scenarios
scenarios <- tibble()

# Generate the 42 scenarios
for(i in 1:7){
  clusterScenario <- expand_grid(theta = theta[i], beta = c(0,beta[i]), m,  phi = phi[i])
  
  clusterScenario$gamma1 <- gamma1[i]
  clusterScenario[clusterScenario$m == 0,]$gamma1 <- 0
  
  clusterScenario$gamma2 <- gamma2[i]
  clusterScenario[clusterScenario$m == 0,]$gamma2 <- 0
  
  clusterScenario$trend <- rep(trend, each = length(m))
  
  scenarios <- bind_rows(scenarios, clusterScenario)
  
}

# Inspect the scenarios and see if they match the paper
scenarios %>%
  select(theta, beta, gamma1, gamma2, phi, m, trend) %>%
  print(n=28)

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

# Construct parameters
simulationsPar <- expand_grid(scenario = 1:nrow(scenarios), rep = 1:nRep, t = 1:n) %>%
  mutate(par = map(scenario, function(x) scenarios[x,]))

# realizations
realizations <- simulationsPar %>%
  mutate(realization = map2_dbl(.x = t, .y = par, .f = function(x, y) simFromNB(t = x, theta = y$theta, beta = y$beta, gamma1 = y$gamma1, gamma2 = y$gamma2, m = y$m, trend = y$trend, phi = y$phi))) %>%
  group_by(scenario, rep) %>%
  nest(baselineData = c(t, realization))

outbreaks <- realizations %>%
  mutate(realization = map(baselineData, .f = function(x) simOutbreak(baselineData = x))) %>%
  ungroup() %>%
  unnest_wider(realization)

con.farrington <- list(
  range = NULL, b = 5, w = 3,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(5,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
  thersholdMethod = "delta"
)

con.noufaily <- list(
  range = NULL, b = 5, w = 3,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(5,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
  thersholdMethod = "Noufaily"
)

compile(file = "../models/PoissonNormal.cpp")
dyn.load(dynlib("../models/PoissonNormal"))
compile(file = "../models/PoissonGamma.cpp")
dyn.load(dynlib("../models/PoissonGamma"))


results <- tibble()
# Loop over all the realizations
# for(i in 1:nrow(realizations)){
for(i in 1:2){
  
  # Unpack the series
  dat <- outbreaks %>%
    filter(row_number() == i) %>%
    unnest(realization) %>%
    select(t, realization)

  # Construct sts for surveillance package
  disease.sts <- sts(observed = dat$realization, frequency = 52)
  
  Farrington <- farringtonFlexible(sts = disease.sts, control = con.farrington)
  
  Noufaily <- farringtonFlexible(sts = disease.sts, control = con.noufaily)
  
  # Construct the design matrix
  designMatrix <- model.matrix(object = realization ~ 1 + t + cos(2*pi*t/52) + sin(2*pi*t/52), data = dat)
  
  # Construct objective function with derivatives based on a compiled C++ template
  PoisG <- MakeADFun(
    data = list(y = dat$realization, n = 1, X = designMatrix),
    parameters = list(beta = rep(0, ncol(designMatrix)),
                      log_phi_u = log(1)),
    hessian = TRUE,
    method = "BFGS",
    DLL = "PoissonGamma",
    silent = TRUE
  )
  
  
  # Optimize the function
  opt <- do.call(what = "optim", args = PoisG)
  
  
  
  
  results <- bind_rows(results,
                       tibble(Farrington = Farrington@alarm[,1],
                              Noufaily = Noufaily@alarm[,1])
                       )
    
}


tmp <- outbreaks %>% 
  unnest_wider(realization) %>%
  unnest(realization) %>%
  mutate(con.farrington = list(con.farrington),
         con.noufaily = list(con.noufaily),
         sts = list(sts(observed = realization, frequency = 52))) 




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




