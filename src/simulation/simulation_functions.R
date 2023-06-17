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
  
  realization <- baselineData %>%
    mutate(outbreak = 0)
  
  # Randomly draw k
  kBaseline <- sample(x = c(2, 3, 5, 10), size = 4, replace = TRUE)
  # Count the number of outbreaks
  nOutbreaks <- 4
  # Calculate the standard deviation of the baseline data
  sdBaseline <- sd(realization$realization)
  # Randomly select a start for the oubreak
  outbreakStarts <- sample(x = refPar$baselineWeeks, size = nOutbreaks, replace = FALSE)
  
  # logs <- tibble()
  # Loop over outbreaks
  for(i in 1:nOutbreaks){
    # Randomly generate a number from a Poisson distribution with mean k
    outbreakSizes <- rpois(n = 1, lambda = kBaseline[i] * sdBaseline)
    # Distribute the cases in the individual outbreaks
    cases <- floor(rlnorm(outbreakSizes, meanlog = 0, sdlog = 0.5)) + outbreakStarts[i]
    # Count up the cases
    tmp <- plyr::count(cases)
    # Add them to the baseline data
    updatedRealization <- realization$realization
    updatedRealization[tmp$x] <-  updatedRealization[tmp$x] + tmp$freq
    updatedOutbreak <- realization$outbreak
    updatedOutbreak[tmp$x] <- updatedOutbreak[tmp$x] + 1
    # Add them back into the tibble
    realization <- realization %>% 
      mutate(realization = updatedRealization, outbreak = updatedOutbreak)
    
  }
  
  # Randomly draw k
  kCurrent <- sample(x = 1:10, size = 1)
  # Count the number of outbreaks
  nOutbreak <- 1
  # Randomly generate a number from a Poisson distribution with mean k
  outbreakSize <- rpois(n = nOutbreak, lambda = kCurrent * sdBaseline)
  # Randomly select a start for the oubreak
  outbreakStart <- sample(x = refPar$curWeeks, size = nOutbreak, replace = FALSE)
  
  # Distribute the cases in the individual outbreaks
  cases <- floor(rlnorm(outbreakSize, meanlog = 0, sdlog = 0.5)) + outbreakStart
  cases[which(cases > max(refPar$curWeeks))] <- max(refPar$curWeeks)
  # Count up the cases
  tmp <- plyr::count(cases)
  # Add them to the baseline data
  updatedRealization <- realization$realization
  updatedRealization[tmp$x] <-  updatedRealization[tmp$x] + tmp$freq
  updatedOutbreak <- realization$outbreak
  updatedOutbreak[tmp$x] <- updatedOutbreak[tmp$x] + 1
  # Add them back into the tibble
  realization <- realization %>% 
    mutate(realization = updatedRealization, outbreak = updatedOutbreak)
  
  return(list(realization=realization,k=kCurrent))
}
