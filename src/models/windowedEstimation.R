
# Import the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Source the negative log-likelihood functions
source(file = "nll_PoisG.R")

# Import the data
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6 agegroups

# Windowed estimation
dat_nest <- dat %>%
  group_by(caseDef, Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  group_by(caseDef) %>%
  nest()

y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

windowedEstimation <- function(data, k = 36, nll, theta, formula, lower, model, excludePastOutbreaks = FALSE){
  
  if(model == "PoissonNormal"){
    library(TMB)
    dyn.load(dynlib("PoissonNormal")) 
    names(theta) <- c(rep("beta", nlevels(data$ageGroup)), "log_sigma_u")
  }else if(model == "PoissonGamma"){
    library(numDeriv)
    names(theta) <- c(rep("beta", nlevels(data$ageGroup)), "phi")
  }else{
    warning("Incorrect specification of 'model'. Should be either 'PoissonNormal' or 'PoissonGamma'", call. = FALSE)
  }
  
  # Extract the months
  Dates <- data %>%
    reframe(Date = unique(Date))
  
  # Count the number of observations
  nObs <- length(Dates$Date)
  
  # Make a placeholder for the 'results' and 'pastOutbreaks'
  results <- tibble()
  if(excludePastOutbreaks == TRUE){
    pastOutbreaks <- tibble()
  }
  # for(i in 1:6){
  for(i in 1:(nObs-k)){
    
    # Compute the dates in window
    datesInWindow <- Dates[i:(i+k-1),]
    #... and the reference date
    refData <- Dates[(i+k-1),]
    
    # Extract the observations in this window
    yWindow <- data %>%
      filter(Date %in% datesInWindow$Date)
    
    
    # Exclude past observations, if they were deemed an outbreak
    # Turned on by 'excludePastOutbreaks = TRUE'
    if(excludePastOutbreaks == TRUE){
      if(nrow(pastOutbreaks) > 0){
        yWindow <- yWindow %>%
          setdiff(pastOutbreaks)
      }
    }
    
    # Construct the design matrix
    designMatrix <- model.matrix(object = formula, data = yWindow)
    
    if(model == "PoissonNormal"){
      
      # Make the function
      PoisLN <- MakeADFun(
        data = list(y = yWindow$y, x = yWindow$n, X = designMatrix),
        parameters = list(u = rep(1, length(yWindow$y)),
                          beta = theta[1:nlevels(yWindow$ageGroup)],
                          log_sigma_u = theta[nlevels(yWindow$ageGroup)+1]),
        random = "u",
        DLL = "PoissonNormal",
        silent = TRUE
      )
      
      # Optimize the function
      opt <- nlminb(start = PoisLN$par, PoisLN$fn, PoisLN$gr, lower = lower)
      
      ## Report on the random effects
      rep <- sdreport(PoisLN)
      # Extract the standard deviation on the fixed effect parameter estimates
      se.theta <- sqrt(diag(rep$cov.fixed))
      
      # Extract the standard deviation from the Gaussian model
      log_sigma <- opt$par[nlevels(yWindow$ageGroup)+1]
      
      # Make inference on random effects and calculate alarms
      ran.ef <- tibble(
        window.date = yWindow$Date,
        ageGroup = yWindow$ageGroup,
        n = yWindow$n,
        y = yWindow$y,
        u = rep$par.random,
        p = pnorm(q = u, mean = 0, sd = exp(log_sigma)),
        alarm = p >= 0.95
      )
      
      par <- tibble(
        Parameter = c(paste0("$\\beta_{", unique(yWindow$ageGroup), "}$"),
                      "$\\log(\\sigma)$"),
        ageGroup = c(unique(yWindow$ageGroup),"All"),
        theta = opt$par,
        se.theta = se.theta
        )
      
      # Combine the results in a tibble
      results <- bind_rows(results,
                           tibble(
                             ref.date = refData$Date,
                             par = list(par),
                             log_sigma = opt$par[nlevels(yWindow$ageGroup)+1],
                             # theta = list(opt$par),
                             # se.theta = list(se.theta),
                             objective = opt$objective,
                             ran.ef = list(ran.ef)
                             )
                           )
      
      # Construct parameter guess for next iteration
      theta <- opt$par
      
    }else if(model == "PoissonGamma"){
      
      # Optimize the parameters 
      opt <- nlminb(start = theta,
                    objective = nll,
                    data = yWindow,
                    designMatrix = designMatrix,
                    # formula = formula,
                    lower = lower)
      
      # Calculate the hessian
      H <- hessian(func = nll,
                   x = opt$par,
                   data = yWindow,
                   designMatrix = designMatrix)
      # ... and the standard deviation on the fixed effect parameter estimates
      se.theta <- sqrt(diag(solve(H)))
      names(se.theta) <- names(opt$par)
      
      # Reconstruct the parameters
      beta <- opt$par[1:nlevels(yWindow$ageGroup)]
      lambda <- c(exp(designMatrix %*% beta - log(yWindow$n)))
      # lambda <- exp(beta - log(yWindow$n))
      phi <- opt$par[nlevels(yWindow$ageGroup)+1]
      
      # Make inference on the random effects
      ran.ef <- tibble(
        window.date = yWindow$Date,
        ageGroup = yWindow$ageGroup,
        n = yWindow$n,
        y = yWindow$y,
        u = ( yWindow$y*phi+1 ) / ( lambda*phi+1 ),
        p = pgamma(q = u, shape = 1/phi, scale = phi),
        alarm = p >= 0.95
        )
      
      # u <- ( yWindow$y*phi+1 ) / ( lambda*phi+1 )
      
      par <- tibble(
        Parameter = c(paste0("$\\beta_{", unique(yWindow$ageGroup), "}$"), "$\\phi$"),
        ageGroup = c(unique(yWindow$ageGroup),"All"),
        theta = opt$par,
        se.theta = se.theta
      )
      
      # Collect the results
      results <- bind_rows(results,
                           tibble(
                             ref.date = refData$Date,
                             par = list(par),
                             phi = opt$par[nlevels(yWindow$ageGroup)+1],
                             # theta = list(opt$par),
                             # se.theta = list(se.theta),
                             objective = opt$objective,
                             ran.ef = list(ran.ef)
                             )
                           )
      
      # Use current estimate as guess for the next iteration
      theta <- c(beta,phi)
    }
    
    if(excludePastOutbreaks == TRUE & sum(ran.ef$alarm) > 0){
      pastOutbreaks <- ran.ef %>%
        filter(alarm == TRUE) %>%
        select(Date = window.date, ageGroup, y, n) %>%
        bind_rows(pastOutbreaks)
    }
    
  }
  
  return(results)
  
}

# theta <- rep(1,7)
# 
# PoissonGamma <- windowedEstimation(data = y,
#                           k = 36,
#                           nll = nll.age6.window,
#                           theta = theta,
#                           formula = y ~ -1 + ageGroup,
#                           lower = rep(1e-6,7),
#                           model = "PoissonGamma",
#                           excludePastOutbreaks = TRUE)
# 
# PoissonNormal <- windowedEstimation(data = y,
#                                     k = 36,
#                                     nll = nll.age6.window,
#                                     theta = theta,
#                                     formula = y ~ -1 + ageGroup,
#                                     lower = c(rep(1e-6,6),-4),
#                                     model = "PoissonNormal",
#                                     excludePastOutbreaks = TRUE)
# 
# tmp <- dat_nest %>% 
#   filter(caseDef %in% c("Legionella","Shiga- og veratoxin producerende E. coli."))
# 
# tmp$data

STEC_res <- dat_nest %>% 
  filter(caseDef %in% c("Shiga- og veratoxin producerende E. coli.")) %>%
  mutate(
    PoissonGamma_excludePastOutbreaks = map(data, function(df){
      windowedEstimation(
      data = df,
      k = 36,
      nll = nll.age6.window,
      theta = rep(1,7),
      formula = y~ -1 + ageGroup,
      lower = rep(1e-6,7),
      model = "PoissonGamma",
      excludePastOutbreaks = TRUE)
    }
    ),
    PoissonGamma = map(data, function(df){
      windowedEstimation(
      data = df,
      k = 36,
      nll = nll.age6.window,
      theta = rep(1,7),
      formula = y~ -1 + ageGroup,
      lower = rep(1e-6,7),
      model = "PoissonGamma",
      excludePastOutbreaks = FALSE)
    }
    ),
    PoissonNormal_excludePastOutbreaks = map(data, function(df) {
      windowedEstimation(
        data = df,
        k = 36,
        nll = nll.age6.window,
        theta = rep(1,7),
        formula = y~ -1 + ageGroup,
        lower = c(rep(1e-6,6),-4),
        model = "PoissonNormal",
        excludePastOutbreaks = TRUE)
      }
      ),
    PoissonNormal = map(data, function(df) {
      windowedEstimation(
        data = df,
        k = 36,
        nll = nll.age6.window,
        theta = rep(1,7),
        formula = y~ -1 + ageGroup,
        lower = c(rep(1e-6,6),-4),
        model = "PoissonNormal",
        excludePastOutbreaks = FALSE)
    }
    )
    )

# STEC_res$PoissonGamma[[1]]$ran.ef
# 
# STEC_res$PoissonNormal[[1]]

write_rds(x = STEC_res, file = "STEC_res.rds")
