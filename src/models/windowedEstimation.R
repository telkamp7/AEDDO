
# Import the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Source the negative log-likelihood functions
source(file = "nll_PoisG.R")



# Windowed estimation
dat_nest <- dat %>%
  group_by(caseDef, Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(monthInYear = as.integer(format(Date, "%m"))) %>%
  group_by(caseDef) %>%
  nest()

y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

y.season <- y %>%
  mutate(monthInYear = as.integer(format(Date, "%m")))

windowedEstimation <- function(data, k = 36, nll, theta, formula, lower, model, excludePastOutbreaks = FALSE){
  
  # Construct the full design matrix
  # designMatrix <- model.matrix(object = formula, data = data)
  
  # Extract the fixed effect parameter names
  # parNames <- colnames(designMatrix)
  
  if(model == "PoissonNormal"){
    library(TMB)
    dyn.load(dynlib("PoissonNormal")) 
    # names(theta) <- c(rep("beta", nlevels(data$ageGroup)), "log_sigma_u")
  }else if(model == "PoissonGamma"){
    library(numDeriv)
    # names(theta) <- c(rep("beta", nlevels(data$ageGroup)), "phi")
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
    refDate <- Dates[(i+k),]
    
    # Extract the observations in this window
    yWindow <- data %>%
      filter(Date %in% datesInWindow$Date)
    #... and the reference data
    refData <- data %>%
      filter(Date %in% refDate)
    
    # Exclude past observations, if they were deemed an outbreak
    # Turned on by 'excludePastOutbreaks = TRUE'
    if(excludePastOutbreaks == TRUE){
      if(nrow(pastOutbreaks) > 0){
        yWindow <- yWindow %>%
          # select(Date, ageGroup, y, n) %>%
          setdiff(pastOutbreaks)
      }
    }
    
    # Construct the design matrix
    designMatrix <- model.matrix(object = formula, data = yWindow)
    
    if(model == "PoissonNormal"){
      
      # Make the function for the window
      PoisLN <- MakeADFun(
        data = list(y = yWindow$y, x = yWindow$n, X = designMatrix),
        parameters = list(u = rep(1, length(yWindow$y)),
                          beta = theta[1:ncol(designMatrix)],
                          log_sigma_u = theta[-(1:ncol(designMatrix))]),
        hessian = TRUE,
        random = "u",
        DLL = "PoissonNormal",
        silent = TRUE
      )
      
      # Pak argumenter ud og brug optim
      # Optimize the function
      opt <- do.call(what = "optim", args = PoisLN)
      
      ## Report on the random effects
      rep <- sdreport(PoisLN)
      
      # Make a design matrix for the reference data
      refDesign <- model.matrix(object = formula, data = refData)
      
      # Make the function for the reference data
      refPoisLN <- MakeADFun(
        data = list(y = refData$y, x = refData$n, X = refDesign),
        parameters = list(u = rep(1, length(refData$y)),
                          beta = theta[1:ncol(refDesign)],
                          log_sigma_u = theta[-(1:ncol(refDesign))]),
        random = "u",
        DLL = "PoissonNormal",
        silent = TRUE
      )
      
      # Calculate the random effects
      refRep <- sdreport(obj = refPoisLN, par.fixed = opt$par, hessian.fixed = opt$hessian)
      
      # Extract the standard deviation on the fixed effect parameter estimates
      se.theta <- sqrt(diag(rep$cov.fixed))
      
      # Extract the standard deviation from the Gaussian model
      log_sigma <- opt$par[ncol(designMatrix)+1]
      
      # Make inference on random effects and calculate alarms
      # ran.ef <- tibble(
      #   window.date = yWindow$Date,
      #   ageGroup = yWindow$ageGroup,
      #   n = yWindow$n,
      #   y = yWindow$y,
      #   u = rep$par.random,
      #   p = pnorm(q = u, mean = 0, sd = exp(log_sigma)),
      #   alarm = p >= 0.95
      # )
      
      # Include lambda at some point
      # c(exp(refDesign %*% opt$par[1:ncol(refDesign)] - log(refData$n)))
      
      ran.ef <- tibble(
        ref.date = refDate$Date,
        ageGroup = refData$ageGroup,
        n = refData$n,
        y = refData$y,
        monthInYear = refData$monthInYear,
        log_sigma = opt$par[-(1:ncol(designMatrix))],
        u = refRep$par.random,
        p = pnorm(q = u, mean = 0, sd = exp(log_sigma)),
        alarm = p >= 0.95
      )
      
      par <- tibble(
        Parameter = c(colnames(designMatrix),"log_sigma"),
        theta = opt$par,
        se.theta = se.theta
        )
      
      # Combine the results in a tibble
      results <- bind_rows(results,
                           tibble(
                             ref.date = refDate$Date,
                             par = list(par),
                             log_sigma = opt$par[-(1:ncol(designMatrix))],
                             # theta = list(opt$par),
                             # se.theta = list(se.theta),
                             objective = opt$objective,
                             window.data = list(yWindow),
                             ran.ef = list(ran.ef), 
                             convergence = opt$convergence,
                             message = opt$message
                             )
                           )
      
      # Construct parameter guess for next iteration
      theta <- opt$par * 0.7
      
    }else if(model == "PoissonGamma"){
      
      # Optimize the parameters 
      opt <- nlminb(start = theta,
                    objective = nll,
                    data = yWindow,
                    designMatrix = designMatrix,
                    lower = lower)
      # Calculate the hessian
      H <- hessian(func = nll,
                   x = opt$par,
                   data = yWindow,
                   designMatrix = designMatrix)
      # ... and the standard deviation on the fixed effect parameter estimates
      se.theta <- sqrt(diag(solve(H)))
      # names(se.theta) <- names(opt$par)
      
      # Reconstruct the parameters
      # beta <- opt$par[1:ncol(designMatrix)]
      # lambda <- c(exp(designMatrix %*% beta - log(yWindow$n)))
      # # lambda <- exp(beta - log(yWindow$n))
      # phi <- opt$par[-(1:ncol(designMatrix))]
      
      # # Make inference on the random effects
      # ran.ef <- tibble(
      #   window.date = yWindow$Date,
      #   ageGroup = yWindow$ageGroup,
      #   n = yWindow$n,
      #   y = yWindow$y,
      #   u = ( yWindow$y*phi+1 ) / ( lambda*phi+1 ),
      #   p = pgamma(q = u, shape = 1/phi, scale = phi),
      #   alarm = p >= 0.95
      #   )
      
      refDesign <- model.matrix(object = formula, data = refData)
      
      ran.ef <- tibble(
        ref.date = refDate$Date,
        ageGroup = refData$ageGroup,
        n = refData$n,
        y = refData$y,
        monthInYear = refData$monthInYear,
        # refDesign = model.matrix(object = formula, data = refData),
        # beta = ,
        phi = exp(opt$par[-(1:ncol(refDesign))]),
        lambda = c(exp(refDesign %*% opt$par[1:ncol(refDesign)] - log(n))),
        u = ( y*phi+1 ) / (lambda * phi + 1),
        p = pgamma(q = u, shape = 1/phi, scale = phi),
        alarm = p >= 0.95
      )
      
      # u <- ( yWindow$y*phi+1 ) / ( lambda*phi+1 )
      
      par <- tibble(
        Parameter = c(colnames(designMatrix),"phi"),
        theta = opt$par,
        se.theta = se.theta
      )
      
      # Collect the results
      results <- bind_rows(results,
                           tibble(
                             ref.date = refDate$Date,
                             par = list(par),
                             phi = opt$par[-(1:ncol(designMatrix))],
                             # theta = list(opt$par),
                             # se.theta = list(se.theta),
                             objective = opt$objective,
                             window.data = list(yWindow),
                             ran.ef = list(ran.ef),
                             convergence = opt$convergence,
                             message = opt$message
                             )
                           )
      
      # Use current estimate as guess for the next iteration
      theta <- opt$par * 0.7
    }
    
    if(excludePastOutbreaks == TRUE & sum(ran.ef$alarm) > 0){
      pastOutbreaks <- ran.ef %>%
        filter(alarm == TRUE) %>%
        select(Date = ref.date, ageGroup, y, n, monthInYear) %>%
        bind_rows(pastOutbreaks)
    }
    
  }
  
  return(results)
  
}

# 
formula <- y ~ -1 + ageGroup + sin(2*pi/12*monthInYear) + cos(2*pi/12*monthInYear)
season.model <- model.matrix(object = formula, data = y.season)

theta <- rep(1,ncol(season.model)+1)
PoissonGamma <- windowedEstimation(data = y.season,
                          k = 36,
                          nll = nll.model,
                          theta = theta,
                          formula = formula,
                          lower = c(rep(1e-6,6),rep(-Inf,2),-8),
                          model = "PoissonGamma",
                          excludePastOutbreaks = TRUE)

tmp <- PoissonGamma %>%
  select(ref.date, par, message) %>%
  unnest(par)


# 
PoissonNormal <- windowedEstimation(data = y.season,
                                    k = 36,
                                    nll = nll.model,
                                    theta = theta,
                                    formula = formula,
                                    lower = c(rep(1e-6,6),rep(-Inf,2),-4),
                                    model = "PoissonNormal",
                                    excludePastOutbreaks = TRUE)
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
      nll = nll.model,
      theta = rep(1,7),
      formula = y~ -1 + ageGroup,
      lower = rep(1e-6,7),
      model = "PoissonGamma",
      excludePastOutbreaks = TRUE)
    }
    ),
    PoissonGamma_seasonal = map(data, function(df){
      windowedEstimation(
      data = df,
      k = 36,
      nll = nll.model,
      theta = rep(1,9),
      formula = y ~ -1 + ageGroup + sin(2*pi/12*monthInYear) + cos(2*pi/12*monthInYear),
      lower = c(rep(1e-6,6),rep(-Inf,2),1e-6),
      model = "PoissonGamma",
      excludePastOutbreaks = TRUE)
    }
    ),
    PoissonNormal_excludePastOutbreaks = map(data, function(df) {
      windowedEstimation(
        data = df,
        k = 36,
        nll = nll.model,
        theta = rep(1,7),
        formula = y~ -1 + ageGroup,
        lower = c(rep(1e-6,6),-4),
        model = "PoissonNormal",
        excludePastOutbreaks = TRUE)
      }
      ),
    PoissonNormal_seasonal = map(data, function(df) {
      windowedEstimation(
        data = df,
        k = 36,
        nll = nll.model,
        theta = rep(1,9),
        formula = y ~ -1 + ageGroup + sin(2*pi/12*monthInYear) + cos(2*pi/12*monthInYear),
        lower = c(rep(1e-6,6),rep(-Inf,2),-4),
        model = "PoissonNormal",
        excludePastOutbreaks = TRUE)
    }
    )
    )

# STEC_res$PoissonGamma[[1]]$ran.ef
# 
# STEC_res$PoissonNormal[[1]]

write_rds(x = STEC_res, file = "STEC_res.rds")
