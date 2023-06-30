# This script tries to replicate the outbreak simulation study, as seen in https://doi.org/10.1002/sim.5595
commandArgs(trailingOnly = TRUE)
args = (commandArgs(trailingOnly = TRUE))

listArgs <- as.list(unlist(lapply(args, function(x) {
  y <- strsplit(x, split = "=")[[1]]
  z <- y[2]
  names(z) <- y[1]
  return(z)
})))


# Converting numbers from character where possible
suppressWarnings(tmpNum <- as.numeric(listArgs))
listArgs[!is.na(tmpNum)] <- as.list(tmpNum[!is.na(tmpNum)])

# Import the libraries
library(tidyr)
library(readr)
library(plyr)
library(dplyr)
library(purrr)
library(doParallel)
library(TMB)
library(surveillance)

# Register number of cores to be used

registerDoParallel(cores = 20)

# Source the simulation functions
source("simulation_functions.R")
source("../models/aeddo.R")

# Determine paramters
refPar <- list(
  theta = c(0.1, -2, 1.5, 0.5, 2.5, 3.75, 5),
  beta = c(0.0025, 0.005, 0.003, 0.002, 0.001, 0.001, 0.0001),
  gamma1 = c(0.6, 0.1, 0.2, 0.5, 1, 0.1, 0.05),
  gamma2 = c(0.6, 0.3, -0.4, 0.5, 0.1, -0.1, 0.01),
  phi = c(1.5, 2, 1, 5, 3, 1.1, 1.2),
  m = 0:1,
  trend = 0:1,
  scenario = 1:28,
  nRep = 100, # Assign the number of replicates
  n = 624, # ... and the size of each scenario (in weeks)
  trainWeeks = 1:312, # Now we assign the weeks used for training the adaptive re-weighting schemes
  baselineWeeks = 313:575, # ... and the weeks for baseline
  curWeeks = 576:624, # ... and the curWeeks
  formula = y ~ 1 + t + cos(2*pi*t/52) + sin(2*pi*t/52),
  alphaStateOfTheArt = 0.005
)

if(length(listArgs) > 0){
  # Check which values to change
  changes <- match(names(listArgs), names(refPar))
  # Stop if trying to change non-existing parameter
  if(any(is.na(changes))){
    stop(paste( paste(names(listArgs)[is.na(changes)], collapse = ", "), 
                "is/are not valid parameter names,"))
  }
  # Change the requested parameters values and return
  for(i in 1:length(changes)){
    refPar[[changes[i]]] <- listArgs[[i]]
  }
}

con.farrington <- list(
  range = NULL, b = 5, w = 3,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = refPar$alphaStateOfTheArt, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded  = 3,
  thersholdMethod = "delta"
)

con.noufaily <- list(
  range = NULL, b = 5, w = 3,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = refPar$alphaStateOfTheArt, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 10, pastWeeksNotIncluded  = 26,
  thersholdMethod = "nbPlugin"
)

arb.dates <- seq.Date(from = as.Date("2000-01-01"), to = as.Date("2000-01-01")+(refPar$n-1)*7, length.out = refPar$n)

# Allocate space for scenarios
scenarios <- tibble()

# Generate the 28 scenarios
for(i in 1:7){
  clusterScenario <- expand_grid(theta = refPar$theta[i], beta = c(0,refPar$beta[i]), m = refPar$m,  phi = refPar$phi[i])
  
  clusterScenario$gamma1 <- refPar$gamma1[i]
  clusterScenario[clusterScenario$m == 0,]$gamma1 <- 0
  
  clusterScenario$gamma2 <- refPar$gamma2[i]
  clusterScenario[clusterScenario$m == 0,]$gamma2 <- 0
  
  clusterScenario$trend <- rep(refPar$trend, each = length(refPar$m))
  
  scenarios <- bind_rows(scenarios, clusterScenario)
  
}

# Inspect the scenarios and see if they match the paper
scenarios %>%
  select(theta, beta, gamma1, gamma2, phi, m, trend) %>%
  print(n=28)

write_rds(x = scenarios, file = "scenarios.Rds")

# Only consider the chosen scenarios
scenariosConsidered <- scenarios[refPar$scenario,]

Data <- foreach(sim = 1:refPar$nRep, .packages = c("dplyr", "tidyr", "purrr", "surveillance", "TMB", "plyr")) %dopar% {
  
  simulationPar <- expand_grid(scenario = refPar$scenario, t = 1:refPar$n) %>%
    mutate(par = map(scenario, function(x) scenarios[x,]))
  
  realizations <- simulationPar %>%
    mutate(realization = map2_dbl(.x = t, .y = par, .f = function(x, y) simFromNB(t = x, theta = y$theta, beta = y$beta, gamma1 = y$gamma1, gamma2 = y$gamma2, m = y$m, trend = y$trend, phi = y$phi))) %>%
    group_by(scenario) %>%
    mutate(Date = arb.dates, n = 1, ageGroup = "All") %>%
    nest(baselineData = c(Date, t, realization, n, ageGroup))
  
  outbreaks <- realizations %>%
    mutate(realization = map(baselineData, .f = function(x) simOutbreak(baselineData = x))) %>%
    ungroup()
  
  j <- outbreaks$scenario
  
  scenarioUnpack <- outbreaks %>%
    filter(scenario == j) %>%
    unnest_wider(realization)
  
  k <- scenarioUnpack %>%
    select(k)
  
  # Calculate baselines FPR
  baseline <- scenarioUnpack %>%
    select(scenario, baselineData) %>%
    unnest(baselineData) %>%
    select(Date, t, y = realization, n, ageGroup)
  
  baseline.sts <- sts(observed = baseline$y, frequency = 52, epoch = baseline$Date, )
  
  farrington_base <- farringtonFlexible(sts = baseline.sts, control = con.farrington)
  noufaily_base <- farringtonFlexible(sts = baseline.sts, control = con.noufaily)
  
  PoisN_base <- aeddo(data = baseline,
                 formula = refPar$formula,
                 trend = TRUE,
                 seasonality = TRUE,
                 theta = c(rep(0, 5)),
                 method = "BFGS",
                 model = "PoissonNormal", 
                 k = 52*3, 
                 sig.level = 0.95,
                 cpp.dir = "../models/",
                 period = "week",
                 excludePastOutbreaks = TRUE)
  
  PoisG_base <- aeddo(data = baseline,
                 formula = refPar$formula,
                 trend = TRUE,
                 seasonality = TRUE,
                 theta = c(rep(0, 5)),
                 method = "BFGS",
                 model = "PoissonGamma", 
                 k = 52*3, 
                 sig.level = 0.95,
                 cpp.dir = "../models/",
                 period = "week",
                 excludePastOutbreaks = TRUE)
  
  
  alarm_Farrington_base <- as_tibble(farrington_base@alarm) %>%
    mutate(Date = as.Date(x = farrington_base@epoch, origin = "1970-01-01")) %>%
    rename(alarm_Farrington = observed1)
  
  alarm_Noufaily_base <- as_tibble(noufaily_base@alarm) %>%
    mutate(Date = as.Date(x = noufaily_base@epoch, origin = "1970-01-01")) %>%
    rename(alarm_Noufaily = observed1)
  
  alarm_PoisN_base <- PoisN_base %>%
    select(ran.ef) %>% 
    unnest(ran.ef) %>%
    select(Date = ref.date, alarm_PoisN = alarm)
  
  alarm_PoisG_base <- PoisG_base %>%
    select(ran.ef) %>% 
    unnest(ran.ef) %>%
    select(Date = ref.date, alarm_PoisG = alarm)
  
  baselineAndFPR <- scenarioUnpack %>%
    select(scenario, baselineData) %>%
    unnest(baselineData) %>%
    select(Date, t, y = realization, n, ageGroup) %>%
    full_join(y = alarm_Farrington_base, by = join_by("Date")) %>%
    full_join(y = alarm_Noufaily_base, by = join_by("Date")) %>%
    full_join(y = alarm_PoisN_base, by = join_by("Date")) %>%
    full_join(y = alarm_PoisG_base, by = join_by("Date")) %>%
    nest() %>%
    rename(baseline = data)
  
  # Calculate outbreaks POD
  series <- scenarioUnpack %>%
    unnest(realization) %>%
    select(Date, t, y = realization, n, ageGroup)
  
  outbreak.sts <- sts(observed = series$y, frequency = 52, epoch = series$Date, )
  
  farrington <- farringtonFlexible(sts = outbreak.sts, control = con.farrington)
  noufaily <- farringtonFlexible(sts = outbreak.sts, control = con.noufaily)

  PoisN <- aeddo(data = series,
                 formula = refPar$formula,
                 trend = TRUE,
                 seasonality = TRUE,
                 theta = c(rep(0, 5)),
                 method = "BFGS",
                 model = "PoissonNormal", 
                 k = 52*3, 
                 sig.level = 0.95,
                 cpp.dir = "../models/",
                 period = "week",
                 excludePastOutbreaks = TRUE)
  
  PoisG <- aeddo(data = series,
                 formula = refPar$formula,
                 trend = TRUE,
                 seasonality = TRUE,
                 theta = c(rep(0, 5)),
                 method = "BFGS",
                 model = "PoissonGamma", 
                 k = 52*3, 
                 sig.level = 0.95,
                 cpp.dir = "../models/",
                 period = "week",
                 excludePastOutbreaks = TRUE)
  
  
  alarm_Farrington <- as_tibble(farrington@alarm) %>%
    mutate(Date = as.Date(x = farrington@epoch, origin = "1970-01-01")) %>%
    rename(alarm_Farrington = observed1)
  
  alarm_Noufaily <- as_tibble(noufaily@alarm) %>%
    mutate(Date = as.Date(x = noufaily@epoch, origin = "1970-01-01")) %>%
    rename(alarm_Noufaily = observed1)
  
  alarm_PoisN <- PoisN %>%
    select(ran.ef) %>% 
    unnest(ran.ef) %>%
    select(Date = ref.date, alarm_PoisN = alarm)
  
  alarm_PoisG <- PoisN %>%
    select(ran.ef) %>% 
    unnest(ran.ef) %>%
    select(Date = ref.date, alarm_PoisG = alarm)
  
  outbreaksAndAlarms <- scenarioUnpack %>%
    unnest(realization) %>%
    select(Date, t, y = realization, n, ageGroup, outbreak, k) %>%
    full_join(y = alarm_Farrington, by = join_by("Date")) %>%
    full_join(y = alarm_Noufaily, by = join_by("Date")) %>%
    full_join(y = alarm_PoisN, by = join_by("Date")) %>%
    full_join(y = alarm_PoisG, by = join_by("Date")) %>%
    mutate(outbreakTF = outbreak > 0) %>% 
    nest() %>%
    rename(outbreaks = data)

  ans <- tibble(baselineAndFPR) %>% 
    bind_cols(outbreaksAndAlarms)
      
  list(ans)
  
}

if(length(listArgs)>0){
  parString <- paste0(names(listArgs),unlist(listArgs), collapse = "_")
} else{
  parString <- "Default"
}

write_rds(x = list(Data=Data, refPar=refPar), file = paste0("~/proj/AEDDO/src/simulation/",parString, ".Rds"))




