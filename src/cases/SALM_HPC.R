# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(surveillance)

# Source the estimation function
source(file = "../models/aeddo.R")

# Load in the data
dat <- read_rds(file = "../../data/processed/dat2.rds")

# Summary statistic of all the data
dat %>%
  group_by(Date, caseDef) %>%
  reframe(y = sum(cases)) %>% 
  group_by(caseDef) %>%
  summarise(meanCases = mean(y), medianCases = median(y))

# Only consider the SALM cases
SALM <- dat %>%
  filter(caseDef == "SALM") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- SALM %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- SALM %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
SALM.sts <- sts(
  observed = observed[,-1],
  epoch = observed$Date,
  epochAsDate = TRUE,
  frequency = 12,
  population = as.matrix(population[,-1])
)


# Farrington ------------------------------------------------------------------------

con.farrington <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "delta"
)

SALM_Farrington <- farringtonFlexible(sts = SALM.sts, con.farrington)
write_rds(x = SALM_Farrington, file = "SALM_Farrington.rds")

# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
)

SALM_Noufaily <- farringtonFlexible(sts = SALM.sts, con.noufaily)
write_rds(x = SALM_Noufaily, file = "SALM_Noufaily.rds")

# Hierarchical Poisson Normal model ---------------------------------------

SALM_PoisN_ageGroup <- aeddo(data = SALM,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(0,7),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,6), -6),
                             upper = rep(1e2, 7),
                             model = "PoissonNormal",
                             k = 36,
                             sig.level = 0.9,
                             cpp.dir = "../models/",
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisN_ageGroup, file = "SALM_PoisN_ageGroup.rds")

start.theta.PoisN <- SALM_PoisN_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SALM_PoisN_ageGroup_trend <- aeddo(data = SALM,
                                   formula = y ~ -1 + t + ageGroup,
                                   trend = TRUE,
                                   theta = c(0, start.theta.PoisN),
                                   method = "L-BFGS-B",
                                   lower = c(-0.5, start.theta.PoisN-6),
                                   upper = c(0.5, start.theta.PoisN+6),
                                   model = "PoissonNormal", 
                                   k = 36, 
                                   sig.level = 0.9,
                                   cpp.dir = "../models/",
                                   CI = TRUE,
                                   excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisN_ageGroup_trend, file = "SALM_PoisN_ageGroup_trend.rds")

SALM_PoisN_ageGroup_seasonality <- aeddo(data = SALM,
                                         formula = y ~ -1 +  ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7]) - 6,
                                         upper = c(start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7]) + 6,
                                         model = "PoissonNormal", 
                                         k = 36, 
                                         sig.level = 0.9,
                                         cpp.dir = "../models/",
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisN_ageGroup_seasonality, file = "SALM_PoisN_ageGroup_seasonality.rds")

SALM_PoisN_ageGroup_trend_seasonality <- aeddo(data = SALM,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7])-6),
                                               upper = c(0.5, c(start.theta.PoisN[1:6], 0,0, start.theta.PoisN[7])+6),
                                               model = "PoissonNormal", 
                                               k = 36, 
                                               sig.level = 0.9,
                                               cpp.dir = "../models/",
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisN_ageGroup_trend_seasonality, file = "SALM_PoisN_ageGroup_trend_seasonality.rds")

# Hierarchical Poisson Gamma model ----------------------------------------

SALM_PoisG_ageGroup <- aeddo(data = SALM,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(1,7),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,7),-6),
                             upper = rep(1e2, 7),
                             model = "PoissonGamma",
                             k = 36,
                             cpp.dir = "../models/",
                             sig.level = 0.9,
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisG_ageGroup, file = "SALM_PoisG_ageGroup.rds")

start.theta.PoisG <- SALM_PoisG_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SALM_PoisG_ageGroup_trend <- aeddo(data = SALM,
                                   formula = y ~ -1 + t + ageGroup,
                                   trend = TRUE,
                                   theta = c(0, start.theta.PoisG),
                                   method = "L-BFGS-B",
                                   lower = c(-0.5, start.theta.PoisG-6),
                                   upper = c(0.5, start.theta.PoisG+6),
                                   model = "PoissonGamma",
                                   k = 36,
                                   cpp.dir = "../models/",
                                   sig.level = 0.9,
                                   CI = TRUE,
                                   excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisG_ageGroup_trend, file = "SALM_PoisG_ageGroup_trend.rds")

SALM_PoisG_ageGroup_seasonality <- aeddo(data = SALM,
                                         formula = y ~ -1 + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]) - 6,
                                         upper = c(start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]) + 6,
                                         model = "PoissonGamma",
                                         k = 36,
                                         cpp.dir = "../models/",
                                         sig.level = 0.9,
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisG_ageGroup_seasonality, file = "SALM_PoisG_ageGroup_seasonality.rds")

SALM_PoisG_ageGroup_trend_seasonality <- aeddo(data = SALM,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]) - 6),
                                               upper = c(0.5, c(start.theta.PoisG[1:6], 0,0, start.theta.PoisG[7]) + 6),
                                               model = "PoissonGamma",
                                               k = 36,
                                               cpp.dir = "../models/",
                                               sig.level = 0.9,
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SALM_PoisG_ageGroup_trend_seasonality, file = "SALM_PoisG_ageGroup_trend_seasonality.rds")




