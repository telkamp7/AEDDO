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

# Only consider the SHIL cases
SHIL <- dat %>%
  filter(caseDef == "SHIL") %>%
  mutate(ageGroup = fct_collapse(ageGroup, `<25 years` = c("<1 year", "1-4 years","5-14 years","15-24 years")),
         ageGroup = fct_collapse(ageGroup, `25+ years` = c("25-64 years","65+ years"))) %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))


# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- SHIL %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- SHIL %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
SHIL.sts <- sts(
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

SHIL_Farrington <- farringtonFlexible(sts = SHIL.sts, con.farrington)
write_rds(x = SHIL_Farrington, file = "SHIL_Farrington.rds")


# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 10, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
)

SHIL_Noufaily <- farringtonFlexible(sts = SHIL.sts, con.noufaily)
write_rds(x = SHIL_Noufaily, file = "SHIL_Noufaily.rds")

# Hierarchical Poisson Normal model ---------------------------------------

SHIL_PoisN_ageGroup <- aeddo(data = SHIL,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(0,3),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,2), -6),
                             upper = rep(1e2, 3),
                             model = "PoissonNormal",
                             k = 36,
                             sig.level = 0.9,
                             cpp.dir = "../models/",
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisN_ageGroup, file = "SHIL_PoisN_ageGroup.rds")

start.theta.PoisN <- SHIL_PoisN_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SHIL_PoisN_ageGroup_trend <- aeddo(data = SHIL,
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

write_rds(x = SHIL_PoisN_ageGroup_trend, file = "SHIL_PoisN_ageGroup_trend.rds")

SHIL_PoisN_ageGroup_seasonality <- aeddo(data = SHIL,
                                         formula = y ~ -1 +  ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) - 6,
                                         upper = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) + 6,
                                         model = "PoissonNormal", 
                                         k = 36, 
                                         sig.level = 0.9,
                                         cpp.dir = "../models/",
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisN_ageGroup_seasonality, file = "SHIL_PoisN_ageGroup_seasonality.rds")

SHIL_PoisN_ageGroup_trend_seasonality <- aeddo(data = SHIL,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])-6),
                                               upper = c(0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])+6),
                                               model = "PoissonNormal", 
                                               k = 36, 
                                               sig.level = 0.9,
                                               cpp.dir = "../models/",
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisN_ageGroup_trend_seasonality, file = "SHIL_PoisN_ageGroup_trend_seasonality.rds")

# Hierarchical Poisson Gamma model ----------------------------------------




SHIL_PoisG_ageGroup <- aeddo(data = SHIL,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(1,3),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,2),-6),
                             upper = rep(1e2, 3),
                             model = "PoissonGamma",
                             k = 36,
                             cpp.dir = "../models/",
                             sig.level = 0.9,
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisG_ageGroup, file = "SHIL_PoisG_ageGroup.rds")

start.theta.PoisG <- SHIL_PoisG_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SHIL_PoisG_ageGroup_trend <- aeddo(data = SHIL,
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

write_rds(x = SHIL_PoisG_ageGroup_trend, file = "SHIL_PoisG_ageGroup_trend.rds")


SHIL_PoisG_ageGroup_seasonality <- aeddo(data = SHIL,
                                         formula = y ~ -1 + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6,
                                         upper = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6,
                                         model = "PoissonGamma",
                                         k = 36,
                                         cpp.dir = "../models/",
                                         sig.level = 0.9,
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisG_ageGroup_seasonality, file = "SHIL_PoisG_ageGroup_seasonality.rds")

SHIL_PoisG_ageGroup_trend_seasonality <- aeddo(data = SHIL,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6),
                                               upper = c(0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6),
                                               model = "PoissonGamma",
                                               k = 36,
                                               cpp.dir = "../models/",
                                               sig.level = 0.9,
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SHIL_PoisG_ageGroup_trend_seasonality, file = "SHIL_PoisG_ageGroup_trend_seasonality.rds")
