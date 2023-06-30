
# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(readr)
library(surveillance)
library(zoo)

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

# Only consider the LIST cases
LIST <- dat %>%
  filter(caseDef == "LIST" & Date <= as.Date("2022-12-31")) %>%
  mutate(ageGroup = fct_collapse(ageGroup, `<65 years` = c("<1 year", "1-4 years", "5-14 years", "15-24 years", "25-64 years"))) %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- LIST %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- LIST %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
LIST.sts <- sts(
  observed = observed[,-1],
  epoch = observed$Date,
  epochAsDate = TRUE,
  frequency = 12,
  population = as.matrix(population[,-1])
)

# Farrington ------------------------------------------------------------------------

con.farrington <- list(
  range = NULL, b = 3, w = 3,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "delta"
)

LIST_Farrington <- farringtonFlexible(sts = LIST.sts, con.farrington)
write_rds(x = LIST_Farrington, file = "LIST_Farrington.rds")


# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 3,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 10, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
)

LIST_Noufaily <- farringtonFlexible(sts = LIST.sts, con.noufaily)
write_rds(x = LIST_Noufaily, file = "LIST_Noufaily.rds")

# Hierarchical Poisson Normal model -------------------------------------------------

LIST_PoisN_ageGroup <- aeddo(data = LIST,
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

write_rds(x = LIST_PoisN_ageGroup, file = "LIST_PoisN_ageGroup.rds")

# Hierarchical Poisson Gamma model ----------------------------------------------------

LIST_PoisG_ageGroup <- aeddo(data = LIST,
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

write_rds(x = LIST_PoisG_ageGroup, file = "LIST_PoisG_ageGroup.rds")
