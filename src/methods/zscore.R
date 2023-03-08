
# Import the libraries
library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(TMB)

# Load in the processed data
dat <- read_rds(file = "../../data/processed/dat.rds")

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

# Dynamically link the C++ template
dyn.load(dynlib(name = "../models/PoissonLognormal"))
# Load the Poisson-normal model
PoisLN <- read_rds(file = "../models/PoissonLognormal.rds")
# ... and generate report
rep <- sdreport(PoisLN, getJointPrecision = TRUE)
rep$par.fixed[12]

PoisLN_res <- y %>%
  mutate(u = rep$par.random,
         log_lambda = rep(rep$par.fixed[1:11], times = 180),
         log_sigma = rep(rep$par.fixed[12:22], times = 180),
         p = pnorm(q = u, mean = 0, sd = exp(log_sigma)),
         alarm = p >= 0.95)

write_rds(x = PoisLN_res, file = "PoisLN_res.rds")

# Load the Poisson-Gamma model
PoisG <- read_rds(file = "../models/PoissonGamma.rds")

# Join the parameters and the results in a tibble
PoisG_res <- PoisG %>%
  keep(names(.) %in% c("par", "results")) %>%
  reduce(.f = inner_join, by = "ageGroup") %>%
  mutate(p = pgamma(q = u, shape = 1/beta, scale = beta), alarm = p >= 0.95) %>%
  select(Date, beta, ageGroup, y:alarm)

write_rds(x = PoisG_res, file = "PoisG_res.rds")

