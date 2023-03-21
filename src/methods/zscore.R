
# Import the libraries
library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(TMB)

# Load in the processed data
# dat <- read_rds(file = "../../data/processed/dat.rds") # 11-agegroups
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6-agegroups

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

# Make months into integers
y.design <- y %>%
  mutate(Month = as.integer(format(Date, "%m")))

# Dynamically link the C++ template
dyn.load(dynlib(name = "../models/PoissonNormal"))
# Load the Poisson-normal model
PoisLN <- read_rds(file = "../models/PoissonNormal.rds")
# ... and generate report
rep <- sdreport(PoisLN, getJointPrecision = TRUE)

# Extract beta
beta <- rep$par.fixed[1:nlevels(y$ageGroup)]

# Construct the design matrix
designMatrix <- model.matrix(y ~ -1 + ageGroup, data = y.design)

PoisLN_res <- y %>%
  mutate(u = rep$par.random,
         log_sigma = rep$par.fixed[nlevels(y.design$ageGroup)+1],
         p = pnorm(q = u, mean = 0, sd = exp(log_sigma)),
         alarm = p >= 0.95)

write_rds(x = PoisLN_res, file = "PoisLN_res.rds")

# Load the Poisson-Gamma model
PoisG <- read_rds(file = "../models/PoissonGamma.rds")

# Join the parameters and the results in a tibble
PoisG_res <- PoisG %>%
  keep(names(.) %in% c("par.inf", "results")) %>%
  reduce(.f = inner_join, by = "ageGroup") %>%
  mutate(p = pgamma(q = u, shape = 1/phi, scale = phi), alarm = p >= 0.95) %>%
  select(Date, phi, ageGroup, y:alarm)

write_rds(x = PoisG_res, file = "PoisG_res.rds")

