
# Import libraries
library(readr)
library(dplyr)
library(TMB)

# Import the data
dat <- read_rds(file = "../../data/processed/dat.rds")

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

compile(file = "PoissonLognormal.cpp")  # Compile the C++ file
# dyn.unload(dynlib("PoissonLognormal"))
dyn.load(dynlib("PoissonLognormal"))    # Dynamically link the C++ code

# Function and derivative
PoisLN <- MakeADFun(
  data = list(y = y$y, ageGroup = y$ageGroup, w = y$n),
  parameters = list(u = rep(1, length(y$y)),
                    log_lambda = rep(log(1), nlevels(y$ageGroup)),
                    log_sigma_u = rep(log(1), nlevels(y$ageGroup))),
  random = "u",
  DLL = "PoissonLognormal"
)

opt <- nlminb(start = PoisLN$par, PoisLN$fn, PoisLN$gr, lower = c(0.01, 0.01))
PoisLN$fn()
PoisLN$gr()


write_rds(x = PoisLN, file = "PoissonLognormal.rds")

opt$par
opt$objective

## report on result
tmp <- sdreport(PoisLN)
tmp$par.random
