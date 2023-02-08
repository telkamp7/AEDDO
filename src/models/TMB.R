
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
  summarize(y = sum(cases))

compile(file = "PoissonLognormal.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonLognormal"))    # Dynamically link the C++ code

# Function and derivative
PoisLN <- MakeADFun(
  data = list(y = y$y, ageGroup = y$ageGroup),
  parameters = list(u = rep(1, length(y$y)),
                    lambda = rep(1, nlevels(y$ageGroup)),
                    log_sigma_u = log(1)),
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
sdreport(PoisLN)

f$par

rap <- sdreport(f,getJointPrecision = TRUE)
summary(rap,"random")
rap$par.random
rap$diag.cov.random
names(rap)

compile(file = "PoissonGamma.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonGamma"))    # Dynamically link the C++ code

# Function and derivative
f <- MakeADFun(
  data = list(y = y$y),
  parameters = list(u = rep(0, length(y$y)),
                    lambda = 1,
                    phi = 1),
  random = "u",
  DLL = "PoissonGamma"
)

f$fn()
f$gr()

opt <- nlminb(start = f$par, f$fn, f$gr, lower = c(0.01, 0.01, 0.01))
TMB:::op_table(f$env$ADFun)
TMB:::tape_print(f$env$ADFun)

