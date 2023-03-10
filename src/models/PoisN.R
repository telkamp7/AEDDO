
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

compile(file = "PoissonNormal.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonNormal"))    # Dynamically link the C++ code

# Make months into integers
y.design <- y %>%
  mutate(Month = as.integer(format(Date, "%m")))

# Construct the design matrix
designMatrix <- model.matrix(y ~ -1 + ageGroup, data = y.design)

# Construct initial parameters
beta <- rep(1, nlevels(y$ageGroup))

designMatrix %*% beta

# Function and derivative
PoisLN <- MakeADFun(
  data = list(y = y$y, x = y$n, X = designMatrix),
  parameters = list(u = rep(1, length(y$y)),
                    beta = beta,
                    log_sigma_u = log(1)),
  random = "u",
  DLL = "PoissonNormal"
)

opt <- nlminb(start = PoisLN$par, PoisLN$fn, PoisLN$gr, lower = c(0.01, 0.01))
PoisLN$fn()
PoisLN$gr()


write_rds(x = PoisLN, file = "PoissonNormal.rds")

opt$par
opt$objective

## report on result
tmp <- sdreport(PoisLN)
tmp$par.random
