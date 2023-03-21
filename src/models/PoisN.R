
# Import libraries
library(readr)
library(dplyr)
library(TMB)

# Import the data
# dat <- read_rds(file = "../../data/processed/dat.rds")
dat <- read_rds(file = "../../data/processed/dat2.rds")

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

# designMatrix %*% beta

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

unique(y.design$Date)

# Windowed estimation
k <- 36

Dates <- y.design %>%
  reframe(Date = unique(Date))

TT <- length(Dates$Date)

results.window <- tibble()

# Construct initial parameters
beta <- rep(1, nlevels(y.window$ageGroup))
sigma <- log(1)

for(i in 1:(TT-k)){
  
  monthsInWindow <- Dates[i:(i+k-1),] 
  
  y.window <- y.design %>%
    filter(Date %in% monthsInWindow$Date)
  
  # Construct the design matrix
  designMatrixWindow <- model.matrix(y ~ -1 + ageGroup, data = y.window)
  
  
  
  # designMatrix %*% beta
  
  # Function and derivative
  PoisLN <- MakeADFun(
    data = list(y = y.window$y, x = y.window$n, X = designMatrixWindow),
    parameters = list(u = rep(1, length(y.window$y)),
                      beta = beta,
                      log_sigma_u = sigma),
    random = "u",
    DLL = "PoissonNormal",
    silent = TRUE
  )
  
  opt <- nlminb(start = PoisLN$par, PoisLN$fn, PoisLN$gr, lower = rep(1e-9,nlevels(y.window$ageGroup)+1))
  
  
  results.window <- bind_rows(results.window,
                              tibble(iter = i+k-1,
                                     optimization = list(opt)))
  
  
  # Construct parameter guess for next iteration
  beta <- opt$par[1:nlevels(y.window$ageGroup)]
  sigma <- opt$par[nlevels(y.window$ageGroup)+1]
  
}

results.window$optimization[[1]]$par



