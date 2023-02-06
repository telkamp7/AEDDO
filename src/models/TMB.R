
# Import libraries
library(readr)
library(dplyr)
library(TMB)
library(ggplot2)

# Import the data
dat <- read_rds(file = "../../data/processed/dat.rds")

dat %>%
  group_by(caseDef, ageLabel) %>%
  summarize(sum(cases)) %>%
  print(n = 100)

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "VTEC") %>%
  group_by(year,maaned) %>%
  summarize(y = sum(cases))

compile(file = "PoissonGamma.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonGamma"))    # Dynamically link the C++ code

alpha <- 1
beta <- 3
x <- rgamma(n = 10, shape = alpha, scale = beta)

dgamma(x = x, shape = alpha, rate = beta)

dgamma(x = x, shape = alpha, rate = beta, log = TRUE)

log(beta^alpha) + log(x^(alpha-1)) - (beta*x) - log(gamma(alpha))

(x^(alpha-1)*exp(-beta*x)*beta^alpha)/gamma(alpha)

# Function and derivative
f <- MakeADFun(
  data = list(y = y$y),
  parameters = list(u = rep(.5, length(y$y)),
                    lambda = 1,
                    phi = 1),
  random = "u",
  DLL = "PoissonGamma"
  )

f$fn()
f$gr()

opt <- nlminb(start = f$par, f$fn, f$gr, lower = c(0.01, 0.01, 0.01))

compile(file = "PoissonGaussian.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonGaussian"))    # Dynamically link the C++ code


# Function and derivative
f <- MakeADFun(
  data = list(y = y$y, ageLabel = y$ageLabel),
  parameters = list(u = rep(1, length(y$y)),
                    lambda = rep(1, nlevels(y$ageLabel)),
                    log_sigma_u = log(1)),
  random = "u",
  DLL = "PoissonGaussian"
)

f$fn()
f$gr()

opt <- nlminb(start = f$par, f$fn, f$gr, lower = c(0.01, 0.01))

opt$par
opt$objective

## report on result
sdreport(f)

rap <- sdreport(f,getJointPrecision = TRUE)
summary(rap,"random")
rap$par.random
rap$diag.cov.random
names(rap)


y <- y %>%
  ungroup() %>%
  mutate(parRandom = rap$par.random)

y %>%
  ggplot(mapping = aes(x = maaned, y = parRandom)) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(ageLabel))


