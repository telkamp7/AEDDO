
# Import libraries
library(readr)
library(dplyr)
library(TMB)
library(ggplot2)

# Import the data
dat <- read_rds(file = "../../data/processed/dat.rds")

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "SYPH") %>%
  group_by(year,maaned, ageLabel) %>%
  summarize(y = sum(cases))

compile(file = "PoissonGamma.cpp")  # Compile the C++ file
dyn.load(dynlib("PoissonGamma"))    # Dynamically link the C++ code

# Function and derivative
f <- MakeADFun(
  data = list(y = y$y),
  parameters = list(u = rep(1, length(y$y)),
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
                    lambda = 1,
                    log_sigma_u = 1),
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


