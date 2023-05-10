
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source(file = "aeddo.R")

# Import the data
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6 agegroups

y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

y.season <- y %>%
  mutate(monthInYear = as.integer(format(Date, "%m")))

formula <- y ~ -1 + ageGroup + sin(2*pi/12*monthInYear) + cos(2*pi/12*monthInYear)
season.model <- model.matrix(object = formula, data = y.season)
theta <- rep(1,ncol(season.model)+1)

PoissonGamma <- aeddo(data = y.season,
                      formula = formula,
                      theta = theta,
                      model = "PoissonGamma")

PoissonGamm2a <- aeddo(data = y.season,
                      formula = formula,
                      theta = theta,
                      method = "L-BFGS-B",
                      lower = c(rep(1e-6,6),rep(-Inf,2),-6),
                      model = "PoissonGamma")

age.formula <- y ~ -1 + ageGroup
age.model <- model.matrix(object = age.formula, data = y.season)
theta <- rep(1,ncol(age.model)+1)
PoissonGamma3 <- aeddo(data = y.season,
                       formula = age.formula,
                       theta = theta,
                       method = "L-BFGS-B",
                       lower = c(rep(1e-6,6),-6),
                       model = "PoissonGamma")

PoissonGamma %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  filter(grepl(pattern = "sin|cos", x = Parameter)) %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(Parameter))



PoissonGamma3 %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  filter(grepl(pattern = "log_phi", x = Parameter)) %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(Parameter))
?MakeADFun


PoissonGamm2a %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, u, phi, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, ref.date, NA)) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup))
