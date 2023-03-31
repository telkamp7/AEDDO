
# Import libraries
library(readr)
library(dplyr)
library(numDeriv)

# Source the negative log-likelihood functions
source(file = "nll_PoisG.R")

# Load in the processed data
# dat <- read_rds(file = "../../data/processed/dat.rds") # 11-agegroups
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6-agegroups

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date) %>%
  reframe(y = sum(cases), n = sum(n))

# Optimize the parameters
opt <- nlminb(start = c(1,1), objective = nll, x = y$y, lower = c(0,0))

# Extract the optimized parameters
lambda <- opt$par[1]
beta <- opt$par[2]

# Reconstruct r and p
r <- 1/beta
p <- 1/(lambda*beta+1)

# Calculate the expectation from the geometric distribution
PoisG.mu <- r*(1-p)/p

u <- (y$y * beta + 1)/(lambda*beta+1)

plot(u)

# Agegroups 11 ---------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

formula <- y ~ -1 + ageGroup

# Optimize the parameters
opt.age <- nlminb(start = c(rep(14,nlevels(y.age$ageGroup)),rep(0.5,nlevels(y.age$ageGroup))),
                  objective = nll.age11,
                  data = y.age,
                  formula = formula,
                  lower = rep(1e-6,nlevels(y.age$ageGroup)*2))

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Define parameters
par.PoisG <- tibble(ageGroup = levels(y.age$ageGroup),
                    beta = opt.age$par[1:n.ageGroup],
                    phi = opt.age$par[(n.ageGroup+1):(n.ageGroup*2)])

# Construct parameters and calculate the expectation for the negative binomial distribution
# par.PoisG %>%
#   # Construct the size and probability
#   mutate(r = 1/phi, p = 1/( beta*phi + 1 )) %>%
#   # Calculate the expectation
#   reframe(r, p, mu = r*(1-p)/p)

# library(ggplot2)
# y.age %>%
#   mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 ) / ( lambda*phi+1 ) ) %>%
#   ggplot(mapping = aes(x = Date, y = u)) +
#   geom_point() +
#   facet_wrap(facets = vars(ageGroup))

# Make inference on individual random effects
y.age <- y.age %>%
  full_join(y = par.PoisG, join_by(ageGroup)) %>%
  mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 )/( lambda * phi + 1 )) %>%
  select(Date:n, u)

# y.age %>%
#   group_by(ageGroup) %>%
#   summarize(mean(y), mean(u))

# Save the results in a list
PoisG <- list(par = par.PoisG, fn = nll.age, results = y.age)

# Write the results
write_rds(x = PoisG, file = "PoissonGamma.rds")


# 6 agegroups -----------------------------------------------------------------------

# Only consider some of the data
y.age <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

formula <- y ~ -1 + ageGroup

# Optimize the parameters
opt.age <- nlminb(start = c(rep(14,nlevels(y.age$ageGroup)),rep(0.5,1)),
                  objective = nll.age6,
                  data = y.age,
                  formula = formula,
                  lower = rep(1e-6,nlevels(y.age$ageGroup)+1))

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Calculate uncertainty directly in R
H <- hessian(nll.age6, opt.age$par, data = y.age, formula = formula)
se.theta <- sqrt(diag(solve(H)))

# Construct parameter table
par.tbl.PoisG <- tibble(Parameter = c(paste0("$\\beta_{", levels(y.age$ageGroup), "}$"), "$\\phi$"),
                        Estimate = opt.age$par,
                        `Std. Error` = se.theta)

# Extract number of agegroups
n.ageGroup <- n_distinct(y.age$ageGroup)

# Define parameters
par.inf.PoisG <- tibble(ageGroup = levels(y.age$ageGroup),
                    beta = opt.age$par[1:n.ageGroup],
                    phi = opt.age$par[-c(1:n.ageGroup)])

# Construct parameters and calculate the expectation for the negative binomial distribution
# par.PoisG %>%
#   # Construct the size and probability
#   mutate(r = 1/phi, p = 1/( beta*phi + 1 )) %>%
#   # Calculate the expectation
#   reframe(r, p, mu = r*(1-p)/p)

# library(ggplot2)
# y.age %>%
#   mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 ) / ( lambda*phi+1 ) ) %>%
#   ggplot(mapping = aes(x = Date, y = u)) +
#   geom_point() +
#   facet_wrap(facets = vars(ageGroup))

# Make inference on individual random effects
results <- y.age %>%
  full_join(y = par.inf.PoisG, join_by(ageGroup)) %>%
  mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 )/( lambda * phi + 1 )) %>%
  select(Date:n, u)

# y.age %>%
#   group_by(ageGroup) %>%
#   summarize(mean(y), mean(u))

# Save the results in a list
PoisG <- list(par.tbl = par.tbl.PoisG,
              par.inf = par.inf.PoisG,
              fn = nll.age6,
              results = results)

# Write the results
write_rds(x = PoisG, file = "PoissonGamma.rds")



# Windowed estimation

# Make months into integers
y.design <- y.age %>%
  mutate(Month = as.integer(format(Date, "%m")))

k <- 36

Dates <- y.design %>%
  reframe(Date = unique(Date))

TT <- length(Dates$Date)

results.window <- tibble()

# Construct initial parameters
beta <- rep(14, nlevels(y.design$ageGroup))
phi <- 0.5

for(i in 1:(TT-k)){
  
  monthsInWindow <- Dates[i:(i+k-1),] 
  
  y.window <- y.design %>%
    filter(Date %in% monthsInWindow$Date)
  
  formula <- y ~ -1 + ageGroup
  
  # Optimize the parameters
  opt.age.window <- nlminb(start = c(beta, phi),
                    objective = nll.age6.window,
                    data = y.window,
                    formula = formula,
                    lower = rep(1e-6,nlevels(y.window$ageGroup)+1))
  
  # Extract number of agegroups
  n.ageGroup <- n_distinct(y.window$ageGroup)
  
  # Calculate uncertainty directly in R
  H <- hessian(nll.age6.window, opt.age.window$par, data = y.window, formula = formula)
  se.theta <- sqrt(diag(solve(H)))
  
  # Construct parameter table
  par.tbl.PoisG <- tibble(Parameter = c(paste0("$\\beta_{", levels(y.window$ageGroup), "}$"), "$\\phi$"),
                          Estimate = opt.age.window$par,
                          `Std. Error` = se.theta)
  
  # Extract number of agegroups
  n.ageGroup <- n_distinct(y.window$ageGroup)
  
  # Define parameters
  par.inf.PoisG <- tibble(ageGroup = levels(y.window$ageGroup),
                          beta = opt.age.window$par[1:n.ageGroup],
                          phi = opt.age.window$par[-c(1:n.ageGroup)])
  
  # Make inference on individual random effects
  results <- y.window %>%
    full_join(y = par.inf.PoisG, join_by(ageGroup)) %>%
    mutate(lambda = exp(beta - log(n)), u = ( y*phi+1 )/( lambda * phi + 1 )) %>%
    select(Date:n, u)
  
  # Combine the results in a tibble
  results.window <- bind_rows(results.window,
                              tibble(iter = i,
                                     par = list(opt.age.window$par),
                                     objective = opt.age.window$objective,
                                     `Random Effects` = list(results$u)))
  
  # Construct parameter guess for next iteration
  beta <- opt.age.window$par[1:nlevels(y.window$ageGroup)]
  phi <- opt.age.window$par[nlevels(y.window$ageGroup)+1]
  
}

write_rds(x = results.window, file = "windowedPoissonGamma.rds")
