
# Import libraries
library(TMB)
library(microbenchmark)
library(ggplot2)

# Compile and load the models
compile(file = "PoissonNormal.cpp")
compile(file = "PoissonGamma.cpp")
dyn.load(dynlib(name = "PoissonNormal"))
dyn.load(dynlib(name = "PoissonGamma"))

# Set a seed
set.seed(42)

# Set number of obs
nObs <- 1e3
theta <- rep(1,2)

sim <- data.frame(y = rpois(n = nObs, lambda = 5), n = 1)

# Construct the design matrix
designMatrix <- model.matrix(object = y ~ 1, data = sim)

# Make the functions
PoisN <- MakeADFun(
  data = list(y = sim$y, n = sim$n, X = designMatrix),
  parameters = list(u = rep(1, length(sim$y)),
                    beta = theta[1:ncol(designMatrix)],
                    log_sigma_u = theta[-(1:ncol(designMatrix))]),
  hessian = TRUE,
  random = "u",
  method = "BFGS",
  DLL = "PoissonNormal",
  silent = TRUE
)
PoisG <- MakeADFun(
  data = list(y = sim$y, n = sim$n, X = designMatrix),
  parameters = list(u = rep(1, length(sim$y)),
                    beta = theta[1:ncol(designMatrix)],
                    log_phi_u = theta[-(1:ncol(designMatrix))]),
  hessian = TRUE,
  method = "BFGS",
  DLL = "PoissonGamma",
  silent = TRUE
)

mbm <- microbenchmark(
  "Poisson Normal" = {
    do.call(what = "optim", args = PoisN)
  },
  "Poisson Gamma" = {
    do.call(what = "optim", args = PoisG)
  }
)

p <- autoplot(mbm)


Compare_Performance <- p +
  stat_ydensity(mapping = aes(fill = expr)) +
  scale_fill_manual(name = "Method", values = c("#990000", "#2F3EEA")) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(filename = "Compare_Performance.png",
       plot = Compare_Performance,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 12,
       units = "in",
       dpi = "print")
