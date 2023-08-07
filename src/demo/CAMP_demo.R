
# Import libraries
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

# DTU colours
dtuPalette <- c("#990000",
                "#2F3EEA",
                "#1FD082", 
                "#030F4F", 
                "#F6D04D",
                "#FC7634",
                "#F7BBB1", 
                "#DADADA", 
                "#E83F48",
                "#008835", 
                "#79238E")

# Set a nice theme
theme_set(
  new = theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 14),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22))
)

# Source the estimation function
source(file = "../models/aeddo.R")

# Load in the data
CAMP <- read_rds(file = "../../data/processed/CAMP.rds") %>%
  filter(Date >= "2015-01-01" & Date <= "2019-12-01") %>%
  rename(y = cases) %>%
  select(Date, ageGroup, y, n)

# Create a stacked bar graph for the number of CAMP cases
p <- CAMP %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
p

# Setup hyperparameters
k <- 36
sig.level <- 0.9
excludePastOutbreaks <- TRUE

# What months are used as reference?
p + annotate(geom = "rect", xmin = CAMP$Date[1]-10, xmax = CAMP$Date[1] %m+% months(k), ymin = -Inf, ymax = Inf, alpha = 0.2)

# Try formulas
# formula <- y ~ -1 + ageGroup ; trend <- FALSE; seasonality <- FALSE; theta <- rep(0,7) ; lower <- c(rep(1e-6, 6), -6); upper <- rep(1e2, 7)
# formula <- y ~ -1 + ageGroup + t ; trend <- TRUE; seasonality <- FALSE; theta <- rep(0,8) ; lower <- c(rep(1e-6, 7), -6); upper <- rep(1e2, 8)
# formula <- y ~ -1 + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear); trend <- FALSE; seasonality <- TRUE; theta <- rep(0,9) ; lower <- c(rep(1e-6, 8), -6); upper <- rep(1e2, 9)
formula <- y ~ -1 + ageGroup + t + sin(pi/6*periodInYear) + cos(pi/6*periodInYear); trend <- TRUE; seasonality <- TRUE; theta <- rep(0,10) ; lower <- c(rep(1e-6, 2), rep(-6,3)); upper <- rep(1e2, 10)

# Try modelling framework
model <- "PoissonNormal"
# model <- "PoissonGamma"

# Estimate the parameters and infer the one-step ahead random effects
results <- aeddo(data = CAMP,
                 formula = formula,
                 trend = trend,
                 seasonality = seasonality,
                 theta = theta,
                 method = "L-BFGS-B",
                 lower = lower,
                 upper = upper,
                 model = model,
                 k = k,
                 sig.level = sig.level,
                 cpp.dir = "../models/",
                 excludePastOutbreaks = excludePastOutbreaks)


# Unnest the results and add the corresponding threshold
if (model == "PoissonNormal"){
  unnest_results <- results %>%
    select(ran.ef) %>%
    unnest(ran.ef) %>%
    mutate(threshold = qnorm(p = sig.level, mean = 0, sd = exp(log_sigma))) 
} else if (model == "PoissonGamma"){
  unnest_results %>%
    select(ran.ef) %>%
    unnest(ran.ef) %>%
    mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi))
}

# Make a nice ggplot of the resulting one-step ahead random effects and thresholds
unnest_results %>%
  ggplot(mapping = aes(x = ref.date, colour = ageGroup)) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_line(mapping = aes(y = threshold), lty = "dashed") +
  facet_wrap(facets = vars(ageGroup), ncol = 2) +
  scale_colour_manual(values = dtuPalette) +
  scale_y_continuous(name = "Estimated one-step ahead random effects") +
  scale_x_date(name = "Month", date_breaks = "years", date_labels = "%Y") +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

