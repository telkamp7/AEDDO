
# Import libraries
library(readr)
library(dplyr)
library(ggplot2)
# library(TMB)

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
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22))
)

# Source the estimation function
source(file = "../models/aeddo.R")

# Load in the data
CAMP <- read_rds(file = "../../data/processed/CAMP.rds")

# Create a stacked bar graph for the number of CAMP cases
CAMP %>%
  ggplot(mapping = aes(x = Date, y = cases       , fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")


# Setup hyperparameters
k <- 36
sig.level <- 0.9
excludePastOutbreaks <- TRUE

# Estimate the parameters and infer the one-step ahead random effects
results <- aeddo(data = CAMP,
                 formula = y ~ -1 + ageGroup,
                 theta = rep(0,3),
                 method = "L-BFGS-B",
                 lower = c(rep(1e-6,2), -6),
                 upper = rep(1e2, 3),
                 model = "PoissonNormal",
                 k = k,
                 sig.level = sig.level,
                 cpp.dir = "../models/",
                 excludePastOutbreaks = excludePastOutbreaks)

