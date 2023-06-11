# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(surveillance)
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
                         
# Set global theme options
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
dat <- read_rds(file = "../../data/processed/dat5.rds")

# Summary statistic of all the data
dat %>%
  group_by(Date, caseDef) %>%
  reframe(y = sum(cases)) %>% 
  group_by(caseDef) %>%
  summarise(meanCases = mean(y), medianCases = median(y))

# Only consider the STEC cases
SHIG <- dat %>%
  filter(caseDef == "SHIG") %>%
  mutate(ageGroup = fct_collapse(ageGroup, `<25 years` = c("<1 year", "1-4 years","5-14 years","15-24 years")),
         ageGroup = fct_collapse(ageGroup, `25+ years` = c("25-64 years","65+ years"))) %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

SHIG_long_plot <- SHIG %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
ggsave(filename = "SHIG_long_plot.png",
       plot = SHIG_long_plot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- SHIG %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- SHIG %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
SHIG.sts <- sts(
  observed = observed[,-1],
  epoch = observed$Date,
  epochAsDate = TRUE,
  frequency = 12,
  population = as.matrix(population[,-1])
)


# Farrington ------------------------------------------------------------------------

con.farrington <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
  thersholdMethod = "delta"
)

SHIG_Farrington <- farringtonFlexible(sts = SHIG.sts, con.farrington)

upperbound_Farrington <- as_tibble(SHIG_Farrington@upperbound) %>%
  mutate(Date = as.Date(x = SHIG_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIG$ageGroup)))

alarm_Farrington <- as_tibble(SHIG_Farrington@alarm) %>%
  mutate(Date = as.Date(x = SHIG_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIG$ageGroup)))

# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
  thersholdMethod = "Noufaily"
)

SHIG_Noufaily <- farringtonFlexible(sts = SHIG.sts, con.noufaily)

upperbound_Noufaily <- as_tibble(SHIG_Noufaily@upperbound) %>%
  mutate(Date = as.Date(x = SHIG_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIG$ageGroup)))

alarm_Noufaily <- as_tibble(SHIG_Noufaily@alarm) %>%
  mutate(Date = as.Date(x = SHIG_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIG$ageGroup)))


# Compare the Farrington and Noufaily method

Compare_stateOfTheArt_SHIG_dat <- SHIG %>%
  left_join(upperbound_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(upperbound_Noufaily, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = c(threshold_Farrington,threshold_Noufaily), names_to = "method", names_prefix = "threshold_", values_to = "threshold") %>%
  pivot_longer(cols = c(alarm_Farrington,alarm_Noufaily), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  filter(method == method2) %>%
  select(-method2) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))
write_rds(x = Compare_stateOfTheArt_SHIG_dat, file = "Compare_stateOfTheArt_SHIG_dat.rds")

Compare_stateOfTheArt_SHIG <- Compare_stateOfTheArt_SHIG_dat %>%
  ggplot(mapping = aes(x = Date, fill = ageGroup)) +
  geom_col(mapping = aes(y = y/n*1e5, alpha = alarm), linewidth = 0.4) +
  geom_line(mapping = aes(x = Date, y = threshold/n*1e5), lty = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
  facet_grid(rows = vars(ageGroup), cols = vars(method), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000") +
  scale_x_date(name = "Month") +
  scale_fill_manual(values = dtuPalette) +
  scale_alpha_manual(values = c(0.3, 1)) +
  guides(fill = "none", alpha = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = as.Date("2011-03-01"), ymin = -Inf, ymax = Inf, alpha = 0.2)
ggsave(filename = "Compare_stateOfTheArt_SHIG.png",
       plot = Compare_stateOfTheArt_SHIG,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Outbreaks invstigated by SSI ------------------------------------------------------

SSI_outbreaks <- tibble(Start = as.Date(x = c("2007-08-06", "2009-04-01", "2020-08-22")),
                        End = as.Date(x = c("2007-08-24", "2009-06-01", "2020-09-9")))

SHIG_SSI_outbreaks <- SSI_outbreaks %>%
  filter(Start > as.Date("2008-01-01")) %>%
  arrange(desc(Start)) %>%
  mutate(outbreak_no = row_number()) %>%
  ggplot() +
  geom_segment(mapping = aes(x = Start, xend = End, y = outbreak_no, yend = outbreak_no), size = 1.2, colour = dtuPalette[2]) +
  geom_point(mapping = aes(x = Start, y = outbreak_no), pch = 17, size = 3,colour = dtuPalette[4]) +
  scale_x_date(name = "Date", limits = c(as.Date(c("2008-01-01", "2022-12-01")))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(filename = "SHIG_SSI_outbreaks.png",
       plot = SHIG_SSI_outbreaks,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  


# Hierarchical Poisson Normal model ---------------------------------------

SHIG_PoisN_ageGroup <- aeddo(data = SHIG,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(0,3),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,2), -6),
                             upper = rep(1e2, 3),
                             model = "PoissonNormal",
                             k = 36,
                             sig.level = 0.9,
                             cpp.dir = "../models/",
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisN_ageGroup, file = "SHIG_PoisN_ageGroup.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisN_ageGroup.rds")

start.theta.PoisN <- SHIG_PoisN_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SHIG_PoisN_ageGroup_trend <- aeddo(data = SHIG,
                                   formula = y ~ -1 + t + ageGroup,
                                   trend = TRUE,
                                   theta = c(0, start.theta.PoisN),
                                   method = "L-BFGS-B",
                                   lower = c(-0.5, start.theta.PoisN-6),
                                   upper = c(0.5, start.theta.PoisN+6),
                                   model = "PoissonNormal", 
                                   k = 36, 
                                   sig.level = 0.9,
                                   cpp.dir = "../models/",
                                   CI = TRUE,
                                   excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisN_ageGroup_trend, file = "SHIG_PoisN_ageGroup_trend.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisN_ageGroup_trend.rds")

SHIG_PoisN_ageGroup_seasonality <- aeddo(data = SHIG,
                                         formula = y ~ -1 +  ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) - 6,
                                         upper = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) + 6,
                                         model = "PoissonNormal", 
                                         k = 36, 
                                         sig.level = 0.9,
                                         cpp.dir = "../models/",
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisN_ageGroup_seasonality, file = "SHIG_PoisN_ageGroup_seasonality.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisN_ageGroup_seasonality.rds")

SHIG_PoisN_ageGroup_trend_seasonality <- aeddo(data = SHIG,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])-6),
                                               upper = c(0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])+6),
                                               model = "PoissonNormal", 
                                               k = 36, 
                                               sig.level = 0.9,
                                               cpp.dir = "../models/",
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisN_ageGroup_trend_seasonality, file = "SHIG_PoisN_ageGroup_trend_seasonality.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisN_ageGroup_trend_seasonality.rds")


# Hierarchical Poisson Gamma model ----------------------------------------




SHIG_PoisG_ageGroup <- aeddo(data = SHIG,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(1,3),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,2),-6),
                             upper = rep(1e2, 3),
                             model = "PoissonGamma",
                             k = 36,
                             cpp.dir = "../models/",
                             sig.level = 0.9,
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisG_ageGroup, file = "SHIG_PoisG_ageGroup.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisG_ageGroup.rds")

start.theta.PoisG <- SHIG_PoisG_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

SHIG_PoisG_ageGroup_trend <- aeddo(data = SHIG,
                                   formula = y ~ -1 + t + ageGroup,
                                   trend = TRUE,
                                   theta = c(0, start.theta.PoisG),
                                   method = "L-BFGS-B",
                                   lower = c(-0.5, start.theta.PoisG-6),
                                   upper = c(0.5, start.theta.PoisG+6),
                                   model = "PoissonGamma",
                                   k = 36,
                                   cpp.dir = "../models/",
                                   sig.level = 0.9,
                                   CI = TRUE,
                                   excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisG_ageGroup_trend, file = "SHIG_PoisG_ageGroup_trend.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisG_ageGroup_trend.rds")


SHIG_PoisG_ageGroup_seasonality <- aeddo(data = SHIG,
                                         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6,
                                         upper = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6,
                                         model = "PoissonGamma",
                                         k = 36,
                                         cpp.dir = "../models/",
                                         sig.level = 0.9,
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisG_ageGroup_seasonality, file = "SHIG_PoisG_ageGroup_seasonality.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisG_ageGroup_trend.rds")

SHIG_PoisG_ageGroup_trend_seasonality <- aeddo(data = SHIG,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0, start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6),
                                               upper = c(0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6),
                                               model = "PoissonGamma",
                                               k = 36,
                                               cpp.dir = "../models/",
                                               sig.level = 0.9,
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = SHIG_PoisG_ageGroup_trend_seasonality, file = "SHIG_PoisG_ageGroup_trend_seasonality.rds")
# SHIG_PoisN_ageGroup <- read_rds(file = "SHIG_PoisG_ageGroup_trend_seasonality.rds")

