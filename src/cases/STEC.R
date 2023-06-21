
# Import libraries
library(dplyr)
library(readr)
library(tidyr)
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
dat <- read_rds(file = "../../data/processed/dat2.rds")

# Summary statistic of all the data
dat %>%
  group_by(Date, caseDef) %>%
  reframe(y = sum(cases)) %>% 
  group_by(caseDef) %>%
  summarise(meanCases = mean(y), medianCases = median(y))


# Only consider the STEC cases
STEC <- dat %>%
  filter(caseDef == "STEC") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) 

# Bar plot
STEC_long_plot <- STEC %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_long_plot.png",
       plot = STEC_long_plot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# STEC %>%
#   ggplot(mapping = aes(x = Date, y = y/n*1e5, fill = ageGroup, group = ageGroup)) +
#     geom_col() +
#     facet_wrap(facets = vars(ageGroup)) +
#     scale_y_continuous(name = "Incidence per 100.000") +
#     scale_x_date(name = "Date") +
#     scale_fill_manual(values = dtuPalette) +
#     guides(fill = "none")

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- STEC %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- STEC %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
STEC.sts <- sts(
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
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "delta"
)

STEC_Farrington <- farringtonFlexible(sts = STEC.sts, con.farrington)

upperbound_Farrington <- as_tibble(STEC_Farrington@upperbound) %>%
  mutate(Date = as.Date(x = STEC_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

alarm_Farrington <- as_tibble(STEC_Farrington@alarm) %>%
  mutate(Date = as.Date(x = STEC_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 2,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
)

STEC_Noufaily <- farringtonFlexible(sts = STEC.sts, con.noufaily)

upperbound_Noufaily <- as_tibble(STEC_Noufaily@upperbound) %>%
  mutate(Date = as.Date(x = STEC_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

alarm_Noufaily <- as_tibble(STEC_Noufaily@alarm) %>%
  mutate(Date = as.Date(x = STEC_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))


# Compare the Farrington and Noufaily method
Compare_stateOfTheArt_STEC_dat <- STEC %>%
  left_join(upperbound_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(upperbound_Noufaily, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = c(threshold_Farrington,threshold_Noufaily), names_to = "method", names_prefix = "threshold_", values_to = "threshold") %>%
  pivot_longer(cols = c(alarm_Farrington,alarm_Noufaily), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  filter(method == method2) %>%
  select(-method2) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))
write_rds(x = Compare_stateOfTheArt_STEC_dat, file = "Compare_stateOfTheArt_STEC_dat.rds")

Compare_stateOfTheArt_STEC <- Compare_stateOfTheArt_STEC_dat %>%
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
ggsave(filename = "Compare_stateOfTheArt_STEC.png",
       plot = Compare_stateOfTheArt_STEC,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Hierarchical Poisson Normal model -------------------------------------------------

STEC_PoisN_ageGroup <- aeddo(data = STEC,
                    formula = y ~ -1 + ageGroup,
                    theta = rep(1,7),
                    method = "L-BFGS-B",
                    lower = c(rep(1e-6,6), -6),
                    upper = rep(1e2, 7),
                    model = "PoissonNormal",
                    k = 36,
                    sig.level = 0.9,
                    cpp.dir = "../models/",
                    CI = TRUE,
                    excludePastOutbreaks = TRUE)

start.theta.PoisN <- STEC_PoisN_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

write_rds(x = STEC_PoisN_ageGroup, file = "STEC_PoisN_ageGroup.rds")
# STEC_PoisN_ageGroup <- read_rds(file = "STEC_PoisN_ageGroup.rds")

STEC_PoisN_ageGroup_trend <- aeddo(data = STEC,
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

write_rds(x = STEC_PoisN_ageGroup_trend, file = "STEC_PoisN_ageGroup_trend.rds")
# STEC_PoisN_ageGroup_trend <- read_rds(file = "STEC_PoisN_ageGroup_trend.rds")

STEC_PoisN_ageGroup_seasonality <- aeddo(data = STEC,
                                         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                         trend = TRUE,
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7])-6,
                                         upper = c(start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7])+6,
                                         model = "PoissonNormal", 
                                         k = 36, 
                                         sig.level = 0.9,
                                         cpp.dir = "../models/",
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = STEC_PoisN_ageGroup_seasonality, file = "STEC_PoisN_ageGroup_seasonality.rds")
# STEC_PoisN_ageGroup_seasonality <- read_rds(file = "STEC_PoisN_ageGroup_seasonality.rds")


STEC_PoisN_ageGroup_trend_seasonality <- aeddo(data = STEC,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0,start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5,c(start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7])-6),
                                               upper = c(0.5,c(start.theta.PoisN[1:6], 0,0,start.theta.PoisN[7])+6),
                                               model = "PoissonNormal", 
                                               k = 36, 
                                               sig.level = 0.9,
                                               cpp.dir = "../models/",
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)


write_rds(x = STEC_PoisN_ageGroup_trend_seasonality, file = "STEC_PoisN_ageGroup_trend_seasonality.rds")
# STEC_PoisN_ageGroup_trend_seasonality <- read_rds(file = "STEC_PoisN_ageGroup_trend_seasonality.rds")


STEC_PoisN_ageGroup_unnested <- STEC_PoisN_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

STEC_PoisN_ageGroup_unnested %>%
  ggplot(mapping = aes(x = Date)) +
  geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
  geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(shape = "none")

STEC_PoisN_ageGroup_trend_unnested <- STEC_PoisN_ageGroup_trend %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

STEC_PoisN_ageGroup_trend_unnested %>%
  ggplot(mapping = aes(x = Date)) +
  geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
  geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(shape = "none")

STEC_PoisN_ageGroup_trend_seasonality_unnested <- STEC_PoisN_ageGroup_trend_seasonality %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

STEC_PoisN_ageGroup_trend_seasonality_par <- STEC_PoisN_ageGroup_trend_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

STEC_PoisN_ageGroup_trend_seasonality_unnested %>%
  ggplot(mapping = aes(x = Date)) +
  geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
  geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(shape = "none")


# Hierarchical Poisson Gamma model ----------------------------------------------------
STEC_PoisG_ageGroup <- aeddo(data = STEC,
                             formula = y ~ -1 + ageGroup,
                             theta = rep(1,7),
                             method = "L-BFGS-B",
                             lower = c(rep(1e-6,6), -6),
                             upper = rep(1e2, 7),
                             model = "PoissonGamma",
                             k = 36,
                             sig.level = 0.9,
                             cpp.dir = "../models/",
                             CI = TRUE,
                             excludePastOutbreaks = TRUE)

write_rds(x = STEC_PoisG_ageGroup, file = "STEC_PoisG_ageGroup.rds")
# STEC_PoisG_ageGroup <- read_rds(file = "STEC_PoisG_ageGroup.rds")

start.theta.PoisG <- STEC_PoisG_ageGroup %>%
  filter(row_number() == 1) %>%
  select(par) %>%
  unnest(par) %>%
  select(theta)%>%
  .$theta

STEC_PoisG_ageGroup_trend <- aeddo(data = STEC,
                                   formula = y ~ -1 + t + ageGroup,
                                   trend = TRUE,
                                   theta = c(0, start.theta.PoisG),
                                   method = "L-BFGS-B",
                                   lower = c(-0.5, start.theta.PoisG-6),
                                   upper = c(0.5, start.theta.PoisG+6),
                                   model = "PoissonGamma", 
                                   k = 36, 
                                   sig.level = 0.9,
                                   cpp.dir = "../models/",
                                   CI = TRUE,
                                   excludePastOutbreaks = TRUE)

write_rds(x = STEC_PoisG_ageGroup_trend, file = "STEC_PoisG_ageGroup_trend.rds")
# STEC_PoisG_ageGroup_trend <- read_rds(file = "STEC_PoisG_ageGroup_trend.rds")

STEC_PoisG_ageGroup_seasonality <- aeddo(data = STEC,
                                         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                         trend = TRUE,
                                         seasonality = TRUE,
                                         theta = c(start.theta.PoisG[1:6],0,0,start.theta.PoisG[7]),
                                         method = "L-BFGS-B",
                                         lower = c(start.theta.PoisG[1:6],0,0,start.theta.PoisG[7])-6,
                                         upper = c(start.theta.PoisG[1:6],0,0,start.theta.PoisG[7])+6,
                                         model = "PoissonGamma", 
                                         k = 36, 
                                         sig.level = 0.9,
                                         cpp.dir = "../models/",
                                         CI = TRUE,
                                         excludePastOutbreaks = TRUE)

write_rds(x = STEC_PoisG_ageGroup_seasonality, file = "STEC_PoisG_ageGroup_seasonality.rds")
#  STEC_PoisG_ageGroup_seasonality <- read_rds(file = "STEC_PoisG_ageGroup_seasonality.rds")


STEC_PoisG_ageGroup_trend_seasonality <- aeddo(data = STEC,
                                               formula = y ~ -1 + t + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                               trend = TRUE,
                                               seasonality = TRUE,
                                               theta = c(0,start.theta.PoisG[1:6],0,0,start.theta.PoisG[7]),
                                               method = "L-BFGS-B",
                                               lower = c(-0.5,c(start.theta.PoisG[1:6],0,0,start.theta.PoisG[7])-6),
                                               upper = c(0.5,c(start.theta.PoisG[1:6],0,0,start.theta.PoisG[7])+6),
                                               model = "PoissonGamma", 
                                               k = 36, 
                                               sig.level = 0.9,
                                               cpp.dir = "../models/",
                                               CI = TRUE,
                                               excludePastOutbreaks = TRUE)

write_rds(x = STEC_PoisG_ageGroup_trend_seasonality, file = "STEC_PoisG_ageGroup_trend_seasonality.rds")
# STEC_PoisG_ageGroup_trend_seasonality <- read_rds(file = "STEC_PoisG_ageGroup_trend_seasonality.rds")

STEC_PoisG_ageGroup_trend_seasonality_unnested <- STEC_PoisG_ageGroup_trend_seasonality %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)

STEC_PoisG_ageGroup_trend_seasonality_par <- STEC_PoisG_ageGroup_trend_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

STEC_PoisG_ageGroup %>%
  summarise(avgLogS = mean(LogS))

STEC_PoisG_ageGroup_seasonality %>%
  summarise(avgLogS = mean(LogS))

STEC_PoisG_ageGroup_unnested <- STEC_PoisG_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.95, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)

STEC_PoisN_ageGroup_tbl <- STEC_PoisN_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup")

STEC_PoisN_ageGroup_trend_tbl <- STEC_PoisN_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend")

STEC_PoisN_ageGroup_seasonality_tbl <- STEC_PoisN_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_seasonality")

STEC_PoisN_ageGroup_trend_seasonality_tbl <- STEC_PoisN_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend_seasonality")

STEC_PoisG_ageGroup_tbl <- STEC_PoisG_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup")

STEC_PoisG_ageGroup_trend_tbl <- STEC_PoisG_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend")

STEC_PoisG_ageGroup_seasonality_tbl <- STEC_PoisG_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_seasonality")

STEC_PoisG_ageGroup_trend_seasonality_tbl <- STEC_PoisG_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend_seasonality")

STEC_novel_tbl <- bind_rows(
  STEC_PoisN_ageGroup_tbl,
  STEC_PoisN_ageGroup_trend_tbl, 
  STEC_PoisN_ageGroup_seasonality_tbl,
  STEC_PoisN_ageGroup_trend_seasonality_tbl,
  STEC_PoisG_ageGroup_tbl,
  STEC_PoisG_ageGroup_trend_tbl,
  STEC_PoisG_ageGroup_seasonality_tbl,
  STEC_PoisG_ageGroup_trend_seasonality_tbl)

STEC_novel_tbl %>% print(n = 68)

write_rds(STEC_novel_tbl, file = "STEC_novel_tbl.rds")
# Compare methods -------------------------------------------------------------------

# Compare the Poisson Normal and Poisson Gamma model
STEC_novel <- full_join(STEC_PoisN_ageGroup_trend_seasonality_unnested,
                        STEC_PoisG_ageGroup_trend_seasonality_unnested,
                        by = join_by(Date, ageGroup))

Compare_novel <- STEC_novel %>%
  pivot_longer(cols = c(`u_Poisson Normal`, `u_Poisson Gamma`), names_to = "method", names_prefix = "u_", values_to = "u") %>%
  pivot_longer(cols = c(`alarm_Poisson Normal`, `alarm_Poisson Gamma`), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  pivot_longer(cols = c(`threshold_Poisson Normal`, `threshold_Poisson Gamma`), names_to = "method3", names_prefix = "threshold_", values_to = "threshold") %>%
  filter(method == method2 & method == method3) %>%
  select(-method2, -method3) %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(y = threshold, group = method), lty = "dashed") +
  facet_grid(rows = vars(ageGroup), cols = vars(method), scales = "free_y") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")
ggsave(filename = "Compare_novel_STEC.png",
       plot = Compare_novel,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  


custom_labeller <- as_labeller(
  c(`ageGroup<1 year`="beta[1~year]", `ageGroup1-4 years`="beta[1-4~years]",
    `ageGroup5-14 years`="beta[5-14~years]",`ageGroup15-24 years`="beta[15-24~years]",
    `ageGroup25-64 years`="beta[25-64~years]", `ageGroup65+ years`="beta[65+~years]",
    `t`="beta[trend]", `sin(pi/6 * monthInYear)` ="beta[sin]",
    `cos(pi/6 * monthInYear)`="beta[cos]", `log_phi`="phi", `log_sigma`="sigma"),
  default = label_parsed
)

STEC_novel_par <- bind_rows(STEC_PoisG_ageGroup_trend_seasonality_par, STEC_PoisN_ageGroup_trend_seasonality_par) 


STEC_novel_par_ageGroup <- STEC_novel_par %>%
  filter(grepl(x = Parameter, pattern = "ageGroup")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_novel_par_ageGroup.png",
       plot = STEC_novel_par_ageGroup,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

STEC_novel_par_trend <- STEC_novel_par %>%
  filter(Parameter == "t") %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_novel_par_trend.png",
       plot = STEC_novel_par_trend,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

STEC_novel_par_seasonality <- STEC_novel_par %>%
  filter(grepl(x = Parameter, pattern = "cos|sin")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_novel_par_seasonality.png",
       plot = STEC_novel_par_seasonality,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 


STEC_novel_par_dispersion <- STEC_novel_par %>%
  filter(grepl(x = Parameter, pattern = "log_sigma|log_phi")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = exp(theta), colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = exp(CI.lwr), colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = exp(CI.upr), colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(Psi)) +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_novel_par_dispersion.png",
       plot = STEC_novel_par_dispersion,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  



# Outbreaks invstigated by SSI ------------------------------------------------------

SSI_outbreaks <- tibble(Start = as.Date(x = c("2007-02-5","2012-09-15", "2018-09-03", "2019-05-06", "2021-12-03")),
                        End = as.Date(x = c("2007-03-31","2012-10-15", "2018-12-02", "2019-06-07", "2022-01-06")))

STEC_SSI_outbreaks <- SSI_outbreaks %>%
  filter(Start > as.Date("2008-01-01")) %>%
  arrange(desc(Start)) %>%
  mutate(outbreak_no = row_number()) %>%
  ggplot() +
  geom_segment(mapping = aes(x = Start, xend = End, y = outbreak_no, yend = outbreak_no), linewidth = 1.2, colour = dtuPalette[6]) +
  geom_point(mapping = aes(x = Start, y = outbreak_no), pch = 17, size = 3, colour = dtuPalette[6]) +
  scale_x_date(name = "Date", limits = c(as.Date(c("2008-01-01", "2022-12-01")))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(filename = "STEC_SSI_outbreaks.png",
       plot = STEC_SSI_outbreaks,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  


# Compare alarms across all the models
SSI_corrected <- SSI_outbreaks %>%
  mutate(method = "SSI", alarm = TRUE) %>%
  select(Date = Start, method:alarm)

STEC_compare <- STEC_novel %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "STEC")

write_rds(x = STEC_compare, file = "STEC_compare.rds")

Compare_alarms <- STEC_compare %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,6)]) +
  scale_y_discrete(limits = rev(levels(STEC_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_STEC.png",
       plot = Compare_alarms,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

