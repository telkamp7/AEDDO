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
source(file = "../src/models/aeddo.R")

# Load in the data
dat <- read_rds(file = "../data/processed/dat2.rds")


# Only consider the STEC cases
STEC <- dat %>%
  filter(caseDef == "STEC") %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(ageGroup = fct_recode(ageGroup, `0-1 year` = "<1 year"))

# Bar plot
STEC_long_plot <- STEC %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
ggsave(filename = "STEC_long_plot.png",
       plot = STEC_long_plot,
       path = "../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

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
  alpha = 0.005, trend = TRUE, pThresholdTrend = 0.05,
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
  alpha = 0.005, trend = TRUE, pThresholdTrend = 1,
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



# Poisson Normal ----------------------------------------------------------

STEC_PoisN_975 <- aeddo(
  data = STEC,
  formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
  trend = TRUE,
  seasonality = TRUE,
  theta = c(0, rep(1, 6), 0, 0, 1),
  method = "L-BFGS-B",
  lower = c(-0.5, rep(1e-6, 6), 0, 0, -6),
  upper = c(0.5, rep(1e2, 6), 0, 0, 1e2),
  model = "PoissonNormal",
  k = 36,
  sig.level = 0.975,
  cpp.dir = "../src/models/",
  excludePastOutbreaks = TRUE)

STEC_PoisN_95 <- aeddo(
  data = STEC,
  formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
  trend = TRUE,
  seasonality = TRUE,
  theta = c(0, rep(1, 6), 0, 0, 1),
  method = "L-BFGS-B",
  lower = c(-0.5, rep(1e-6, 6), 0, 0, -6),
  upper = c(0.5, rep(1e2, 6), 0, 0, 1e2),
  model = "PoissonNormal",
  k = 36,
  sig.level = 0.95,
  cpp.dir = "../src/models/",
  excludePastOutbreaks = TRUE)

# Poisson Gamma -----------------------------------------------------------

STEC_PoisG_975 <- aeddo(
  data = STEC,
  formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
  trend = TRUE,
  seasonality = TRUE,
  theta = c(0,rep(1,6),0,0,1),
  method = "L-BFGS-B",
  lower = c(-0.5,rep(1e-6,6),0,0,-6),
  upper = c(0.5,rep(1e2,6),0,0,1e2),
  model = "PoissonGamma",
  k = 36,
  sig.level = 0.975,
  cpp.dir = "../src/models/",
  excludePastOutbreaks = TRUE)

STEC_PoisG_95 <- aeddo(
  data = STEC,
  formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
  trend = TRUE,
  seasonality = TRUE,
  theta = c(0,rep(1,6),0,0,1),
  method = "L-BFGS-B",
  lower = c(-0.5,rep(1e-6,6),0,0,-6),
  upper = c(0.5,rep(1e2,6),0,0,1e2),
  model = "PoissonGamma",
  k = 36,
  sig.level = 0.95,
  cpp.dir = "../src/models/",
  excludePastOutbreaks = TRUE)


# SSI outbreaks -----------------------------------------------------------

SSI_outbreaks <- tibble(Start = as.Date(x = c("2007-02-5","2012-09-15", "2018-09-03", "2019-05-06", "2021-12-03")),
                        End = as.Date(x = c("2007-03-31","2012-10-15", "2018-12-02", "2019-06-07", "2022-01-06")))

SSI_corrected <- SSI_outbreaks %>%
  mutate(method = "SSI", alarm = TRUE) %>%
  select(Date = Start, method:alarm)

STEC_PoisN_975_unnested <- STEC_PoisN_975 %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

STEC_PoisN_95_unnested <- STEC_PoisN_95 %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)


STEC_PoisG_975_unnested <- STEC_PoisG_975 %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.95, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)

STEC_PoisG_95_unnested <- STEC_PoisG_95 %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.95, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)


STEC_novel_975 <- full_join(STEC_PoisN_975_unnested,
                        STEC_PoisG_975_unnested,
                        by = join_by(Date, ageGroup))

STEC_novel_95 <- full_join(STEC_PoisN_95_unnested,
                            STEC_PoisG_95_unnested,
                            by = join_by(Date, ageGroup))

STEC_compare_975 <- STEC_novel_975 %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "STEC")

STEC_compare_95 <- STEC_novel_95 %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "STEC")

Compare_alarms_975 <- STEC_compare_975 %>%
  filter(Date >= as.Date("2011-01-01")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,6)]) +
  scale_y_discrete(limits = rev(levels(STEC_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_975.png",
       plot = Compare_alarms_975,
       path = "../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

Compare_alarms_95 <- STEC_compare_95 %>%
  filter(Date >= as.Date("2011-01-01")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,6)]) +
  scale_y_discrete(limits = rev(levels(STEC_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_95.png",
       plot = Compare_alarms_95,
       path = "../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")
