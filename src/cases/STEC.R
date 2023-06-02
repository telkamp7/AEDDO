
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
dat <- read_rds(file = "../../data/processed/dat4.rds")

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
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(monthInYear = as.integer(format(Date, "%m")))

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


STEC %>%
  ggplot(mapping = aes(x = Date, y = y/n*1e5, fill = ageGroup, group = ageGroup)) +
    geom_col() +
    facet_wrap(facets = vars(ageGroup)) +
    scale_y_continuous(name = "Incidence per 100.000") +
    scale_x_date(name = "Date") +
    scale_fill_manual(values = dtuPalette) +
    guides(fill = "none")



# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- STEC %>%
  select(-n, -monthInYear) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- STEC %>%
  select(-y, -monthInYear) %>%
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
  limit54 = c(5,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
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
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(5,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWooksNotIncluded = NULL,
  thersholdMethod = "Noufaily"
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



# Hierarchical Poisson Normal model -------------------------------------------------

STEC_PoisN <- aeddo(data = STEC,
                    formula = y ~ -1 + ageGroup,
                    theta = rep(1,7),
                    method = "L-BFGS-B",
                    lower = c(rep(1e-6,6), -6),
                    upper = rep(1e2, 7),
                    model = "PoissonNormal", k = 36, cpp.dir = "../models/",
                    excludePastOutbreaks = TRUE)


STEC_PoisN_unnested <- STEC_PoisN %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)
  

# Hierarchical Poisson Gamma model ----------------------------------------------------

STEC_PoisG_ageGroup <- aeddo(data = STEC,
                    formula = y ~ -1 + ageGroup,
                    theta = rep(1,7),
                    method = "L-BFGS-B",
                    lower = c(rep(1e-6,6), -6),
                    upper = rep(1e2, 7),
                    model = "PoissonGamma", k = 36, cpp.dir = "../models/",
                    excludePastOutbreaks = TRUE)

STEC_PoisG_ageGroup %>%
  summarise(avgLogS = mean(logS))

STEC_PoisG_ageGroupSeasonality <- aeddo(data = STEC,
                                        formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
                                        theta = rep(1,9),
                                        method = "L-BFGS-B", 
                                        lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
                                        upper = rep(1e2, 9),
                                        model = "PoissonGamma", k = 36, cpp.dir = "../models/",
                                        excludePastOutbreaks = TRUE)

STEC_PoisG_ageGroupSeasonality %>%
  summarise(avgLogS = mean(logS))
  

STEC_PoisG_ageGroup_unnested <- STEC_PoisG_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.95, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)


# Compare methods -------------------------------------------------------------------

# Compare the Farrington and Noufaily method
Compare_stateOfTheArt <- STEC %>%
  left_join(upperbound_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(upperbound_Noufaily, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = c(threshold_Farrington,threshold_Noufaily), names_to = "method", names_prefix = "threshold_", values_to = "threshold") %>%
  pivot_longer(cols = c(alarm_Farrington,alarm_Noufaily), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  filter(method == method2) %>%
  select(-method2) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA)) %>%
  ggplot(mapping = aes(x = Date, fill = ageGroup)) +
  geom_col(mapping = aes(y = y/n*1e5, alpha = alarm), size = 0.4) +
  geom_line(mapping = aes(x = Date, y = threshold/n*1e5), lty = "dashed", size = 0.4, inherit.aes = FALSE) +
  facet_grid(rows = vars(ageGroup), cols = vars(method), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000") +
  scale_x_date(name = "Date") +
  scale_fill_manual(values = dtuPalette) +
  scale_alpha_manual(values = c(0.3, 1)) +
  guides(fill = "none", alpha = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = as.Date("2011-03-01"), ymin = -Inf, ymax = Inf, alpha = 0.2)
ggsave(filename = "Compare_stateOfTheArt.png",
       plot = Compare_stateOfTheArt,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Compare the Poisson Normal and Poisson Gamma model
STEC_novel <- full_join(STEC_PoisN_unnested, STEC_PoisG_ageGroup_unnested, by = join_by(Date, ageGroup))

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
ggsave(filename = "Compare_novel.png",
       plot = Compare_novel,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  



# Compare alarms across all the models
STEC_compare <- STEC_novel %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method))

Compare_alarms <- STEC_compare %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11)]) +
  scale_y_discrete(limits = rev(levels(STEC_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms.png",
       plot = Compare_alarms,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  


