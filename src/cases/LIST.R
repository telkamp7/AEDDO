
# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(forcats)
library(readr)
library(surveillance)
library(ggplot2)
library(zoo)

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

# Only consider the LIST cases
LIST <- dat %>%
  filter(caseDef == "LIST" & Date <= as.Date("2022-12-31")) %>%
  mutate(ageGroup = fct_collapse(ageGroup, `<65 years` = c("<1 year", "1-4 years", "5-14 years", "15-24 years", "25-64 years"))) %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

LIST_meanAndStandardDeviation <- LIST %>%
  group_by(ageGroup) %>%
  mutate(Mean = across(y, ~ rollapply(., 36, mean, fill = NA)),
         `Variance` = across(y, ~ rollapply(., 36, var, fill = NA))) %>%
  ungroup() %>%
  pivot_longer(cols = Mean:`Variance`, names_to = "Statistic", values_to = "value") %>%
  ggplot(mapping = aes(x = Date, y = value$y, colour = Statistic, group = Statistic)) +
  geom_line(linewidth=1.2) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_x_date(name = "Month")+
  scale_y_continuous(name = "Value") +
  scale_color_manual(name = "Statistic", values = dtuPalette[10:11])
ggsave(filename = "LIST_meanAndStandardDeviation.png",
       plot = LIST_meanAndStandardDeviation,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

LIST_long_plot <- LIST %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
ggsave(filename = "LIST_long_plot.png",
       plot = LIST_long_plot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- LIST %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- LIST %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
LIST.sts <- sts(
  observed = observed[,-1],
  epoch = observed$Date,
  epochAsDate = TRUE,
  frequency = 12,
  population = as.matrix(population[,-1])
)

# Farrington ------------------------------------------------------------------------

con.farrington <- list(
  range = NULL, b = 3, w = 3,
  reweight = TRUE, weightsThreshold = 1,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "delta"
)

LIST_Farrington <- farringtonFlexible(sts = LIST.sts, con.farrington)

upperbound_Farrington <- as_tibble(LIST_Farrington@upperbound) %>%
  mutate(Date = as.Date(x = LIST_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(LIST$ageGroup)))

alarm_Farrington <- as_tibble(LIST_Farrington@alarm) %>%
  mutate(Date = as.Date(x = LIST_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(LIST$ageGroup)))

# Noufaily --------------------------------------------------------------------------

con.noufaily <- list(
  range = NULL, b = 3, w = 3,
  reweight = TRUE, weightsThreshold = 2.58,
  verbose = TRUE, glmWarnings = TRUE,
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 10, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
)

LIST_Noufaily <- farringtonFlexible(sts = LIST.sts, con.noufaily)

upperbound_Noufaily <- as_tibble(LIST_Noufaily@upperbound) %>%
  mutate(Date = as.Date(x = LIST_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(LIST$ageGroup)))

alarm_Noufaily <- as_tibble(LIST_Noufaily@alarm) %>%
  mutate(Date = as.Date(x = LIST_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(LIST$ageGroup)))


# Compare the Farrington and Noufaily method

Compare_stateOfTheArt_LIST_dat <- LIST %>%
  left_join(upperbound_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(upperbound_Noufaily, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = c(threshold_Farrington,threshold_Noufaily), names_to = "method", names_prefix = "threshold_", values_to = "threshold") %>%
  pivot_longer(cols = c(alarm_Farrington,alarm_Noufaily), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  filter(method == method2) %>%
  select(-method2) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))
write_rds(x = Compare_stateOfTheArt_LIST_dat, file = "Compare_stateOfTheArt_LIST_dat.rds")


# Compare_stateOfTheArt_LIST_dat %>%
#   ggplot(mapping = aes(x = Date, fill = ageGroup)) +
#   geom_col(mapping = aes(y = y, alpha = alarm, group = ageGroup)) +
#   facet_wrap(facets = vars(method)) +
#   scale_y_continuous(name = "Incidence per 100.000") +
#   scale_x_date(name = "Month") +
#   scale_fill_manual(name = "Age group", values = dtuPalette) +
#   scale_alpha_manual(values = c(0.3, 1)) +
#   guides(alpha = "none") +
#   theme(panel.spacing.y = unit(1, "lines"), 
#         # axis.ticks.x = element_blank(),
#         axis.text.x = element_text(vjust = -1.2)) +
#   annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = as.Date("2011-03-01"), ymin = -Inf, ymax = Inf, alpha = 0.2)

Compare_stateOfTheArt_LIST <- Compare_stateOfTheArt_LIST_dat %>%
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
ggsave(filename = "Compare_stateOfTheArt_LIST.png",
       plot = Compare_stateOfTheArt_LIST,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")


# Hierarchical Poisson Normal model -------------------------------------------------

LIST_PoisN_ageGroup <- aeddo(data = LIST,
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

write_rds(x = LIST_PoisN_ageGroup, file = "LIST_PoisN_ageGroup.rds")
# LIST_PoisN_ageGroup <- read_rds(file = "LIST_PoisN_ageGroup.rds")

# start.theta.PoisN <- LIST_PoisN_ageGroup %>%
#   filter(row_number() == 1) %>%
#   select(par) %>%
#   unnest(par) %>%
#   select(theta)%>%
#   .$theta
# 
# LIST_PoisN_ageGroup_trend <- aeddo(data = LIST,
#                                    formula = y ~ -1 + t + ageGroup,
#                                    trend = TRUE,
#                                    theta = c(0, start.theta.PoisN),
#                                    method = "L-BFGS-B",
#                                    lower = c(-0.5, start.theta.PoisN-6),
#                                    upper = c(0.5, start.theta.PoisN+6),
#                                    model = "PoissonNormal", 
#                                    k = 36, 
#                                    sig.level = 0.9,
#                                    cpp.dir = "../models/",
#                                    CI = TRUE,
#                                    excludePastOutbreaks = TRUE)
# 
# write_rds(x = LIST_PoisN_ageGroup_trend, file = "LIST_PoisN_ageGroup_trend.rds")
# # LIST_PoisN_ageGroup_trend <- read_rds(file = "LIST_PoisN_ageGroup_trend.rds")
# 
# LIST_PoisN_ageGroup_seasonality <- aeddo(data = LIST,
#                                          formula = y ~ -1 +  ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
#                                          seasonality = TRUE,
#                                          theta = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
#                                          method = "L-BFGS-B",
#                                          lower = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) - 6,
#                                          upper = c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]) + 6,
#                                          model = "PoissonNormal", 
#                                          k = 36, 
#                                          sig.level = 0.9,
#                                          cpp.dir = "../models/",
#                                          CI = TRUE,
#                                          excludePastOutbreaks = TRUE)
# 
# write_rds(x = LIST_PoisN_ageGroup_seasonality, file = "LIST_PoisN_ageGroup_seasonality.rds")
# LIST_PoisN_ageGroup_seasonality <- read_rds(file = "LIST_PoisN_ageGroup_seasonality.rds")

# LIST_PoisN_ageGroup_trend_seasonality <- aeddo(data = LIST,
#                                                formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
#                                                trend = TRUE,
#                                                seasonality = TRUE,
#                                                theta = c(0, start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3]),
#                                                method = "L-BFGS-B",
#                                                lower = c(-0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])-6),
#                                                upper = c(0.5, c(start.theta.PoisN[1:2], 0,0, start.theta.PoisN[3])+6),
#                                                model = "PoissonNormal", 
#                                                k = 36, 
#                                                sig.level = 0.9,
#                                                cpp.dir = "../models/",
#                                                CI = TRUE,
#                                                excludePastOutbreaks = TRUE)
# 
# write_rds(x = LIST_PoisN_ageGroup_trend_seasonality, file = "LIST_PoisN_ageGroup_trend_seasonality.rds")
# LIST_PoisN_ageGroup_trend_seasonality <- read_rds(file = "LIST_PoisN_ageGroup_trend_seasonality.rds")


LIST_PoisN_ageGroup_par <- LIST_PoisN_ageGroup  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

LIST_PoisN_ageGroup_tbl <- LIST_PoisN_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup")

LIST_PoisN_ageGroup_unnested <- LIST_PoisN_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

LIST_PoisN_ageGroup_par_unnested <- LIST_PoisN_ageGroup %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

# Inspect the periods with relatively bad convergence from the Poisson Normal model
LIST_PoisN_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, log_sigma) %>%
  filter(log_sigma < -4) %>%
  print(n = 28)


LIST_PoisN_ageGroup_unnested %>%
  ggplot(mapping = aes(x = Date)) +
  geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
  geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Month") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(shape = "none")




# LIST_PoisN_ageGroup_trend_tbl <- LIST_PoisN_ageGroup_trend %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisN_ageGroup_trend")
# 
# LIST_PoisN_ageGroup_trend_unnested <- LIST_PoisN_ageGroup_trend %>% 
#   select(ran.ef) %>%
#   unnest(ran.ef) %>%
#   mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
#   select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)
# 
# LIST_PoisN_ageGroup_trend_unnested %>%
#   ggplot(mapping = aes(x = Date)) +
#   geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
#   geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_y_continuous(name = expression(u[t[1]])) +
#   scale_x_date(name = "Month") +
#   scale_colour_manual(values = dtuPalette) +
#   scale_shape_manual(values = c(1,19)) +
#   guides(shape = "none")




# LIST_PoisN_ageGroup_seasonality_tbl <- LIST_PoisN_ageGroup_seasonality %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisN_ageGroup_seasonality")
# 
# LIST_PoisN_ageGroup_seasonality_unnested <- LIST_PoisN_ageGroup_seasonality %>% 
#   select(ran.ef) %>%
#   unnest(ran.ef) %>%
#   mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
#   select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)
# 
# LIST_PoisN_ageGroup_seasonality_unnested %>%
#   ggplot(mapping = aes(x = Date)) +
#   geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
#   geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_y_continuous(name = expression(u[t[1]])) +
#   scale_x_date(name = "Month") +
#   scale_colour_manual(values = dtuPalette) +
#   scale_shape_manual(values = c(1,19)) +
#   guides(shape = "none")



# LIST_PoisN_ageGroup_trend_seasonality_tbl <- LIST_PoisN_ageGroup_trend_seasonality %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisN_ageGroup_trend_seasonality")

# LIST_PoisN_ageGroup_trend_seasonality_unnested <- LIST_PoisN_ageGroup_trend_seasonality %>% 
#   select(ran.ef) %>%
#   unnest(ran.ef) %>%
#   mutate(threshold = qnorm(p = 0.90, mean = 0, sd = exp(log_sigma))) %>%
#   select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)
# 
# LIST_PoisN_ageGroup_trend_seasonality_unnested %>%
#   ggplot(mapping = aes(x = Date)) +
#   geom_point(mapping = aes(y = `u_Poisson Normal`, colour = ageGroup, group = ageGroup, shape = `alarm_Poisson Normal`), size = 2) +
#   geom_line(mapping = aes(y = `threshold_Poisson Normal`), linewidth = 0.4) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_y_continuous(name = expression(u[t[1]])) +
#   scale_x_date(name = "Month") +
#   scale_colour_manual(values = dtuPalette) +
#   scale_shape_manual(values = c(1,19)) +
#   guides(shape = "none")


# Hierarchical Poisson Gamma model ----------------------------------------------------

# LIST_PoisG_ageGroup <- aeddo(data = LIST,
#                              formula = y ~ -1 + ageGroup,
#                              theta = rep(1,3),
#                              method = "L-BFGS-B",
#                              lower = c(rep(1e-6,2),-6),
#                              upper = rep(1e2, 3),
#                              model = "PoissonGamma", k = 36, cpp.dir = "../models/",
#                              excludePastOutbreaks = TRUE)

LIST_PoisG_ageGroup <- aeddo(data = LIST,
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

write_rds(x = LIST_PoisG_ageGroup, file = "LIST_PoisG_ageGroup.rds")
# LIST_PoisG_ageGroup <- read_rds(file = "LIST_PoisG_ageGroup.rds")

# start.theta.PoisG <- LIST_PoisG_ageGroup %>%
#   filter(row_number() == 1) %>%
#   select(par) %>%
#   unnest(par) %>%
#   select(theta)%>%
#   .$theta
# 
# LIST_PoisG_ageGroup_trend <- aeddo(data = LIST,
#                              formula = y ~ -1 + t + ageGroup,
#                              trend = TRUE,
#                              theta = c(0, start.theta.PoisG),
#                              method = "L-BFGS-B",
#                              lower = c(-0.5, start.theta.PoisG-6),
#                              upper = c(0.5, start.theta.PoisG+6),
#                              model = "PoissonGamma",
#                              k = 36,
#                              cpp.dir = "../models/",
#                              sig.level = 0.9,
#                              CI = TRUE,
#                              excludePastOutbreaks = TRUE)
# 
# write_rds(x = LIST_PoisG_ageGroup_trend, file = "LIST_PoisG_ageGroup_trend.rds")
# # LIST_PoisG_ageGroup_trend <- read_rds(file = "LIST_PoisG_ageGroup_trend.rds")
# 
# 
# LIST_PoisG_ageGroup_seasonality <- aeddo(data = LIST,
#                                    formula = y ~ -1 + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
#                                    seasonality = TRUE,
#                                    theta = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
#                                    method = "L-BFGS-B",
#                                    lower = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6,
#                                    upper = c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6,
#                                    model = "PoissonGamma",
#                                    k = 36,
#                                    cpp.dir = "../models/",
#                                    sig.level = 0.9,
#                                    CI = TRUE,
#                                    excludePastOutbreaks = TRUE)
# 
# write_rds(x = LIST_PoisG_ageGroup_seasonality, file = "LIST_PoisG_ageGroup_seasonality.rds")
# LIST_PoisG_ageGroup_seasonality <- read_rds(file = "LIST_PoisG_ageGroup_seasonality.rds")

# LIST_PoisG_ageGroup_trend_seasonality <- aeddo(data = LIST,
#                                          formula = y ~ -1 + t + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
#                                          trend = TRUE,
#                                          seasonality = TRUE,
#                                          theta = c(0, start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]),
#                                          method = "L-BFGS-B",
#                                          lower = c(-0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) - 6),
#                                          upper = c(0.5, c(start.theta.PoisG[1:2], 0,0, start.theta.PoisG[3]) + 6),
#                                          model = "PoissonGamma",
#                                          k = 36,
#                                          cpp.dir = "../models/",
#                                          sig.level = 0.9,
#                                          CI = TRUE,
#                                          excludePastOutbreaks = TRUE)

# write_rds(x = LIST_PoisG_ageGroup_trend_seasonality, file = "LIST_PoisG_ageGroup_trend_seasonality.rds")
# LIST_PoisG_ageGroup_trend_seasonality <- read_rds(file = "LIST_PoisG_ageGroup_trend_seasonality.rds")


LIST_PoisG_ageGroup %>%
  summarise(avgLogS = mean(LogS))

LIST_PoisG_ageGroup_par <- LIST_PoisG_ageGroup  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

LIST_PoisG_ageGroup_unnested <- LIST_PoisG_ageGroup %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)


LIST_PoisG_ageGroup_tbl <- LIST_PoisG_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup")

# LIST_PoisG_ageGroup_trend_tbl <- LIST_PoisG_ageGroup_trend %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisG_ageGroup")
# 
# LIST_PoisG_ageGroup_seasonality_tbl <- LIST_PoisG_ageGroup_seasonality %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisG_ageGroup")
# 
# LIST_PoisG_ageGroup_seasonality <- aeddo(data = LIST,
#                                          formula = y ~ -1 + ageGroup + sin(pi/6*periodInYear) + cos(pi/6*periodInYear),
#                                          seasonality = TRUE,
#                                          theta = rep(1, 5),
#                                          method = "L-BFGS-B",
#                                          lower = c(rep(1e-6,3),rep(-Inf,2), -6),
#                                          upper = rep(1e2, 5),
#                                          model = "PoissonGamma",
#                                          k = 36,
#                                          sig.level = 0.9,
#                                          cpp.dir = "../models/",
#                                          excludePastOutbreaks = TRUE)
# 
# LIST_PoisG_ageGroup_seasonality_tbl <- LIST_PoisG_ageGroup_seasonality %>%
#   select(ref.date, par, LogS) %>%
#   mutate(avgLogS = mean(LogS)) %>%
#   filter(row_number() == n()) %>%
#   select(-LogS) %>%
#   unnest(par) %>%
#   mutate(method = "PoisG_ageGroup_seasonality")
# 
# LIST_PoisG_ageGroup_unnested <- LIST_PoisG_ageGroup %>% 
#   select(ran.ef) %>%
#   unnest(ran.ef) %>%
#   mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi)) %>%
#   select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)
# 
# LIST_PoisG_ageGroup_par_unnested <- LIST_PoisG_ageGroup %>% 
#   select(ref.date, par) %>%
#   unnest(par) %>% 
#   mutate(Method = "Poisson Gamma")
# 
# LIST_PoisG_ageGroup %>%
#   select(ref.date , par) %>%
#   unnest(par) %>%
#   filter(Parameter == "log_phi") %>%
#   print(n = 144)


# Compare the Poisson Normal and Poisson Gamma model
LIST_novel_tbl <- bind_rows(
  LIST_PoisN_ageGroup_tbl,
  # LIST_PoisN_ageGroup_trend_tbl, 
  # LIST_PoisN_ageGroup_seasonality_tbl,
  LIST_PoisG_ageGroup_tbl,
  # LIST_PoisG_ageGroup_trend_tbl,
  # LIST_PoisG_ageGroup_seasonality_tbl
  )



LIST_novel_tbl %>% 
  print(n = 24)

write_rds(LIST_novel_tbl, file = "LIST_novel_tbl.rds")


LIST_novel <- full_join(LIST_PoisN_ageGroup_unnested, LIST_PoisG_ageGroup_unnested, by = join_by(Date, ageGroup))

Compare_novel_dat <- LIST_novel %>%
  pivot_longer(cols = c(`u_Poisson Normal`, `u_Poisson Gamma`), names_to = "method", names_prefix = "u_", values_to = "u") %>%
  pivot_longer(cols = c(`alarm_Poisson Normal`, `alarm_Poisson Gamma`), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  pivot_longer(cols = c(`threshold_Poisson Normal`, `threshold_Poisson Gamma`), names_to = "method3", names_prefix = "threshold_", values_to = "threshold") %>%
  filter(method == method2 & method == method3) %>%
  select(-method2, -method3) %>%
  mutate(method = factor(method, levels = c("Poisson Normal", "Poisson Gamma")))
write_rds(x = Compare_novel_dat, file = "Compare_novel_dat.rds")

Compare_novel_LIST <- Compare_novel_dat %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(y = threshold, group = method), lty = "dashed") +
  facet_grid(rows = vars(method), cols = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = expression(u[t[1]])) +
  scale_x_date(name = "Month") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")
ggsave(filename = "Compare_novel_LIST.png",
       plot = Compare_novel_LIST,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  


LIST_novel_par <- bind_rows(LIST_PoisN_ageGroup_par, LIST_PoisG_ageGroup_par)

custom_labeller <- as_labeller(
  c(`ageGroup<65 years`="beta['<65 years']",
    `ageGroup65+ years`="beta[65+~years]",
    `log_phi`="phi", `log_sigma`="sigma"),
  default = label_parsed
)

LIST_novel_par_plot <- LIST_novel_par %>%
  mutate_at(vars(theta:CI.upr), list(~case_when(Parameter == "log_sigma" ~ exp(.),
                                               Parameter == "log_phi" ~ exp(.),
                                               TRUE ~ .))) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), scales = "free_y", labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(widehat(theta))) +
  scale_x_date(name = "Month")
ggsave(filename = "LIST_novel_par_plot.png",
       plot = LIST_novel_par_plot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

LIST_novel_par_ageGroup <- LIST_novel_par %>%
  filter(grepl(x = Parameter, pattern = "ageGroup")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter)) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "LIST_novel_par_ageGroup.png",
       plot = LIST_novel_par_ageGroup,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

LIST_novel_par_dispersion <- LIST_novel_par %>%
  filter(grepl(x = Parameter, pattern = "log_sigma|log_phi")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter)) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "LIST_novel_par_dispersion.png",
       plot = LIST_novel_par_dispersion,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

# Compare the Poisson Normal and Poisson Gamma model
LIST_novel <- full_join(LIST_PoisN_ageGroup_unnested, LIST_PoisG_ageGroup_unnested, by = join_by(Date, ageGroup))

Compare_novel <- LIST_novel %>%
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
  scale_x_date(name = "Month") +
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


# Outbreaks invstigated by SSI ------------------------------------------------------

SSI_caseDistribution <- tribble(
  ~FUD, ~Date, ~y,
  "FUD1373", as.Date("2014-03-31"), 1,
  "FUD1373", as.Date("2014-05-26"), 1,
  "FUD1373", as.Date("2014-06-02"), 1,
  "FUD1373", as.Date("2014-06-16"), 2,
  "FUD1373", as.Date("2014-07-07"), 2,
  "FUD1373", as.Date("2014-07-14"), 2,
  "FUD1373", as.Date("2014-07-21"), 4,
  "FUD1373", as.Date("2014-07-28"), 4,
  "FUD1373", as.Date("2014-08-04"), 5,
  "FUD1373", as.Date("2014-08-11"), 4,
  "FUD1373", as.Date("2014-08-18"), 7,
  "FUD1373", as.Date("2014-08-25"), 1,
  "FUD1373", as.Date("2014-09-15"), 1,
  "FUD1373", as.Date("2014-09-22"), 1,
  "FUD1373", as.Date("2014-10-06"), 1,
  "FUD2080", as.Date("2022-05-09"), 2,
  "FUD2080", as.Date("2022-05-16"), 3,
  "FUD2080", as.Date("2022-05-23"), 3,
  "FUD2080", as.Date("2022-06-06"), 1,
  "FUD2074", as.Date("2022-10-01"), 1,
  "FUD2074", as.Date("2022-04-01"), 2,
  "FUD2074", as.Date("2022-05-01"), 1,
  "FUD2074", as.Date("2022-06-01"), 4,
  "FUD2127", as.Date("2022-08-15"), 1,
  "FUD2127", as.Date("2022-09-19"), 2,
  "FUD2127", as.Date("2022-09-26"), 1,
  "FUD2127", as.Date("2022-10-10"), 1,
  "FUD2127", as.Date("2022-10-24"), 3,
  "FUD2127", as.Date("2022-10-31"), 2,
  "FUD2127", as.Date("2022-11-21"), 1,
  "FUD887", as.Date("2009-04-06"), 7,
  "FUD887", as.Date("2009-05-11"), 1,
  "FUD1939", as.Date("2020-10-01"), 1,
  "FUD1939", as.Date("2020-11-01"), 1,
  "FUD1939", as.Date("2021-02-01"), 1,
  "FUD1939", as.Date("2021-04-01"), 1,
  "FUD1939", as.Date("2021-05-01"), 1,
  "FUD1939", as.Date("2021-06-01"), 1,
  "FUD1939", as.Date("2021-07-01"), 2,
  "FUD1939", as.Date("2021-08-01"), 1,
  "FUD1939", as.Date("2021-12-01"), 1,
  "FUD1939", as.Date("2022-06-01"), 1,
  "FUD1559", as.Date("2016-12-01"), 1,
  "FUD1559", as.Date("2017-01-01"), 2,
  "FUD1559", as.Date("2017-04-01"), 1,
  "FUD1559", as.Date("2017-05-01"), 1,
  "FUD1559", as.Date("2017-11-01"), 1,
  "FUD1559", as.Date("2018-11-01"), 1,
  "FUD1559", as.Date("2019-01-01"), 1,
  "FUD1559", as.Date("2019-02-01"), 1,
  "FUD1597", as.Date("2017-05-01"), 1,
  "FUD1597", as.Date("2017-08-01"), 4)

SSI_caseDistribution %>%
  ggplot(mapping = aes(x = Date, y = y)) +
  geom_col() +
  facet_wrap(facets = vars(FUD), ncol = 1)


SSI_outbreaks <- tibble(FUD = c("FUD887", "FUD1373", "FUD1597", "FUD1559",NA_character_, "FUD1970","FUD1939","FUD2074","FUD2080","FUD2127"),
                        Start = as.Date(x = c("2009-05-06","2013-09-01", "2017-05-01", "2016-12-01", "2016-01-01", "2018-12-01", "2020-10-01", "2021-10-01", "2022-05-13", "2022-08-15")),
                        Detected = as.Date(x = c(NA_character_, "2014-06-26", "2017-08-23", rep(NA_character_, 7))),
                        End = as.Date(x = c("2009-05-11","2014-10-01", "2017-09-01", "2019-02-01", "2019-09-20", "2021-11-01", "2022-05-01", "2022-06-01", "2022-06-06", "2022-11-27")))


tmp <- full_join(SSI_caseDistribution, SSI_outbreaks, by = join_by(FUD)) %>%
  print(n=25)



LIST_SSI_outbreaks <- SSI_outbreaks %>%
  arrange(desc(Start)) %>%
  mutate(outbreak_no = row_number()) %>%
  ggplot() +
  geom_segment(mapping = aes(x = Start, xend = End, y = outbreak_no, yend = outbreak_no), linewidth = 1.2, colour = dtuPalette[2]) +
  geom_point(mapping = aes(x = Start, y = outbreak_no), pch = 17, size = 3,colour = dtuPalette[2]) +
  geom_point(mapping = aes(x = Detected, y = outbreak_no), pch = 19, size = 3, colour = dtuPalette[2]) +
  scale_x_date(name = "Date", limits = c(as.Date(c("2008-01-01", "2022-12-01")))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(filename = "LIST_SSI_outbreaks.png",
       plot = LIST_SSI_outbreaks,
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

LIST_compare <- LIST_novel %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "LIST")

write_rds(x = LIST_compare, file = "LIST_compare.rds")
# LIST_compare <- read_rds(file = "LIST_compare.rds")

Compare_alarms <- LIST_compare %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  scale_y_discrete(limits = rev(levels(LIST_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_LIST.png",
       plot = Compare_alarms,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

