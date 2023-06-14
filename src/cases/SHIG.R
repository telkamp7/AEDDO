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
  noPeriods = 1, pastWeeksNotIncluded = NULL,
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
  alpha = 0.05, trend = TRUE, pThresholdTrend = 1,
  limit54 = c(0,4), powertrans = "2/3",
  fitFun = "algo.farrington.fitGLM.flexible",
  populationOffset = TRUE,
  noPeriods = 1, pastWeeksNotIncluded = NULL,
  thersholdMethod = "nbPlugin"
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
  geom_segment(mapping = aes(x = Start, xend = End, y = outbreak_no, yend = outbreak_no), linewidth = 1.2, colour = dtuPalette[4]) +
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
# SHIG_PoisN_ageGroup_trend <- read_rds(file = "SHIG_PoisN_ageGroup_trend.rds")

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
# SHIG_PoisN_ageGroup_seasonality <- read_rds(file = "SHIG_PoisN_ageGroup_seasonality.rds")

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
# SHIG_PoisN_ageGroup_trend_seasonality <- read_rds(file = "SHIG_PoisN_ageGroup_trend_seasonality.rds")


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
# SHIG_PoisG_ageGroup <- read_rds(file = "SHIG_PoisG_ageGroup.rds")

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
# SHIG_PoisG_ageGroup_trend <- read_rds(file = "SHIG_PoisG_ageGroup_trend.rds")


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
# SHIG_PoisG_ageGroup_seasonality <- read_rds(file = "SHIG_PoisG_ageGroup_seasonality.rds")

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
# SHIG_PoisG_ageGroup_trend_seasonality <- read_rds(file = "SHIG_PoisG_ageGroup_trend_seasonality.rds")





SHIG_PoisN_ageGroup_tbl <- SHIG_PoisN_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup")

SHIG_PoisN_ageGroup_trend_tbl <- SHIG_PoisN_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend")

SHIG_PoisN_ageGroup_seasonality_tbl <- SHIG_PoisN_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_seasonality")

SHIG_PoisN_ageGroup_trend_seasonality_tbl <- SHIG_PoisN_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend_seasonality")

SHIG_PoisG_ageGroup_tbl <- SHIG_PoisG_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup")

SHIG_PoisG_ageGroup_trend_tbl <- SHIG_PoisG_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend")

SHIG_PoisG_ageGroup_seasonality_tbl <- SHIG_PoisG_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_seasonality")

SHIG_PoisG_ageGroup_trend_seasonality_tbl <- SHIG_PoisG_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend_seasonality")

SHIG_novel_tbl <- bind_rows(
  SHIG_PoisN_ageGroup_tbl,
  SHIG_PoisN_ageGroup_trend_tbl, 
  SHIG_PoisN_ageGroup_seasonality_tbl,
  SHIG_PoisN_ageGroup_trend_seasonality_tbl,
  SHIG_PoisG_ageGroup_tbl,
  SHIG_PoisG_ageGroup_trend_tbl,
  SHIG_PoisG_ageGroup_seasonality_tbl,
  SHIG_PoisG_ageGroup_trend_seasonality_tbl)

SHIG_novel_tbl %>% print(n = 68)

write_rds(SHIG_novel_tbl, file = "SHIG_novel_tbl.rds")
# Compare methods -------------------------------------------------------------------

SHIG_PoisN_ageGroup_trend_unnested <- SHIG_PoisN_ageGroup_trend  %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

SHIG_PoisG_ageGroup_trend_unnested <- SHIG_PoisG_ageGroup_trend %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)


# Compare the Poisson Normal and Poisson Gamma model
SHIG_novel <- full_join(SHIG_PoisN_ageGroup_trend_unnested,
                        SHIG_PoisG_ageGroup_trend_unnested,
                        by = join_by(Date, ageGroup))

Compare_novel <- SHIG_novel %>%
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
ggsave(filename = "Compare_novel_SHIG.png",
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
    `ageGroup<25 years`="beta['<25~years']",`ageGroup25+ years`="beta[25+~years]",
    `t`="beta[trend]", `sin(pi/6 * monthInYear)` ="beta[sin]",
    `cos(pi/6 * monthInYear)`="beta[cos]", `log_phi`="phi", `log_sigma`="sigma"),
  default = label_parsed
)

SHIG_PoisG_ageGroup_trend_par <- SHIG_PoisG_ageGroup_trend  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

SHIG_PoisN_ageGroup_trend_par <- SHIG_PoisN_ageGroup_trend  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

SHIG_novel_par <- bind_rows(SHIG_PoisG_ageGroup_trend_par, SHIG_PoisN_ageGroup_trend_par) 


SHIG_novel_par_ageGroup <- SHIG_novel_par %>%
  filter(grepl(x = Parameter, pattern = "ageGroup")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIG_novel_par_ageGroup.png",
       plot = SHIG_novel_par_ageGroup,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

SHIG_novel_par_trend <- SHIG_novel_par %>%
  filter(Parameter == "t") %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIG_novel_par_trend.png",
       plot = SHIG_novel_par_trend,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

# SHIG_novel_par_seasonality <- SHIG_novel_par %>%
#   filter(grepl(x = Parameter, pattern = "cos|sin")) %>%
#   ggplot(mapping = aes(x = ref.date)) +
#   geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
#   geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
#   geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
#   facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
#   scale_color_manual(values = dtuPalette) +
#   scale_y_continuous(name = expression(beta[i])) +
#   scale_x_date(name = "Month")
# ggsave(filename = "SHIG_novel_par_seasonality.png",
#        plot = SHIG_novel_par_seasonality,
#        path = "../../figures/",
#        device = png,
#        width = 16,
#        height = 8,
#        units = "in",
#        dpi = "print") 

SHIG_novel_par_dispersion <- SHIG_novel_par %>%
  filter(grepl(x = Parameter, pattern = "log_sigma|log_phi")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = exp(theta), colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = exp(CI.lwr), colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = exp(CI.upr), colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(Psi)) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIG_novel_par_dispersion.png",
       plot = SHIG_novel_par_dispersion,
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

SHIG_compare <- SHIG_novel %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "SHIG")

write_rds(x = SHIG_compare, file = "SHIG_compare.rds")

Compare_alarms <- SHIG_compare %>%
  filter(Date >= as.Date("2011-01-01")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,4)]) +
  scale_y_discrete(limits = rev(levels(SHIG_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_SHIG.png",
       plot = Compare_alarms,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  
