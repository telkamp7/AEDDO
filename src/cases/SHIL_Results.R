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
dat <- read_rds(file = "../../data/processed/dat2.rds")

# Only consider the SHIL cases
SHIL <- dat %>%
  filter(caseDef == "SHIL") %>%
  mutate(ageGroup = fct_collapse(ageGroup, `<25 years` = c("<1 year", "1-4 years","5-14 years","15-24 years")),
         ageGroup = fct_collapse(ageGroup, `25+ years` = c("25-64 years","65+ years"))) %>%
  group_by(Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))

SHIL_long_plot <- SHIL %>%
  ggplot(mapping = aes(x = Date, y = y, fill = ageGroup)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(name = "Age group", values = dtuPalette) +
  scale_y_continuous(name = "Number of cases") +
  scale_x_date(name = "Month")
ggsave(filename = "SHIL_long_plot.png",
       plot = SHIL_long_plot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Prepare to use surveillance package -----------------------------------------------

# Widen observations into a matrix format
observed <- SHIL %>%
  select(-n) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = y) %>%
  arrange(Date)

# Widen population sizes into a matrix format
population <- SHIL %>%
  select(-y) %>%
  pivot_wider(
    names_from = c(ageGroup),
    names_sep = ".",
    values_from = n) %>%
  arrange(Date)

# Convert observations into an 'sts' class
SHIL.sts <- sts(
  observed = observed[,-1],
  epoch = observed$Date,
  epochAsDate = TRUE,
  frequency = 12,
  population = as.matrix(population[,-1])
)

# Farrington ------------------------------------------------------------------------

SHIL_Farrington <- read_rds(file = "SHIL_Farrington.rds")

upperbound_Farrington <- as_tibble(SHIL_Farrington@upperbound) %>%
  mutate(Date = as.Date(x = SHIL_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIL$ageGroup)))

alarm_Farrington <- as_tibble(SHIL_Farrington@alarm) %>%
  mutate(Date = as.Date(x = SHIL_Farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Farrington") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIL$ageGroup)))

# Noufaily --------------------------------------------------------------------------

SHIL_Noufaily <- read_rds(file = "SHIL_Noufaily.rds")

upperbound_Noufaily <- as_tibble(SHIL_Noufaily@upperbound) %>%
  mutate(Date = as.Date(x = SHIL_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIL$ageGroup)))

alarm_Noufaily <- as_tibble(SHIL_Noufaily@alarm) %>%
  mutate(Date = as.Date(x = SHIL_Noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm_Noufaily") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(SHIL$ageGroup)))


# Compare the Farrington and Noufaily method

Compare_stateOfTheArt_SHIL_dat <- SHIL %>%
  left_join(upperbound_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  left_join(upperbound_Noufaily, by = join_by(Date, ageGroup)) %>%
  left_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = c(threshold_Farrington,threshold_Noufaily), names_to = "method", names_prefix = "threshold_", values_to = "threshold") %>%
  pivot_longer(cols = c(alarm_Farrington,alarm_Noufaily), names_to = "method2", names_prefix = "alarm_", values_to = "alarm") %>%
  filter(method == method2) %>%
  select(-method2) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))
write_rds(x = Compare_stateOfTheArt_SHIL_dat, file = "Compare_stateOfTheArt_SHIL_dat.rds")

Compare_stateOfTheArt_SHIL <- Compare_stateOfTheArt_SHIL_dat %>%
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
ggsave(filename = "Compare_stateOfTheArt_SHIL.png",
       plot = Compare_stateOfTheArt_SHIL,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Outbreaks invstigated by SSI ------------------------------------------------------

SSI_outbreaks <- tibble(Start = as.Date(x = c("2007-08-06", "2009-04-01", "2020-08-22")),
                        End = as.Date(x = c("2007-08-24", "2009-06-01", "2020-09-9")))

SHIL_SSI_outbreaks <- SSI_outbreaks %>%
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
ggsave(filename = "SHIL_SSI_outbreaks.png",
       plot = SHIL_SSI_outbreaks,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

# Hierarchical Poisson Normal model ---------------------------------------

SHIL_PoisN_ageGroup <- read_rds(file = "SHIL_PoisN_ageGroup.rds")
SHIL_PoisN_ageGroup_trend <- read_rds(file = "SHIL_PoisN_ageGroup_trend.rds")
SHIL_PoisN_ageGroup_seasonality <- read_rds(file = "SHIL_PoisN_ageGroup_seasonality.rds")
SHIL_PoisN_ageGroup_trend_seasonality <- read_rds(file = "SHIL_PoisN_ageGroup_trend_seasonality.rds")


# Hierarchical Poisson Gamma model ----------------------------------------

SHIL_PoisG_ageGroup <- read_rds(file = "SHIL_PoisG_ageGroup.rds")
SHIL_PoisG_ageGroup_trend <- read_rds(file = "SHIL_PoisG_ageGroup_trend.rds")
SHIL_PoisG_ageGroup_seasonality <- read_rds(file = "SHIL_PoisG_ageGroup_seasonality.rds")
SHIL_PoisG_ageGroup_trend_seasonality <- read_rds(file = "SHIL_PoisG_ageGroup_trend_seasonality.rds")



SHIL_PoisN_ageGroup_tbl <- SHIL_PoisN_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup")

SHIL_PoisN_ageGroup_trend_tbl <- SHIL_PoisN_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend")

SHIL_PoisN_ageGroup_seasonality_tbl <- SHIL_PoisN_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_seasonality")

SHIL_PoisN_ageGroup_trend_seasonality_tbl <- SHIL_PoisN_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisN_ageGroup_trend_seasonality")

SHIL_PoisG_ageGroup_tbl <- SHIL_PoisG_ageGroup %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup")

SHIL_PoisG_ageGroup_trend_tbl <- SHIL_PoisG_ageGroup_trend %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend")

SHIL_PoisG_ageGroup_seasonality_tbl <- SHIL_PoisG_ageGroup_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_seasonality")

SHIL_PoisG_ageGroup_trend_seasonality_tbl <- SHIL_PoisG_ageGroup_trend_seasonality %>%
  select(ref.date, par, LogS) %>%
  mutate(avgLogS = mean(LogS)) %>%
  filter(row_number() == n()) %>%
  select(-LogS) %>%
  unnest(par) %>%
  mutate(method = "PoisG_ageGroup_trend_seasonality")

SHIL_novel_tbl <- bind_rows(
  SHIL_PoisN_ageGroup_tbl,
  SHIL_PoisN_ageGroup_trend_tbl, 
  SHIL_PoisN_ageGroup_seasonality_tbl,
  SHIL_PoisN_ageGroup_trend_seasonality_tbl,
  SHIL_PoisG_ageGroup_tbl,
  SHIL_PoisG_ageGroup_trend_tbl,
  SHIL_PoisG_ageGroup_seasonality_tbl,
  SHIL_PoisG_ageGroup_trend_seasonality_tbl)

SHIL_novel_tbl %>% print(n = 68)

write_rds(SHIL_novel_tbl, file = "SHIL_novel_tbl.rds")
# Compare methods -------------------------------------------------------------------

SHIL_PoisN_ageGroup_trend_unnested <- SHIL_PoisN_ageGroup_trend  %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qnorm(p = 0.9, mean = 0, sd = exp(log_sigma))) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Normal` = u, `alarm_Poisson Normal` = alarm, `threshold_Poisson Normal` = threshold)

SHIL_PoisG_ageGroup_trend_unnested <- SHIL_PoisG_ageGroup_trend %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(threshold = qgamma(p = 0.9, shape = 1/phi, scale = phi)) %>%
  select(Date = ref.date, ageGroup, `u_Poisson Gamma` = u, `alarm_Poisson Gamma` = alarm, `threshold_Poisson Gamma` = threshold)


# Compare the Poisson Normal and Poisson Gamma model
SHIL_novel <- full_join(SHIL_PoisN_ageGroup_trend_unnested,
                        SHIL_PoisG_ageGroup_trend_unnested,
                        by = join_by(Date, ageGroup))
write_rds(x = SHIL_novel, file = "SHIL_novel.rds")

Compare_novel <- SHIL_novel %>%
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
ggsave(filename = "Compare_novel_SHIL.png",
       plot = Compare_novel,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  

SHIL_PoisN_ageGroup_par <- SHIL_PoisN_ageGroup  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

SHIL_PoisN_ageGroup_trend_par <- SHIL_PoisN_ageGroup_trend  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

SHIL_PoisN_ageGroup_seasonality_par <- SHIL_PoisN_ageGroup_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

SHIL_PoisN_ageGroup_trend_seasonality_par <- SHIL_PoisN_ageGroup_trend_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Normal")

SHIL_PoisG_ageGroup_par <- SHIL_PoisG_ageGroup  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

SHIL_PoisG_ageGroup_trend_par <- SHIL_PoisG_ageGroup_trend  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

SHIL_PoisG_ageGroup_seasonality_par <- SHIL_PoisG_ageGroup_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

SHIL_PoisG_ageGroup_trend_seasonality_par <- SHIL_PoisG_ageGroup_trend_seasonality  %>% 
  select(ref.date, par) %>%
  unnest(par) %>% 
  mutate(Method = "Poisson Gamma")

custom_labeller <- as_labeller(
  c(`ageGroup<25 years`="beta['<25~years']",`ageGroup25+ years`="beta[25+~years]",
    `t`="beta[trend]", `sin(pi/6 * periodInYear)` ="beta[sin]",
    `cos(pi/6 * periodInYear)`="beta[cos]", `log_phi`="phi", `log_sigma`="sigma",
    `Constant`="Constant", `Trend`="Trend", `Seasonality`="Seasonality", `Combined`="Combined"),
  default = label_parsed
)

SHIL_novel_par <- bind_rows(
  SHIL_PoisN_ageGroup_par %>% mutate(Model = "Constant"),
  SHIL_PoisG_ageGroup_par %>% mutate(Model = "Constant"),
  SHIL_PoisN_ageGroup_trend_par %>% mutate(Model = "Trend"),
  SHIL_PoisG_ageGroup_trend_par %>% mutate(Model = "Trend"),
  SHIL_PoisN_ageGroup_seasonality_par %>% mutate(Model = "Seasonality"),
  SHIL_PoisG_ageGroup_seasonality_par %>% mutate(Model = "Seasonality"),
  SHIL_PoisN_ageGroup_trend_seasonality_par %>% mutate(Model = "Combined"),
  SHIL_PoisG_ageGroup_trend_seasonality_par %>% mutate(Model = "Combined")) %>%
  mutate(Model = factor(Model, levels = c("Constant", "Trend", "Seasonality", "Combined")),
         Parameter = factor(Parameter, levels = c("ageGroup<25 years", "ageGroup25+ years","t",
                                                  "sin(pi/6 * periodInYear)","cos(pi/6 * periodInYear)",
                                                  "log_sigma", "log_phi")))


SHIL_novel_par_plot <- SHIL_novel_par %>%
  mutate_at(vars(theta:CI.upr), list(~case_when(Parameter == "log_sigma" ~ exp(.),
                                                Parameter == "log_phi" ~ exp(.),
                                                TRUE ~ .))) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_grid(rows = vars(Parameter), cols = vars(Model), scales = "free_y", labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(widehat(theta))) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIL_novel_par_plot.png",
       plot = SHIL_novel_par_plot,
       path = "../../figures/",
       device = png,
       width =22,
       height = 28,
       units = "in",
       dpi = "print")  



SHIL_novel_par <- bind_rows(SHIL_PoisG_ageGroup_trend_par, SHIL_PoisN_ageGroup_trend_par) 

SHIL_novel_par_ageGroup <- SHIL_novel_par %>%
  filter(grepl(x = Parameter, pattern = "ageGroup")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIL_novel_par_ageGroup.png",
       plot = SHIL_novel_par_ageGroup,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

SHIL_novel_par_trend <- SHIL_novel_par %>%
  filter(Parameter == "t") %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = theta, colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = CI.lwr, colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = CI.upr, colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(beta[i])) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIL_novel_par_trend.png",
       plot = SHIL_novel_par_trend,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print") 

SHIL_novel_par_dispersion <- SHIL_novel_par %>%
  filter(grepl(x = Parameter, pattern = "log_sigma|log_phi")) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = exp(theta), colour = Method), linewidth = 1) +
  geom_line(mapping = aes(y = exp(CI.lwr), colour = Method), lty = "dashed") +
  geom_line(mapping = aes(y = exp(CI.upr), colour = Method), lty = "dashed") +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_y_continuous(name = expression(Psi)) +
  scale_x_date(name = "Month")
ggsave(filename = "SHIL_novel_par_dispersion.png",
       plot = SHIL_novel_par_dispersion,
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

SHIL_compare <- SHIL_novel %>%
  select(Date, ageGroup, `alarm_Poisson Normal`, `alarm_Poisson Gamma`) %>%
  full_join(alarm_Farrington, by = join_by(Date, ageGroup)) %>%
  full_join(alarm_Noufaily, by = join_by(Date, ageGroup)) %>%
  pivot_longer(cols = `alarm_Poisson Normal`:alarm_Noufaily, names_to = "method", names_prefix = "alarm_", values_to = "alarm") %>%
  group_by(Date, method) %>%
  reframe(alarm = any(alarm)) %>%
  bind_rows(SSI_corrected) %>%
  mutate(alarmDate = if_else(alarm, Date, NA), method = factor(method), disease = "SHIL")

write_rds(x = SHIL_compare, file = "SHIL_compare.rds")

Compare_alarms <- SHIL_compare %>%
  filter(Date >= as.Date("2011-01-01")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17, size = 4) +
  scale_color_manual(values = dtuPalette[c(7,9:11,4)]) +
  scale_y_discrete(limits = rev(levels(SHIL_compare$method))) +
  scale_x_date(name = "Date") +
  guides(colour = "none") +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))
ggsave(filename = "Compare_alarms_SHIL.png",
       plot = Compare_alarms,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")  
