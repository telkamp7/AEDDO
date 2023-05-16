

# Import libraries
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
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
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18))
)

# Source the aeddo function
source(file = "aeddo.R")

# Import the data
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6 agegroups

# Nest the data
dat_nest <- dat %>%
  group_by(caseDef, Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(monthInYear = as.integer(format(Date, "%m"))) %>%
  group_by(caseDef) %>%
  nest()

PoissonGamma <- dat_nest %>%
  filter(caseDef %in% c("Shiga- og veratoxin producerende E. coli.")) %>%
  mutate(
    M1 = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup,
        theta = rep(1,7),
        method = "L-BFGS-B", 
        lower = c(rep(1e-6,6),-6), 
        model = "PoissonGamma",
        k = 36, 
        excludePastOutbreaks = TRUE)
    }
    ),
    M2 = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
        theta = rep(1,9),
        method = "L-BFGS-B", 
        lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
        model = "PoissonGamma",
        k = 36, 
        excludePastOutbreaks = TRUE)
    }
    )
  )

PoissonGamma %>%
  ungroup() %>%
  select(M1, M2) %>%
  unnest(c(M1, M2), names_sep = "_") %>%
  select(ref.date = M1_ref.date, M1 = M1_value, M2 = M2_value) %>%
  pivot_longer(cols = M1:M2, names_to = "Model") %>%
  ggplot(mapping = aes(x = ref.date, y = value, colour = Model)) +
  geom_line() +
  scale_x_date(name = "Month") +
  scale_y_continuous(name = "Objective") +
  scale_colour_manual(values = dtuPalette)
  

PoissonGamma %>%
  ungroup() %>%
  select(M1, M2) %>%
  unnest(c(M1, M2), names_sep = "_") %>%
  select(ref.date = M1_ref.date, M1 = M1_ran.ef, M2 = M2_ran.ef) %>%
  unnest(c(M1, M2), names_sep = "_") %>%
  select(ref.date, ageGroup = M1_ageGroup, n = M1_n, y = M1_y, monthInYear = M1_monthInYear, M1_phi:M1_alarm, M2_phi:M2_alarm) %>%
  rename(M1 = M1_u, M2 = M2_u) %>%
  pivot_longer(cols = c(M1, M2), names_to = "Model_u", values_to = "u") %>%
  rename(M1 = M1_phi, M2 = M2_phi) %>%
  pivot_longer(cols = c(M1, M2), names_to = "Model_phi", values_to = "phi") %>%
  rename(M1 = M1_alarm, M2 = M2_alarm) %>%
  pivot_longer(cols = c(M1, M2), names_to = "Model_alarm", values_to = "Alarm") %>%
  filter(Model_u == Model_phi & Model_u == Model_alarm) %>%
  rename(Model = Model_u) %>%
  select(ref.date:monthInYear, Model, u, phi, Alarm) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = Alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_grid(rows = vars(Model), cols = vars(ageGroup)) +
  scale_x_date(name = "Month") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none",shape = "none")



PoissonGamma %>%
  ungroup() %>%
  select(M1, M2) %>%
  unnest(c(M1, M2), names_sep = "_") %>%
  select(ref.date = M1_ref.date, M1 = M1_ran.ef, M2 = M2_ran.ef) %>%
  unnest(c(M1, M2), names_sep = "_") %>%
  select(ref.date, ageGroup = M1_ageGroup, M1_alarm, M2_alarm) %>%
  rename(M1 = M1_alarm, M2 = M2_alarm) %>%
  pivot_longer(cols = c(M1, M2), names_to = "Model", values_to = "Alarm") %>%
  mutate(Model = factor(Model, levels = c("M1", "M2")),
         alarmDate = if_else(Alarm, ref.date, as.Date(NA_character_))) %>%
  ggplot(mapping = aes(x = alarmDate, y = Model, colour = Model)) +
  geom_point(shape = 17) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_y_discrete(limits = rev(c("M1", "M2"))) +
  guides(colour = "none")
