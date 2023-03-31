
# Import libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(TMB)

# Set locale
Sys.setlocale(category = "LC_ALL", locale = "en")

# Dynamically link the C++ template
# dyn.load(dynlib(name = "../models/PoissonLognormal"))
dyn.load(dynlib(name = "../models/PoissonNormal"))
# Set global theme options
theme_set(
  new = theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 18))
)

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

# Load data
# dat <- read_rds(file = "../../data/processed/dat.rds") # 11-agegroups
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6 agegroups

# # Only consider some of the data
# y <- dat %>%
#   filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
#   group_by(Date, ageGroup) %>%
#   reframe(y = sum(cases), n = sum(n), month = factor(format(x = Date, "%b")), year = factor(format(x = Date, "%Y")))


monthLevels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  mutate(month = factor(format(x = Date, "%b"), levels = monthLevels),
         year = factor(format(x = Date, "%Y"))) %>%
  group_by(month, year, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n))
y$month

# Epidemiological plots
y %>%
  filter(year %in% as.character(2012:2022)) %>%
  group_by(year, month) %>%
  reframe(y = sum(y), n = sum(n)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name = "Number of cases per 100.000") +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.")
ggsave(filename = "EpixSTEC.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

# Some diseases
diseases <- c("Gonoré",
              "Shiga- og veratoxin producerende E. coli.",
              "Legionella",
              "Leptospirosis")

dat %>%
  filter(caseDef == "Gonoré") %>%
  summarise(age = n_distinct(ageGroup), ld = n_distinct(landsdel))

tmp <- dat %>%
  filter(caseDef == "Gonoré") %>%
  group_by(Date) %>%
  reframe(y = sum(cases), n = sum(n))

dat %>%
  group_by(caseDef) %>%
  summarize(nAgeGroup = n_distinct(ageGroup), nLandsdel = n_distinct(landsdel)) %>%
  print(n = 26)

# Epidemiological plots
dat %>%
  filter(Date >= as.Date("2012-01-01") & caseDef %in% diseases) %>%
  mutate(month = factor(format(x = Date, "%b")), year = factor(format(x = Date, "%Y"))) %>%
  group_by(year, month, caseDef) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(caseDef = factor(caseDef, labels = diseases))


# Epidemiological plots
dat %>%
  filter(Date >= as.Date("2012-01-01") & caseDef %in% diseases) %>%
  mutate(month = factor(format(x = Date, "%b")), year = factor(format(x = Date, "%Y"))) %>%
  group_by(year, month, caseDef) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(caseDef = factor(caseDef, labels = diseases)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(caseDef), scales = "free_y") +
  scale_y_continuous(name = "Number of cases per 100.000") +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1))

dat %>%
  group_by(caseDef) %>%
  summarize(miny = min(cases), maxy = max(cases), meany = mean(cases), sumy = sum(cases),
            minn = min(n), maxn = max(n), meann = mean(n)) %>% 
  arrange(sumy) %>%
  print(n = 26)

dat %>%
  filter(n == 236)


# Extract case definitions
caseDef <- unique(dat$caseDef)

# Loop over case definitions and create overview plots
for(c in caseDef){
  
  # caseDefxLandsdel
  dat %>% 
    filter(caseDef == c) %>% 
    group_by(Date, landsdel, ageGroup) %>%
    reframe(cases = sum(cases), n = sum(n)) %>%
    ggplot(mapping = aes(x = Date, y = cases/n * 1e5, colour = ageGroup)) +
    geom_point() +
    facet_wrap(facets = vars(landsdel)) +
    scale_y_continuous(name = "Number of cases per 100.000") +
    scale_colour_manual(name = "Age group", values = dtuPalette) +
    guides(colour = guide_legend(nrow = 1)) +
    ggtitle(label = c)
  ggsave(filename = paste0(
    str_replace_all(string = c, pattern = " |-|/|\\.", ""),"xLandsdel.png"
    ),
    path = "../../figures/",
    device = png,
    width = 16,
    height = 8,
    units = "in",
    dpi = "print")
  
  # caseDefxAgeGroup
  dat %>% 
    filter(caseDef == c) %>% 
    group_by(Date, ageGroup) %>%
    reframe(cases = sum(cases), n = sum(n)) %>%
    ggplot(mapping = aes(x = Date, y = cases/n * 1e5, colour = ageGroup)) +
    geom_point() +
    facet_wrap(facets = vars(ageGroup), scales = "free_y") +
    scale_y_continuous(name = "Number of cases per 100.000") +
    scale_colour_manual(values = dtuPalette) +
    guides(colour = "none") +
    ggtitle(label = c)
  ggsave(filename = paste0(
    str_replace_all(string = c, pattern = " |-|/|\\.", ""),"xAgeGroup.png"
  ),
  path = "../../figures/",
  device = png,
  width = 16,
  height = 8,
  units = "in",
  dpi = "print")
}

# Extract case definitions
landsdel <- unique(dat$landsdel)

# Loop over 'landsdele' and create overview plots
for(l in landsdel){
  dat %>%
    filter(landsdel == l) %>%
    ggplot(mapping = aes(x = Date, y = cases/n * 1e5, colour = ageGroup)) +
    geom_point() +
    facet_wrap(facets = vars(caseDef)) +
    scale_y_continuous(name = "Number of cases per 100.000") +
    scale_colour_manual(name = "Age group", values = dtuPalette) +
    guides(colour = guide_legend(nrow = 1)) +
    ggtitle(label = l)
  ggsave(filename = paste0(l,"xCaseDef.png"), path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
}

# Load the Poisson-normal model
PoisLN_res <- read_rds(file = "../methods/PoisLN_res.rds")

PoisLN_res %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA)) %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_point(mapping = aes(y = u)) +
  geom_line(mapping = aes(x = Date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierarchical Poisson Normal model")
ggsave(filename = "PoisLNxSTEC.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

### Deprecated
# PoisLN <- read_rds(file = "../models/PoissonLognormal.rds")
# # ... and generate report
# rep <- sdreport(PoisLN, getJointPrecision = TRUE)
# 
# PoisLN.res <- y %>%
#   mutate(`Random effects` = rep$par.random,
#          `Std. Error` = sqrt(rep$diag.cov.random))
# 
# PoisLN.res %>%
#   ggplot(mapping = aes(x = Date, y = `Random effects`, colour = ageGroup)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
#   scale_colour_manual(values = dtuPalette) +
#   guides(colour = "none") +
#   ggtitle(label = "Shiga- og veratoxin producerende E. coli.")

# Load the Poisson-Gamma model
PoisG_res <- read_rds(file = "../methods/PoisG_res.rds")

PoisG_res %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA)) %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_point(mapping = aes(y = u)) +
  geom_line(mapping = aes(x = Date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "PoisGxSTEC.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

### Deprecated
# PoisG <- read_rds(file = "../models/PoissonGamma.rds")
# 
# PoisG$results %>%
#   rename(`Random effects` = u) %>%
#   ggplot(mapping = aes(x = Date, y = `Random effects`, colour = ageGroup)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_y_continuous(name = expression(paste("Random effect, ", u[t]^a))) +
#   scale_colour_manual(values = dtuPalette) +
#   guides(colour = "none") +
#   ggtitle(label = "Shiga- og veratoxin producerende E. coli.")

# Load the Farrington-method
STEC_farrington <- read_rds(file = "../models/STEC_farrington.rds")

observed <- as_tibble(STEC_farrington@observed) %>%
  mutate(Date = as.Date(x = STEC_farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "y") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

upperbound <- as_tibble(STEC_farrington@upperbound) %>%
  mutate(Date = as.Date(x = STEC_farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

alarm <- as_tibble(STEC_farrington@alarm) %>%
  mutate(Date = as.Date(x = STEC_farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

population <- as_tibble(STEC_farrington@populationFrac) %>%
  mutate(Date = as.Date(x = STEC_farrington@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "n") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

STEC_farrington_tbl <- observed %>% 
  full_join(y = upperbound, join_by(Date, ageGroup)) %>%
  full_join(y = alarm, join_by(Date, ageGroup)) %>%
  full_join(y = population, join_by(Date, ageGroup)) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))

STEC_farrington_tbl %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_point(mapping = aes(y = y/n*1e5)) +
  geom_line(mapping = aes(y = threshold/n*1e5), lty = "dotted", size = 0.5) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL), outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Number of cases per 100.000") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Farrington method")
ggsave(filename = "STEC_farrington.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

# Load the Farrington-method
STEC_noufaily <- read_rds(file = "../models/STEC_noufaily.rds")

observed <- as_tibble(STEC_noufaily@observed) %>%
  mutate(Date = as.Date(x = STEC_noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "y") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

upperbound <- as_tibble(STEC_noufaily@upperbound) %>%
  mutate(Date = as.Date(x = STEC_noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "threshold") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

alarm <- as_tibble(STEC_noufaily@alarm) %>%
  mutate(Date = as.Date(x = STEC_noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "alarm") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

population <- as_tibble(STEC_noufaily@populationFrac) %>%
  mutate(Date = as.Date(x = STEC_noufaily@epoch, origin = "1970-01-01")) %>%
  pivot_longer(cols = -Date, 
               names_to = "ageGroup",
               values_to = "n") %>%
  mutate(ageGroup = factor(x = ageGroup, levels = levels(dat$ageGroup)))

STEC_noufaily_tbl <- observed %>% 
  full_join(y = upperbound, join_by(Date, ageGroup)) %>%
  full_join(y = alarm, join_by(Date, ageGroup)) %>%
  full_join(y = population, join_by(Date, ageGroup)) %>%
  mutate(dateOfAlarm = if_else(alarm, Date, NA))

STEC_noufaily_tbl %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_point(mapping = aes(y = y/n*1e5)) +
  geom_line(mapping = aes(y = threshold/n*1e5), lty = "dotted", size = 0.5) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL), outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Number of cases per 100.000") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Noufaily method")
ggsave(filename = "STEC_noufaily.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")



# STEC Result -----------------------------------------------------------------------

STEC_res <- read_rds(file = "../models/STEC_res.rds")

STEC_res %>%
  select(PoissonGamma) %>%
  unnest(PoissonGamma) %>% 
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line() +
  facet_wrap(facets = vars(Parameter), scales = "free_y")

tmp <- STEC_res %>%
  select(PoissonGamma) %>%
  unnest(PoissonGamma) %>% 
  ungroup() %>%
  unnest(ran.ef) %>%
  select(ref.date, window.date, u,  phi, ageGroup, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, window.date, NA))

STEC_res %>%
  select(PoissonGamma) %>%
  unnest(PoissonGamma) %>% 
  ungroup() %>%
  unnest(ran.ef) %>%
  select(ref.date, window.date, u,  phi, ageGroup, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, window.date, NA)) %>%
  filter(ref.date == window.date) %>%
  ggplot(mapping = aes(x = window.date, y = u, colour = ageGroup)) +
  geom_line() +
  geom_point() +
  geom_line(mapping = aes(x = window.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "windowedSTECPoisG.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonGamma_excludePastOutbreaks) %>%
  unnest(PoissonGamma_excludePastOutbreaks) %>% 
  ungroup() %>%
  unnest(ran.ef) %>%
  select(ref.date, window.date, u,  phi, ageGroup, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, window.date, NA)) %>%
  filter(ref.date == window.date) %>%
  ggplot(mapping = aes(x = window.date, y = u, colour = ageGroup)) +
  geom_line() +
  geom_point() +
  geom_line(mapping = aes(x = window.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "windowedSTECPoisGExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonNormal) %>%
  unnest(PoissonNormal) %>% 
  ungroup() %>%
  unnest(ran.ef) %>%
  select(ref.date, window.date, u,  log_sigma, ageGroup, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, window.date, NA)) %>%
  filter(ref.date == window.date) %>%
  ggplot(mapping = aes(x = window.date, y = u, colour = ageGroup)) +
  geom_line() +
  geom_point() +
  geom_line(mapping = aes(x = window.date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierarchical Poisson Normal model")
ggsave(filename = "windowedSTEDPoisN.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>% 
  ungroup() %>%
  unnest(ran.ef) %>%
  select(ref.date, window.date, u,  log_sigma, ageGroup, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, window.date, NA)) %>%
  filter(ref.date == window.date) %>%
  ggplot(mapping = aes(x = window.date, y = u, colour = ageGroup)) +
  geom_line() +
  geom_point() +
  geom_line(mapping = aes(x = window.date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierarchical Poisson Normal model")
ggsave(filename = "windowedSTEDPoisNExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>% 
  ungroup() %>%
  unnest(par) %>%
  select(ref.date, Parameter, ageGroup, theta, se.theta) %>%
  filter(ageGroup != "All") %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = ageGroup)) +
  geom_line() +
  facet_wrap(facets = vars(Parameter)) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")

STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>% 
  ungroup() %>%
  unnest(par) %>%
  select(ref.date, Parameter, ageGroup, theta, se.theta) %>%
  filter(ageGroup == "All") %>%
  ggplot(mapping = aes(x = ref.date, y = exp(theta), colour = ageGroup)) +
  geom_ribbon(mapping = aes(x = ref.date, ymin = exp(theta - 2 * se.theta), ymax = exp(theta + 2 * se.theta)), alpha = 0.4, fill = "grey70", inherit.aes = FALSE) +
  geom_line() +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
  


## STEC Agegroup epidimiological
dat %>% 
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>% 
  mutate(Year = format(Date, "%Y"), Month = factor(format(x = Date, "%b"), levels = monthLevels)) %>%
  group_by(Year, Month, ageGroup) %>%
  reframe(cases = sum(cases), n = sum(n)) %>%
  filter(Year %in% as.character(2012:2022)) %>%
  ggplot(mapping = aes(x = Month, y = cases/n * 1e5, colour = Year, group = Year)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Number of cases per 100.000") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.")
ggsave(filename = "STECxEpixAgeGroup.png",
path = "../../figures/",
device = png,
width = 16,
height = 8,
units = "in",
dpi = "print")
  

