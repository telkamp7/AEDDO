
# Import libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(TMB)
library(scales)

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
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 35))
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
                

tibble(z = seq(0,5,length.out = 101)) %>%
  mutate(density = dlnorm(z, meanlog = 0, sdlog = 0.5)) %>%
  mutate(discreteDensity = case_when(z < 1 ~ cumsum(z),
                                     z > 1 & z < 2 ~ cumsum(z))) %>% 
  print(n = 101)

plot(diff(cumsum(dlnorm(x = seq(0, 5, length.out = 101), meanlog = 0, sdlog = 0.5))))

PDFLogNormal <- tibble(z = c(0:4),
       density = c(
  integrate(f = dlnorm, lower = 0, upper = 1, meanlog = 0, sdlog = 0.5)$value,
  integrate(f = dlnorm, lower = 1, upper = 2, meanlog = 0, sdlog = 0.5)$value,
  integrate(f = dlnorm, lower = 2, upper = 3, meanlog = 0, sdlog = 0.5)$value,
  integrate(f = dlnorm, lower = 3, upper = 4, meanlog = 0, sdlog = 0.5)$value,
  integrate(f = dlnorm, lower = 4, upper = 5, meanlog = 0, sdlog = 0.5)$value)) %>%
  ggplot(mapping = aes(x = z, y = density)) +
  geom_step(colour = "#1FD082", linewidth = 2) +
  scale_y_continuous(name = "PMF")
ggsave(filename = "PDFLogNormal.png", plot = PDFLogNormal, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
# 
# PDFLogNormal <- ggplot(data = data.frame(z = c(0,3)), aes(z)) +
#   stat_function(fun = dlnorm, n = 101, args = list(mean = 0, sd = 0.5), colour = "#1FD082", linewidth = 1.2) +
#   scale_y_continuous(name = "PDF", breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0,1)) +
#   scale_x_continuous(name = )
# ggsave(filename = "PDFLogNormal.png", plot = PDFLogNormal, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")


# Load data
# dat <- read_rds(file = "../../data/processed/dat.rds") # 11-agegroups
dat <- read_rds(file = "../../data/processed/dat2.rds") # 6 agegroups

monthLevels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
finalDat <- dat %>%
  group_by(Date, caseDef) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(month = factor(format(x = Date, "%b"), levels = monthLevels),
         year = factor(format(x = Date, "%Y")))

# LIST_long_plot <- finalDat %>%
#   filter(caseDef == "LIST") %>%
#   ggplot(mapping = aes(x = Date, y = y/n*1e5, fill = "All")) +
#   geom_col() +
#   scale_y_continuous(name = "Incidence per 100.000") +
#   scale_x_date(name = "Date") +
#   scale_fill_manual(values = dtuPalette) +
#   guides(fill = "none")
# LIST_long_plot
# ggsave(filename = "LIST_long_plot.png", plot = LIST_long_plot, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
# 
# 
# STEC_long_plot <- dat %>%
#   filter(caseDef == "STEC") %>%
#   ggplot(mapping = aes(x = Date, y = cases/n*1e5, fill = ageGroup, group = ageGroup)) +
#   geom_col() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_y_continuous(name = "Incidence per 100.000") +
#   scale_x_date(name = "Date") +
#   scale_fill_manual(values = dtuPalette) +
#   guides(fill = "none")
# STEC_long_plot
# ggsave(filename = "STEC_long_plot.png", plot = STEC_long_plot, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
# 
# 
# SHIG_long_plot <- dat %>%
#   filter(caseDef == "SHIG") %>%
#   ggplot(mapping = aes(x = Date, y = cases/n*1e5, fill = ageGroup, group = ageGroup)) +
#   geom_col() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_y_continuous(name = "Incidence per 100.000") +
#   scale_x_date(name = "Date") +
#   scale_fill_manual(values = dtuPalette) +
#   guides(fill = "none")
# SHIG_long_plot
# ggsave(filename = "SHIG_long_plot.png", plot = SHIG_long_plot, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")


# Epidemiological plots
stecEpiPlot <- finalDat %>%
  filter(caseDef == "STEC" & year %in% as.character(2012:2022)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.1)) +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(override.aes = list(size = 3, linewidth = 2), nrow = 1)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.title = element_text(size = 19),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 25))
yearLegend <- get_legend(stecEpiPlot)

shigellaEpiPlot <- finalDat %>%
  filter(caseDef == "SHIL" & year %in% as.character(2012:2022)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.1)) +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.title = element_text(size = 19))


listEpiPlot <- finalDat %>%
  filter(caseDef == "LIST" & year %in% as.character(2012:2022)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.1)) +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.title = element_text(size = 19))

salmEpiPlot <- finalDat %>%
  filter(caseDef == "SALM" & year %in% as.character(2012:2022)) %>%
  ggplot(mapping = aes(x = month, y = y/n * 1e5, colour = year, group = year)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.1)) +
  scale_x_discrete(name = "Month") +
  scale_colour_manual(name = "Year", values = dtuPalette) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
          axis.title = element_text(size = 19))

diseasePlots <- plot_grid(listEpiPlot + guides(colour = "none"),
                          shigellaEpiPlot + guides(colour = "none"),
                          stecEpiPlot + guides(colour = "none"),
                          salmEpiPlot + guides(colour = "none"),
                          labels = c("(a)", "(b)", "(c)", "(d)"),
                          label_size = 18,
                          align = "hv")

finalEpiPlot <- plot_grid(yearLegend, diseasePlots, ncol = 1, rel_heights = c(0.05,0.9))
finalEpiPlot
ggsave(filename = "EpiPlot.png", plot = finalEpiPlot, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

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
  scale_y_continuous(name = "Incidence per 100.000") +
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
    scale_y_continuous(name = "Incidence per 100.000") +
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
    geom_line(size = 0.4) +
    geom_point(size = 2) +
    facet_wrap(facets = vars(ageGroup), scales = "free_y") +
    scale_y_continuous(name = "Incidence per 100.000") +
    scale_x_date(name = "Time [Months]") +
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
    scale_y_continuous(name = "Incidence per 100.000") +
    scale_colour_manual(name = "Age group", values = dtuPalette) +
    guides(colour = guide_legend(nrow = 1)) +
    ggtitle(label = l)
  ggsave(filename = paste0(l,"xCaseDef.png"), path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
}



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


tmp <- dat %>% 
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>% 
  group_by(Date, ageGroup) %>%
  reframe(cases = sum(cases), n = sum(n)) %>%
  left_join(STEC_farrington_tbl)

dat %>% 
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>% 
  group_by(Date, ageGroup) %>%
  reframe(cases = sum(cases), n = sum(n)) %>%
  left_join(STEC_farrington_tbl) %>%
  mutate(alarm = if_else(is.na(alarm), FALSE, alarm)) %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_line(mapping = aes(y = cases/n*1e5), size = 0.4) +
  geom_point(mapping = aes(y = cases/n*1e5, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = Date, y = threshold/n*1e5), lty = "dashed", size = 0.4, inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000") +
  scale_x_date(name = "Time [Months]") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Farrington method") +
  annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = as.Date("2011-03-01"), ymin = -Inf, ymax = Inf, alpha = 0.2)
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

dat %>% 
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>% 
  group_by(Date, ageGroup) %>%
  reframe(cases = sum(cases), n = sum(n)) %>%
  left_join(STEC_noufaily_tbl) %>%
  mutate(alarm = if_else(is.na(alarm), FALSE, alarm)) %>%
  ggplot(mapping = aes(x = Date, colour = ageGroup)) +
  geom_line(mapping = aes(y = cases/n*1e5), size = 0.4) +
  geom_point(mapping = aes(y = cases/n*1e5, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = Date, y = threshold/n*1e5), lty = "dashed", size = 0.4, inherit.aes = FALSE) +
  # geom_rug(mapping = aes(x = dateOfAlarm, y = NULL), outside = TRUE, sides = "b", inherit.aes = FALSE) +
  # coord_cartesian(clip = "off") +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000") +
  scale_x_date(name = "Time [Months]") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        # axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Noufaily method") +
  annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = as.Date("2011-03-01"), ymin = -Inf, ymax = Inf, alpha = 0.2)
ggsave(filename = "STEC_noufaily.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")



# STEC Result -----------------------------------------------------------------------

STEC_res <- read_rds(file = "../models/STEC_res.rds")

STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(Parameter == "log_sigma") %>%
  ggplot(mapping = aes(x = ref.date, y = exp(theta))) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = exp(theta - 2*se.theta), ymax = exp(theta + 2*se.theta)), inherit.aes = FALSE, alpha = 0.2) +
  ylab(label = expression(sigma)) +
  xlab(label = "Time [Months]") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
ggsave(filename = "phiSTECPoisNExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")


# PoissonGamma %>%
#   select(ref.date, par) %>%
#   unnest(par) %>%
#   mutate(Parameter = factor(Parameter)) %>%
#   filter(Parameter == "phi") %>%
#   ggplot(mapping = aes(x = ref.date, y = exp(theta))) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(mapping = aes(x= ref.date, ymin = exp(theta) / exp(2*se.theta), ymax = exp(theta) * exp(2*se.theta)), inherit.aes = FALSE, alpha = 0.2) +
#   ylab(label = expression(sigma)) +
#   xlab(label = "Time [Months]") +
#   ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
# 
# PoissonNormal %>%select(ref.date, par) %>%
#   unnest(par) %>%
#   mutate(Parameter = factor(Parameter)) %>%
#   filter(Parameter == "log_sigma") %>%
#   print(n = 144)

STEC_res %>%
  select(PoissonNormal_seasonal) %>%
  unnest(PoissonNormal_seasonal) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(Parameter == "log_sigma") %>%
  ggplot(mapping = aes(x = ref.date, y = exp(theta))) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = exp(theta) / exp(2*se.theta), ymax = exp(theta) * exp(2*se.theta)), inherit.aes = FALSE, alpha = 0.2) +
  ylab(label = expression(sigma)) +
  xlab(label = "Time [Months]") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
ggsave(filename = "phiSTECPoisNSeasonal.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

custom_labeller <- as_labeller(
  c(`ageGroup<1 year`="beta[1~year]", `ageGroup1-4 years`="beta[1-4~years]",
    `ageGroup5-14 years`="beta[5-14~years]",`ageGroup15-24 years`="beta[15-24~years]",
    `ageGroup25-64 years`="beta[25-64~years]", `ageGroup65+ years`="beta[65+~years]"),
  default = label_parsed
  )

STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(Parameter != "log_sigma") %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_colour_manual(values = dtuPalette) +
  xlab(label = "Time [Months]") +
  ylab(label = "") +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
ggsave(filename = "thetaSTECPoisNExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

custom_labeller_seasonal <- as_labeller(
  c(`ageGroup<1 year`="beta[1~year]", `ageGroup1-4 years`="beta[1-4~years]",
    `ageGroup5-14 years`="beta[5-14~years]",`ageGroup15-24 years`="beta[15-24~years]",
    `ageGroup25-64 years`="beta[25-64~years]", `ageGroup65+ years`="beta[65+~years]",
    `sin(2 * pi/12 * monthInYear)`="beta[sin]", `cos(2 * pi/12 * monthInYear)`="beta[cos]"),
  default = label_parsed
  )
  
STEC_res %>%
    select(PoissonNormal_seasonal) %>%
    unnest(PoissonNormal_seasonal) %>%
    ungroup() %>%
    select(ref.date, par) %>%
    unnest(par) %>%
    mutate(Parameter = factor(Parameter)) %>%
    filter(grepl(x = Parameter, patter = "ageGroup")) %>%
    ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
    geom_line(linewidth = 1) +
    geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
    facet_wrap(facets = vars(Parameter), labeller = custom_labeller_seasonal) +
    scale_colour_manual(values = dtuPalette) +
    xlab(label = "Time [Months]") +
    ylab(label = "") +
    guides(colour = "none") +
    ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
ggsave(filename = "thetaSTECPoisNSeasonal.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonNormal_seasonal) %>%
  unnest(PoissonNormal_seasonal) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(grepl(x = Parameter, patter = "cos|sin")) %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller_seasonal) +
  scale_colour_manual(values = dtuPalette) +
  xlab(label = "Time [Months]") +
  ylab(label = "") +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierachical Poisson Normal model")
ggsave(filename = "thetaSTECPoisNSeasonalCosSin.png",
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
  select(ref.date, par) %>%
  unnest(par) %>%
  filter(Parameter == "phi") %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  ylab(label = expression(phi)) +
  xlab(label = "Time [Months]") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "phiSTECPoisGExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonGamma_seasonal) %>%
  unnest(PoissonGamma_seasonal) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  filter(Parameter == "phi") %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  ylab(label = expression(phi)) +
  xlab(label = "Time [Months]") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "phiSTECPoisGSeasonal.png",
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
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(Parameter != "phi") %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller) +
  scale_colour_manual(values = dtuPalette) +
  xlab(label = "Time [Months]") +
  ylab(label = "") +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "thetaSTECPoisGExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonGamma_seasonal) %>%
  unnest(PoissonGamma_seasonal) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(grepl(x = Parameter, pattern = "ageGroup")) %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller_seasonal) +
  scale_colour_manual(values = dtuPalette) +
  xlab(label = "Time [Months]") +
  ylab(label = "") +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "thetaSTECPoisGSeasonal.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

STEC_res %>%
  select(PoissonGamma_seasonal) %>%
  unnest(PoissonGamma_seasonal) %>%
  ungroup() %>%
  select(ref.date, par) %>%
  unnest(par) %>%
  mutate(Parameter = factor(Parameter)) %>%
  filter(grepl(x = Parameter, pattern = "cos|sin")) %>%
  ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
  geom_line(linewidth = 1) +
  geom_ribbon(mapping = aes(x= ref.date, ymin = theta - 2*se.theta, ymax = theta + 2*se.theta), inherit.aes = FALSE, alpha = 0.2) +
  facet_wrap(facets = vars(Parameter), labeller = custom_labeller_seasonal) +
  scale_colour_manual(values = dtuPalette) +
  xlab(label = "Time [Months]") +
  ylab(label = "") +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compound Poisson Gamma model")
ggsave(filename = "thetaSTECPoisGSeasonalCosSin.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

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
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, u, phi, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, ref.date, NA)) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_x_date(name = "Time [Months]") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
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
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, u, phi, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, ref.date, NA)) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_x_date(name = "Time [Months]") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
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
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, u, log_sigma, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, ref.date, NA)) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_x_date(name = "Time [Months]") +
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none",shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
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
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(ref.date, ageGroup, u, log_sigma, alarm) %>%
  mutate(dateOfAlarm = if_else(alarm, ref.date, NA)) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup)) +
  geom_line(mapping = aes(y = u), linewidth = 0.4) +
  geom_point(mapping = aes(y = u, shape = alarm), size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  geom_rug(mapping = aes(x = dateOfAlarm, y = NULL),
           outside = TRUE, sides = "b", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[it]))) +
  scale_x_date(name = "Time [Months]") + 
  scale_colour_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none",shape = "none") +
  theme(panel.spacing.y = unit(1, "lines"), 
        axis.text.x = element_text(vjust = -1.2)) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Hierarchical Poisson Normal model")
ggsave(filename = "windowedSTEDPoisNExclude.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# STEC_res %>%
#   select(PoissonNormal_excludePastOutbreaks) %>%
#   unnest(PoissonNormal_excludePastOutbreaks) %>% 
#   ungroup() %>%
#   unnest(par) %>%
#   select(ref.date, Parameter, theta, se.theta) %>%
#   filter(Parameter != "log_sigma") %>%
#   ggplot(mapping = aes(x = ref.date, y = theta, colour = Parameter)) +
#   geom_line() +
#   facet_wrap(facets = vars(Parameter)) +
#   scale_colour_manual(values = dtuPalette) +
#   guides(colour = "none")

# STEC_res %>%
#   select(PoissonNormal_excludePastOutbreaks) %>%
#   unnest(PoissonNormal_excludePastOutbreaks) %>% 
#   ungroup() %>%
#   unnest(par) %>%
#   select(ref.date, Parameter, ageGroup, theta, se.theta) %>%
#   filter(ageGroup == "All") %>%
#   ggplot(mapping = aes(x = ref.date, y = exp(theta), colour = ageGroup)) +
#   geom_ribbon(mapping = aes(x = ref.date, ymin = exp(theta - 2 * se.theta), ymax = exp(theta + 2 * se.theta)), alpha = 0.4, fill = "grey70", inherit.aes = FALSE) +
#   geom_line() +
#   scale_colour_manual(values = dtuPalette) +
#   guides(colour = "none")
  


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
  scale_y_continuous(name = "Incidence per 100.000") +
  scale_x_discrete(name = "Month") +
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
  


STEC_farrington_comp <- STEC_farrington_tbl %>%
  select(Date, ageGroup, alarm) %>%
  mutate(method = "Farrington")
STEC_noufaily_comp <- STEC_noufaily_tbl %>%
  select(Date, ageGroup, alarm)%>%
  mutate(method = "Noufaily")
STEC_PoisN_comp <- STEC_res %>%
  select(PoissonNormal_excludePastOutbreaks) %>%
  unnest(PoissonNormal_excludePastOutbreaks) %>% 
  ungroup() %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(Date = ref.date, ageGroup, alarm) %>%
  mutate(method = "Poisson Normal")
STEC_PoisG_comp <- STEC_res %>%
  select(PoissonGamma_excludePastOutbreaks) %>%
  unnest(PoissonGamma_excludePastOutbreaks) %>% 
  ungroup() %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  select(Date = ref.date, ageGroup, alarm) %>%
  mutate(method = "Poisson Gamma")

compare <- bind_rows(
  STEC_PoisG_comp,
  STEC_PoisN_comp,
  STEC_noufaily_comp, 
  STEC_farrington_comp) %>%
  mutate(method = factor(method,
                         levels = c("Farrington",
                                    "Noufaily",
                                    "Poisson Normal",
                                    "Poisson Gamma")))

comp_dat <- dat %>% 
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>% 
  group_by(Date, ageGroup) %>%
  reframe(cases = sum(cases), n = sum(n)) 


compare %>% 
  mutate(alarmDate = if_else(alarm, Date, as.Date(NA_character_))) %>%
  ggplot(mapping = aes(x = method, y = alarmDate, colour = "one")) +
  geom_point(shape = 17) +
  coord_flip() +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = "#E83F48") +
  scale_x_discrete(limits = rev(levels(compare$method)), position = "top") +
  scale_y_date(name = "Time [Months]") + 
  guides(colour = "none") +
  theme(axis.title.y = element_blank()) +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.", subtitle = "Compare methods")
ggsave(filename = "compareMethods.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

library(cowplot)

scaleX <- scale_x_date(
  name = "Time [Months]",
  limits = as.Date(c("2010-11-01", "2023-02-01")),
  expand = expansion(add = c(0)))

alarm1 <- compare %>% 
  mutate(alarmDate = if_else(alarm, Date, as.Date(NA_character_))) %>%
  filter(ageGroup %in% c("<1 year", "1-4 years", "5-14 years")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(name = "Alarms:", values = dtuPalette[c(7,9:11)]) +
  scale_y_discrete(limits = rev(levels(compare$method))) +
  scaleX + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))

alarm2 <- compare %>% 
  mutate(alarmDate = if_else(alarm, Date, as.Date(NA_character_))) %>%
  filter(ageGroup %in% c("15-24 years", "25-64 years", "65+ years")) %>%
  ggplot(mapping = aes(x = alarmDate, y = method, colour = method)) +
  geom_point(shape = 17) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette[c(7,9:11)]) +
  scale_y_discrete(limits = rev(levels(compare$method))) +
  scaleX +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(2.67, "lines"))

# alarm1 <- compare %>% 
#   mutate(alarmDate = if_else(alarm, Date, as.Date(NA_character_))) %>%
#   filter(ageGroup %in% c("<1 year", "1-4 years", "5-14 years")) %>%
#   ggplot(mapping = aes(x = method, y = alarmDate, colour = method)) +
#   geom_point(shape = 17) +
#   coord_flip() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_color_manual(name = "Alarms", values = dtuPalette[7:10]) +
#   scale_x_discrete(limits = rev(levels(compare$method)), position = "top") +
#   scale_y_date(name = "Time [Months]") + 
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# alarm2 <- compare %>% 
#   mutate(alarmDate = if_else(alarm, Date, as.Date(NA_character_))) %>%
#   filter(ageGroup %in% c("15-24 years", "25-64 years", "65+ years")) %>%
#   ggplot(mapping = aes(x = method, y = alarmDate, colour = method)) +
#   geom_point(shape = 17) +
#   coord_flip() +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_color_manual(values = dtuPalette[7:10]) +
#   scale_x_discrete(limits = rev(levels(compare$method)), position = "top") +
#   scale_y_date(name = "Time [Months]") + 
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())

library(scales)

cases1 <- compare %>%
  left_join(comp_dat, by = c("Date", "ageGroup")) %>%
  filter(ageGroup %in% c("<1 year", "1-4 years", "5-14 years")) %>%
  ggplot(mapping = aes(x = Date, y = cases/n*1e5, colour = ageGroup, shape = alarm, group = ageGroup)) +
  geom_line(size = 0.4) +
  geom_point(size = 2) +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.01)) +
  scaleX +
  scale_colour_manual(values = dtuPalette[1:3]) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  # theme(axis.title.y = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       panel.spacing.x = unit(2, "lines"))
  theme(axis.title.y = element_text(size = 12),
        panel.spacing.x = unit(.8, "lines"))
cases2 <- compare %>%
  left_join(comp_dat, by = c("Date", "ageGroup")) %>%
  filter(ageGroup %in% c("15-24 years", "25-64 years", "65+ years")) %>%
  ggplot(mapping = aes(x = Date, y = cases/n*1e5, colour = ageGroup, shape = alarm, group = ageGroup)) +
  geom_line(size = 0.4) +
  geom_point(size = 2) +
  facet_wrap(facets = vars(ageGroup), scales = "free_y") +
  scale_y_continuous(name = "Incidence per 100.000", labels = label_number(accuracy = 0.01)) +
  scaleX +
  scale_colour_manual(values = dtuPalette[4:6]) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none") +
  # theme(axis.title.y = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank(),
  #       panel.spacing.x = unit(2, "lines"))
  theme(axis.title.y = element_text(size = 12),
        panel.spacing.x = unit(.8, "lines"))

legend <- get_legend(alarm1)

aligned_plots <- align_plots(
  cases1 + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank()),
  alarm1 + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 strip.text = element_blank(),
                 legend.position = "none"),
  cases2 + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank()),
  alarm2 + theme(strip.text = element_blank(),
                 legend.position = "none"),
  align = "hv",
  axis = "tblr")



plot_grid(legend,
          aligned_plots[[1]], 
          NULL,
          aligned_plots[[2]], 
          NULL,
          aligned_plots[[3]], 
          NULL,
          aligned_plots[[4]], 
          ncol = 1,rel_heights = c(0.05,0.32,-.045,0.155,-.02,0.32,-.045,0.155))
ggsave(filename = "compareMethodsDeluxe.png",
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")




