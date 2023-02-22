
# Import libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(TMB)

# Dynamically link the C++ template
dyn.load(dynlib(name = "../models/PoissonLognormal"))

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
dat <- read_rds(file = "../../data/processed/dat.rds")

# Only consider some of the data
y <- dat %>%
  filter(caseDef == "Shiga- og veratoxin producerende E. coli.") %>%
  group_by(Date, ageGroup) %>%
  mutate(y = sum(cases)) %>%
  select(Date, ageGroup, y, n)

# Extract case definitions
caseDef <- unique(dat$caseDef)

# Loop over case definitions and create overview plots
for(c in caseDef){
  
  # caseDefxLandsdel
  dat %>%
    filter(caseDef == c) %>%
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

# Load the Poisson-lognormal model
PoisLN <- read_rds(file = "../models/PoissonLognormal.rds")
# ... and generate report
rep <- sdreport(PoisLN, getJointPrecision = TRUE)

PoisLN.res <- y %>%
  ungroup() %>%
  mutate(`Random effects` = rep$par.random,
         `Std. Error` = sqrt(rep$diag.cov.random))

PoisLN.res %>%
  ggplot(mapping = aes(x = Date, y = `Random effects`, colour = ageGroup)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = vars(ageGroup)) +
  scale_y_continuous(name = expression(paste("Random effect, ", u[t]^a))) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none") +
  ggtitle(label = "Shiga- og veratoxin producerende E. coli.")
ggsave(filename = "VTECxRandomEffects.png", path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")


