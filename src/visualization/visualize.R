
# Import libraries
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# Set global theme options
theme_set(
  new = theme_bw() +
    theme(legend.position = "top")
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
    scale_x_date() +
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
    scale_x_date() +
    scale_colour_manual(name = "Age group", values = dtuPalette) +
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
    scale_x_date() +
    scale_colour_manual(name = "Age group", values = dtuPalette) +
    guides(colour = guide_legend(nrow = 1)) +
    ggtitle(label = l)
  ggsave(filename = paste0(l,"xCaseDef.png"), path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
}
