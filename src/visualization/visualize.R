
# Import libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load data
dat <- read_rds(file = "../../data/processed/dat.rds")

# Extract case definitions
caseDef <- unique(dat$caseDef)

# Loop over cases and create overview plots
for(c in caseDef){
  dat %>%
    filter(caseDef == c) %>%
    ggplot(mapping = aes(x = maaned, y = cases/n * 1e5, colour = year)) +
    geom_point() +
    facet_grid(rows = vars(ageLabel), cols = vars(landsdel)) +
    scale_y_continuous(name = "Number of cases per 100.000") +
    scale_x_discrete(name = "Month") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top") +
    ggtitle(label = c)
  ggsave(filename = c, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")
}

# An example plot of 'VTEC' in 'Nordsjælland'

dat %>%
  filter(landsdel == "Nordsjælland") %>%
  ggplot(mapping = aes(x = maaned, y = cases/n * 1e5, colour = caseDef)) +
  geom_point() +
  facet_grid(rows = vars(ageLabel), cols = vars(year))
