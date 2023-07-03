
# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set global theme options
theme_set(
  new = theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 20),
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


# Find simulations files
filesInDir <- list.files()
indexFiles <- grepl(pattern = "scenario...Rds|scenario..Rds", x = filesInDir)
scenarioFiles <- filesInDir[indexFiles]
scenarioFiles <- scenarioFiles[scenarioFiles!="scenarios.Rds"]

# Load in simulation files and unpack
results <- tibble()
for(i in 1:length(scenarioFiles)){
  scenario_no <- parse_number(scenarioFiles[i])
  iterData <- read_rds(file = scenarioFiles[i])
  for(j in 1:length(iterData$Data)){
    
    unpackData <- iterData$Data[[j]][[1]] %>%
      mutate(scenario = scenario_no, sim = j)
    
    results <- bind_rows(results,
                         unpackData) 
  }
}

Realizations <- results %>%
  filter(scenario %in% c(8,12,13,20) & sim == 1) %>%
  select(scenario, outbreaks) %>%
  unnest(outbreaks) %>%
  mutate(outbreakDate = if_else(outbreakTF, t, NA_integer_)) %>%
  ggplot(mapping = aes(x = t, y = y)) +
  geom_line(alpha = 0.5) +
  geom_point(mapping = aes(x = outbreakDate)) +
  facet_wrap(facets = vars(scenario), scales = "free_y") +
  scale_x_continuous(name = "Week") +
  scale_y_continuous(name = "Number of cases") +
  annotate(geom = "rect", xmin = 575, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2)
ggsave(filename = "Realizations.png",
       plot = Realizations,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

results %>%
  arrange(scenario) %>%
  filter(row_number() == 1)

# False positive rate
FPR <- results %>%
  arrange(scenario) %>%
  select(sim, scenario, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(sim, scenario) %>%
  slice_tail(n = 49) %>%
  reframe(FP_Farrington = sum(alarm_Farrington == TRUE  & outbreakTF == FALSE),
          TN_Farrington = sum(alarm_Farrington == FALSE & outbreakTF == FALSE),
          FP_Noufaily   = sum(alarm_Noufaily == TRUE    & outbreakTF == FALSE),
          TN_Noufaily   = sum(alarm_Noufaily == FALSE   & outbreakTF == FALSE),
          FP_PoisN      = sum(alarm_PoisN == TRUE       & outbreakTF == FALSE),
          TN_PoisN      = sum(alarm_PoisN == FALSE      & outbreakTF == FALSE),
          FP_PoisG      = sum(alarm_PoisG == TRUE       & outbreakTF == FALSE),
          TN_PoisG      = sum(alarm_PoisG == FALSE      & outbreakTF == FALSE)) %>%
  pivot_longer(cols = FP_Farrington:TN_PoisG,
               names_to = c("Statistic", "Method"),
               names_pattern = "(\\w+)_(\\w+)",
               values_to = "value") %>%
  pivot_wider(names_from = Statistic, values_from = value) %>%
  group_by(sim, scenario, Method) %>%
  reframe(FPR = FP/(FP+TN)) %>%
  mutate(Method = factor(Method,levels = c("Farrington",
                                           "Noufaily",
                                           "PoisN",
                                           "PoisG"),
                         labels = c("Farrington",
                                    "Noufaily", 
                                    "Poisson Normal", 
                                    "Poisson Gamma")),
         scenario = factor(scenario, levels = 1:28,
                           labels = 1:28)) 
write_rds(x = FPR, file = "FPR.Rds")

FPRPlot <- FPR %>%
  ggplot(mapping = aes(x = scenario, y = FPR, group = scenario, fill = Method)) +
  geom_boxplot(fatten = 3, alpha = 0.7) +
  facet_wrap(facets = vars(Method)) +
  scale_x_discrete(name = "Scenario") +
  scale_fill_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(size = 16))
ggsave(filename = "FPRPlot.png",
       plot = FPRPlot,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

badPerformanceScenarios <- results %>%
  arrange(scenario) %>%
  select(sim, scenario, outbreaks) %>%
  unnest(outbreaks) %>% 
  filter(t %in% 576:624 & k == 10) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  mutate(Method = factor(Method,
                          levels = c("alarm_Farrington",
                                     "alarm_Noufaily",
                                     "alarm_PoisN",
                                     "alarm_PoisG"),
                          labels = c("Farrington",
                                     "Noufaily", 
                                     "Poisson Normal", 
                                     "Poisson Gamma"))) %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(scenario, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  group_by(scenario) %>%
  mutate(badScenario = any(POD < 0.6)) %>% 
  filter(badScenario)

# Probability of detection
POD <- results %>%
  arrange(scenario) %>%
  select(sim, scenario, outbreaks) %>%
  unnest(outbreaks) %>% 
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(scenario, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  group_by(Method, k) %>%
  mutate(medianPOD = median(POD), Method = factor(Method,
                                                  levels = c("alarm_Farrington",
                                                             "alarm_Noufaily",
                                                             "alarm_PoisN",
                                                             "alarm_PoisG"),
                                                  labels = c("Farrington",
                                                             "Noufaily", 
                                                             "Poisson Normal", 
                                                             "Poisson Gamma")))
write_rds(x = POD, file = "POD.Rds")

PropDetect <- POD %>% 
  ggplot(mapping = aes(x = k, colour = Method)) +
  geom_line(mapping = aes(y = POD, group = scenario), alpha = 0.6) +
  geom_line(mapping = aes(y = medianPOD), linewidth=2) +
  # geom_text(data = badPerformanceScenarios, mapping = aes(x = k, y = POD, label = scenario), inherit.aes = FALSE) +
  facet_wrap(facets = vars(Method)) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(color = "none") 
ggsave(filename = "PropDetect.png",
       plot = PropDetect,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")



POD %>% 
  filter(scenario %in% 13:16) %>%
  ggplot(mapping = aes(x = k, colour = Method)) +
  geom_line(mapping = aes(y = POD, group = scenario), alpha = 0.6) +
  # geom_line(mapping = aes(y = medianPOD), linewidth=2) +
  # geom_text(data = badPerformanceScenarios, mapping = aes(x = k, y = POD, label = scenario), inherit.aes = FALSE) +
  facet_wrap(facets = vars(Method), labeller = label_parsed) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(color = "none") 








