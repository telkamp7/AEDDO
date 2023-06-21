
# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

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


# Find simulations files
filesInDir <- list.files()
indexFiles <- grepl(pattern = "scenario...Rds|scenario..Rds", x = filesInDir)
scenarioFiles <- filesInDir[indexFiles]
scenarioFiles <- scenarioFiles[scenarioFiles!="scenarios.Rds"]

# Load in simulation files and unpack
results <- tibble()
for(i in 1:length(scenarioFiles)){
  iterData <- read_rds(file = scenarioFiles[i])
  for(j in 1:length(iterData$Data)){
    
    unpackData <- iterData$Data[[j]][[1]]
    
    results <- bind_rows(results,
                         unpackData)
  }
}

Realizations <- results %>%
  filter(scenario %in% c(5, 7,  12, 28) & sim == 1) %>%
  select(scenario, data) %>%
  unnest(data) %>%
  ggplot(mapping = aes(x = t, y = y)) +
  geom_line() +
  facet_wrap(facets = vars(scenario), scales = "free_y") +
  scale_x_continuous(name = "week") +
  scale_y_continuous(name = "count")
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
  select(sim, scenario,
         FP_Farrington, FP_Noufaily, FP_PoisN, FP_PoisG,
         TN_Farrington, TN_Noufaily, TN_PoisN, TN_PoisG) %>%
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
                           labels = 1:28)) %>%
  ggplot(mapping = aes(x = scenario, y = FPR, group = scenario, fill = Method)) +
  geom_boxplot(fatten = 3, alpha = 0.7) +
  facet_wrap(facets = vars(Method)) +
  scale_x_discrete(name = "Scenario") +
  scale_fill_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(fill = "none")
ggsave(filename = "FPR.png",
       plot = FPR,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

badPerformanceScenarios <- results %>%
  arrange(scenario) %>%
  select(sim, scenario, data) %>%
  unnest(data) %>% 
  filter(t %in% 576:624 & k == 10) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  mutate(Method = factor(Method,
                          levels = c("alarm_Farrington",
                                     "alarm_Noufaily",
                                     "alarm_PoisN",
                                     "alarm_PoisG"),
                          labels = c("Farrington",
                                     "Noufaily", 
                                     "'Poisson Normal'", 
                                     "'Poisson Gamma'"))) %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(scenario, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  group_by(scenario) %>%
  mutate(badScenario = any(POD < 0.6)) %>% 
  filter(badScenario)

# Probability of detection
PropDetect <- results %>%
  arrange(scenario) %>%
  select(sim, scenario, data) %>%
  unnest(data) %>% 
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
                                                             "'Poisson Normal'", 
                                                             "'Poisson Gamma'"))) %>%
  ggplot(mapping = aes(x = k, colour = Method)) +
  geom_line(mapping = aes(y = POD, group = scenario), alpha = 0.6) +
  geom_line(mapping = aes(y = medianPOD), linewidth=1.2) +
  geom_text(data = badPerformanceScenarios, mapping = aes(x = k, y = POD, label = scenario), inherit.aes = FALSE) +
  facet_wrap(facets = vars(Method), labeller = label_parsed) +
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



results %>%
  arrange(scenario) %>%
  select(sim, scenario, data) %>%
  unnest(data) %>% 
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(scenario, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  filter(Method == "alarm_Noufaily" & k == 10) %>%
  print(n = 28)

results %>%
  arrange(scenario) %>%
  select(sim, scenario, data) %>%
  unnest(data) %>% 
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(scenario, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  filter(scenario == 7 & k == 10)

results %>%
  filter(scenario == 7) %>%
  unnest(data) %>%
  filter(t %in% 576:624) %>% 
  select(sim, t, y) %>%
  ggplot(mapping = aes(x = t, y = y)) + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.05), fun.max = function(x) quantile(x, 0.95),alpha = 0.1, fill ="#2F3EEA") + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.1), fun.max = function(x) quantile(x, 0.9), alpha = 0.15, fill ="#2F3EEA") +
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.15), fun.max = function(x) quantile(x, 0.85), alpha = 0.2, fill ="#2F3EEA") +
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.2), fun.max = function(x) quantile(x, 0.8), alpha = 0.25, fill ="#2F3EEA") +
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.25), fun.max = function(x) quantile(x, 0.75), alpha = 0.3, fill ="#2F3EEA") + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.3), fun.max = function(x) quantile(x, 0.7), alpha = 0.35, fill ="#2F3EEA") + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.35), fun.max = function(x) quantile(x, 0.65), alpha = 0.4, fill ="#2F3EEA") + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.4), fun.max = function(x) quantile(x, 0.6), alpha = 0.45, fill ="#2F3EEA") + 
  stat_summary(geom="ribbon", fun.min = function(x) quantile(x, 0.45), fun.max = function(x) quantile(x, 0.55), alpha = 0.5, fill ="#2F3EEA") + 
  stat_summary(geom="line", fun=median, linewidth=1.2, colour = "#030F4F")

results %>%
  filter(scenario == 7) %>%
  unnest(data) %>%
  filter(t %in% 576:624) %>% 
  select(sim, t, y) %>%
  ggplot(mapping = aes(x = t, y = y)) +
  geom_bin2d(binwidth = c(1,1))

results %>%
  filter(scenario == 7) %>%
  unnest(data) %>%
  filter(t %in% 576:624) %>% 
  select(sim, t, y) %>%
  ggplot(mapping = aes(x = y, y = after_stat(density))) +
  geom_density(alpha = 0.25)

results %>%
  arrange(scenario) %>%
  select(sim, scenario, data) %>%
  unnest(data) %>% 
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(sim, scenario, Method, k) %>%
  reframe(Detected = outbreakTF == TRUE & Alarms == TRUE) %>%
  mutate(POD = Detected / (Detected + Undetected)) %>%
  ggplot(mapping = aes(x = k, y = POD, group = sim)) +
  geom_line() +
  facet_wrap(facets = vars(Method))
