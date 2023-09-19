
# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
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
indexFiles <- grepl(pattern = "alphaNovel", x = filesInDir)
scenarioFiles <- filesInDir[indexFiles]

# Load in simulation files and unpack
results <- tibble()
for(i in 1:length(scenarioFiles)){
  alpha_num <- as.numeric(str_extract(scenarioFiles[i], "\\d+\\.\\d+"))
  iterData <- read_rds(file = scenarioFiles[i])
  for(j in 1:length(iterData$Data)){
    
    unpackData <- iterData$Data[[j]][[1]] %>%
      mutate(alpha = alpha_num, sim = j)
    
    results <- bind_rows(results,
                         unpackData) 
  }
}

# False positive rate
FPR <- results %>%
  arrange(alpha) %>%
  select(sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(sim, alpha) %>%
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
  group_by(sim, alpha, Method) %>%
  reframe(FPR = FP/(FP+TN)) %>%
  mutate(Method = factor(Method,levels = c("Farrington",
                                           "Noufaily",
                                           "PoisN",
                                           "PoisG"),
                         labels = c("Farrington",
                                    "Noufaily", 
                                    "Poisson Normal", 
                                    "Poisson Gamma"))) 
write_rds(x = FPR, file = "FPRalpha.Rds")

# Probability of detection
POD <- results %>%
  arrange(alpha) %>%
  select(sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>% 
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(sim, alpha, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>% 
  group_by(alpha, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  group_by(alpha, Method, k) %>%
  mutate(medianPOD = median(POD), Method = factor(Method,
                                                  levels = c("alarm_Farrington",
                                                             "alarm_Noufaily",
                                                             "alarm_PoisN",
                                                             "alarm_PoisG"),
                                                  labels = c("Farrington",
                                                             "Noufaily", 
                                                             "Poisson Normal", 
                                                             "Poisson Gamma")))
write_rds(x = POD, file = "POD_alpha.Rds")

# Making beautiful plot
profilePODxFPR_num <- POD %>% 
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  filter(k == 5) %>%
  ggplot(mapping = aes(x = meanFPR, y = medianPOD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_label(aes(label = alpha)) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) 
ggsave(filename = "profilePODxFPR_num.png",
       plot = profilePODxFPR_num,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

profilePODxFPR_shape <- POD %>% 
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  filter(k == 5) %>%
  ggplot(mapping = aes(x = meanFPR, y = medianPOD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_point(mapping = aes(x = meanFPR, y = medianPOD, shape = factor(alpha)), size = 3) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  scale_shape_discrete(name = "alpha") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
ggsave(filename = "profilePODxFPR_shape.png",
       plot = profilePODxFPR_shape,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

profilePODxFPR_facet <- POD %>% 
  filter(k > 2) %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  ggplot(mapping = aes(x = meanFPR, y = medianPOD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_point(mapping = aes(x = meanFPR, y = medianPOD, shape = factor(alpha)), size = 3) +
  facet_wrap(facets = vars(k), ncol = 2) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  scale_shape_discrete(name = "alpha") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
ggsave(filename = "profilePODxFPR_facet.png",
       plot = profilePODxFPR_facet,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

scenarioIllustration <- tibble(x = 1:100, constant = 6, trend = 4 + 0.03*x, seasonality = 5 + sin(2*pi*x/50) + cos(2*pi*x/50), combined = 5 + 0.03*x + sin(2*pi*x/50) + cos(2*pi*x/50)) %>%
  pivot_longer(cols = constant:combined) %>%
  mutate(name = factor(name, 
                       levels = c("constant", "trend", "seasonality", "combined"),
                       labels = c("constant", "trend", "seasonality", "combined"))) %>%
  ggplot(mapping = aes(x = x, y = value)) +
  geom_line(size = 2) +
  facet_wrap(facets = vars(name)) +
  theme_bw() +
  theme(strip.text = element_text(size = 26),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26))
ggsave(filename = "scenarioIllustration.png",
       plot = scenarioIllustration,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")
