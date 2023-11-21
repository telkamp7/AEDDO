
# Setup -------------------------------------------------------------------


# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(broom)

# Set working directory
setwd("src/simulation/")

# Set global theme options
theme_set(
  new = theme_bw() +
    theme(legend.position = "top",
          strip.text = element_text(size = 20),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 30),
          legend.key.size = unit(x = 2, units = "cm"))
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
indexFiles <- grepl(pattern = "scenario._alphaState", x = filesInDir)
indexFiles2 <- grepl(pattern = "scenario.._alphaState", x = filesInDir)
scenarioFiles <- filesInDir[indexFiles|indexFiles2]


# Load in data ------------------------------------------------------------

# Load in simulation files and unpack
results <- tibble()
for(i in 1:length(scenarioFiles)){
  alpha_num <- as.numeric(str_extract(scenarioFiles[i], "\\d+\\.\\d+"))
  scenario_num <- as.numeric(str_extract(scenarioFiles[i], "\\d+"))
  iterData <- read_rds(file = scenarioFiles[i])
  for(j in 1:length(iterData$Data)){

    unpackData <- iterData$Data[[j]][[1]] %>%
      mutate(scenario = scenario_num,
             alpha = alpha_num,
             sim = j)

    results <- bind_rows(results,
                         unpackData)
  }
}

# Unnest the data ---------------------------------------------------------

baseline_results <- results %>%
  arrange(scenario, alpha, sim) %>%
  select(scenario, sim, alpha, baseline) %>%
  unnest(baseline) %>%
  group_by(scenario, sim, alpha) %>%
  slice_tail(n = 49) %>%
  pivot_longer(
    cols = alarm_Farrington:alarm_PoisG,
    names_to = "method",
    values_to = "alarm",
    names_prefix = "alarm_") %>%
  group_by(scenario, alpha, method, sim) %>%
  reframe(
    FP = sum(alarm == TRUE),
    TN = sum(alarm == FALSE),
    specificity = TN/(TN+FP),
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Farrington",
                 "Noufaily",
                 "PoisN",
                 "PoisG"),
      labels = c("Farrington",
                 "Noufaily",
                 "Poisson Normal",
                 "Poisson Gamma")),
    alpha = factor(
      alpha,
      levels = c(0.005, 0.01, 0.025, 0.05, 0.1)
    ),
    scenario = factor(
      scenario,
      levels = 1:28
    ))

POD_results <- results %>%
  arrange(scenario, alpha, sim) %>%
  select(scenario, sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(scenario, sim, alpha, k) %>%
  slice_tail(n = 49) %>%
  pivot_longer(
    cols = alarm_Farrington:alarm_PoisG,
    names_to = "method",
    values_to = "alarm",
    names_prefix = "alarm_") %>%
  group_by(scenario, sim, alpha, method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & alarm == TRUE)) %>%
  group_by(scenario, alpha, method, k) %>%
  reframe(POD = mean(Detected)) %>%
  mutate(
    method = factor(
      method,
      levels = c("Farrington",
                 "Noufaily",
                 "PoisN",
                 "PoisG"),
      labels = c("Farrington",
                 "Noufaily",
                 "Poisson Normal",
                 "Poisson Gamma")),
    alpha = factor(
      alpha,
      levels = c(0.005, 0.01, 0.025, 0.05, 0.1)
    ),
    k = factor(
      k,
      levels = 1:10
    ),
    scenario = factor(
      scenario,
      levels = 1:28
    ))




outbreak_results <- results %>%
  arrange(scenario, alpha, sim) %>%
  select(scenario, sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(scenario, sim, alpha, k) %>%
  slice_tail(n = 49) %>%
  pivot_longer(
    cols = alarm_Farrington:alarm_PoisG,
    names_to = "method",
    values_to = "alarm",
    names_prefix = "alarm_") %>%
  group_by(scenario, alpha, k, method) %>%
  reframe(
    TP = sum(outbreakTF == TRUE & alarm == TRUE),
    FP = sum(outbreakTF == FALSE & alarm == TRUE),
    TN = sum(outbreakTF == FALSE & alarm == FALSE),
    FN = sum(outbreakTF == TRUE & alarm == FALSE),
    sensitivity = TP/(TP+FN),
    specificity = TN/(TN+FP),
    `LR+` = sensitivity/(1-specificity),
    `LR-` = (1-specificity)/sensitivity,
    DOR = `LR+`/`LR-`,
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("Farrington",
                 "Noufaily",
                 "PoisN",
                 "PoisG"),
      labels = c("Farrington",
                 "Noufaily",
                 "Poisson Normal",
                 "Poisson Gamma")),
    alpha = factor(
      alpha,
      levels = c(0.005, 0.01, 0.025, 0.05, 0.1)
    ),
    k = factor(
      k,
      levels = 1:10
    ),
    scenario = factor(
      scenario,
      levels = 1:28
    ))


# DOR ---------------------------------------------------------------------

outbreak_results_DOR <- outbreak_results %>%
  mutate(TP = TP + 0.5, FP = FP + 0.5, TN = TN + 0.5, FN = FN + 0.5,
         sensitivity = TP/(TP+FN),
         specificity = TN/(TN+FP),
         `LR+` = sensitivity/(1-specificity),
         `LR-` = (1-specificity)/sensitivity,
         DOR = `LR+`/`LR-`,
         logDOR = log(DOR))

# k = 3
medianDORk3 <- outbreak_results_DOR %>%
  filter(k == 3) %>%
  group_by(alpha, k, method) %>%
  reframe(medianlogDOR = median(logDOR, na.rm = TRUE))

logDORk3 <- outbreak_results_DOR %>%
  filter(k == 3) %>%
  ggplot(mapping = aes(x = alpha)) +
  geom_line(mapping = aes(y = logDOR, group = scenario, colour = method), alpha = 0.5) +
  geom_line(data = medianDORk3, mapping = aes(y = medianlogDOR, colour = method, group = method), linewidth = 2) +
  facet_wrap(facets = vars(method)) +
  scale_y_continuous(name = expression("log(DOR) (k = 3)")) +
  scale_x_discrete(name = expression(alpha)) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(colour = "none")
ggsave(filename = "logDORk3.png",
       plot = logDORk3,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# k = 5
medianDORk5 <- outbreak_results_DOR %>%
  filter(k == 5) %>%
  group_by(alpha, k, method) %>%
  reframe(medianlogDOR = median(log(DOR), na.rm = TRUE))

logDORk5 <- outbreak_results_DOR %>%
  filter(k == 5) %>%
  ggplot(mapping = aes(x = alpha)) +
  geom_line(mapping = aes(y = logDOR, group = scenario, colour = method), alpha = 0.5) +
  geom_line(data = medianDORk5, mapping = aes(y = medianlogDOR, colour = method, group = method), linewidth = 2) +
  facet_wrap(facets = vars(method)) +
  scale_y_continuous(name = expression("log(DOR) (k = 5)")) +
  scale_x_discrete(name = expression(alpha)) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(colour = "none")
ggsave(filename = "logDORk5.png",
       plot = logDORk5,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# k = 7
medianDORk7 <- outbreak_results_DOR %>%
  filter(k == 7) %>%
  group_by(alpha, k, method) %>%
  reframe(medianlogDOR = median(log(DOR), na.rm = TRUE))

logDORk7 <- outbreak_results_DOR %>%
  filter(k == 7) %>%
  ggplot(mapping = aes(x = alpha)) +
  geom_line(mapping = aes(y = logDOR, group = scenario, colour = method), alpha = 0.5) +
  geom_line(data = medianDORk7, mapping = aes(y = medianlogDOR, colour = method, group = method), linewidth = 2) +
  facet_wrap(facets = vars(method)) +
  scale_y_continuous(name = expression("log(DOR) (k = 7)")) +
  scale_x_discrete(name = expression(alpha)) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(colour = "none")
ggsave(filename = "logDORk7.png",
       plot = logDORk7,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")


custom_labeller <- as_labeller(
  c(`0.005` = "alpha[0.005]",
    `0.01` = "alpha[0.01]",
    `0.025` = "alpha[0.025]",
    `0.05` = "alpha[0.05]",
    `0.1` = "alpha[0.1]"),
  default = label_parsed
)

logDOR_alpha <- outbreak_results_DOR %>%
  group_by(alpha, k, method) %>%
  reframe(medianlogDOR = median(logDOR)) %>%
  ggplot(
    mapping = aes(
      x = k,
      y = medianlogDOR,
      colour = method,
      group = method)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(facets = vars(alpha), labeller = custom_labeller) +
  scale_y_continuous(name = expression("log(DOR)")) +
  scale_x_discrete(name = "k") +
  scale_color_manual(
    name = "",
    values = dtuPalette[c(7,9:11,5)])
ggsave(filename = "logDOR_alpha.png",
       plot = logDOR_alpha,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")



# FPR ---------------------------------------------------------------------

boxplot_labeller <- as_labeller(
  c(`0.005` = "alpha[0.005]",
    `0.01` = "alpha[0.01]",
    `0.025` = "alpha[0.025]",
    `0.05` = "alpha[0.05]",
    `0.1` = "alpha[0.1]",
    `Farrington` = "Farrington",
    `Noufaily` = "Noufaily",
    `Poisson Normal` = "'Poisson Normal'",
    `Poisson Gamma` = "'Poisson Gamma'"),
  default = label_parsed
)

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

nominal_values <- expand_grid(
  method = c("Farrington",
             "Noufaily",
             "Poisson Normal",
             "Poisson Gamma"),
  alpha_num = c(0.005, 0.01, 0.025, 0.05, 0.1)) %>%
  mutate(alpha = factor(alpha_num))

FPR_alpha_methods <- baseline_results %>%
  mutate(FPR = FP/(FP+TN)) %>%
  ggplot(mapping = aes(x = scenario, y = FPR, fill = method)) +
  geom_boxplot() +
  geom_hline(
    data = nominal_values,
    mapping = aes(yintercept = alpha_num),
    linetype = "dashed",
    linewidth = 0.6) +
  facet_grid(
    rows = vars(alpha),
    cols = vars(method),
    labeller = boxplot_labeller) +
  scale_x_discrete(
    name = "Scenario",
    breaks = every_nth(n = 4)) +
  scale_fill_manual(values = dtuPalette[c(7,9:11,5)]) +
  guides(fill = "none")
ggsave(filename = "FPR_alpha_methods.png",
       plot = FPR_alpha_methods,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

FPR_alpha_methods_median <- baseline_results %>%
  mutate(FPR = FP/(FP+TN)) %>%
  group_by(scenario, alpha, method) %>%
  reframe(medianFPR = median(FPR)) %>%
  ggplot(
    mapping = aes(
      x = scenario,
      y = medianFPR,
      colour = method,
      group = method)) +
  geom_point(shape = "-", size = 12) +
  geom_hline(
    data = nominal_values,
    mapping = aes(yintercept = alpha_num),
    linetype = "dashed",
    linewidth = 0.6) +
  facet_wrap(
    facets = vars(alpha),
    labeller = boxplot_labeller) +
  scale_x_discrete(
    name = "Scenario",
    breaks = every_nth(n = 4)) +
  scale_colour_manual(
    name = "",
    values = dtuPalette[c(7,9:11,5)]) +
  scale_y_continuous(name = "median(FPR)") +
  guides(colour = guide_legend(override.aes = list(size = 30)))
ggsave(filename = "FPR_alpha_methods_median.png",
       plot = FPR_alpha_methods_median,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

FPR_alpha_methods_mean <- baseline_results %>%
  mutate(FPR = FP/(FP+TN)) %>%
  group_by(scenario, alpha, method) %>%
  reframe(meanFPR = mean(FPR)) %>%
  ggplot(
    mapping = aes(
      x = scenario,
      y = meanFPR,
      colour = method,
      group = method)) +
  geom_point(shape = "-", size = 12) +
  geom_hline(
    data = nominal_values,
    mapping = aes(yintercept = alpha_num),
    linetype = "dashed",
    linewidth = 0.6) +
  facet_wrap(
    facets = vars(alpha),
    labeller = boxplot_labeller) +
  scale_x_discrete(
    name = "Scenario",
    breaks = every_nth(n = 4)) +
  scale_colour_manual(
    name = "",
    values = dtuPalette[c(7,9:11,5)]) +
  scale_y_continuous(name = "mean(FPR)") +
  guides(colour = guide_legend(override.aes = list(size = 30)))
ggsave(filename = "FPR_alpha_methods_mean.png",
       plot = FPR_alpha_methods_mean,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")


# POD ---------------------------------------------------------------------

POD_labeller <- as_labeller(
  c(`0.005` = "alpha[0.005]",
    `0.01` = "alpha[0.01]",
    `0.025` = "alpha[0.025]",
    `0.05` = "alpha[0.05]",
    `0.1` = "alpha[0.1]",
    `Farrington` = "Farrington",
    `Noufaily` = "Noufaily",
    `Poisson Normal` = "'Poisson Normal'",
    `Poisson Gamma` = "'Poisson Gamma'"),
  default = label_parsed
)

meadian_POD_results <- POD_results %>%
  group_by(alpha, method, k) %>%
  reframe(medianPOD = median(POD))

POD_alpha_methods <- POD_results %>%
  ggplot(mapping = aes(x = k, y = POD, colour = method)) +
  geom_line(mapping = aes(group = scenario)) +
  geom_line(data = meadian_POD_results, mapping = aes(y = medianPOD, group = method), linewidth = 2) +
  facet_grid(rows = vars(alpha), cols = vars(method), labeller = POD_labeller) +
  scale_colour_manual(values = dtuPalette[c(7, 9:11, 5)]) +
  guides(colour = "none") +
  theme(panel.spacing.y = unit(1, "lines"))
ggsave(filename = "POD_alpha_methods.png",
       plot = POD_alpha_methods,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

POD_alpha_methods_median <- meadian_POD_results %>%
  ggplot(mapping = aes(x = k, y = medianPOD, colour = method, group = method)) +
  geom_line(linewidth = 2) +
  facet_wrap(facets = vars(alpha), labeller = POD_labeller) +
  scale_colour_manual(
    name = "",
    values = dtuPalette[c(7, 9:11, 5)]) +
  scale_y_continuous(name = "median(POD)")
ggsave(filename = "POD_alpha_methods_median.png",
       plot = POD_alpha_methods_median,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# old ---------------------------------------------------------------------


# Specificity
tmp <- results %>%
  arrange(scenario, alpha) %>%
  select(scenario, sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(scenario, sim, alpha) %>%
  slice_tail(n = 49) %>%
  reframe(TP_Farrington = sum(alarm_Farrington == TRUE  & outbreakTF == TRUE),
          FN_Farrington = sum(alarm_Farrington == FALSE & outbreakTF == TRUE),
          FP_Farrington = sum(alarm_Farrington == TRUE  & outbreakTF == FALSE),
          TN_Farrington = sum(alarm_Farrington == FALSE & outbreakTF == FALSE),
          TP_Noufaily   = sum(alarm_Noufaily == TRUE    & outbreakTF == TRUE),
          FN_Noufaily   = sum(alarm_Noufaily == FALSE   & outbreakTF == TRUE),
          FP_Noufaily   = sum(alarm_Noufaily == TRUE    & outbreakTF == FALSE),
          TN_Noufaily   = sum(alarm_Noufaily == FALSE   & outbreakTF == FALSE),
          TP_PoisN      = sum(alarm_PoisN == TRUE       & outbreakTF == TRUE),
          FN_PoisN      = sum(alarm_PoisN == FALSE      & outbreakTF == TRUE),
          FP_PoisN      = sum(alarm_PoisN == TRUE       & outbreakTF == FALSE),
          TN_PoisN      = sum(alarm_PoisN == FALSE      & outbreakTF == FALSE),
          TP_PoisG      = sum(alarm_PoisG == TRUE       & outbreakTF == TRUE),
          FN_PoisG      = sum(alarm_PoisG == FALSE      & outbreakTF == TRUE),
          FP_PoisG      = sum(alarm_PoisG == TRUE       & outbreakTF == FALSE),
          TN_PoisG      = sum(alarm_PoisG == FALSE      & outbreakTF == FALSE)) %>%
  mutate(sensitivity_Farrington = TP_Farrington/(TP_Farrington+FN_Farrington),
         specificity_Farrington = TN_Farrington/(TN_Farrington+FP_Farrington),
         sensitivity_Noufaily = TP_Noufaily/(TP_Noufaily+FN_Noufaily),
         specificity_Noufaily = TN_Noufaily/(TN_Noufaily+FP_Noufaily),
         sensitivity_PoisN = TP_PoisN/(TP_PoisN+FN_PoisN),
         specificity_PoisN = TN_PoisN/(TN_PoisN+FP_PoisN),
         sensitivity_PoisG = TP_PoisG/(TP_PoisG+FN_PoisG),
         specificity_PoisG = TN_PoisG/(TN_PoisG+FP_PoisG),
         `LR+_Farrington` = sensitivity_Farrington/(1-specificity_Farrington),
         `LR+_Noufaily` = sensitivity_Noufaily/(1-specificity_Noufaily),
         `LR+_PoisN` = sensitivity_PoisN/(1-specificity_PoisN),
         `LR+_PoisG` = sensitivity_PoisG/(1-specificity_PoisG),
         `LR-_Farrington` = (1-specificity_Farrington)/sensitivity_Farrington,
         `LR-_Noufaily` = (1-specificity_Noufaily)/sensitivity_Noufaily,
         `LR-_PoisN` = (1-specificity_PoisN)/sensitivity_PoisN,
         `LR-_PoisG` = (1-specificity_PoisG)/sensitivity_PoisG,
         DOR_Farrington = `LR+_Farrington`/`LR-_Farrington`,
         DOR_Noufaily = `LR+_Noufaily`/`LR-_Noufaily`,
         DOR_PoisN = `LR+_PoisN`/`LR-_PoisN`,
         DOR_PoisG = `LR+_PoisG`/`LR-_PoisG`) %>%
  select(scenario, alpha, sim, k, `LR+_Farrington`:DOR_PoisG) %>%
  pivot_longer(cols = DOR_Farrington:DOR_PoisG) %>%
  ggplot(mapping = aes(x = factor(alpha), y = value, colour = factor(k))) +
  geom_boxplot() +
  facet_wrap(facets = vars(name))




# Calculate FPR across scenarios
FPR_scernario_alpha <- results %>%
  arrange(scenario, alpha) %>%
  select(scenario, sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  group_by(scenario, sim, alpha) %>%
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
  group_by(scenario, sim, alpha, Method) %>%
  reframe(FPR = FP/(FP+TN)) %>%
  mutate(Method = factor(Method,levels = c("Farrington",
                                           "Noufaily",
                                           "PoisN",
                                           "PoisG"),
                         labels = c("Farrington",
                                    "Noufaily",
                                    "Poisson Normal",
                                    "Poisson Gamma")))
write_rds(x = FPR_scernario_alpha, file = "FPR_scernario_alpha.Rds")

FPR_scenario_alpha_mean <- FPR_scernario_alpha %>%
  group_by(alpha, Method) %>%
  reframe(meanFPR = mean(FPR))

POD_scenario_alpha <- results %>%
  arrange(scenario, alpha) %>%
  select(scenario, sim, alpha, outbreaks) %>%
  unnest(outbreaks) %>%
  filter(t %in% 576:624) %>%
  pivot_longer(cols = alarm_Farrington:alarm_PoisG, names_to = "Method", values_to = "Alarms") %>%
  group_by(scenario, sim, alpha, Method, k) %>%
  reframe(Detected = any(outbreakTF == TRUE & Alarms == TRUE)) %>%
  group_by(scenario, alpha, Method, k) %>%
  reframe(POD = mean(Detected)) %>%
  group_by(scenario, alpha, Method, k) %>%
  mutate(Method = factor(Method,
                         levels = c("alarm_Farrington",
                                    "alarm_Noufaily",
                                    "alarm_PoisN",
                                    "alarm_PoisG"),
                         labels = c("Farrington",
                                    "Noufaily",
                                    "Poisson Normal",
                                    "Poisson Gamma")))
write_rds(x = POD_scenario_alpha, file = "POD_scenario_alpha.Rds")

POD_scenario_alpha %>%
  full_join(FPR_scenario_alpha_mean, by = join_by(alpha, Method)) %>%
  reframe(`POD/meanFPR` = POD/meanFPR) %>%
  group_by(alpha, Method, k) %>%
  reframe(meanPODFPR = mean(`POD/meanFPR`)) %>%
  filter(k %in% 2:7) %>%
  ggplot(mapping = aes(x = factor(alpha), y = meanPODFPR, colour = factor(k), group = k)) +
  geom_line() +
  facet_wrap(facets = vars(Method))



POD_scenario_alpha %>%
  full_join(FPR_scenario_alpha_mean, by = join_by(alpha, Method)) %>%
  reframe(`ODDS(POD)/ODDS(meanFPR)` = (POD/(1-POD))/(meanFPR/(1-meanFPR))) %>%
  group_by(alpha, Method, k) %>%
  reframe(meanODDS = mean(`ODDS(POD)/ODDS(meanFPR)`)) %>%
  filter(k %in% 2:7) %>%
  ggplot(mapping = aes(x = factor(alpha), y = meanODDS, colour = factor(k), group = k)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette)

# Making beautiful plot
profilePODxFPR_num <- POD_scenario_alpha %>%
  full_join(FPR_scenario_alpha_mean, by = join_by(alpha, Method)) %>%
  filter(k == 5) %>%
  ungroup(scenario) %>%
  reframe(POD = median(POD), meanFPR = median(meanFPR)) %>%
  ggplot(mapping = aes(x = meanFPR, y = POD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_label(aes(label = alpha)) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)])

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

FPR_mean <- FPR %>%
  group_by(alpha, Method) %>%
  reframe(meanFPR = mean(FPR))

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
  mutate(Method = factor(Method,
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
  ggplot(mapping = aes(x = meanFPR, y = POD), nudge_x = 1, nudge_y = 1, size = 0.1) +
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
  ggplot(mapping = aes(x = meanFPR, y = POD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_point(mapping = aes(x = meanFPR, y = POD, shape = factor(alpha)), size = 3) +
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
  ggplot(mapping = aes(x = meanFPR, y = POD), nudge_x = 1, nudge_y = 1, size = 0.1) +
  geom_line(mapping = aes(colour = Method), linewidth = 1) +
  geom_point(mapping = aes(x = meanFPR, y = POD, shape = factor(alpha)), size = 3) +
  facet_wrap(facets = vars(k), ncol = 2) +
  scale_color_manual(values = dtuPalette[c(7,9:11,5)]) +
  scale_shape_discrete(name = "alpha") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        axis.text = element_text(size = 16), strip.text = element_text(size = 16))
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
  geom_line(linewidth = 2) +
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

# POD/FPR with k = 5
POD_FPR_k5 <- POD %>%
  filter(k == 5) %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  reframe(`POD/meanFPR` = POD/meanFPR) %>%
  ggplot(mapping = aes(x = factor(alpha), y = `POD/meanFPR`, colour = factor(Method), group = Method)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette[c(7,8:10)]) +
  guides(colour = "none")
ggsave(filename = "POD_FPR_k5.png",
       plot = POD_FPR_k5,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

ODDS_POD_FPR_k5 <- POD %>%
  filter(k == 5) %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  reframe(`ODDS(POD)/ODDS(meanFPR)` = (POD/(1-POD))/(meanFPR/(1-meanFPR))) %>%
  ggplot(mapping = aes(x = factor(alpha), y = `ODDS(POD)/ODDS(meanFPR)`, colour = factor(Method), group = Method)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette[c(7,8:10)]) +
  guides(colour = "none")
ggsave(filename = "ODDS_POD_FPR_k5.png",
       plot = ODDS_POD_FPR_k5,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")

# Subset of k
POD_FPR <- POD %>%
  filter(k %in% c(2:6)) %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  reframe(`POD/meanFPR` = POD/meanFPR) %>%
  ggplot(mapping = aes(x = factor(alpha), y = `POD/meanFPR`, colour = factor(k), group = k)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette)
ggsave(filename = "POD_FPR.png",
       plot = POD_FPR,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")




ODDS_POD_FPR <- POD %>%
  filter(k %in% c(2:6)) %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  reframe(`ODDS(POD)/ODDS(meanFPR)` = (POD/(1-POD))/(meanFPR/(1-meanFPR))) %>%
  ggplot(mapping = aes(x = factor(alpha), y = `ODDS(POD)/ODDS(meanFPR)`, colour = factor(k), group = k)) +
  geom_line(linewidth = 1) +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette)
ggsave(filename = "ODDS_POD_FPR.png",
       plot = ODDS_POD_FPR,
       path = "../../figures/",
       device = png,
       width = 16,
       height = 8,
       units = "in",
       dpi = "print")


POD %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  reframe(p1 = POD/(1-POD),
          p2 = meanFPR/(1-meanFPR),
          `Odds ratio` = p1/p2) %>%
  ggplot(mapping = aes(x = alpha, y = `Odds ratio`, colour = factor(k), group = factor(k))) +
  geom_line() +
  facet_wrap(facets = vars(Method)) +
  scale_color_manual(name = "k", values = dtuPalette)

binary_logistic_regression_POD <- function(df) {
  glm(POD ~ alpha + k, family = binomial(link = "logit"), data = df)
}
binary_logistic_regression_FPR <- function(df) {
  glm(meanFPR ~ alpha + k, family = binomial(link = "logit"), data = df)
}

binary_logistic_regression_tbl <- POD %>%
  full_join(FPR_mean, by = join_by(alpha, Method)) %>%
  group_by(Method) %>%
  nest() %>%
  mutate(fit_POD = map(data, binary_logistic_regression_POD),
         fit_FPR = map(data, binary_logistic_regression_FPR),
         tidied_POD = map(fit_POD, tidy),
         tidied_FPR = map(fit_FPR, tidy))

binary_logistic_regression_tbl %>%
  unnest(tidied_POD)
binary_logistic_regression_tbl %>%
  unnest(tidied_FPR)

