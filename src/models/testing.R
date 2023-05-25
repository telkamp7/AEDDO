
# Import libraries
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(lubridate)
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
                
custom_labeller <- as_labeller(
  c(`ageGroup<1 year`="beta['<1'~year]", `ageGroup1-4 years`="beta[1-4~years]",
    `ageGroup5-14 years`="beta[5-14~years]",`ageGroup15-24 years`="beta[15-24~years]",
    `ageGroup25-64 years`="beta[25-64~years]", `ageGroup65+ years`="beta[65+~years]",
    `M_ageGroup`="M[ageGroup]", `M_ageGroupSeasonality`="M[ageGroupSeasonality]",
    `sin(pi/6 * monthInYear)`="beta[sin]", `cos(pi/6 * monthInYear)`="beta[cos]",
    `log_phi`="log(phi)", `log_sigma`="log(sigma)", `<1 year`="'<1'~year", `1-4 years`="1-4~years",
    `5-14 years`="5-14~years", `15-24 years`="15-24~years", `25-64 years`="25-64~years",
    `65+ years`="65+~years"),
  default = label_parsed
)
                
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
dat <- read_rds(file = "../../data/processed/dat4.rds") # 6 agegroups

ggplot(data = dat, mapping = aes(x = Date, y = cases/n *1e5)) +
  geom_line() +
  facet_grid(rows = vars(caseDef), cols = vars(ageGroup), scales = "free_y")

# Nest the data
dat_nest <- dat %>%
  group_by(caseDef, Date, ageGroup) %>%
  reframe(y = sum(cases), n = sum(n)) %>%
  mutate(monthInYear = as.integer(format(Date, "%m"))) %>%
  group_by(caseDef) %>%
  nest()

results <- dat_nest %>% 
  filter(caseDef %in% c("STEC")) %>%
  mutate(
    M_ageGroup_PoisG = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup,
        theta = rep(1,7),
        method = "L-BFGS-B",
        lower = c(rep(1e-6,6), -6),
        upper = rep(1e2, 7),
        model = "PoissonGamma",
        k = 36,
        excludePastOutbreaks = TRUE
      )
    }),
    M_ageGroupSeasonality_PoisG = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
        theta = rep(1,9),
        method = "L-BFGS-B", 
        lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
        upper = rep(1e2, 9),
        model = "PoissonGamma",
        k = 36,
        excludePastOutbreaks = TRUE
      )
    }),
    M_ageGroup_PoisN = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup,
        theta = rep(1,7),
        method = "L-BFGS-B",
        lower = c(rep(10,6), -6),
        upper = c(rep(20,6),3),
        model = "PoissonNormal",
        k = 36,
        excludePastOutbreaks = TRUE  
      )
    }),
    M_ageGroupSeasonality_PoisN = map(data, function(df){
      aeddo(
        data = df,
        formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
        theta = rep(1,9),
        method = "L-BFGS-B",
        lower = c(rep(10,6),rep(-10,2),-6),
        upper = c(rep(20, 6),rep(10,2),3),
        model = "PoissonNormal",
        k = 36,
        excludePastOutbreaks = TRUE
      )
    }
    )
  )
# 
# results1 <- dat_nest %>% 
#   filter(caseDef %in% c("STEC")) %>%
#   mutate(
#     M_ageGroupSeasonality_PoisN = map(data, function(df){
#       aeddo(
#         data = df,
#         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
#         theta = rep(1,9),
#         method = "L-BFGS-B",
#         lower = c(rep(10,6),rep(-10,2),-6),
#         upper = c(rep(20, 6),rep(10,2),10),
#         model = "PoissonNormal",
#         k = 36,
#         excludePastOutbreaks = TRUE
#         )
#       }
#       ) 
#   )

# results_PoisN <- dat_nest %>% 
#   filter(caseDef %in% c("STEC")) %>%
#   mutate(
#     M_ageGroup_PoisN = map(data, function(df){
#       aeddo(
#         data = df,
#         formula = y ~ -1 + ageGroup,
#         theta = rep(1,7),
#         method = "BFGS",
#         model = "PoissonNormal",
#         k = 36,
#         excludePastOutbreaks = TRUE  
#         )
#       }),
#     M_ageGroupSeasonality_PoisN = map(data, function(df){
#       aeddo(
#         data = df,
#         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
#         theta = rep(1,9),
#         method = "L-BFGS-B", 
#         model = "PoissonNormal",
#         k = 36,
#         excludePastOutbreaks = TRUE
#         )
#       })
#   )
# 
# STEC <- dat %>%
#   filter(caseDef == "STEC") %>%
#   rename(y = cases) %>%
#   mutate(monthInYear = as.integer(format(Date, "%m")))

# 
# aeddo(data = STEC,
#       formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
#       theta = rep(1,9),
#       method = "L-BFGS-B", 
#       lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
#       model = "PoissonNormal",
#       k = 36,
#       excludePastOutbreaks = TRUE)

lab1 <- c(expression(M[ageGroup]),
          expression(M[ageGroupSeasonality]))

AICxSTEC_PoisG <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisG, M_ageGroupSeasonality_PoisG) %>%
  unnest(c(M_ageGroup_PoisG, M_ageGroupSeasonality_PoisG), names_sep = "_") %>%
  select(ref.date = M_ageGroup_PoisG_ref.date, M_ageGroup_PoisG_value, M_ageGroupSeasonality_PoisG_value,
         M_ageGroup_PoisG_par, M_ageGroupSeasonality_PoisG_par, M_ageGroup_PoisG_window.data, 
         M_ageGroupSeasonality_PoisG_window.data) %>%
  rename(M_ageGroup_PoisG = M_ageGroup_PoisG_value, M_ageGroupSeasonality_PoisG = M_ageGroupSeasonality_PoisG_value) %>%
  pivot_longer(cols = M_ageGroup_PoisG:M_ageGroupSeasonality_PoisG, names_to = "M_value") %>%
  rename(M_ageGroup_PoisG = M_ageGroup_PoisG_par, M_ageGroupSeasonality_PoisG = M_ageGroupSeasonality_PoisG_par) %>%
  pivot_longer(cols = M_ageGroup_PoisG:M_ageGroupSeasonality_PoisG, names_to = "M_par", values_to = "par") %>%
  rename(M_ageGroup_PoisG = M_ageGroup_PoisG_window.data, M_ageGroupSeasonality_PoisG = M_ageGroupSeasonality_PoisG_window.data) %>%
  pivot_longer(cols = M_ageGroup_PoisG:M_ageGroupSeasonality_PoisG, names_to = "M_data", values_to = "data") %>%
  filter(M_value == M_par & M_value == M_data) %>%
  select(M = M_value, ref.date, value, par, data) %>%
  mutate(nPar = map_int(par, nrow), nData = map_int(data, nrow), AIC = 2*nPar + 2*value, M = gsub(pattern = "_PoisG", replacement = "", x = M)) %>%
  ggplot(mapping = aes(x = ref.date, colour = M, group = M)) +
  geom_line(mapping = aes(y = value)) +
  geom_line(mapping = aes(y = nData*2), lty = "dashed") +
  scale_y_continuous(name = "AIC", sec.axis = sec_axis(~ ./2, name = "Number of observations")) +
  scale_x_date(name = "Date") +
  scale_colour_manual(name = "Model", values = dtuPalette, labels = lab1)
ggsave(filename = "AICxSTEC_PoisG.png", plot = AICxSTEC_PoisG, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

M_ageGroup_PoisG_par_tbl <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisG) %>%
  unnest(c(M_ageGroup_PoisG), names_sep = "_") %>%
  select(ref.date = M_ageGroup_PoisG_ref.date, M_ageGroup_PoisG_par) %>%
  unnest(c(M_ageGroup_PoisG_par), names_sep = "_") %>%
  rename(M_par = M_ageGroup_PoisG_par_Parameter, M_theta = M_ageGroup_PoisG_par_theta, M_se.theta = M_ageGroup_PoisG_par_se.theta) %>%
  mutate(M = "M_ageGroup")

M_ageGroupSeasonality_PoisG_par_tbl <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroupSeasonality_PoisG) %>%
  unnest(c(M_ageGroupSeasonality_PoisG), names_sep = "_") %>%
  select(ref.date = M_ageGroupSeasonality_PoisG_ref.date, M_ageGroupSeasonality_PoisG_par) %>%
  unnest(c(M_ageGroupSeasonality_PoisG_par), names_sep = "_") %>%
  rename(M_par = M_ageGroupSeasonality_PoisG_par_Parameter, M_theta = M_ageGroupSeasonality_PoisG_par_theta, M_se.theta = M_ageGroupSeasonality_PoisG_par_se.theta) %>%
  mutate(M = "M_ageGroupSeasonality")

M_PoisG_par_tbl <- bind_rows(M_ageGroup_PoisG_par_tbl, M_ageGroupSeasonality_PoisG_par_tbl) %>%
  mutate(M = factor(M, levels = c("M_ageGroup", "M_ageGroupSeasonality")),
         M_par = factor(M_par, levels = c(
           paste0("ageGroup",
                  c("<1 year",
                    "1-4 years",
                    "5-14 years",
                    "15-24 years",
                    "25-64 years",
                    "65+ years")),
           "log_phi",
           "sin(pi/6 * monthInYear)",
           "cos(pi/6 * monthInYear)")))

ageGroupParxSTEC_PoisG <- M_PoisG_par_tbl %>% 
  filter(grepl(pattern = "ageGroup", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_grid(rows = vars(M), cols = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "ageGroupParxSTEC_PoisG.png", plot = ageGroupParxSTEC_PoisG, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

SeasonalityParxSTEC_PoisG <- M_ageGroupSeasonality_PoisG_par_tbl %>%
  filter(grepl(pattern = "sin|cos", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_wrap(facets = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "SeasonalityParxSTEC_PoisG.png", plot = SeasonalityParxSTEC_PoisG, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

log_phixSTEC_PoisG <- M_PoisG_par_tbl %>% 
  filter(grepl(pattern = "log_phi", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_grid(rows = vars(M), cols = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "log_phixSTEC_PoisG.png", plot = log_phixSTEC_PoisG, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

M_ageGroup_OutbreakDetection_PoisG <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisG) %>%
  unnest(M_ageGroup_PoisG) %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(M = "M_ageGroup")
  
M_ageGroupSeasonality_PoisG_OutbreakDetection <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroupSeasonality_PoisG) %>%
  unnest(M_ageGroupSeasonality_PoisG) %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(M = "M_ageGroupSeasonality")

M_OutbreakDetection_PoisG <- bind_rows(M_ageGroup_OutbreakDetection_PoisG, M_ageGroupSeasonality_PoisG_OutbreakDetection)

OutbreakDetectionxSTEC_PoisG <- M_OutbreakDetection_PoisG %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_grid(rows = vars(M), cols = vars(ageGroup), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  scale_x_date(name = "Date") +
  guides(colour = "none", shape = "none")
ggsave(filename = "OutbreakDetectionxSTEC_PoisG.png", plot = OutbreakDetectionxSTEC_PoisG, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print") 



AICxSTEC_PoisN <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisN, M_ageGroupSeasonality_PoisN) %>%
  unnest(c(M_ageGroup_PoisN, M_ageGroupSeasonality_PoisN), names_sep = "_") %>%
  select(ref.date = M_ageGroup_PoisN_ref.date, M_ageGroup_PoisN_value, M_ageGroupSeasonality_PoisN_value,
         M_ageGroup_PoisN_par, M_ageGroupSeasonality_PoisN_par, M_ageGroup_PoisN_window.data, 
         M_ageGroupSeasonality_PoisN_window.data) %>%
  rename(M_ageGroup_PoisN = M_ageGroup_PoisN_value, M_ageGroupSeasonality_PoisN = M_ageGroupSeasonality_PoisN_value) %>%
  pivot_longer(cols = M_ageGroup_PoisN:M_ageGroupSeasonality_PoisN, names_to = "M_value") %>%
  rename(M_ageGroup_PoisN = M_ageGroup_PoisN_par, M_ageGroupSeasonality_PoisN = M_ageGroupSeasonality_PoisN_par) %>%
  pivot_longer(cols = M_ageGroup_PoisN:M_ageGroupSeasonality_PoisN, names_to = "M_par", values_to = "par") %>%
  rename(M_ageGroup_PoisN = M_ageGroup_PoisN_window.data, M_ageGroupSeasonality_PoisN = M_ageGroupSeasonality_PoisN_window.data) %>%
  pivot_longer(cols = M_ageGroup_PoisN:M_ageGroupSeasonality_PoisN, names_to = "M_data", values_to = "data") %>%
  filter(M_value == M_par & M_value == M_data) %>%
  select(M = M_value, ref.date, value, par, data) %>%
  mutate(nPar = map_int(par, nrow), nData = map_int(data, nrow), AIC = 2*nPar + 2*value, M = gsub(pattern = "_PoisN", replacement = "", x = M)) %>%
  ggplot(mapping = aes(x = ref.date, colour = M, group = M)) +
  geom_line(mapping = aes(y = value)) +
  geom_line(mapping = aes(y = nData*2), lty = "dashed") +
  scale_y_continuous(name = "AIC", sec.axis = sec_axis(~ ./2, name = "Number of observations")) +
  scale_x_date(name = "Date") +
  scale_colour_manual(name = "Model", values = dtuPalette, labels = lab1)
ggsave(filename = "AICxSTEC_PoisN.png", plot = AICxSTEC_PoisN, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

M_ageGroup_PoisN_par_tbl <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisN) %>%
  unnest(c(M_ageGroup_PoisN), names_sep = "_") %>%
  select(ref.date = M_ageGroup_PoisN_ref.date, M_ageGroup_PoisN_par) %>%
  unnest(c(M_ageGroup_PoisN_par), names_sep = "_") %>%
  rename(M_par = M_ageGroup_PoisN_par_Parameter, M_theta = M_ageGroup_PoisN_par_theta, M_se.theta = M_ageGroup_PoisN_par_se.theta) %>%
  mutate(M = "M_ageGroup")

M_ageGroupSeasonality_PoisN_par_tbl <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroupSeasonality_PoisN) %>%
  unnest(c(M_ageGroupSeasonality_PoisN), names_sep = "_") %>%
  select(ref.date = M_ageGroupSeasonality_PoisN_ref.date, M_ageGroupSeasonality_PoisN_par) %>%
  unnest(c(M_ageGroupSeasonality_PoisN_par), names_sep = "_") %>%
  rename(M_par = M_ageGroupSeasonality_PoisN_par_Parameter, M_theta = M_ageGroupSeasonality_PoisN_par_theta, M_se.theta = M_ageGroupSeasonality_PoisN_par_se.theta) %>%
  mutate(M = "M_ageGroupSeasonality")

M_PoisN_par_tbl <- bind_rows(M_ageGroup_PoisN_par_tbl, M_ageGroupSeasonality_PoisN_par_tbl) %>%
  mutate(M = factor(M, levels = c("M_ageGroup", "M_ageGroupSeasonality")),
         M_par = factor(M_par, levels = c(
           paste0("ageGroup",
                  c("<1 year",
                    "1-4 years",
                    "5-14 years",
                    "15-24 years",
                    "25-64 years",
                    "65+ years")),
           "log_sigma",
           "sin(pi/6 * monthInYear)",
           "cos(pi/6 * monthInYear)")))


ageGroupParxSTEC_PoisN <- M_PoisN_par_tbl %>% 
  filter(grepl(pattern = "ageGroup", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_grid(rows = vars(M), cols = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "ageGroupParxSTEC_PoisN.png", plot = ageGroupParxSTEC_PoisN, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

SeasonalityParxSTEC_PoisN <- M_ageGroupSeasonality_PoisN_par_tbl %>%
  filter(grepl(pattern = "sin|cos", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_wrap(facets = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "SeasonalityParxSTEC_PoisN.png", plot = SeasonalityParxSTEC_PoisN, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")

log_sigmaxSTEC_PoisN <- M_PoisN_par_tbl %>% 
  filter(grepl(pattern = "log_sigma", x = M_par)) %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_grid(rows = vars(M), cols = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")
ggsave(filename = "log_sigmaxSTEC_PoisN.png", plot = log_sigmaxSTEC_PoisN, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print")


log_sigmaxSTEC_PoisN <- M_PoisN_par_tbl %>% 
  filter(grepl(pattern = "log_sigma", x = M_par) & M == "M_ageGroupSeasonality") %>%
  ggplot(mapping = aes(x = ref.date)) +
  geom_line(mapping = aes(y = M_theta, colour = M_par)) +
  geom_line(mapping = aes(y = M_theta - 2*M_se.theta), linetype = "dashed") +
  geom_line(mapping = aes(y = M_theta + 2*M_se.theta), linetype = "dashed") +
  facet_grid(rows = vars(M), cols = vars(M_par), labeller = custom_labeller) +
  scale_y_continuous(name = "Parameter value") +
  scale_x_date(name = "Date") +
  coord_cartesian(ylim = c(-6, 3)) +
  scale_colour_manual(values = dtuPalette) +
  guides(colour = "none")

M_ageGroup_OutbreakDetection_PoisN <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroup_PoisN) %>%
  unnest(M_ageGroup_PoisN) %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(M = "M_ageGroup")

M_ageGroupSeasonality_PoisN_OutbreakDetection <- results %>%
  ungroup() %>%
  filter(caseDef == "STEC") %>%
  select(M_ageGroupSeasonality_PoisN) %>%
  unnest(M_ageGroupSeasonality_PoisN) %>% 
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  mutate(M = "M_ageGroupSeasonality")

M_OutbreakDetection_PoisN <- bind_rows(M_ageGroup_OutbreakDetection_PoisN, M_ageGroupSeasonality_PoisN_OutbreakDetection)

OutbreakDetectionxSTEC_PoisN <- M_OutbreakDetection_PoisN %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qnorm(p = 0.95, mean = 0, sd = exp(log_sigma))),
            lty = "dashed", inherit.aes = FALSE) +
  facet_grid(rows = vars(M), cols = vars(ageGroup), labeller = custom_labeller) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  scale_x_date(name = "Date") +
  guides(colour = "none", shape = "none")
ggsave(filename = "OutbreakDetectionxSTEC_PoisN.png", plot = OutbreakDetectionxSTEC_PoisN, path = "../../figures/", device = png, width = 16, height = 8, units = "in", dpi = "print") 
















 

















LIST_dat <- dat %>%
  filter(caseDef == "LIST") %>%
  group_by(Date) %>%
  reframe(y = sum(cases), n = sum(n), ageGroup = "All") %>%
  mutate(monthInYear = as.integer(format(Date, "%m")))

LIST_dat %>%
  ggplot(mapping = aes(x = Date, y = y)) +
  geom_line() +
  geom_point()

results_LIST <- aeddo(data = LIST_dat,
                      formula = y ~ 1,
                      theta = c(1,1),
                      method = "L-BFGS-B",
                      lower = c(1e-6,-6),
                      model = "PoissonGamma",
                      k = 24,
                      excludePastOutbreaks = TRUE)








# PoissonGamma %>%
#   ungroup() %>%
#   select(M1, M2) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date = M1_ref.date, M1 = M1_value, M2 = M2_value) %>%
#   pivot_longer(cols = M1:M2, names_to = "Model") %>%
#   ggplot(mapping = aes(x = ref.date, y = value, colour = Model)) +
#   geom_line() +
#   scale_x_date(name = "Month") +
#   scale_y_continuous(name = "Objective") +
#   scale_colour_manual(values = dtuPalette)

# formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
#         theta = rep(1,9),
#         method = "L-BFGS-B", 
#         lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
#         model = "PoissonGamma",

results %>%
  filter(caseDef == "STEC") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  filter(alarm) %>%
  distinct(ref.date)

results %>%
  filter(caseDef == "STEC") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")


LIST_dat %>%
  ggplot(mapping = aes(x = Date, y = y)) +
  geom_line() +
  geom_point() +
  annotate(geom = "rect", xmin = as.Date("2008-01-01")-10, xmax = ymd(as.Date("2008-01-01")) %m+% months(24), ymin = -Inf, ymax = Inf, alpha = 0.2)

results_LIST %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = y, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")

results_LIST %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")

results_LIST %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  filter(alarm) %>%
  distinct(ref.date)

results_LIST %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")

results_LIST %>%
  unnest(par) %>%
  filter(Parameter == "(Intercept)") %>%
  ggplot(mapping = aes(x = ref.date, y = theta)) +
  geom_line()


results %>%
  filter(caseDef == "SHIG") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = y, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line() +
  geom_point() +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")

results %>%
  filter(caseDef == "SHIG") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")

results %>%
  filter(caseDef == "SHIG") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  filter(alarm) %>%
  distinct(ref.date)

results %>%
  filter(caseDef == "SHIG") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")


results %>%
  filter(caseDef == "SALM") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  filter(alarm) %>%
  distinct(ref.date)

results %>%
  filter(caseDef == "SALM") %>%
  unnest(M_ageGroup) %>%
  select(ran.ef) %>%
  unnest(ran.ef) %>%
  ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = alarm)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line(mapping = aes(x = ref.date,
                          y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
            lty = "dashed", inherit.aes = FALSE) +
  facet_wrap(facets = vars(ageGroup)) +
  scale_color_manual(values = dtuPalette) +
  scale_shape_manual(values = c(1,19)) +
  guides(colour = "none", shape = "none")



# PoissonGamma <- dat_nest %>%
#   filter(caseDef %in% c("Shiga- og veratoxin producerende E. coli.")) %>%
#   mutate(
#     M1 = map(data, function(df){
#       aeddo(
#         data = df,
#         formula = y ~ -1 + ageGroup,
#         theta = rep(1,7),
#         method = "L-BFGS-B", 
#         lower = c(rep(1e-6,6),-6), 
#         model = "PoissonGamma",
#         k = 36, 
#         excludePastOutbreaks = TRUE)
#     }
#     ),
#     M2 = map(data, function(df){
#       aeddo(
#         data = df,
#         formula = y ~ -1 + ageGroup + sin(pi/6*monthInYear) + cos(pi/6*monthInYear),
#         theta = rep(1,9),
#         method = "L-BFGS-B", 
#         lower = c(rep(1e-6,6),rep(-Inf,2),-6), 
#         model = "PoissonGamma",
#         k = 36, 
#         excludePastOutbreaks = TRUE)
#     }
#     )
#   )
# 
# PoissonGamma %>%
#   ungroup() %>%
#   select(M1, M2) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date = M1_ref.date, M1 = M1_value, M2 = M2_value) %>%
#   pivot_longer(cols = M1:M2, names_to = "Model") %>%
#   ggplot(mapping = aes(x = ref.date, y = value, colour = Model)) +
#   geom_line() +
#   scale_x_date(name = "Month") +
#   scale_y_continuous(name = "Objective") +
#   scale_colour_manual(values = dtuPalette)
#   
# 
# PoissonGamma %>%
#   ungroup() %>%
#   select(M1, M2) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date = M1_ref.date, M1 = M1_ran.ef, M2 = M2_ran.ef) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date, ageGroup = M1_ageGroup, n = M1_n, y = M1_y, monthInYear = M1_monthInYear, M1_phi:M1_alarm, M2_phi:M2_alarm) %>%
#   rename(M1 = M1_u, M2 = M2_u) %>%
#   pivot_longer(cols = c(M1, M2), names_to = "Model_u", values_to = "u") %>%
#   rename(M1 = M1_phi, M2 = M2_phi) %>%
#   pivot_longer(cols = c(M1, M2), names_to = "Model_phi", values_to = "phi") %>%
#   rename(M1 = M1_alarm, M2 = M2_alarm) %>%
#   pivot_longer(cols = c(M1, M2), names_to = "Model_alarm", values_to = "Alarm") %>%
#   filter(Model_u == Model_phi & Model_u == Model_alarm) %>%
#   rename(Model = Model_u) %>%
#   select(ref.date:monthInYear, Model, u, phi, Alarm) %>%
#   ggplot(mapping = aes(x = ref.date, y = u, colour = ageGroup, group = ageGroup, shape = Alarm)) +
#   geom_line(linewidth = 0.4) +
#   geom_point(size = 2) +
#   geom_line(mapping = aes(x = ref.date,
#                           y = qgamma(p = 0.95, shape = 1/phi, scale = phi)),
#             lty = "dashed", inherit.aes = FALSE) +
#   facet_grid(rows = vars(Model), cols = vars(ageGroup)) +
#   scale_x_date(name = "Month") +
#   scale_colour_manual(values = dtuPalette) +
#   scale_shape_manual(values = c(1,19)) +
#   guides(colour = "none",shape = "none")
# 
# 
# 
# PoissonGamma %>%
#   ungroup() %>%
#   select(M1, M2) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date = M1_ref.date, M1 = M1_ran.ef, M2 = M2_ran.ef) %>%
#   unnest(c(M1, M2), names_sep = "_") %>%
#   select(ref.date, ageGroup = M1_ageGroup, M1_alarm, M2_alarm) %>%
#   rename(M1 = M1_alarm, M2 = M2_alarm) %>%
#   pivot_longer(cols = c(M1, M2), names_to = "Model", values_to = "Alarm") %>%
#   mutate(Model = factor(Model, levels = c("M1", "M2")),
#          alarmDate = if_else(Alarm, ref.date, as.Date(NA_character_))) %>%
#   ggplot(mapping = aes(x = alarmDate, y = Model, colour = Model)) +
#   geom_point(shape = 17) +
#   facet_wrap(facets = vars(ageGroup)) +
#   scale_color_manual(values = dtuPalette) +
#   scale_y_discrete(limits = rev(c("M1", "M2"))) +
#   guides(colour = "none")
