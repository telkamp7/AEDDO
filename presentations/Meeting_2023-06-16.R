

library(readr)
library(dplyr)


LIST_PoisN_ageGroup <- read_rds(file = "C:/GIT/AEDDO/src/cases/LIST_PoisN_ageGroup.rds") %>%
  mutate(method = "PoisN")
LIST_PoisG_ageGroup <- read_rds(file = "C:/GIT/AEDDO/src/cases/LIST_PoisG_ageGroup.rds")%>%
  mutate(method = "PoisG")

STEC_PoisN_ageGroup_trend_seasonality <- read_rds(file = "C:/GIT/AEDDO/src/cases/STEC_PoisN_ageGroup_trend_seasonality.rds")%>%
  mutate(method = "PoisN")
STEC_PoisG_ageGroup_trend_seasonality <- read_rds(file = "C:/GIT/AEDDO/src/cases/STEC_PoisG_ageGroup_trend_seasonality.rds") %>%
  mutate(method = "PoisG")


LIST <- bind_rows(LIST_PoisN_ageGroup, LIST_PoisG_ageGroup)
STEC <- bind_rows(STEC_PoisN_ageGroup_trend_seasonality, STEC_PoisG_ageGroup_trend_seasonality)


LIST %>%
  select(ran.ef, method) %>%
  unnest(ran.ef) %>%
  select(Date = ref.date, ageGroup, y, method, u, p, alarm) %>%
  filter(alarm)



STEC %>%
  select(ran.ef, method) %>%
  unnest(ran.ef) %>%
  select(Date = ref.date, ageGroup, y, method, u, p, alarm) %>%
  filter(alarm)
