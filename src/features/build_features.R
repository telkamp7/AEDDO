
# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(lubridate)
# library(openxlsx)
library(readxl)

locale(date_names = "da")
# Set locale
Sys.setlocale(category = "LC_ALL", locale = "da_DK.utf8")

# Import the data
FOLK1A <- read_csv2("../../data/raw/FOLK1A.csv")
# diseaseData <- read.xlsx(xlsxFile = "../../data/raw/disease_data_raw.xlsx")
diseaseData <- read_xlsx(path = "../../data/raw/disease_data_raw.xlsx")
NUTS <- read_csv2("../../data/raw/NUTS_V1_2007.csv")

# Check for 'Uoplyst'
diseaseData %>%
  group_by(ageGroup) %>%
  summarize(sum(cases)) %>%print(n=144)

# tmp <- diseaseData %>%
#   group_by(caseDef, year) %>%
#   summarize(combi = n_distinct(ageGroup, landsdel), nAge = n_distinct(ageGroup), nLandsdel = n_distinct(landsdel))

# unique(diseaseData$ageGroup)
# unique(diseaseData$landsdel)

# Change '?r' to 'years'
diseaseData <- diseaseData %>%
  mutate(ageGroup = if_else(
    ageGroup == "<1 år",
    str_replace(ageGroup, pattern = "år", replacement = "year"),
    str_replace(ageGroup, pattern = "år", replacement = "years")
    ))

# Make reference table for 'landsdele' and 'kommuner'
# NutsCor <- NUTS %>%
#   select(KODE, NIVEAU, TITEL) %>%
#   pivot_wider(id_cols = KODE, names_from = NIVEAU, values_from = TITEL) %>%
#   rename(Region = `1`, Landsdel = `2`, Kommune = `3`) %>%
#   fill(Region, Landsdel, .direction = "down") %>%
#   mutate(Kommune = replace(Kommune, Kommune=="K?benhavn", "Copenhagen")) %>%
#   mutate(Landsdel = str_replace(Landsdel, "Landsdel ", "")) %>%
#   mutate(Landsdel = str_replace(Landsdel, "Byen København", "København by")) %>%
#   drop_na()
# ... Note: 'Byen København' is turned into 'København by', to match diseaseData

# # Join FOLK1A and with reference table, to obtain 'landsdele'
# FOLK1AxNutsCor <- inner_join(x = FOLK1A,
#                              y = NutsCor,
#                              by = c("OMRÅDE" = "Kommune"))
# 
# # # Calculate FOLK1A per 'landsdel'
# FOLK1AxLandsdel <- FOLK1AxNutsCor %>%
#   group_by(ALDER, TID, Landsdel) %>%
#   reframe(n = sum(INDHOLD))

FOLK1AxAgegroup <- FOLK1A %>%
  group_by(ALDER, TID) %>%
  reframe(n = sum(INDHOLD))

# FOLK1AxAgegroup <- FOLK1A %>%
#   group_by(ALDER, TID) %>%
#   reframe(n = sum(n))

# 
# unique(FOLK1AxLandsdel$ALDER)

# FOLK1AxLandsdel <- FOLK1AxNutsCor %>%
#   group_by(ALDER, TID, Landsdel) %>%
#   summarise(n = sum(INDHOLD))
# 
# FOLK1AxLandsdel %>% 
#   summarize(combi = n_distinct(ALDER, Landsdel), nAge = n_distinct(ALDER), nLandsdel = n_distinct(Landsdel))
# 
# # Calculate FOLK1A per 'Region'
# FOLK1AxRegion <- FOLK1AxNutsCor %>%
#   group_by(ALDER, TID, Region) %>%
#   reframe(n = sum(INDHOLD))

# Finalize the data set
# datTmp <- diseaseData %>%
#   filter(month != "Uoplyst") %>%
#   mutate(Date = as.Date(
#     paste0(year,"-",month,"-1"), format = "%Y-%B-%e")
#     ) %>%
#   mutate(QUARTER = quarter(Date)) %>%
#   mutate(TID = paste0(year,"Q",QUARTER)) %>%
#   inner_join(y = FOLK1AxLandsdel, by = c("ageGroup" = "ALDER",
#                                          "landsdel" = "Landsdel",
#                                          "TID" = "TID")) %>%
#   select(Date, ageGroup, landsdel, caseDef, cases, n)

datTmp <- diseaseData %>%
  filter(month != "Uoplyst") %>%
  mutate(Date = as.Date(
    paste0(year,"-",month,"-1"), format = "%Y-%B-%e")
  ) %>%
  mutate(QUARTER = quarter(Date)) %>%
  mutate(TID = paste0(year,"Q",QUARTER)) %>%
  inner_join(y = FOLK1AxAgegroup, by = c("ageGroup" = "ALDER",
                                         "TID" = "TID")) %>%
  select(Date, ageGroup, caseDef, cases, n)
# ... Note: maaned and year is converted to quarters, in order to add 'n', which
# is the number of individuals from FOLK1A in each 'landsdel' and age group
# Hereafter, the data is combined in 'datTmp', which contains the final features
# ... Note: Some combination of 'ageGroup' and 'landsdele' are missing in the 
# data.

# Extract levels
# maanedLevels <- unique(datTmp$maaned)
# yearLevels <- unique(datTmp$year)
ageGroupLevels <- unique(datTmp$ageGroup)
caseDefLevels <- unique(datTmp$caseDef)
# caseDefLabels <- c(
#   "AIDS",
#   "Botulisme",
#   "Gonoré",
#   "Hepatitis A",
#   "Hepatitis B",
#   "Hepatitis C",
#   "Hæmophilus influenza meningitis",
#   "HIV infektion",
#   "Legionella",
#   "Leptospirosis",
#   "Mæslinger",
#   "Andre meningitis",
#   "Meningokoksygdom",
#   "MPOX",
#   "Fåresyge",
#   "Neuroborreliose",
#   "Ornitose",
#   "Kighoste",
#   "Pneumokok meningitis",
#   "Røde hunde",
#   "Shigella",
#   "Syfilis",
#   "Stivkrampe",
#   "Tubercolosis",
#   "Tyfus / paratyfus",
#   "Shiga- og veratoxin producerende E. coli."
# )

caseDefLabels <- c("STEC",
                   "SHIG",
                   "LIST",
                   "SALM")

# Finalize 'dat'
dat <- datTmp %>%
  mutate(ageGroup = factor(x = ageGroup, levels = ageGroupLevels),
         caseDef = factor(x = caseDef, level = caseDefLevels, labels = caseDefLabels),
         cases = as.integer(cases),
         n = as.integer(n)) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "<1 year", "1-4 years")) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "5-14 years", after = 2))

# Save the processed data, 'dat'
write_rds(x = dat, file = "../../data/processed/dat.rds")


dat2 <- dat %>%
  mutate(ageGroup = fct_collapse(ageGroup, `25-64 years` = c("25-34 years", "35-44 years", "45-54 years", "55-64 years"))) %>%
  mutate(ageGroup = fct_collapse(ageGroup,`65+ years` = c("65-74 years", "75-84 years", "85+ years"))) %>%
  group_by(Date, ageGroup, caseDef) %>%
  reframe(cases = sum(cases), n = sum(n))


dat2 <- dat2 %>%
  filter(caseDef %in% c("LIST", "SALM") & Date <= as.Date("2021-12-01") | caseDef %in% c("SHIG", "STEC"))


# dat2 <- dat %>%
#   mutate(ageGroup = fct_recode(ageGroup,
#                                `25-64 years` = "25-34 years",
#                                `25-64 years` = "35-44 years",
#                                `25-64 years` = "45-54 years",
#                                `25-64 years` = "55-64 years",
#                                `65+ years` = "65-74 years",
#                                `65+ years` = "75-84 years",
#                                `65+ years` = "85+ years")) %>%
#   group_by(Date, ageGroup, caseDef) %>%
#   reframe(cases = sum(cases), n = sum(n))

dat3 <- dat %>%
  filter(Date <= as.Date("2022-12-01")) %>%
  mutate(ageGroup = fct_collapse(ageGroup, `25-64 years` = c("25-34 years", "35-44 years", "45-54 years", "55-64 years"))) %>%
  mutate(ageGroup = fct_collapse(ageGroup,`65+ years` = c("65-74 years", "75-84 years", "85+ years"))) %>%
  group_by(Date, ageGroup, caseDef) %>%
  reframe(cases = sum(cases), n = sum(n))


write_rds(x = dat3, file = "../../data/processed/dat5.rds")



