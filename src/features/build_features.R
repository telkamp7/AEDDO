
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
# FOLK1Akon <- read_csv2("../../data/raw/FOLK1Akon.csv")
diseaseData <- read_rds(file = "../../data/raw/disease_data_raw_new.Rds")
NUTS <- read_csv2("../../data/raw/NUTS_V1_2007.csv")

# # # Calculate FOLK1A per 'agegroup'
FOLK1AxAgegroup <- FOLK1A %>%
  group_by(ALDER, TID) %>%
  reframe(n = sum(INDHOLD))

# Check for 'Uoplyst'
diseaseData %>%
  group_by(ageGroup) %>%
  summarize(sum(cases)) %>%print(n=144)


# Change 'år' to 'years'
diseaseData <- diseaseData %>%
  mutate(ageGroup = if_else(
    ageGroup == "<1 år",
    str_replace(ageGroup, pattern = "år", replacement = "year"),
    str_replace(ageGroup, pattern = "år", replacement = "years")
  ))

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
#   "GonorÃ©",
#   "Hepatitis A",
#   "Hepatitis B",
#   "Hepatitis C",
#   "HÃ¦mophilus influenza meningitis",
#   "HIV infektion",
#   "Legionella",
#   "Leptospirosis",
#   "MÃ¦slinger",
#   "Andre meningitis",
#   "Meningokoksygdom",
#   "MPOX",
#   "FÃ¥resyge",
#   "Neuroborreliose",
#   "Ornitose",
#   "Kighoste",
#   "Pneumokok meningitis",
#   "RÃ¸de hunde",
#   "Shigella",
#   "Syfilis",
#   "Stivkrampe",
#   "Tubercolosis",
#   "Tyfus / paratyfus",
#   "Shiga- og veratoxin producerende E. coli."
# )

caseDefLabels <- c("STEC",
                   "SHIL",
                   "VTEC",
                   "LIST",
                   "SALM")

# Finalize 'dat'
dat <- datTmp %>%
  mutate(ageGroup = factor(x = ageGroup, levels = ageGroupLevels),
         caseDef = factor(x = caseDef, level = caseDefLevels, labels = caseDefLabels),
         cases = as.integer(cases),
         n = as.integer(n)) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "<1 year", "1-4 years")) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "5-14 years", after = 2)) %>%
  ungroup()

# Save the processed data, 'dat'
write_rds(x = dat, file = "../../data/processed/dat.rds")


dat2 <- dat %>%
  mutate(ageGroup = fct_collapse(ageGroup, `25-64 years` = c("25-34 years", "35-44 years", "45-54 years", "55-64 years"))) %>%
  mutate(ageGroup = fct_collapse(ageGroup,`65+ years` = c("65-74 years", "75-84 years", "85+ years"))) %>%
  group_by(Date, ageGroup, caseDef) %>%
  reframe(cases = sum(cases), n = sum(n))

write_rds(x = dat2, file = "../../data/processed/dat2.rds")
