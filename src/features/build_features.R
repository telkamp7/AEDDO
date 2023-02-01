
# Import libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(lubridate)
library(openxlsx)

locale(date_names = "da")
# Set locale
Sys.setlocale(category = "LC_ALL", locale = "da_DK.utf8")

# Import the data
FOLK1A <- read_csv2("../../data/raw/FOLK1A.csv")
diseaseData <- read.xlsx(xlsxFile = "../../data/raw/disease_data_raw.xlsx")
NUTS <- read_csv2("../../data/raw/NUTS_V1_2007.csv")

# Make reference table for 'landsdele' and 'kommuner'
NutsCor <- NUTS %>%
  select(KODE, NIVEAU, TITEL) %>%
  pivot_wider(id_cols = KODE, names_from = NIVEAU, values_from = TITEL) %>%
  rename(Region = `1`, Landsdel = `2`, Kommune = `3`) %>%
  fill(Region, Landsdel, .direction = "down") %>%
  mutate(Kommune = replace(Kommune, Kommune=="København", "Copenhagen")) %>%
  mutate(Landsdel = str_replace(Landsdel, "Landsdel ", "")) %>%
  mutate(Landsdel = str_replace(Landsdel, "Byen København", "København by")) %>%
  drop_na()
# ... Note: 'Byen København' is turned into 'København by', to match diseaseData

# Join FOLK1A and with reference table, to obtain 'landsdele'
FOLK1AxNutsCor <- inner_join(x = FOLK1A,
                             y = NutsCor,
                             by = c("OMRÅDE" = "Kommune"))

# Calculate FOLK1A per 'landsdel'
FOLK1AxLandsdel <- FOLK1AxNutsCor %>%
  group_by(ALDER, TID, Landsdel) %>%
  summarise(n = sum(INDHOLD))

# Calculate FOLK1A per 'Region'
FOLK1AxRegion <- FOLK1AxNutsCor %>%
  group_by(ALDER, TID, Region) %>%
  summarise(n = sum(INDHOLD))

# Finalize the data set
datTmp <- diseaseData %>%
  filter(maaned != "Uoplyst") %>%
  mutate(TIME = as.Date(paste0(year,"-",maaned,"-1"), format = "%Y-%B-%e")) %>%
  mutate(QUARTER = quarter(TIME)) %>%
  mutate(TID = paste0(year,"Q",QUARTER)) %>%
  inner_join(y = FOLK1AxLandsdel, by = c("age_label" = "ALDER",
                                         "landsdel_navn" = "Landsdel",
                                         "TID" = "TID")) %>%
  select(maaned, year, ageLabel = age_label,
         landsdel = landsdel_navn, caseDef = casedef, cases = value, n)
# ... Note: maaned and year is converted to quarters, in order to add 'n', which
# is the number of individuals from FOLK1A in each 'landsdel' and age group
# Hereafter, the data is combined in 'datTmp', which contains the final features

# Extract levels
maanedLevels <- unique(datTmp$maaned)
yearLevels <- unique(datTmp$year)
ageLabelLevels <- unique(datTmp$ageLabel)
landsdelLevels <- unique(datTmp$landsdel)
caseDefLevels <- unique(datTmp$caseDef)

# Finalize 'dat'
dat <- datTmp %>%
  mutate(maaned = factor(x = maaned, levels = maanedLevels),
         year = factor(x = year, levels = yearLevels),
         ageLabel = factor(x = ageLabel, levels = ageLabelLevels),
         landsdel = factor(x = landsdel, level = landsdelLevels),
         caseDef = factor(x = caseDef, levels = caseDefLevels),
         cases = as.integer(cases),
         n = as.integer(n)) %>%
  mutate(ageLabel = fct_relevel(ageLabel, "<1 år", "1-4 år")) %>%
  mutate(ageLabel = fct_relevel(ageLabel, "5-14 år", after = 2))

# Save the processed data, 'dat'
write_rds(x = dat, file = "../../data/processed/dat.rds")





