# Import relevant libraries
library(openxlsx)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(forcats)


diseaseCode <- "CAMP"

diseaseData <- tibble()

for (agegrp in 1:12){
  tt <- read_csv2(file = paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",agegrp,"&sygdomskode=",diseaseCode,"&xaxis=Aar&yaxis=Maaned&aar=2008|2023&show=Table&lang=DA"))
  # setDT(tt)
  # print(c(ldel, agegrp, dim(tt)))
  print(c(agegrp, dim(tt)))
  # Check if there have been any registered cases in the group
  # -- There is a bug in the API where it returns only a single year for a group, if
  # -- there havent been any cases.
  if(dim(tt)[2] == 2){
    tt2 <- expand_grid(month = c("Januar", "Februar", "Marts", "April", "Maj", "Juni",
                                 "Juli", "August", "September", "Oktober", "November", "December", "Uoplyst"),
                       year = as.character(2008:2022),
                       cases = 0,
                       ageGroup = agegrp)
  }else{
    tt2 <- tt %>%
      rename("month" = ...1) %>%
      pivot_longer(cols = -month, names_to = "year", values_to = "cases") %>%
      mutate(
        # landsdel = ldel,
        ageGroup = agegrp
      )
  }
  diseaseData <- bind_rows(diseaseData, tt2)
}

diseaseData <- diseaseData %>%
  mutate(caseDef = diseaseCode) %>%
  group_by(caseDef)

tt <-read_csv2(file = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=CAMP&xaxis=Landsdelkode&yaxis=Aldersgruppe&aar=1994|2023&show=Table&lang=DA", locale=locale(encoding="latin1"))
# tt3 <- read_csv2(file = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SHIL&xaxis=Aar&yaxis=Kon&aar=2001|2023&show=Table&lang=DA", locale=locale(encoding="latin1"))


diseaseData <- diseaseData %>%
  mutate(ageGroup = tt$...1[ageGroup])

# Import the data
FOLK1A <- read_csv2("../../data/raw/FOLK1A.csv")

# # # Calculate FOLK1A per 'agegroup'
FOLK1AxAgegroup <- FOLK1A %>%
  group_by(ALDER, TID) %>%
  reframe(n = sum(INDHOLD))


# Check for 'Uoplyst'
diseaseData %>%
  group_by(ageGroup) %>%
  summarize(sum(cases)) %>% print(n=144)

# Change '?r' to 'years'
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

caseDefLabels <- c("CAMP")

# Finalize 'dat'
dat <- datTmp %>%
  filter(Date <= as.Date("2022-12-31")) %>%
  mutate(ageGroup = factor(x = ageGroup, levels = ageGroupLevels),
         caseDef = factor(x = caseDef, level = caseDefLevels, labels = caseDefLabels),
         cases = as.integer(cases),
         n = as.integer(n)) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "<1 year", "1-4 years")) %>%
  mutate(ageGroup = fct_relevel(ageGroup, "5-14 years", after = 2)) %>%
  ungroup()

dat2 <- dat %>%
  mutate(ageGroup = fct_collapse(ageGroup, `25-64 years` = c("25-34 years", "35-44 years", "45-54 years", "55-64 years"))) %>%
  mutate(ageGroup = fct_collapse(ageGroup,`65+ years` = c("65-74 years", "75-84 years", "85+ years"))) %>%
  group_by(Date, ageGroup, caseDef) %>%
  reframe(cases = sum(cases), n = sum(n))
write_rds(x = dat2, file = "../../data/processed/CAMP.rds")  
