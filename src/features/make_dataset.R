# Adapted by: Kasper Schou Telkamp
# Originally made by: Lasse Engbo Christiansen

# Import relevant libraries
library(data.table)
library(openxlsx)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)


# Disease codes ---------------------------------------------------------------------

diseaseCodes <- c("STEC",
                  "SHIL",
                  "VTEC",
                  "LIST",
                  "SALM")

# Scrape data from multiple diseases from SSI API - and gather them in a data.table
# diseaseData <- data.table()
diseaseData <- tibble()

for (cdf in diseaseCodes){
  dat <- tibble()
  # dat <- data.table()
  print(cdf)
  # for (ldel in 1:12){
  for (agegrp in 1:12){
    # tt <- read.table(file=paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",agegrp,"&sygdomskode=",cdf,"&xaxis=Aar&yaxis=Maaned&aar=2008|2023&show=Table&lang=DA"), sep=";", header=TRUE)
    tt <- read_csv2(file = paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",agegrp,"&sygdomskode=",cdf,"&xaxis=Aar&yaxis=Maaned&aar=2008|2023&show=Table&lang=DA"))
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
                         # landsdel = ldel,
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
    dat <- bind_rows(dat, tt2)

  }
  
  dat <- dat %>%
    mutate(caseDef = cdf) %>%
    group_by(caseDef)
  
  diseaseData <- rbind(diseaseData, dat)
}

# Så mangler vi bare de rigtige køn og aldersgrupper
tt <-read_csv2(file = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SALM&xaxis=Landsdelkode&yaxis=Aldersgruppe&aar=1994|2023&show=Table&lang=DA", locale=locale(encoding="latin1"))
# tt3 <- read_csv2(file = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SHIL&xaxis=Aar&yaxis=Kon&aar=2001|2023&show=Table&lang=DA", locale=locale(encoding="latin1"))


diseaseData <- diseaseData %>%
  mutate(ageGroup = tt$...1[ageGroup])

# diseaseData <- diseaseData %>%
#   unnest(data) %>%
#   mutate(ageGroup = tt$...1[ageGroup],
#          sex = tt3$...1[sex+1]
#          ) %>%
#   nest()

write_rds(x = diseaseData, file = "../../data/raw/disease_data_raw_new.Rds")

# Population data -------------------------------------------------------------------

# Define endpoint for metadata
metadataEndpoint <- "https://api.statbank.dk/v1/tableinfo"

# Define table to look up
tableId <- "folk1a"

# Define english as preferedd language
language <- "en"

# Define query for metadata
callBodyMetadata <- list(
  lang = language,
  table = tableId
)

# GET meta data 
result <- POST(url = metadataEndpoint,
               body = callBodyMetadata,
               encode = "json")

# Extract contents of the result
fullResult <- fromJSON(content(result))

# Define endpoint for the data
dataEndpoint <- "https://api.statbank.dk/v1/data"

# Find index for all Municipalities
muniIdx <- grep(pattern = c("All|Region"),
                x = fullResult$variables$values[[1]]$text,
                invert = TRUE)

# Find all municipality IDs
muniId <- as.integer(fullResult$variables$values[[1]]$id[muniIdx])

# Define agegroups
agegroups <- c("sum(<1 year=0)",
               "sum(1-4 years=1;2;3;4)",
               "sum(5-14 years=5;6;7;8;9;10;11;12;13;14)",
               "sum(15-24 years=15;16;17;18;19;20;21;22;23;24)",
               "sum(25-34 years=25;26;27;28;29;30;31;32;33;34)",
               "sum(35-44 years=35;36;37;38;39;40;41;42;43;44)",
               "sum(45-54 years=45;46;47;48;49;50;51;52;53;54)",
               "sum(55-64 years=55;56;57;58;59;60;61;62;63;64)",
               "sum(65-74 years=65;66;67;68;69;70;71;72;73;74)",
               "sum(75-84 years=75;76;77;78;79;80;81;82;83;84)",
               paste0("sum(85+ years=",paste(85:125,collapse = ";"),")")
)

FOLK1A <- tibble()
for(age in agegroups){
  
  variables <- list(
    list(code = "område", values = muniId),
    list(code = "køn", values = "TOT"),
    list(code = "alder", values = list(age)),
    list(code = "tid", values = list("*"))
  )
  
  # Define query for metadata
  callBodyData <- list(
    table = tableId,
    lang = language,
    format = "CSV",
    variables = variables
  )
  
  # GET data 
  result <- POST(url = dataEndpoint,
                 body = callBodyData,
                 encode = "json")
  # Read CSV file
  age.iter <- read_csv2(content(result, type = "text", encoding = "UTF-8"))
  
  FOLK1A <- bind_rows(FOLK1A, age.iter)
  
}

# Save FOLK1A data
write_csv2(x = FOLK1A, file = "../../data/raw/FOLK1A.csv")

variables <- list(
  list(code = "område", values = "000"),
  list(code = "køn", values = 1:2),
  list(code = "alder", values = paste0("sum(All=0",paste(0:125,collapse = ";"),")")),
  list(code = "tid", values = list("*"))
)

# Define query for metadata
callBodyData <- list(
  table = tableId,
  lang = language,
  format = "CSV",
  variables = variables
)
# GET data 
result <- POST(url = dataEndpoint,
               body = callBodyData,
               encode = "json")

# Read CSV file
FOLK1Akon <- read_csv2(content(result, type = "text", encoding = "UTF-8"))

# Save FOLK1A data
write_csv2(x = FOLK1Akon, file = "../../data/raw/FOLK1Akon.csv")

# Download CSV containing specification of 'landsdele' and 'kommuner'
urlNutsV12007 <- "https://www.dst.dk/klassifikationsbilag/5d18d1e0-400b-4505-92ad-6782915980a3csv_da"

# Provide destination for file
destFile <- "../../data/raw/NUTS_V1_2007.csv"

# Apply download.file function in 
download.file(urlNutsV12007, destFile)

