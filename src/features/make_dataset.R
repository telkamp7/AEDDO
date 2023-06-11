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

diseaseCodes <- c("VTEC",
                  "SHIG",
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

      # tt2 <- melt(tt, id.vars = 1, variable.name = "year")
      # tt2[, landsdel := ldel][, age_group := age_group]
      # dat <- rbind(dat, tt2)
    }
  # }

  dat <- dat %>%
    mutate(caseDef = cdf)

  diseaseData <- rbind(diseaseData, dat)
}


# for (cdf in diseaseCodes){
#   dat <- tibble()
#   print(cdf)
#   for (agegrp in 1:12){
#     tt <- read_csv2(file = paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",agegrp,"&sygdomskode=",cdf,"&xaxis=Aar&yaxis=Maaned&aar=2001|2023&show=Table&lang=DA"))
#     print(c(agegrp, dim(tt)))
#     if(dim(tt)[2] == 2){
#       tt2 <- expand_grid(month = c("Januar", "Februar", "Marts", "April", "Maj", "Juni",
#                                    "Juli", "August", "September", "Oktober", "November", "December", "Uoplyst"),
#                          year = as.character(1994:2022),
#                          cases = 0,
#                          ageGroup = agegrp)
#     }else{
#       tt2 <- tt %>%
#         rename("month" = ...1) %>%
#         pivot_longer(cols = -month, names_to = "year", values_to = "cases") %>%
#         mutate(ageGroup = agegrp)
#     }
#     dat <- bind_rows(dat, tt2)
#     
#   }
#   dat <- dat %>%
#     mutate(caseDef = cdf)
#   
#   diseaseData <- rbind(diseaseData, dat)
# }


# # Disease codes
# diseaseCodes <- c(
#   "AIDS",  # AIDS
#   "BOTU",  # Botulisme
#   "GONO",  # Gonoré
#   "HEPA",  # Hepatitis A
#   "HEPB",  # Hepatitis B
#   "HEPC",  # Hepatitis C
#   "HIB",   # Hæmophilus influenza meningitis
#   "HIV",   # HIV infektion
#   "LEGI",  # Legionella
#   "LEPT",  # Leptospirosis
#   "MEAS",  # Mæslinger
#   "MENAN", # Andre meningitis
#   "MENI",  # Meningokoksygdom
#   "MPOX",  # MPOX
#   "MUMP",  # Fåresyge
#   "NEBO",  # Neuroborreliose
#   "ORNI",  # Ornitose
#   "PERT",  # Kighoste
#   "PNEU",  # Pneumokik meningitis
#   "RUBE",  # Røde hunde
#   "SHIG",  # Shigella
#   "SYPH",  # Syfilis
#   "TETA",  # Stivkrampe
#   "TUBE",  # Tubercolosis
#   "TYPH",  # Tyfus / paratyfus
#   "VTEC"   # Shigatoxin producerende / veratoxin producerende E. coli.
# )

# for (cdf in diseaseCodes){
#   dat <- tibble()
#   # dat <- data.table()
#   print(cdf)
#   for (ldel in 1:12){
#     for (agegrp in 1:12){
#       # tt <- read.table(file=paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",age_group,"&landsdelkode=",ldel,"&sygdomskode=",casedef,"&xaxis=Aar&yaxis=Maaned&aar=1994|2023&show=Table&lang=DA"), sep=";", header=TRUE)
#       tt <- read_csv2(file = paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",agegrp,"&landsdelkode=",ldel,"&sygdomskode=",cdf,"&xaxis=Aar&yaxis=Maaned&aar=1994|2022&show=Table&lang=DA"))
#       # setDT(tt)
#       print(c(ldel, agegrp, dim(tt)))
#       # Check if there have been any registered cases in the group
#       # -- There is a bug in the API where it returns only a single year for a group, if
#       # -- there havent been any cases.
#       if(dim(tt)[2] == 2){
#         tt2 <- expand_grid(month = c("Januar", "Februar", "Marts", "April", "Maj", "Juni",
#                                      "Juli", "August", "September", "Oktober", "November", "December", "Uoplyst"),
#                            year = as.character(1994:2022),
#                            cases = 0,
#                            landsdel = ldel,
#                            ageGroup = agegrp)
#       }else{
#       tt2 <- tt %>%
#         rename("month" = ...1) %>%
#         pivot_longer(cols = -month, names_to = "year", values_to = "cases") %>%
#         mutate(landsdel = ldel, ageGroup = agegrp)
#       }
#       dat <- bind_rows(dat, tt2)
#       
#       # tt2 <- melt(tt, id.vars = 1, variable.name = "year")
#       # tt2[, landsdel := ldel][, age_group := age_group]
#       # dat <- rbind(dat, tt2)
#     }
#   }
#   
#   dat <- dat %>%
#     mutate(caseDef = cdf)
# 
#   diseaseData <- rbind(diseaseData, dat)
#   
# }

# setnames(dat, "X", "maaned")
# dat[, year := as.integer(substr(year, 2,5))]

# Så mangler vi bare de rigtige landsdelsnavne og aldersgrupper
tt <-read_csv2(file = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SALM&xaxis=Landsdelkode&yaxis=Aldersgruppe&aar=1994|2023&show=Table&lang=DA", locale=locale(encoding="latin1"))
# diseaseData <- diseaseData %>%
#   mutate(ageGroup = tt$...1[ageGroup],
#          landsdel = colnames(tt)[-1][landsdel])
diseaseData <- diseaseData %>%
  mutate(ageGroup = tt$...1[ageGroup])

# tt <- fread(input = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SYPH&xaxis=Landsdelkode&yaxis=Aldersgruppe&aar=1994|2023&show=Table&lang=DA", encoding = "Latin-1")
# dat[data.table(age_group=1:12, age_label = tt$V1), on = "age_group", age_label := i.age_label]
# dat[data.table(landsdel=1:12, landsdel_navn = names(tt)[-1]), on = "landsdel", landsdel_navn := i.landsdel_navn]
# dat[, casedef := casedef]

# Convert into tibble, to allow for easy management in the tidyverse
write.xlsx(diseaseData, file="../../data/raw/disease_data_raw.xlsx")

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

variables <- list(
  list(code = "område", values = muniId),
  list(code = "køn", values = "TOT"),
  list(code = "alder", values = 0:125),
  list(code = "tid", values = list("*"))
)

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

# Download CSV containing specification of 'landsdele' and 'kommuner'
urlNutsV12007 <- "https://www.dst.dk/klassifikationsbilag/5d18d1e0-400b-4505-92ad-6782915980a3csv_da"

# Provide destination for file
destFile <- "../../data/raw/NUTS_V1_2007.csv"

# Apply download.file function in 
download.file(urlNutsV12007, destFile)

