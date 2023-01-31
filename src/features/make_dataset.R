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
library(ggplot2)


# Disease codes ---------------------------------------------------------------------

# Disease codes
# SYPH - Syphilis
# LEGI - Legionella
# SHIG - Shigella
# TUBE - Tubercolosis
# LEPT - Leptospirosis
# VTEC - Shigatoxin producing Escherichia coli.
diseaseCodes <- c("SYPH", "LEGI", "SHIG", "TUBE", "LEPT", "VTEC")

# Scrape data from multiple diseases from SSI API - and gather them in a data.table
dat_samlet <- data.table()
for (casedef in diseaseCodes){
  dat <- data.table()
  print(casedef)
  for (ldel in 1:12){
    for (age_group in 1:12){
      tt <- read.table(file=paste0("https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?aldersgruppe=",age_group,"&landsdelkode=",ldel,"&sygdomskode=",casedef,"&xaxis=Aar&yaxis=Maaned&aar=1994|2023&show=Table&lang=DA"), sep=";", header=TRUE)
      setDT(tt)
      print(c(ldel, age_group, dim(tt)))
      tt2 <- melt(tt, id.vars = 1, variable.name = "year")
      tt2[, landsdel := ldel][, age_group := age_group]
      dat <- rbind(dat, tt2)
    }
  }
  
  setnames(dat, "X", "maaned")
  dat[, year := as.integer(substr(year, 2,5))]
  
  # Så mangler vi bare de rigtige landsdelsnavne og aldersgrupper
  tt <- fread(input = "https://statistik.ssi.dk/api/ssi/surveillance/DiseaseCsv?sygdomskode=SYPH&xaxis=Landsdelkode&yaxis=Aldersgruppe&aar=1994|2023&show=Table&lang=DA", encoding = "Latin-1")
  dat[data.table(age_group=1:12, age_label = tt$V1), on = "age_group", age_label := i.age_label]
  dat[data.table(landsdel=1:12, landsdel_navn = names(tt)[-1]), on = "landsdel", landsdel_navn := i.landsdel_navn]
  dat[, casedef := casedef]
  
  dat_samlet <- rbind(dat_samlet, dat)
  
}

# ggplot(dat_samlet[year >= 2020], aes(x = maaned, y = value, colour = factor(year), shape = age_label)) +
#   geom_point() +
#   facet_wrap(facets = vars(landsdel_navn))


write.xlsx(dat_samlet, file="../../data/raw/disease_data_raw.xlsx")

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
agegroups <- c("sum(<1 år=0)",
               "sum(1-4 år=1;2;3;4)",
               "sum(5-14 år=5;6;7;8;9;10;11;12;13;14)",
               "sum(15-24 år=15;16;17;18;19;20;21;22;23;24)",
               "sum(25-34 år=25;26;27;28;29;30;31;32;33;34)",
               "sum(35-44 år=35;36;37;38;39;40;41;42;43;44)",
               "sum(45-54 år=45;46;47;48;49;50;51;52;53;54)",
               "sum(55-64 år=55;56;57;58;59;60;61;62;63;64)",
               "sum(65-74 år=65;66;67;68;69;70;71;72;73;74)",
               "sum(75-84 år=75;76;77;78;79;80;81;82;83;84)",
               paste0("sum(85+ år=",paste(85:125,collapse = ";"),")")
)

FOLK1A <- tibble()
for(age in agegroups){
  
  variables2 <- list(
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
    variables = variables2
  )
  
  # GET data 
  result <- POST(url = dataEndpoint,
                 body = callBodyData,
                 encode = "json")
  # Read CSV file
  age.iter <- read_csv2(content(result, type = "text", encoding = "UTF-8"))
  
  FOLK1A <- bind_rows(FOLK1A, age.iter)
  
}

write_csv2(x = FOLK1A, file = "../../data/raw/FOLK1A.csv")

# FOLK1A$ALDER
# FOLK1A %>%
#   ggplot(mapping = aes(x = TID, y = INDHOLD)) +
#   geom_point() +
#   facet_wrap(facets = vars(ALDER))

