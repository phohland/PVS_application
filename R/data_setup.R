#'
#' PVS
#'
#' R file to setup raw data
#'


# Libraries ---------------------------------------------------------------

required_packages <- c("PVS", "tidyverse", "janitor", "readxl", "writexl")

for (package_name in required_packages) {
  if (!require(package_name,
               character.only = TRUE,
               quietly = TRUE)) {
    install.packages(package_name)
    library(package_name, character.only = TRUE)
  } else {
    library(package_name, character.only = TRUE)
  }
}

# Data setup --------------------------------------------------------------

# Blood pressure
raw_bp <-
  read.csv(
    "data/Blood_pressure_invasive.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_bp <- as_tibble(raw_bp)
names(raw_bp)[1] <- "pmid"
names(raw_bp)[2] <- "time"
names(raw_bp)[3] <- "systolicbp"
names(raw_bp)[4] <- "diastolicbp"
names(raw_bp)[5] <- "map"
raw_bp <- raw_bp %>%
  clean_names %>%
  select("pmid", "time", "systolicbp", "diastolicbp", "map") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    systolicbp = as.double(systolicbp),
    diastolicbp = as.double(diastolicbp),
    map = as.double(map),
  )

# Body temperature
raw_temp <-
  read.csv(
    "data/Body_temperature.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_temp <- as_tibble(raw_temp)
names(raw_temp)[1] <- "pmid"
names(raw_temp)[2] <- "time"
names(raw_temp)[3] <- "temperature"
raw_temp <- raw_temp %>%
  clean_names %>%
  select("pmid", "time", "temperature") %>%
  mutate(time = as.POSIXct(time, tz = ""),
         temperature = as.double(temperature))

# Heart rate
raw_hr <-
  read.csv(
    "data/Heart_rate.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_hr <- as_tibble(raw_hr)
names(raw_hr)[1] <- "pmid"
names(raw_hr)[2] <- "time"
names(raw_hr)[3] <- "hr"
raw_hr <- raw_hr %>%
  clean_names %>%
  select("pmid", "time", "hr") %>%
  mutate(time = as.POSIXct(time, tz = ""),
         hr = as.double(hr))

# Patient information
raw_master <-
  read.csv(
    "data/Patient_master_data.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_master <- as_tibble(raw_master)
raw_master <- raw_master %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    date_of_birth = as.Date(date_of_birth),
    sex = as.factor(sex)
  )

# Patient information with sepsis labels
raw_master_label <-
  read.csv(
    "data/Patient_master_data_label.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_master_label <- as_tibble(raw_master_label)
raw_master_label <- raw_master_label %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    discharge = as.Date(discharge),
    admission = as.Date(admission),
    label = as.factor(label)
  )

# Pulse
raw_pulse <-
  read.csv(
    "data/Pulse.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_pulse <- as_tibble(raw_pulse)
names(raw_pulse)[2] <- "time"
names(raw_pulse)[3] <- "pulse"
raw_pulse <- raw_pulse %>%
  clean_names %>%
  select("pmid", "time", "pulse") %>%
  mutate(time = as.POSIXct(time, tz = ""),
         pulse = as.double(pulse))

# Respiratory Rate
raw_resp <-
  read.csv(
    "data/Respiratory_rate.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_resp <- as_tibble(raw_resp)
names(raw_resp)[2] <- "time"
names(raw_resp)[3] <- "rr"
raw_resp <- raw_resp %>%
  clean_names %>%
  select("pmid", "time", "rr") %>%
  mutate(time = as.POSIXct(time, tz = ""),
         rr = as.double(rr))

# Oxygen saturation
raw_tom <-
  read.csv(
    "data/Transcutaneous_oxygen_measurement.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
raw_tom <- as_tibble(raw_tom)
names(raw_tom)[2] <- "time"
names(raw_tom)[3] <- "sao2"
raw_tom <- raw_tom %>%
  clean_names %>%
  select("pmid", "time", "sao2") %>%
  mutate(time = as.POSIXct(time, tz = ""),
         sao2 = as.double(sao2))

# Merge all data sets into one final data set
raw_list <-
  list(raw_hr, raw_pulse, raw_tom, raw_resp, raw_bp, raw_temp)
raw <- reduce(raw_list, left_join, by = c("pmid", "time"))

# Delete time duplicates
raw <- raw %>% distinct(time, .keep_all = TRUE)

# Deletion of unrealistic entries
raw <- raw %>%
  mutate(
    temperature = ifelse(temperature > 34, temperature, NA)
  )

# Drop NAs
raw <- na.omit(raw)

# Rename column names in raw data set
names(raw)[1] <- "id"
names(raw)[2] <- "timepoint"

# Merge patient information
master_list <- list(raw_master, raw_master_label)
master <- reduce(master_list, left_join, by = "pmid")
master <- master[master$pmid %in% raw$id, ]

# Add the sex as "condition"
raw$condition <-
  as.factor(master$sex[match(as.character(raw$id), master$pmid)])

# File creation -----------------------------------------------------------

# raw
write_xlsx(raw, "output/raw.xlsx")

# master
write_xlsx(master, "output/master.xlsx")
