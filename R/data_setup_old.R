#'
#' Quantifying patient vital status: a translational composite score
#'
#' R file to setup raw data
#'

# Libraries ---------------------------------------------------------------

# Prior a installation of each package is needed
install.packages(c("hRELSA",
                   "tidyverse",
                   "janitor",
                   "readxl"))

library(hRELSA)

library(tidyverse)
library(janitor)
library(readxl)

# Data setup -------------------------------------------------------------

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
raw <- NULL
raw <-
  merge(
    raw_hr,
    raw_pulse,
    by = c("pmid", "time"),
    all.x = TRUE,
    all.y = TRUE
  )
raw <-
  merge(
    raw,
    raw_tom,
    by = c("pmid", "time"),
    all.x = TRUE,
    all.y = TRUE
  )
raw <-
  merge(
    raw,
    raw_resp,
    by = c("pmid", "time"),
    all.x = TRUE,
    all.y = TRUE
  )
raw <-
  merge(
    raw,
    raw_bp,
    by = c("pmid", "time"),
    all.x = TRUE,
    all.y = TRUE
  )
raw <-
  merge(
    raw,
    raw_temp,
    by = c("pmid", "time"),
    all.x = TRUE,
    all.y = TRUE
  )
raw <- as_tibble(raw)

# Delete duplicates
raw <- raw %>% distinct(time, .keep_all = TRUE)

# Drop NAs
raw <- na.omit(raw)

# Rename column names in raw data set
names(raw)[1] <- "id"
names(raw)[2] <- "timepoint"

# Merge patient information and sepsis labels
raw_master_full <-
  merge(
    raw_master,
    raw_master_label,
    by = c("pmid"),
    all.x = TRUE,
    all.y = TRUE
  )
raw_master_full <-
  raw_master_full[raw_master_full$pmid %in% raw$id, ]
levels(raw_master_full$label) <-
  c(levels(raw_master_full$label), "NoLabel")
raw_master_full[is.na(raw_master_full)] <- "NoLabel"
raw_master_full <- as_tibble(raw_master_full)

# Add the sepsis label as "treatment"
raw$treatment <- NA
raw$treatment <- as.factor(raw$treatment)
levels(raw$treatment) <- levels(raw_master_full$label)
for (l in 1:nrow(raw)) {
  label <-
    raw_master_full$label[raw_master_full$pmid == as.character(raw$id[l])]
  raw$treatment[l] <- as.factor(label)
}

# Add the sex as "condition"
raw$condition <- NA
raw$condition <- as.factor(raw$condition)
levels(raw$condition) <- levels(raw_master_full$sex)
for (l in 1:nrow(raw)) {
  sex <-
    raw_master_full$sex[raw_master_full$pmid == as.character(raw$id[l])]
  raw$condition[l] <- as.factor(sex)
}

# Deletion of unrealistic entries
raw <- raw %>%
  mutate(
    systolicbp = ifelse(systolicbp > 0, systolicbp, NA),
    map = ifelse(map > 0, map, NA),
    rr = ifelse(rr > 0, rr, NA),
    temperature = ifelse(temperature > 34, temperature, NA),
    map = ifelse(map < systolicbp, map, NA)
  )

# File creation -----------------------------------------------------------

# raw
write.csv(raw, file = "output/raw.csv")

#raw_master_full
write.csv(raw_master_full, file = "output/raw_master_full.csv")
