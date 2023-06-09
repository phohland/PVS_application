#'
#' hRELSA:
#' A translational approach to quantify a patient’s disease severity in real-time
#'
#' R file to reproduce the data and figures
#'

# Libraries ---------------------------------------------------------------

# Prior a installation of each package is needed:
install.packages(
  c(
    "hRELSA",
    "tidyverse",
    "janitor",
    "readxl",
    "ggplot2",
    "patchwork",
    "ggsignif",
    "psych",
    "car",
    "lsr",
    "BSDA"
  )
)

library(hRELSA)

library(tidyverse)
library(janitor)
library(readxl)

library(ggplot2)
library(patchwork)
library(ggsignif)

library(psych) # for describeBy and cor
library(car) # for leveneTest
library(lsr) # for cohensD
library(BSDA) # for z.test

# Data setup -------------------------------------------------------------

# The used raw data can be found in the output folder. The data_setup.R shows,
# how it was created. Here they just get fetched.

raw <- as_tibble(read.csv("output/raw.csv")[,-1])
double_vars <-
  c("hr",
    "pulse",
    "sao2",
    "rr",
    "systolicbp",
    "diastolicbp",
    "map",
    "temperature")
raw <- raw %>% mutate(
  timepoint = as.POSIXct(timepoint, tz = ""),
  treatment = as.factor(treatment),
  condition = as.factor(condition)
) %>%
  mutate(across(all_of(double_vars), as.double))

raw_master_full <-
  as_tibble(read.csv("output/raw_master_full.csv")[,-1])
raw_master_full <- raw_master_full %>% mutate(
  pmid = as.factor(pmid),
  sex = as.factor(sex),
  date_of_birth = as.Date(date_of_birth),
  #date_of_die = as.Date(date_of_die),
  admission = as.Date(admission),
  discharge = as.Date(discharge),
  label = as.factor(label)
)

# hRELSA  -----------------------------------------------------------------

# Select the variables to use the hRELSA with
vars <- c("hr", "sao2", "rr", "map", "temperature")
turnvars <- NULL
ambivars <- c("hr", "rr", "map", "temperature")
zvars <- NULL
dropvars <- NULL

# Get reference values
raw_norm <- read_excel("data/reference_values.xlsx", sheet = 1)

# Deletion of pulse
raw_norm <- raw_norm %>% select("age", all_of(vars))

# Generate day column
raw <-
  hrelsa_days(raw,
              format = "timecode",
              formthis = "timepoint",
              newdayone = TRUE)

# Format the data
dat <-
  hrelsa_format(
    raw,
    id = "id",
    time = "time",
    treatment = "treatment",
    condition = "condition",
    vars = vars,
    included_realtime = "timepoint"
  )

# Write dat in a csv
write.csv(dat, file = "output/dat.csv")

# Fetch data for maximal severity evaluation
reference_dat <- dat %>% filter(treatment == "NoSepsis")
bsl <-
  hrelsa_adaptive_baselines(
    dat,
    reference_dat,
    vars = vars,
    turnvars = turnvars,
    ambivars = ambivars,
    realtime = "timepoint",
    dob_dat = raw_master_full,
    dob_data_id_col = 1,
    dob_data_dob_col = 3,
    norm_dat = raw_norm,
    norm_dat_names = names(raw_norm)
  )
pre <- bsl$pre
age_pre <- bsl$age_pre

# Write pre in a csv
write.csv(pre, file = "output/pre.csv")

# Generate final data
final <-
  hrelsa_final(
    pre,
    bsl,
    drop = dropvars,
    turnvars = turnvars,
    ambivars = ambivars,
    zvars = zvars
  )

# Generate some analysis
analysis <- hrelsa_analysis(final)

# Time label setup --------------------------------------------------------

final <- arrange(final, id, time)
dat <- arrange(dat, id, time)

final$timepoint <- dat$timepoint

# Ventilation (one time point)
ventilation <-
  read.csv(
    "data/Ventilation.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
ventilation <- as_tibble(ventilation)
names(ventilation)[1] <- "pmid"
names(ventilation)[4] <- "time"
names(ventilation)[5] <- "ventilation"
ventilation <- ventilation %>%
  select("pmid", "time", "ventilation") %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    time = as.POSIXct(time, tz = ""),
    ventilation = as.factor(ventilation)
  )

# Temperature Regulations
temp_reg <-
  read.csv(
    "data/Temperature_regulation.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
regulated_ids <- unique(temp_reg$PMID)

# SIRS labels (duration)
sirs <-
  read.csv(
    "data/SIRS_label.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
sirs <- as_tibble(sirs)
names(sirs)[1] <- "pmid"
names(sirs)[3] <- "start"
names(sirs)[4] <- "end"
sirs <- sirs %>%
  select("pmid", "start", "end") %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    start = as.POSIXct(start, tz = ""),
    end = as.POSIXct(end, tz = ""),
  )

# Sepsis labels (duration)
sepsis <-
  read.csv(
    "data/Sepsis_label.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
sepsis <- as_tibble(sepsis)
names(sepsis)[1] <- "pmid"
names(sepsis)[3] <- "start"
names(sepsis)[4] <- "end"
names(sepsis)[5] <- "sepsis"
sepsis <- sepsis %>%
  select("pmid", "start", "end", "sepsis") %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    start = as.POSIXct(start, tz = ""),
    end = as.POSIXct(end, tz = ""),
    sepsis = as.factor(sepsis)
  )

# Catecholamine labels (duration)
cat <-
  read.csv(
    "data/Medication_catecholamines.csv",
    sep = ";",
    header = TRUE,
    dec = ".",
    row.names = NULL
  )
cat <- as_tibble(cat)
names(cat)[1] <- "pmid"
names(cat)[3] <- "start"
names(cat)[4] <- "end"
cat <- cat %>%
  select("pmid", "start", "end") %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    start = as.POSIXct(start, tz = ""),
    end = as.POSIXct(end, tz = "")
  )

# Adds SIRS TRUE or FALSE to each entry in final
final$sirs <- FALSE

for (i in 1:nrow(final)) {
  if (final$id[i] %in% sirs$pmid) {
    for (o in 1:nrow(sirs)) {
      if (sirs$pmid[o] == final$id[i]) {
        if (final$timepoint[i] >= sirs$start[o] & final$timepoint[i] <= sirs$end[o]) {
          final$sirs[i] <- TRUE
        }
      }
    }
  }
}

# Adds sepsis variations TRUE or FALSE to each entry in final
final$proven_sepsis <- FALSE
final$sus_sepsis <- FALSE
final$infection <- FALSE

for (i in 1:nrow(final)) {
  if (final$id[i] %in% sepsis$pmid) {
    for (o in 1:nrow(sepsis)) {
      if (sepsis$pmid[o] == final$id[i]) {
        if (final$timepoint[i] >= sepsis$start[o] & final$timepoint[i] <= sepsis$end[o]) {
          if (sepsis$sepsis[o] == "Proven Sepsis") {
            final$proven_sepsis[i] <- TRUE
          } else if (sepsis$sepsis[o] == "Suspected Sepsis") {
            final$sus_sepsis[i] <- TRUE
          } else if (sepsis$sepsis[o] == "Infection") {
            final$infection[i] <- TRUE
          }
        }
      }
    }
  }
}

# Adds ventilation TRUE or FALSE to each entry in final
final$ventilation <- FALSE

matching_indices <- which(final$id %in% ventilation$pmid)
final_matching <- final[matching_indices, ]
ventilation_matching <- ventilation[match(final_matching$id, ventilation$pmid), ]

matching_rows <- final_matching$timepoint == ventilation_matching$time

final$ventilation[matching_indices] <- matching_rows

# Adds catecholamines TRUE or FALSE to each entry in final
final$medication <- FALSE

matching_indices <- which(final$id %in% cat$pmid)
final_matching <- final[matching_indices, ]
cat_matching <- cat[match(final_matching$id, cat$pmid), ]

matching_rows <- final_matching$timepoint >= cat_matching$start & final_matching$timepoint <= cat_matching$end

final$medication[matching_indices] <- matching_rows

# Descriptive patient data analysis ---------------------------------------

desc <- raw_master_full
desc$stay <- desc$discharge - desc$admission

length(unique(desc$pmid))
summary(desc$sex)
summary(desc$date_of_birth)
summary(desc$label)

mean(desc$stay)
median(desc$stay)
min(desc$stay)
max(desc$stay)

sepsis_desc <- desc %>% filter(label == "Sepsis")
mean(sepsis_desc$stay)
median(sepsis_desc$stay)
min(sepsis_desc$stay)
max(sepsis_desc$stay)

nosepsis_desc <- desc %>% filter(label == "NoSepsis")
mean(nosepsis_desc$stay)
median(nosepsis_desc$stay)
min(nosepsis_desc$stay)
max(nosepsis_desc$stay)


sd(nosepsis_desc$stay)
sd(sepsis_desc$stay)
leveneTest(as.numeric(desc$stay), desc$label)
t.test(
  as.numeric(nosepsis_desc$stay),
  as.numeric(sepsis_desc$stay),
  var.equal = FALSE,
  alternative = "two.sided"
)
cohensD(as.numeric(nosepsis_desc$stay),
        as.numeric(sepsis_desc$stay))


age_pre <- arrange(age_pre, age)
age_pre <- arrange(age_pre, desc(age))

sepsis_age <- age_pre %>% filter(treatment == "Sepsis")
mean(sepsis_age$age)
nosepsis_age <- age_pre %>% filter(treatment == "NoSepsis")
mean(nosepsis_age$age)

entries <- as.data.frame(table(final$id))[, 2]
summary(entries)

# Variables used for composite scoring ------------------------------------

# Please note:
# Before generating the correlation matrix include every variable in hRELSA calculation

# prep_final <- (final %>% select(all_of(vars)))
# prep_final <- prep_final %>% select("hr", "pulse", "sao2", "rr", "systolicbp", "map", "temperature")
# export <- round(cor(prep_final, method = "pearson", use = "complete.obs"),2)
# write.csv(export, file = "figs/table2.csv")

# The highest severity values of each variable and what this means --------

analysis

# The Patient with the highest severity value  ---------------------------

analysis # To find out the patient with the highest severity

table(final$id) # check entries

which_id <- "101137"
which_time <- "498480"

dat %>% filter(id == which_id, time == which_time)
age_pre %>% filter(id == which_id, time == which_time)
final  %>% filter(id == which_id, time == which_time) %>% select("sirs", "proven_sepsis", "sus_sepsis", "infection", "ventilation", "medication")

xx <- final %>% filter(id == which_id)
arrange(xx, rms$rms)

# Create Plot with Use Case Creator at the bottom of this script

# The Patient with the lowest average disease severity  ------------------

find_low <- final %>%
  group_by(id, treatment) %>%
  summarise(mean = mean(as.numeric(unlist(rms)), na.rm = TRUE))

arrange(find_low, mean)

table(final$id) # check entries

# Create Plot with Use Case Creator at the bottom of this script

# Does sex influence the disease severity?  -------------------------------

curve <- final %>%
  group_by(id, condition) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE)) %>%
  mutate(condition = factor(
    condition,
    levels = c("f", "m"),
    labels = c("female", "male")
  ))

p4 <- curve %>%
  ggplot(aes(x = condition, y = max, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    comparisons = list(c("female", "male")),
    map_signif_level = TRUE,
    color = c("#000000"),
    textsize = 8,
    annotation = "ns."
  ) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))
p4

ggsave(
  "figs/figure4.tiff",
  plot = p4,
  dpi = 300,
  compression = "lzw"
)

# Statistics
describeBy(curve$max, curve$condition)
leveneTest(curve$max, curve$condition)  # Same VAR
z.test(
  x = curve$max[curve$condition == "female"],
  y = curve$max[curve$condition == "male"],
  sigma.x = sd(curve$max[curve$condition == "female"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$condition == "male"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ condition, data = curve)

# Comparison of different sepsis labels in terms of highest severity --------

# Proven Sepsis ###
proven_sepsis_ids <- unique(final$id[final$proven_sepsis == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoProvenSepsis"
curve$treatment[curve$id %in% proven_sepsis_ids] <- "ProvenSepsis"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoProvenSepsis", "ProvenSepsis"),
                     labels = c("No Proven Sepsis", "Proven Sepsis"))
)

curve <- curve %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

# Figure
px1 <- curve %>%
  ggplot(aes(x = treatment, y = max, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    comparisons = list(c("No Proven Sepsis", "Proven Sepsis")),
    map_signif_level = TRUE,
    color = c("#000000"),
    textsize = 5,
    annotation = "****"
  ) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))
px1


# Statistics
describeBy(curve$max, curve$treatment)
leveneTest(curve$max, curve$treatment)  # Same VAR
z.test(
  x = curve$max[curve$treatment == "No Proven Sepsis"],
  y = curve$max[curve$treatment == "Proven Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No Proven Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Proven Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

# SIRS ###
sirs_ids <- unique(final$id[final$sirs == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoSIRS"
curve$treatment[curve$id %in% sirs_ids] <- "SIRS"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoSIRS", "SIRS"),
                     labels = c("No SIRS", "SIRS"))
)

curve <- curve %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

# Figure
px2 <- curve %>%
  ggplot(aes(x = treatment, y = max, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    comparisons = list(c("No SIRS", "SIRS")),
    map_signif_level = TRUE,
    color = c("#000000"),
    textsize = 5,
    annotation = "****"
  ) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))
px2


# Statistics
describeBy(curve$max, curve$treatment)
leveneTest(curve$max, curve$treatment)  # Same VAR
z.test(
  x = curve$max[curve$treatment == "No SIRS"],
  y = curve$max[curve$treatment == "SIRS"],
  sigma.x = sd(curve$max[curve$treatment == "No SIRS"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "SIRS"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

# Sus Sepsis ###
sus_ids <- unique(final$id[final$sus_sepsis == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoSuspectedSepsis"
curve$treatment[curve$id %in% sus_ids] <- "SuspectedSepsis"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoSuspectedSepsis", "SuspectedSepsis"),
                     labels = c("No Suspected Sepsis", "Suspected Sepsis"))
)

curve <- curve %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

# Figure
px3 <- curve %>%
  ggplot(aes(x = treatment, y = max, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    comparisons = list(c("No Suspected Sepsis", "Suspected Sepsis")),
    map_signif_level = TRUE,
    color = c("#000000"),
    textsize = 5,
    annotation = "****"
  ) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))
px3


# Statistics
describeBy(curve$max, curve$treatment)
leveneTest(curve$max, curve$treatment)  # Same VAR
z.test(
  x = curve$max[curve$treatment == "No Suspected Sepsis"],
  y = curve$max[curve$treatment == "Suspected Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No Suspected Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Suspected Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

# Panel Plot ###
labelcomparison <- px2 +  px1 + px3
labelcomparison <- labelcomparison + plot_layout(ncol = 3)
labelcomparison

ggsave(
  "figs/figure2.tiff",
  plot = labelcomparison,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)

# Figure 6
sirs_ids <- unique(final$id[final$sirs == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoSIRS"
curve$treatment[curve$id %in% sirs_ids] <- "SIRS"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoSIRS", "SIRS"),
                     labels = c("No SIRS", "SIRS"))
)

curve <- curve %>%
  group_by(id, treatment)

sorted_IDs <- reorder(curve$id, curve$rms$rms, FUN = max)

p6 <- curve %>%
  ggplot(aes(x = sorted_IDs, y = rms$rms, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs (x = "", y = "PVS") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 18)
  ) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))
p6

ggsave(
  "figs/figure6.tiff",
  plot = p6,
  dpi = 300,
  compression = "lzw"
)

# hRELSA in comparison with the POPS --------

# Calculate the POPS score with oxygen saturation, heart rate, resp. rate and temperature
dat$pops <- 0
for (s in 1:nrow(dat)) {
  if (anyNA(dat[s, ])) {
    pops <- NA
  } else {
    pops <- 0
    if (dat$sao2[s] < 94 & dat$sao2[s] >= 90) {
      pops <- pops + 1
    } else if (dat$sao2[s] < 90) {
      pops <- pops + 2
    }
    
    if ((dat$hr[s] < 110 &
         dat$hr[s] >= 90) | (dat$hr[s] > 160 & dat$hr[s] <= 180)) {
      pops <- pops + 1
    } else if (dat$hr[s] < 90 | dat$hr[s] > 180) {
      pops <- pops + 2
    }
    
    if ((dat$rr[s] < 30 &
         dat$rr[s] >= 25) | (dat$rr[s] > 40 & dat$rr[s] <= 50)) {
      pops <- pops + 1
    } else if (dat$rr[s] < 25 | dat$rr[s] > 50) {
      pops <- pops + 2
    }
    
    if ((dat$temperature[s] < 36 &
         dat$temperature[s] >= 35) |
        (dat$temperature[s] > 37.5 & dat$temperature[s] <= 39)) {
      pops <- pops + 1
    } else if (dat$temperature[s] < 35 | dat$temperature[s] > 39) {
      pops <- pops + 2
    }
  }
  
  dat$pops[s] <- pops
}

# Figure 8
curve0 <- final %>%
  mutate(
    hours = round(time / (60 * 60), digits = 0),
    minutes = round(time / 60, digits = 0),
    days = round(time / ((60 * 60) * 24))
  ) %>%
  group_by(treatment, id) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

curve00 <- dat %>%
  mutate(
    hours = round(time / (60 * 60), digits = 0),
    minutes = round(time / 60, digits = 0),
    days = round(time / ((60 * 60) * 24))
  ) %>%
  group_by(treatment, id) %>%
  summarise(max = max(as.numeric(pops), na.rm = TRUE))

popstest <- lm(curve0$max ~ curve00$max)
summary(popstest)

p8 <-
  ggplot(fortify(popstest), aes(x = curve00$max, y = curve0$max)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", color = "#1F77B4") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylim(0, 1) +
  xlim(0, 8) +
  labs (x = "highest POPS",
        y = expression('hRELSA'[max]),
        color = "") +
  scale_color_manual(values = c("#1F77B4")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18)
  )
p8

ggsave(
  "figs/figure8.tiff",
  plot = p8,
  dpi = 300,
  compression = "lzw"
)

# Mean increase before and after an event ---------------------------------

### Before and after Proven Sepsis (max out of next 4)
proven_sepsis_final <- final[final$id %in% proven_sepsis_ids, ]
proven_sepsis_befaf <- tibble(col1 = factor(), col2 = double(), col3 = double())
names(proven_sepsis_befaf) <- c("id", "before", "after")

proven_sepsis_final <- arrange(proven_sepsis_final, id, time)
false_before <- TRUE
after <- NA
before <- NA
for (i in 1:nrow(proven_sepsis_final)) {
  if(proven_sepsis_final$proven_sepsis[i] == "TRUE") {
    
    if(false_before) {
      after <- max(proven_sepsis_final$rms$rms[i:(i+3)])
      proven_sepsis_befaf <- proven_sepsis_befaf %>% add_row(id = proven_sepsis_final$id[i], before = before, after = after)
    }
    
    false_before <- FALSE
    
  } else {
    false_before <- TRUE
  }
  
  before <- proven_sepsis_final$rms$rms[i]
}

proven_sepsis_befaf$diff <- proven_sepsis_befaf$after - proven_sepsis_befaf$before
mean(proven_sepsis_befaf$diff)

### Before and after SIRS (max out of next 4)
sirs_final <- final[final$id %in% sirs_ids, ]
sirs_befaf <- tibble(col1 = factor(), col2 = double(), col3 = double())
names(sirs_befaf) <- c("id", "before", "after")

sirs_final <- arrange(sirs_final, id, time)
false_before <- TRUE
after <- NA
before <- NA
for (i in 1:nrow(sirs_final)) {
  if(sirs_final$sirs[i] == "TRUE") {
    
    if(false_before) {
      after <- max(sirs_final$rms$rms[i:(i+3)])
      sirs_befaf <- sirs_befaf %>% add_row(id = sirs_final$id[i], before = before, after = after)
    }
    
    false_before <- FALSE
    
  } else {
    false_before <- TRUE
  }
  
  before <- sirs_final$rms$rms[i]
}

sirs_befaf$diff <- sirs_befaf$after - sirs_befaf$before
mean(sirs_befaf$diff)

### Before and after Suspected Sepsis (max out of next 4)
sus_sepsis_final <- final[final$id %in% sus_ids, ]
sus_sepsis_befaf <- tibble(col1 = factor(), col2 = double(), col3 = double())
names(sus_sepsis_befaf) <- c("id", "before", "after")

sus_sepsis_final <- arrange(sus_sepsis_final, id, time)
false_before <- TRUE
after <- NA
before <- NA
for (i in 1:nrow(sus_sepsis_final)) {
  if(sus_sepsis_final$sus_sepsis[i] == "TRUE") {
    
    if(false_before) {
      after <- max(sus_sepsis_final$rms$rms[i:(i+3)])
      sus_sepsis_befaf <- sus_sepsis_befaf %>% add_row(id = sus_sepsis_final$id[i], before = before, after = after)
    }
    
    false_before <- FALSE
    
  } else {
    false_before <- TRUE
  }
  
  before <- sus_sepsis_final$rms$rms[i]
}

sus_sepsis_befaf$diff <- sus_sepsis_befaf$after - sus_sepsis_befaf$before
mean(sus_sepsis_befaf$diff, na.rm = TRUE)

# Comparison between max and before discharge -----------------------------

## SIRS
sirs_ids <- unique(final$id[final$sirs == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoSIRS"
curve$treatment[curve$id %in% sirs_ids] <- "SIRS"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoSIRS", "SIRS"),
                     labels = c("No SIRS", "SIRS"))
)

max_pvs <- aggregate(rms$rms ~ id, curve, max)
last_pvs <- aggregate(rms$rms ~ id, curve, tail, n = 1)

curve0 <- merge(max_pvs, last_pvs, by = "id")
curve0$treatment <- curve$treatment[match(curve0$id, curve$id)]
colnames(curve0) <- c("id", "max_pvs", "last_pvs", "treatment")

curve00 <- data.frame(
  id = rep(curve0$id, times = 2),
  pvs = c(curve0$max_pvs, curve0$last_pvs),
  type = rep(c("max", "last"), each = nrow(curve0)),
  treatment = rep(curve0$treatment, times = 2)
)
curve00 <- as_tibble(curve00)
curve00$type <- factor(curve00$type, levels = c("max", "last"))

# Figure
pm1 <- curve00 %>%
  ggplot(aes(x = treatment, y = pvs, color = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # geom_signif(
  #   comparisons = list(c("max", "last")),
  #   map_signif_level = TRUE,
  #   color = c("#000000"),
  #   textsize = 5,
  #   annotation = "****"
  # ) +
  labs (x = "", y = "PVS", color = "PVS type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("max" = "#FF7F0E", "last" = "#1F77B4"))
pm1

### Statistics
#SIRS
curve00_sirs <- curve00[curve00$treatment == "SIRS", ]
describeBy(curve00_sirs$pvs, curve00_sirs$type)
leveneTest(curve00_sirs$pvs, curve00_sirs$type)  # Same VAR
z.test(
  x = curve00_sirs$pvs[curve00_sirs$type == "max"],
  y = curve00_sirs$pvs[curve00_sirs$type == "last"],
  sigma.x = sd(curve00_sirs$pvs[curve00_sirs$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_sirs$pvs[curve00_sirs$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_sirs)

# No SIRS
curve00_nosirs <- curve00[curve00$treatment == "No SIRS", ]
describeBy(curve00_nosirs$pvs, curve00_nosirs$type)
leveneTest(curve00_nosirs$pvs, curve00_nosirs$type)  # Same VAR
z.test(
  x = curve00_nosirs$pvs[curve00_nosirs$type == "max"],
  y = curve00_nosirs$pvs[curve00_nosirs$type == "last"],
  sigma.x = sd(curve00_nosirs$pvs[curve00_nosirs$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_nosirs$pvs[curve00_nosirs$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_nosirs)

## Proven Sepsis
proven_sepsis_ids <- unique(final$id[final$proven_sepsis == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoProvenSepsis"
curve$treatment[curve$id %in% proven_sepsis_ids] <- "ProvenSepsis"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoProvenSepsis", "ProvenSepsis"),
                     labels = c("No Proven Sepsis", "Proven Sepsis"))
)

max_pvs <- aggregate(rms$rms ~ id, curve, max)
last_pvs <- aggregate(rms$rms ~ id, curve, tail, n = 1)

curve0 <- merge(max_pvs, last_pvs, by = "id")
curve0$treatment <- curve$treatment[match(curve0$id, curve$id)]
colnames(curve0) <- c("id", "max_pvs", "last_pvs", "treatment")

curve00 <- data.frame(
  id = rep(curve0$id, times = 2),
  pvs = c(curve0$max_pvs, curve0$last_pvs),
  type = rep(c("max", "last"), each = nrow(curve0)),
  treatment = rep(curve0$treatment, times = 2)
)
curve00 <- as_tibble(curve00)
curve00$type <- factor(curve00$type, levels = c("max", "last"))

# Figure
pm2 <- curve00 %>%
  ggplot(aes(x = treatment, y = pvs, color = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # geom_signif(
  #   comparisons = list(c("max", "last")),
  #   map_signif_level = TRUE,
  #   color = c("#000000"),
  #   textsize = 5,
  #   annotation = "****"
  # ) +
  labs (x = "", y = "PVS", color = "PVS type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("max" = "#FF7F0E", "last" = "#1F77B4"))
pm2

### Statistics
#Proven Spesis
curve00_prov <- curve00[curve00$treatment == "Proven Sepsis", ]
describeBy(curve00_prov$pvs, curve00_prov$type)
leveneTest(curve00_prov$pvs, curve00_prov$type)  # Same VAR
z.test(
  x = curve00_prov$pvs[curve00_prov$type == "max"],
  y = curve00_prov$pvs[curve00_prov$type == "last"],
  sigma.x = sd(curve00_prov$pvs[curve00_prov$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_prov$pvs[curve00_prov$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_prov)

# No Proven Sepsis
curve00_noprov <- curve00[curve00$treatment == "No Proven Sepsis", ]
describeBy(curve00_noprov$pvs, curve00_noprov$type)
leveneTest(curve00_noprov$pvs, curve00_noprov$type)  # No Same VAR
z.test(
  x = curve00_noprov$pvs[curve00_noprov$type == "max"],
  y = curve00_noprov$pvs[curve00_noprov$type == "last"],
  sigma.x = sd(curve00_noprov$pvs[curve00_noprov$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_noprov$pvs[curve00_noprov$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_noprov)

## Suspected Sepsis
sus_ids <- unique(final$id[final$sus_sepsis == TRUE])
curve <- final
curve$treatment <- NULL
curve$treatment <- "NoSuspectedSepsis"
curve$treatment[curve$id %in% sus_ids] <- "SuspectedSepsis"
curve <- curve %>% mutate(
  treatment = factor(treatment,
                     levels = c("NoSuspectedSepsis", "SuspectedSepsis"),
                     labels = c("No Suspected Sepsis", "Suspected Sepsis"))
)

max_pvs <- aggregate(rms$rms ~ id, curve, max)
last_pvs <- aggregate(rms$rms ~ id, curve, tail, n = 1)

curve0 <- merge(max_pvs, last_pvs, by = "id")
curve0$treatment <- curve$treatment[match(curve0$id, curve$id)]
colnames(curve0) <- c("id", "max_pvs", "last_pvs", "treatment")

curve00 <- data.frame(
  id = rep(curve0$id, times = 2),
  pvs = c(curve0$max_pvs, curve0$last_pvs),
  type = rep(c("max", "last"), each = nrow(curve0)),
  treatment = rep(curve0$treatment, times = 2)
)
curve00 <- as_tibble(curve00)
curve00$type <- factor(curve00$type, levels = c("max", "last"))

# Figure
pm3 <- curve00 %>%
  ggplot(aes(x = treatment, y = pvs, color = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # geom_signif(
  #   comparisons = list(c("max", "last")),
  #   map_signif_level = TRUE,
  #   color = c("#000000"),
  #   textsize = 5,
  #   annotation = "****"
  # ) +
  labs (x = "", y = "PVS", color = "PVS type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c("max" = "#FF7F0E", "last" = "#1F77B4"))
pm3

### Statistics
#Suspected Spesis
curve00_sus <- curve00[curve00$treatment == "Suspected Sepsis", ]
describeBy(curve00_sus$pvs, curve00_sus$type)
leveneTest(curve00_sus$pvs, curve00_sus$type)  # Same VAR
z.test(
  x = curve00_sus$pvs[curve00_sus$type == "max"],
  y = curve00_sus$pvs[curve00_sus$type == "last"],
  sigma.x = sd(curve00_sus$pvs[curve00_sus$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_sus$pvs[curve00_sus$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_sus)

# No Suspected Sepsis
curve00_nosus <- curve00[curve00$treatment == "No Suspected Sepsis", ]
describeBy(curve00_nosus$pvs, curve00_nosus$type)
leveneTest(curve00_nosus$pvs, curve00_nosus$type)  # No Same VAR
z.test(
  x = curve00_nosus$pvs[curve00_nosus$type == "max"],
  y = curve00_nosus$pvs[curve00_nosus$type == "last"],
  sigma.x = sd(curve00_nosus$pvs[curve00_nosus$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_nosus$pvs[curve00_nosus$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_nosus)

# Panel Plot ###
maxlast <- pm1 +  pm2 + pm3
maxlast <- maxlast + plot_layout(ncol = 3)
maxlast

ggsave(
  "figs/figure2.tiff",
  plot = maxlast,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)




# Use Case Generator ------------------------------------------------------

#Settings
which_id <- "100419"
plot_names <- c("A", "B", "C", "D", "E", "F")

line_size = 1.5

# Plot
curve0 <- final %>%
  filter(id == which_id) %>%
  mutate(count = row_number())

curve00 <- dat %>%
  filter(id == which_id) %>%
  mutate(count = row_number())

line_size = 1.5

p01 <- ggplot() +
  ggtitle(plot_names[1]) +
  geom_line(
    data = curve0,
    aes(x = count, y = rms$rms, color = "PVS"),
    size = line_size
  ) +
  labs (x = "time", y = "PVS", colour = "") +
  ylim(0, 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_manual(values = c("#FF5251")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
        )

p02 <- ggplot() +
  ggtitle(plot_names[2]) +
  geom_line(
    data = curve00,
    aes(x = count, y = hr),
    size = line_size
  ) +
  labs (x = "time", y = "heart rate") +
  ylim(min(dat$hr, na.rm = TRUE), max(dat$hr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 135, 145), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p03 <- ggplot() +
  ggtitle(plot_names[3]) +
  geom_line(
    data = curve00,
    aes(x = count, y = sao2),
    size = line_size
  ) +
  labs (x = "time", y = "oxygen saturation") +
  ylim(min(dat$sao2, na.rm = TRUE), max(dat$sao2, na.rm = TRUE)) +
  geom_hline(yintercept = 100, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p04 <- ggplot() +
  ggtitle(plot_names[4]) +
  geom_line(
    data = curve00,
    aes(x = count, y = rr),
    size = line_size
  ) +
  labs (x = "time", y = "respiratory rate") +
  ylim(min(dat$rr, na.rm = TRUE), max(dat$rr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 35, 42.5), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p05 <- ggplot() +
  ggtitle(plot_names[5]) +
  geom_line(
    data = curve00,
    aes(x = count, y = map),
    size = line_size
  ) +
  labs (x = "time", y = "mean arterial pressure") +
  ylim(min(dat$map, na.rm = TRUE), max(dat$map, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 70, 55), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p06 <- ggplot() +
  ggtitle(plot_names[6]) +
  geom_line(
    data = curve00,
    aes(x = count, y = temperature),
    size = line_size
  ) +
  labs (x = "time", y = "temperature") +
  ylim(min(dat$temperature, na.rm = TRUE), max(dat$temperature, na.rm = TRUE)) +
  geom_hline(yintercept = 37.3, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

usecase <- p01 +  p02 + p03 + p04 + p05 + p06
usecase <- usecase + plot_layout(ncol = 3)
usecase

# ggsave(
#   "figs/figure2.tiff",
#   plot = usecase,
#   dpi = 300,
#   width = 14,
#   height = 7,
#   units = "in",
#   compression = "lzw"
# )