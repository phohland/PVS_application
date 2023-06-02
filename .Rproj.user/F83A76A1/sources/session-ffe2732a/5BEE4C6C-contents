#'
#' hRELSA:
#' A translational approach to quantify a patient’s disease severity in real-time
#'
#' R file to reproduce the data and figures
#'

# Libraries ---------------------------------------------------------------

# Prior a installation of each package is needed

library(hRELSA)

library(tidyverse)
library(janitor)
library(readxl)

library(ggplot2)
library(patchwork)
library(ggsignif)

library(lme4)
library(lmerTest)
library(psych)
library(lubridate)
library(car)

# Data setup -------------------------------------------------------------

# Blood pressure
raw_bp <- read.csv("data/Blood_pressure_invasive.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL)
raw_bp <- as_tibble(raw_bp)
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
raw_temp <- read.csv("data/Body_temperature.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_temp <- as_tibble(raw_temp)
names(raw_temp)[2] <- "time"
names(raw_temp)[3] <- "temperature"
raw_temp <- raw_temp %>%
  clean_names %>%
  select("pmid", "time", "temperature") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    temperature = as.double(temperature)
    
  )

# Heart rate
raw_hr <- read.csv("data/Heart_rate.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_hr <- as_tibble(raw_hr)
names(raw_hr)[2] <- "time"
names(raw_hr)[3] <- "hr"
raw_hr <- raw_hr %>%
  clean_names %>%
  select("pmid", "time", "hr") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    hr = as.double(hr)
  )

# Patient information
raw_master <- read.csv("data/Patient_master_data.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_master <- as_tibble(raw_master)
raw_master <- raw_master %>%
  clean_names %>%
  mutate(
    pmid = as.factor(pmid),
    date_of_birth = as.Date(date_of_birth),
    sex = as.factor(sex)
  )

# Patient information with sepsis labels
raw_master_label <- read.csv("data/Patient_master_data_label.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL)
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
raw_pulse <- read.csv("data/Pulse.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_pulse <- as_tibble(raw_pulse)
names(raw_pulse)[2] <- "time"
names(raw_pulse)[3] <- "pulse"
raw_pulse <- raw_pulse %>%
  clean_names %>%
  select("pmid", "time", "pulse") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    pulse = as.double(pulse)
  )

# Respiratory Rate
raw_resp <- read.csv("data/Respiratory_rate.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_resp <- as_tibble(raw_resp)
names(raw_resp)[2] <- "time"
names(raw_resp)[3] <- "rr"
raw_resp <- raw_resp %>%
  clean_names %>%
  select("pmid", "time", "rr") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    rr = as.double(rr)
  )

# Oxygen saturation
raw_tom <- read.csv("data/Transcutaneous_oxygen_measurement.csv", sep = ";", header = TRUE, dec = ".", row.names = NULL) 
raw_tom <- as_tibble(raw_tom)
names(raw_tom)[2] <- "time"
names(raw_tom)[3] <- "sao2"
raw_tom <- raw_tom %>%
  clean_names %>%
  select("pmid", "time", "sao2") %>%
  mutate(
    time = as.POSIXct(time, tz = ""),
    sao2 = as.double(sao2)
  )

# Merge all data sets into one final data set
raw <- NULL
raw <- merge(raw_hr, raw_pulse, by = c("pmid", "time"), all.x = TRUE, all.y = TRUE)
raw <- merge(raw, raw_tom, by = c("pmid", "time"), all.x = TRUE, all.y = TRUE)
raw <- merge(raw, raw_resp, by = c("pmid", "time"), all.x = TRUE, all.y = TRUE)
raw <- merge(raw, raw_bp, by = c("pmid", "time"), all.x = TRUE, all.y = TRUE)
raw <- merge(raw, raw_temp, by = c("pmid", "time"), all.x = TRUE, all.y = TRUE)
raw <- as_tibble(raw)

# Delete duplicates
raw <- raw %>% distinct(time, .keep_all = TRUE)

# Drop NAs
raw <- na.omit(raw)

# Rename column names in raw data set
names(raw)[1] <- "id"
names(raw)[2] <- "timepoint"

# Write raw in a csv
write.csv(raw, file = "data/raw.csv")

# Merge patient information and sepsis labels
raw_master_full <- merge(raw_master, raw_master_label, by = c("pmid"), all.x = TRUE, all.y = TRUE)
raw_master_full <- raw_master_full[raw_master_full$pmid %in% raw$id, ]
levels(raw_master_full$label) <- c(levels(raw_master_full$label), "NoLabel")
raw_master_full[is.na(raw_master_full)] <- "NoLabel"
raw_master_full <- as_tibble(raw_master_full)

# Add the sepsis label as "treatment"
raw$treatment <- NA
raw$treatment <- as.factor(raw$treatment)
levels(raw$treatment) <- levels(raw_master_full$label)
for (l in 1:nrow(raw)) {
  label <- raw_master_full$label[raw_master_full$pmid == as.character(raw$id[l])]
  raw$treatment[l] <- as.factor(label)
}

# Add the sex as "condition"
raw$condition <- NA
raw$condition <- as.factor(raw$condition)
levels(raw$condition) <- levels(raw_master_full$sex)
for (l in 1:nrow(raw)) {
  sex <- raw_master_full$sex[raw_master_full$pmid == as.character(raw$id[l])]
  raw$condition[l] <- as.factor(sex)
}

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
raw <- hrelsa_days(raw, format = "timecode", formthis = "timepoint", newdayone = TRUE)

# Deletion of unrealistic entries
raw <- raw %>%
  filter(systolicbp > 0 | is.na(systolicbp)) %>%
  filter(map > 0 | is.na(map)) %>%
  filter(temperature > 30 | is.na(temperature)) %>%
  filter(map < systolicbp)

# Format the data
dat <- hrelsa_format(raw, id = "id", time = "time", treatment = "treatment", condition = "condition", vars = vars, included_realtime = "timepoint")

# Fetch data for maximal severity evaluation
reference_dat <- dat %>% filter(treatment == "NoSepsis")
bsl <- hrelsa_adaptive_baselines(dat, reference_dat, vars = vars, turnvars = turnvars, ambivars = ambivars, realtime = "timepoint",
           dob_dat = raw_master_full, dob_data_id_col = 1, dob_data_dob_col = 3,
           norm_dat = raw_norm, norm_dat_names = names(raw_norm))
pre <- bsl$pre
age_pre <- bsl$age_pre

# Generate final data
final <- hrelsa_final(pre, bsl, drop = dropvars, turnvars = turnvars, ambivars = ambivars, zvars = zvars)

# Generate some analysis
analysis <- hrelsa_analysis(final)

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
t.test(as.numeric(nosepsis_desc$stay), as.numeric(sepsis_desc$stay), var.equal = FALSE, alternative = "two.sided")
cohensD(as.numeric(nosepsis_desc$stay), as.numeric(sepsis_desc$stay))


age_pre <- arrange(age_pre, age)
age_pre <- arrange(age_pre, desc(age))

sepsis_age <- age_pre %>% filter(treatment == "Sepsis")
mean(sepsis_age$age)
nosepsis_age <- age_pre %>% filter(treatment == "NoSepsis")
mean(nosepsis_age$age)

entries <- as.data.frame(table(final$id))[,2]
summary(entries)

# Variables used for composite scoring ------------------------------------

# Please note:
# Before generating the correlation matrix include every variable in hRELSA calculation

prep_final <- (final %>% select(all_of(vars)))
prep_final <- prep_final %>% select("hr", "pulse", "sao2", "rr", "systolicbp", "map", "temperature")
export <- round(cor(prep_final, method = "pearson", use = "complete.obs"),2)
write.csv(export, file = "figs/table2.csv")

# The highest severity values of each variable and what this means --------

analysis

# The Patient with the highest severity value  ---------------------------

analysis # To find out the patient with the highest severity

table(final$id) # check entries

which_id <- "101137"
which_patient <- 70
which_time <- "498480"
time_in_days <- (as.numeric(which_time))/((60*60)*24)

dat %>% filter(id == which_id, time == which_time)
age_pre %>% filter(id == which_id, time == which_time)
final  %>% filter(id == which_id, time == which_time)

xx <- final %>% filter(id == which_id)
arrange(xx, rms$rms)

# Plot
which_patient_id <- as.character(unlist(final[final$id2 == which_patient, 1][1, ]))

summary(dat$hr) # To find out min and max, replace hr with regarding variable

curve0 <- final %>% 
  filter( id2 == which_patient) %>%
  mutate(hours = time/(60*60),
         minutes = time/60,
         days = time/((60*60)*24))

curve00 <- dat %>%
  filter( id == which_patient_id) %>%
  mutate(hours = time/(60*60),
         minutes = time/60,
         days = time/((60*60)*24))

point_size = 1
line_size = 0.5
alpha = 1

p01 <-ggplot() +
  ggtitle("A") +
  geom_point(data = curve0, aes(x = days, y = rms$rms, color = "hRELSA"), size = point_size, alpha = alpha) +
  geom_line(data = curve0, aes(x = days, y = rms$rms, color = "hRELSA"), size = line_size, alpha = alpha) +
  labs (x = "days", y = "hRELSA", colour = "") +
  ylim(0,1) +
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 1,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  scale_color_manual(values = c("#FF5251")) +
  theme_classic() +
  theme(legend.position = "none")

p02 <-ggplot() +
  ggtitle("B") +
  geom_point(data = curve00, aes(x = days, y = hr, color = "heart rate"), size = point_size, alpha = alpha_raw) +
  geom_line(data = curve00, aes(x = days, y = hr, color = "heart rate"), size = line_size, alpha = alpha_raw) +
  labs (x = "days", y = "heart rate (/min)", colour = "") +
  ylim(52,232) + # Min and Max Heart Rate of whole Data set
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 135,linetype = 2) + # Reference value
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none")

p03 <-ggplot() +
  ggtitle("C") +
  geom_point(data = curve00, aes(x = days, y = sao2, color = "oxygen saturation"), size = point_size, alpha = alpha_raw) +
  geom_line(data = curve00, aes(x = days, y = sao2, color = "oxygen saturation"), size = line_size, alpha = alpha_raw) +
  labs (x = "days", y = "oxygen saturation (%)", colour = "") +
  ylim(14,100) + # Min and Max Heart Rate of whole Data set
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 100,linetype = 2) + # Reference value
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none")

p04 <-ggplot() +
  ggtitle("D") +
  geom_point(data = curve00, aes(x = days, y = rr, color = "respiratory rate"), size = point_size, alpha = alpha_raw) +
  geom_line(data = curve00, aes(x = days, y = rr, color = "respiratory rate"), size = line_size, alpha = alpha_raw) +
  labs (x = "days", y = "respiratory rate (/min)", colour = "") +
  ylim(0, 158) + # Min and Max Heart Rate of whole Data set
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 35,linetype = 2) + # Reference value
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none")

p05 <-ggplot() +
  ggtitle("E") +
  geom_point(data = curve00, aes(x = days, y = map, color = "mean arterial blood pressure"), size = point_size, alpha = alpha_raw) +
  geom_line(data = curve00, aes(x = days, y = map, color = "mean arterial blood pressure"), size = line_size, alpha = alpha_raw) +
  labs (x = "days", y = "mean arterial blood pressure (mmHg)", colour = "") +
  ylim(24, 196) + # Min and Max Heart Rate of whole Data set
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 70,linetype = 2) + # Reference value
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none")

p06 <-ggplot() +
  ggtitle("F") +
  geom_point(data = curve00, aes(x = days, y = temperature, color = "temperature"), size = point_size, alpha = alpha_raw) +
  geom_line(data = curve00, aes(x = days, y = temperature, color = "temperature"), size = line_size, alpha = alpha_raw) +
  labs (x = "days", y = "temperature (°C)", colour = "") +
  ylim(31.1,40.8) + # Min and Max Heart Rate of whole Data set
  xlim(time_in_days-2,time_in_days+2) +
  geom_hline(yintercept = 37.3,linetype = 2) + # Reference value
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(legend.position = "none")

plot2 <- p01 +  p02 + p03 + p04 + p05 + p06
plot2 <- plot2 + plot_layout(ncol = 3)
ggsave("figs/figure2.tiff", plot = plot2, dpi = 300, width = 14, height = 7, units = "in", compression = "lzw")


# The Patient with the lowest average disease severity  ------------------

find_low <- final %>% 
  group_by(id, treatment) %>%
  summarise(mean = mean(as.numeric(unlist(rms)), na.rm = TRUE))

arrange(find_low, mean)

table(final$id) # check entries



# Does sex influence the disease severity?  -------------------------------

curve <- final %>% 
  group_by(id, condition) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE)) %>%
  mutate(condition = factor(condition, levels = c("f","m"), labels = c("female", "male")))

p4 <- curve %>%
  ggplot(aes(x = condition, y = max, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 3, alpha = 0.75) +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_signif(comparisons = list(c("female", "male")), 
              map_signif_level = TRUE,
              color = c("#000000"),
              textsize = 5,
              annotation = "ns.") +
  labs (x = "", y = expression('hRELSA'[max])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18, face = "bold")) +
  scale_color_manual(values = c("#1F77B4","#FF7F0E"))
p4

ggsave("figs/figure4.tiff", plot = p4, dpi = 300, compression = "lzw")

# Statistics
describeBy(curve$max, curve$condition)
leveneTest(curve$max, curve$condition)  # Same VAR
t.test(max ~ condition, data = curve, var.equal = TRUE, alternative = "two.sided")

cohensD(max ~ condition, data = curve)

# Comparison of non-septic and septic patients in terms of highest severity --------

curve <- final %>% 
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE)) %>%
  mutate(treatment = factor(treatment, levels = c("NoSepsis","Sepsis"), labels = c("No Sepsis", "Sepsis")))

# Figure 5
p5 <- curve %>%
  ggplot(aes(x = treatment, y = max, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 3, alpha = 0.75) +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_signif(comparisons = list(c("No Sepsis", "Sepsis")), 
              map_signif_level = TRUE,
              color = c("#000000"),
              textsize = 5,
              annotation = "****") +
  labs (x = "", y = expression('hRELSA'[max])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18, face = "bold")) +
  scale_color_manual(values = c("#1F77B4","#FF7F0E"))
p5

ggsave("figs/figure5.tiff", plot = p5, dpi = 300, compression = "lzw")

# Statistics
describeBy(curve$max, curve$treatment)
leveneTest(curve$max, curve$treatment)  # Same VAR
t.test(max ~ treatment, data = curve, var.equal = TRUE, alternative = "two.sided")

levels(curve$treatment) <- c("NoSepsis", "Sepsis", NA)
cohensD(max ~ treatment, data = curve)

# Figure 6
curve <- final %>% 
  group_by(id, treatment)

sorted_IDs <- reorder(final$id, final$rms$rms, FUN = max)

p6 <- curve %>%
  ggplot(aes(x = sorted_IDs, y = rms$rms, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept=1, linetype="dashed") +
  labs (x = "", y = "hRELSA") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 18)) +
  scale_color_manual(values = c("#1F77B4","#FF7F0E"))
p6

ggsave("figs/figure6.tiff", plot = p6, dpi = 300, compression = "lzw")

#Figure 7
curve <- final %>% 
  group_by(id, treatment) %>%
  summarise(mean = mean(as.numeric(unlist(rms)), na.rm = TRUE)) %>%
  mutate(treatment = factor(treatment, levels = c("NoSepsis","Sepsis"), labels = c("No Sepsis", "Sepsis")))

p7 <- curve %>%
  ggplot(aes(x = treatment, y = mean, color = treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 3, alpha = 0.75) +
  geom_hline(yintercept=1, linetype="dashed") +
  geom_signif(comparisons = list(c("No Sepsis", "Sepsis")), 
              map_signif_level = TRUE,
              color = c("#000000"),
              textsize = 5,
              annotation = "*") +
  labs (x = "", y = expression('hRELSA'[mean])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18, face = "bold")) +
  scale_color_manual(values = c("#1F77B4","#FF7F0E"))
p7

ggsave("figs/figure7.tiff", plot = p7, dpi = 300, compression = "lzw")

describeBy(curve$mean, curve$treatment)
leveneTest(curve$mean, curve$treatment)  # Same VAR
t.test(mean ~ treatment, data = curve, var.equal = TRUE, alternative = "two.sided")

levels(curve$treatment) <- c("NoSepsis", "Sepsis", NA)
cohensD(mean ~ treatment, data = curve)

# hRELSA in comparison with the POPS --------

# Calculate the POPS score with oxygen saturation, heart rate, resp. rate and temperature
dat$pops <- 0
for (s in 1:nrow(dat)) {
  pops <- 0
  if (dat$sao2[s] < 94 & dat$sao2[s] >= 90) {
    pops <- pops + 1
  } else if (dat$sao2[s] < 90) {
    pops <- pops + 2
  }
  
  if ((dat$hr[s] < 110 & dat$hr[s] >= 90) | (dat$hr[s] > 160 & dat$hr[s] <= 180)) {
    pops <- pops + 1
  } else if (dat$hr[s] < 90 | dat$hr[s] > 180) {
    pops <- pops + 2
  }
  
  if ((dat$rr[s] < 30 & dat$rr[s] >= 25) | (dat$rr[s] > 40 & dat$rr[s] <= 50)) {
    pops <- pops + 1
  } else if (dat$rr[s] < 25 | dat$rr[s] > 50) {
    pops <- pops + 2
  }
  
  if ((dat$temperature[s] < 36 & dat$temperature[s] >= 35) | (dat$temperature[s] > 37.5 & dat$temperature[s] <= 39)) {
    pops <- pops + 1
  } else if (dat$temperature[s] < 35 | dat$temperature[s] > 39) {
    pops <- pops + 2
  }
  
  dat$pops[s] <- pops
}

# Figure 8
curve0 <- final %>%
  mutate(hours = round(time/(60*60), digits = 0),
         minutes = round(time/60, digits = 0),
         days = round(time/((60*60)*24))) %>%
  group_by(treatment, id) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

curve00 <- dat %>%
  mutate(hours = round(time/(60*60), digits = 0),
         minutes = round(time/60, digits = 0),
         days = round(time/((60*60)*24))) %>%
  group_by(treatment, id) %>%
  summarise(max = max(as.numeric(pops), na.rm = TRUE))

popstest <- lm(curve0$max ~ curve00$max)
summary(popstest)

p8 <- ggplot(fortify(popstest), aes(x = curve00$max, y = curve0$max)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_smooth(method = "lm", color = "#1F77B4") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0,1) +
  xlim(0,8) +
  labs (x = "highest POPS", y = expression('hRELSA'[max]), color = "") +
  scale_color_manual(values = c( "#1F77B4")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18))
p8

ggsave("figs/figure8.tiff", plot = p8, dpi = 300, compression = "lzw")
