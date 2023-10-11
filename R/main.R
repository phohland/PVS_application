#'
#' PVS
#'
#' R file to reproduce the data and figures
#'

# Preparation: Libraries ---------------------------------------------------------------

# Required packages
required_packages <- c(
  "hRELSA",
  "tidyverse",
  "janitor",
  "readxl",
  "ggplot2",
  "patchwork",
  "ggsignif",
  "psych",
  # for describeBy and cor
  "lsr",
  # for cohensD
  "BSDA",
  # for z.test
  "zoo",
  # for rollmean
  "writexl",
  "survival",
  "rstatix", "coin"
  # for wilcoxon effect sizes
)

# Install required packages if not already and load them
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

# Preparation: Data Setup -------------------------------------------------------------

# The R/data_setup.R is used to format the raw data.

# raw data
raw <- read_xlsx("output/raw.xlsx")
double_vars <-
  c("hr",
    "pulse",
    "sao2",
    "rr",
    "systolicbp",
    "diastolicbp",
    "map",
    "temperature")
raw <- raw %>% mutate(timepoint = as.POSIXct(timepoint, tz = ""),
                      condition = as.factor(condition)) %>%
  mutate(across(all_of(double_vars), as.double))

master <- read_xlsx("output/master.xlsx")
master <- master %>% mutate(
  pmid = as.factor(pmid),
  sex = as.factor(sex),
  date_of_birth = as.Date(date_of_birth),
  admission = as.Date(admission),
  discharge = as.Date(discharge),
  label = as.factor(label)
)

# SIRS labels
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

# Adds SIRS TRUE or FALSE to each entry in row
raw$sirs <- FALSE

for (i in 1:nrow(raw)) {
  if (raw$id[i] %in% sirs$pmid) {
    sirs_subset <- sirs[sirs$pmid == raw$id[i], ]
    raw$sirs[i] <-
      any(raw$timepoint[i] >= sirs_subset$start &
            raw$timepoint[i] <= sirs_subset$end)
  }
}

# SIRS label setup (for raw)
sirs_ids <- unique(raw$id[raw$sirs == TRUE])
raw$treatment <- factor("NoSIRS", levels = c("NoSIRS", "SIRS"))
raw$treatment[raw$id %in% sirs_ids] <- factor("SIRS")

# Sepsis label setup (for raw)
raw$proven_sepsis <- FALSE
raw$sus_sepsis <- FALSE
raw$infection <- FALSE

for (i in 1:nrow(raw)) {
  if (raw$id[i] %in% sepsis$pmid) {
    for (o in 1:nrow(sepsis)) {
      if (sepsis$pmid[o] == raw$id[i]) {
        if (raw$timepoint[i] >= sepsis$start[o] & raw$timepoint[i] <= sepsis$end[o]) {
          if (sepsis$sepsis[o] == "Proven Sepsis") {
            raw$proven_sepsis[i] <- TRUE
          } else if (sepsis$sepsis[o] == "Suspected Sepsis") {
            raw$sus_sepsis[i] <- TRUE
          } else if (sepsis$sepsis[o] == "Infection") {
            raw$infection[i] <- TRUE
          }
        }
      }
    }
  }
}

proven_sepsis_ids <- unique(raw$id[raw$proven_sepsis == TRUE])
sus_ids <- unique(raw$id[raw$sus_sepsis == TRUE])

# Preparation: PVS Calculation  -----------------------------------------------------------------

# Select variables
vars <- c("hr", "sao2", "rr", "map", "temperature")
turnvars <- NULL
ambivars <- c("hr", "rr",  "map", "temperature")
zvars <- NULL
dropvars <- NULL

# Get reference values
raw_norm <- read_excel("data/reference_values.xlsx", sheet = 1)
raw_norm <- raw_norm %>% select("age", all_of(vars))

# Generate time column (in seconds)
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

# Rolling means
dat_before_rolling_means <- dat
dat <- dat %>% mutate(across(all_of(vars), ~rollapply(., width = 5, FUN = mean, align = "center", fill = NA)))
dat <- dat %>% mutate(window_id = ceiling(row_number()/5))
dat <- dat %>% group_by(window_id) %>%
  filter(row_number() == 3) %>%
  ungroup()
dat <- dat %>% select(-window_id)

# Output dat.xlsv
write_xlsx(dat,"output/dat.xlsx")

# Fetch data for maximal severity evaluation
reference_dat <- dat %>% filter(treatment == "NoSIRS")

# Baseline calculation (maximum severity evalulation)
bsl <-
  hrelsa_adaptive_baselines(
    dat,
    reference_dat,
    vars = vars,
    turnvars = turnvars,
    ambivars = ambivars,
    realtime = "timepoint",
    dob_dat = raw_master,
    dob_data_id_col = 1,
    dob_data_dob_col = 3,
    norm_dat = raw_norm,
    norm_dat_names = names(raw_norm)
  )
pre <- bsl$pre
age_pre <- bsl$age_pre

# Output pre.xlsx
write_xlsx(pre, "output/pre.xlsx")

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

# Preparation: Sepsis Label Setup --------------------------------------------------------

final <- arrange(final, id, time)
dat <- arrange(dat, id, time)

final$timepoint <- dat$timepoint

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

#Output final.xlsx
output_final <- final
output_final$rms <- final$rms$rms

write_xlsx(output_final, "output/final.xlsx")


# 4.1	Cohort characterization  ---------------------------------------

desc <- master
desc$stay <- desc$discharge - desc$admission
desc$age_at_admission <- desc$admission - desc$date_of_birth
#desc$age_at_discharge <- desc$discharge - desc$date_of_birth
desc$SIRS <- ifelse(desc$pmid %in% sirs_ids, TRUE, FALSE)
desc$SuspectedSepsis <- ifelse(desc$pmid %in% sus_ids, TRUE, FALSE)
desc$ProvenSepsis <- ifelse(desc$pmid %in% proven_sepsis_ids, TRUE, FALSE)
desc <- desc %>% select(-date_of_die, - date_of_birth, -admission, -discharge)

length(unique(desc$pmid))

mean(desc$stay)
mean(desc$age_at_admission)

female_desc <- desc %>% filter(sex == "f")
length(female_desc$pmid)
mean(female_desc$stay)
mean(female_desc$age_at_admission)

male_desc <- desc %>% filter(sex == "m")
length(male_desc$pmid)
mean(male_desc$stay)
mean(male_desc$age_at_admission)

nosirs_desc <- desc %>% filter(SIRS == FALSE)
length(nosirs_desc$pmid)
mean(nosirs_desc$stay)
mean(nosirs_desc$age_at_admission)

sirs_desc <- desc %>% filter(SIRS == TRUE)
length(sirs_desc$pmid)
mean(sirs_desc$stay)
mean(sirs_desc$age_at_admission)

sus_desc <- desc %>% filter(SuspectedSepsis == TRUE)
length(sus_desc$pmid)
mean(sus_desc$stay)
mean(sus_desc$age_at_admission)

prov_desc <- desc %>% filter(ProvenSepsis == TRUE)
length(prov_desc$pmid)
mean(prov_desc$stay)
mean(prov_desc$age_at_admission)

# p-values
shapiro_nosirs <- shapiro.test(as.numeric(nosirs_desc$stay))
shapiro_sirs <- shapiro.test(as.numeric(sirs_desc$stay))

stay_values <- c(as.numeric(sirs_desc$stay), as.numeric(nosirs_desc$stay))
condition_values <- c(rep("SIRS", length(sirs_desc$stay)), rep("No SIRS", length(nosirs_desc$stay)))
data_frame <- data.frame(stay = stay_values, condition = factor(condition_values))

wilcox <- wilcox_test(data_frame, stay ~ condition)
eff <- wilcox_effsize(data_frame, stay ~ condition)

z.test(
  x = nosirs_desc$stay,
  y = sirs_desc$stay,
  sigma.x = sd(nosirs_desc$stay, na.rm = TRUE),
  sigma.y = sd(sirs_desc$stay, na.rm = TRUE),
  alternative = "two.sided"
)


age_pre <- arrange(age_pre, age)
age_pre <- arrange(age_pre, desc(age))

sirs_age <- age_pre %>% filter(treatment == "SIRS")
mean(sirs_age$age)
nosirs_age <- age_pre %>% filter(treatment == "NoSIRS")
mean(nosirs_age$age)

entries <- as.data.frame(table(final$id))[, 2]
summary(entries)

#Labels (First run further calculations to generate the variables)
# proven_sepsis_ids
# sus_ids
#desc %>% filter (!(pmid %in% c(sirs_ids, sus_ids, proven_sepsis_ids)))

#Variables
dat_analysis <- dat #%>% filter(treatment == "SIRS")
summary(dat_analysis$hr)
summary(dat_analysis$sao2)
summary(dat_analysis$rr)
summary(dat_analysis$map)
summary(dat_analysis$temperature)


# 4.2 Variable selection for PVS composite scoring + S2 Between  sex comparison of PVSmax --------

# Please note:
# Before generating the correlation matrix include every variable in PVS calculation

# prep_final <- (final %>% select(all_of(vars)))
# prep_final <- prep_final %>% select("hr", "pulse", "sao2", "rr", "systolicbp", "map", "temperature")
# export <- round(cor(prep_final, method = "pearson", use = "complete.obs"),2)
# write_xlsx(data.frame(export), "figs/cor_table.xlsx")
# 
# cortest <- lm(final$hr ~ final$systolicbp)
# summary(cortest)

# Comparison between Sexes
curve <- final %>%
  group_by(id, condition) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE)) %>%
  mutate(condition = factor(
    condition,
    levels = c("f", "m"),
    labels = c("Female", "Male")
  ))

ps3 <- curve %>%
  ggplot(aes(x = condition, y = max, color = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 1),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    comparisons = list(c("Female", "Male")),
    map_signif_level = TRUE,
    color = c("#000000"),
    textsize = 6,
    annotation = "ns."
  ) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18, face = "bold")
  ) +
  scale_color_manual(values = c( "#B30000", "#001F5C"))
ps3

ggsave(
  "figs/supp_2.tiff",
  plot = ps3,
  dpi = 300,
  compression = "lzw"
)

# Statistics
describeBy(curve$max, curve$condition)
IQR(curve$max[curve$condition == "Female"])
IQR(curve$max[curve$condition == "Male"])

shapiro_female <- shapiro.test(curve$max[curve$condition == "Female"])
shapiro_male <- shapiro.test(curve$max[curve$condition == "Male"])
wilcox <- wilcox_test(as.data.frame(curve), max ~ condition)
eff <- wilcox_effsize(as.data.frame(curve), max ~ condition)

z.test(
  x = curve$max[curve$condition == "female"],
  y = curve$max[curve$condition == "male"],
  sigma.x = sd(curve$max[curve$condition == "female"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$condition == "male"], na.rm = TRUE),
  alternative = "two.sided"
)




# 4.3 Disease severity-related vital sign analysis in SIRS patients -----------

raw$hr_scaled <- scale(raw$hr)
raw$sao2_scaled <- scale(raw$sao2)
raw$rr_scaled <- scale(raw$rr)
raw$map_scaled <- scale(raw$map)
raw$temperature_scaled <- scale(raw$temperature)

# Fit the GLMM model
glm_model <- glmer(sirs ~ hr_scaled + sao2_scaled + rr_scaled + map_scaled + temperature_scaled + (1|id), data = raw, family = binomial(link = "logit"),
                   control = glmerControl(optimizer = "Nelder_Mead"))

# Likelihood ratio test for the significance of random effects
null_model <- glmer(sirs ~ (1| id), family = binomial(link = "logit"), data = raw)

lr_test <- anova(null_model, glm_model)
print(lr_test)

conf_intervals            <- confint(glm_model)
odds_ratio_intervals <- exp(conf_intervals)

# 4.4	Using the PVS in disease severity monitoring ---------------------------

## Patient 1 (101137)
analysis # To find out the patient with the highest severity

table(final$id) # check entries: 179

which_id <- "101137"
which_time <- "1"

dat %>% filter(id == which_id, time == which_time)
age_pre %>% filter(id == which_id, time == which_time)
final  %>% filter(id == which_id, time == which_time) #%>% select("sirs", "sus_sepsis")
final  %>% filter(id == which_id, time >= 498180, time <= 499620)

xx <- final %>% filter(id == which_id)
mean(xx$rms$rms, na.rm = TRUE)
arrange(xx, time)

# Comparison Suspected Sepsis PVS and no suspected Sepsis PVS
describeBy(xx$rms$rms, xx$sus_sepsis)
shapiro_sus <- shapiro.test(xx$rms$rms[xx$sus_sepsis == "TRUE"])
shapiro_nosus <- shapiro.test(xx$rms$rms[xx$sus_sepsis == "FALSE"])

xx$nrms <- xx$rms$rms
wilcox <- wilcox_test(as.data.frame(xx), nrms ~ sus_sepsis)
eff <- wilcox_effsize(as.data.frame(xx), nrms ~ sus_sepsis)

wilcox <- wilcox.test(xx$rms$rms[xx$sus_sepsis == "TRUE"], xx$rms$rms[xx$sus_sepsis == "FALSE"], exact = FALSE)
cohensD(rms$rms ~ sus_sepsis, data = xx)

z.test(
  x = xx$rms$rms[xx$sus_sepsis == "TRUE"],
  y = xx$rms$rms[xx$sus_sepsis == "FALSE"],
  sigma.x = sd(xx$rms$rms[xx$sus_sepsis == "TRUE"], na.rm = TRUE),
  sigma.y = sd(xx$rms$rms[xx$sus_sepsis == "FALSE"], na.rm = TRUE),
  alternative = "two.sided"
)


## Patient 2 (100750)
analysis # To find out the patient with the highest severity

table(final$id) # check entries: 1144

which_id <- "100750"
which_time <- "631020"

dat %>% filter(id == which_id, time == which_time)
age_pre %>% filter(id == which_id, time == which_time)
final  %>% filter(id == which_id, time == which_time) #%>% select("sirs", "sus_sepsis")
final  %>% filter(id == which_id, time >= 498180, time <= 499620)

xx <- final %>% filter(id == which_id)
mean(xx$rms$rms, na.rm = TRUE)
arrange(xx, rms$rms)

# Comparison Suspected Sepsis PVS and no suspected Sepsis PVS
describeBy(xx$rms$rms, xx$proven_sepsis)
shapiro_prov <- shapiro.test(xx$rms$rms[xx$proven_sepsis == "TRUE"])
shapiro_noprov <- shapiro.test(xx$rms$rms[xx$proven_sepsis == "FALSE"])
IQR(xx$rms$rms[xx$proven_sepsis == "TRUE"])
IQR(xx$rms$rms[xx$proven_sepsis == "FALSE"])

xx$nrms <- xx$rms$rms
wilcox <- wilcox_test(as.data.frame(xx), nrms ~ proven_sepsis)
eff <- wilcox_effsize(as.data.frame(xx), nrms ~ proven_sepsis)

wilcox <- wilcox.test(xx$rms$rms[xx$proven_sepsis == "TRUE"], xx$rms$rms[xx$proven_sepsis == "FALSE"], exact = FALSE)
cohensD(rms$rms ~ proven_sepsis, data = xx)

z.test(
  x = xx$rms$rms[xx$proven_sepsis == "TRUE"],
  y = xx$rms$rms[xx$proven_sepsis == "FALSE"],
  sigma.x = sd(xx$rms$rms[xx$proven_sepsis == "TRUE"], na.rm = TRUE),
  sigma.y = sd(xx$rms$rms[xx$proven_sepsis == "FALSE"], na.rm = TRUE),
  alternative = "two.sided"
)

# Slope calculation
filtered_data <- final %>%
  filter(id %in% c("100750")) %>% # 101137, 100750
  group_by(id) %>%
  slice(seq_len(which.max(rms$rms)))

arrange(filtered_data, desc(rms$rms))
arrange(filtered_data, time)

(1.09-0.19)/(166-1) # Patient 1 0.0055
(0.82-0.22)/(945-1) # Patient 2 0.00064

model <- lm(rms$rms ~ time, data = filtered_data)
summary(model)



# 4.5	Larger PVS values corresponded to higher disease severity compared to values at patient discharge --------

# Proven Sepsis ###
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

# Statistics
describeBy(curve$max, curve$treatment)
z.test(
  x = curve$max[curve$treatment == "No Proven Sepsis"],
  y = curve$max[curve$treatment == "Proven Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No Proven Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Proven Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No Proven Sepsis"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "Proven Sepsis"])
wilcox <- wilcox.test(curve$max[curve$treatment == "No Proven Sepsis"], curve$max[curve$treatment == "Proven Sepsis"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve), max ~ treatment)

# SIRS ###
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

# Statistics
describeBy(curve$max, curve$treatment)
z.test(
  x = curve$max[curve$treatment == "No SIRS"],
  y = curve$max[curve$treatment == "SIRS"],
  sigma.x = sd(curve$max[curve$treatment == "No SIRS"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "SIRS"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No SIRS"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "SIRS"])
wilcox <- wilcox.test(curve$max[curve$treatment == "No SIRS"], curve$max[curve$treatment == "SIRS"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve), max ~ treatment)

# Sus Sepsis ###
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

# Statistics
describeBy(curve$max, curve$treatment)
z.test(
  x = curve$max[curve$treatment == "No Suspected Sepsis"],
  y = curve$max[curve$treatment == "Suspected Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No Suspected Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Suspected Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(max ~ treatment, data = curve)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No Suspected Sepsis"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "Suspected Sepsis"])
wilcox <- wilcox.test(curve$max[curve$treatment == "No Suspected Sepsis"], curve$max[curve$treatment == "Suspected Sepsis"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve), max ~ treatment)

###

set.seed(123)
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
  ggtitle("A") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    y_position = c(1.05, 1.05, 1.15), xmin = c(0.8, 1.8, 0.8), xmax = c(1.2, 2.2, 1.8),
    annotation = c("****", "****", "****"), tip_length = 0.03, color = "black", textsize = 6) +
  labs (x = "", y = "PVS", color = "PVS Type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18)
  ) +
  scale_color_manual(values = c("max" = "#B30000", "last" = "#001F5C"))

pm1

### Statistics
#SIRS
curve00_sirs <- curve00[curve00$treatment == "SIRS", ]
describeBy(curve00_sirs$pvs, curve00_sirs$type)

z.test(
  x = curve00_sirs$pvs[curve00_sirs$type == "max"],
  y = curve00_sirs$pvs[curve00_sirs$type == "last"],
  sigma.x = sd(curve00_sirs$pvs[curve00_sirs$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_sirs$pvs[curve00_sirs$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)
cohensD(pvs ~ type, data = curve00_sirs)

shapiro_max <- shapiro.test(curve00_sirs$pvs[curve00_sirs$type == "max"])
shapiro_last <- shapiro.test(curve00_sirs$pvs[curve00_sirs$type == "last"])
wilcox <- wilcox.test(curve00_sirs$pvs[curve00_sirs$type == "max"], curve00_sirs$pvs[curve00_sirs$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_sirs), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_sirs), pvs ~ type)

# No SIRS
curve00_nosirs <- curve00[curve00$treatment == "No SIRS", ]
describeBy(curve00_nosirs$pvs, curve00_nosirs$type)
z.test(
  x = curve00_nosirs$pvs[curve00_nosirs$type == "max"],
  y = curve00_nosirs$pvs[curve00_nosirs$type == "last"],
  sigma.x = sd(curve00_nosirs$pvs[curve00_nosirs$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_nosirs$pvs[curve00_nosirs$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_nosirs)

shapiro_max <- shapiro.test(curve00_nosirs$pvs[curve00_nosirs$type == "max"])
shapiro_last <- shapiro.test(curve00_nosirs$pvs[curve00_nosirs$type == "last"])
wilcox <- wilcox.test(curve00_nosirs$pvs[curve00_nosirs$type == "max"], curve00_nosirs$pvs[curve00_nosirs$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_nosirs), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_nosirs), pvs ~ type)

# Last comparison 
curve00_last <- curve00[curve00$type == "last", ]
describeBy(curve00_last$pvs, curve00_last$treatment)
z.test(
  x = curve00_last$pvs[curve00_last$treatment == "SIRS"],
  y = curve00_last$pvs[curve00_last$treatment == "No SIRS"],
  sigma.x = sd(curve00_last$pvs[curve00_last$treatment == "SIRS"], na.rm = TRUE),
  sigma.y = sd(curve00_last$pvs[curve00_last$treatment == "No SIRS"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ treatment, data = curve00_last)

shapiro_with <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "SIRS"])
shapiro_without <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "No SIRS"])
wilcox <- wilcox.test(curve00_last$pvs[curve00_last$treatment == "SIRS"], curve00_last$pvs[curve00_last$treatment == "No SIRS"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_last), pvs ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve00_last), pvs ~ treatment)

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
  ggtitle("C") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    y_position = c(1.05, 1.05, 1.15), xmin = c(0.8, 1.8, 0.8), xmax = c(1.2, 2.2, 1.8),
    annotation = c("****", "**", "**"), tip_length = 0.03, color = "black", textsize = 6) +
  labs (x = "", y = "PVS", color = "PVS Type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18)
  ) +
  scale_color_manual(values = c("max" = "#B30000", "last" = "#001F5C"))
pm2

### Statistics
#Proven Spesis
curve00_prov <- curve00[curve00$treatment == "Proven Sepsis", ]
describeBy(curve00_prov$pvs, curve00_prov$type)
z.test(
  x = curve00_prov$pvs[curve00_prov$type == "max"],
  y = curve00_prov$pvs[curve00_prov$type == "last"],
  sigma.x = sd(curve00_prov$pvs[curve00_prov$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_prov$pvs[curve00_prov$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_prov)

shapiro_max <- shapiro.test(curve00_prov$pvs[curve00_prov$type == "max"])
shapiro_last <- shapiro.test(curve00_prov$pvs[curve00_prov$type == "last"])
wilcox <- wilcox.test(curve00_prov$pvs[curve00_prov$type == "max"], curve00_prov$pvs[curve00_prov$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_prov), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_prov), pvs ~ type)

# No Proven Sepsis
curve00_noprov <- curve00[curve00$treatment == "No Proven Sepsis", ]
describeBy(curve00_noprov$pvs, curve00_noprov$type)
z.test(
  x = curve00_noprov$pvs[curve00_noprov$type == "max"],
  y = curve00_noprov$pvs[curve00_noprov$type == "last"],
  sigma.x = sd(curve00_noprov$pvs[curve00_noprov$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_noprov$pvs[curve00_noprov$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_noprov)

shapiro_max <- shapiro.test(curve00_noprov$pvs[curve00_noprov$type == "max"])
shapiro_last <- shapiro.test(curve00_noprov$pvs[curve00_noprov$type == "last"])
wilcox <- wilcox.test(curve00_noprov$pvs[curve00_noprov$type == "max"], curve00_noprov$pvs[curve00_noprov$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_noprov), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_noprov), pvs ~ type)

# Last comparison 
curve00_last <- curve00[curve00$type == "last", ]
describeBy(curve00_last$pvs, curve00_last$treatment)
z.test(
  x = curve00_last$pvs[curve00_last$treatment == "Proven Sepsis"],
  y = curve00_last$pvs[curve00_last$treatment == "No Proven Sepsis"],
  sigma.x = sd(curve00_last$pvs[curve00_last$treatment == "Proven Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve00_last$pvs[curve00_last$treatment == "No Proven Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ treatment, data = curve00_last)

shapiro_with <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "Proven Sepsis"])
shapiro_without <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "No Proven Sepsis"])
wilcox <- wilcox.test(curve00_last$pvs[curve00_last$treatment == "Proven Sepsis"], curve00_last$pvs[curve00_last$treatment == "No Proven Sepsis"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_last), pvs ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve00_last), pvs ~ treatment)

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
  ggtitle("B") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.35),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    y_position = c(1.05, 1.05, 1.15), xmin = c(0.8, 1.8, 0.8), xmax = c(1.2, 2.2, 1.8),
    annotation = c("****", "****", "***"), tip_length = 0.03, color = "black", textsize = 6) +
  labs (x = "", y = "PVS", color = "PVS Type") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18)
  ) +
  scale_color_manual(values = c("max" = "#B30000", "last" = "#001F5C"))
pm3

### Statistics
#Suspected Spesis
curve00_sus <- curve00[curve00$treatment == "Suspected Sepsis", ]
describeBy(curve00_sus$pvs, curve00_sus$type)
z.test(
  x = curve00_sus$pvs[curve00_sus$type == "max"],
  y = curve00_sus$pvs[curve00_sus$type == "last"],
  sigma.x = sd(curve00_sus$pvs[curve00_sus$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_sus$pvs[curve00_sus$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_sus)

shapiro_max <- shapiro.test(curve00_sus$pvs[curve00_sus$type == "max"])
shapiro_last <- shapiro.test(curve00_sus$pvs[curve00_sus$type == "last"])
wilcox <- wilcox.test(curve00_sus$pvs[curve00_sus$type == "max"], curve00_sus$pvs[curve00_sus$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_sus), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_sus), pvs ~ type)

# No Suspected Sepsis
curve00_nosus <- curve00[curve00$treatment == "No Suspected Sepsis", ]
describeBy(curve00_nosus$pvs, curve00_nosus$type)
z.test(
  x = curve00_nosus$pvs[curve00_nosus$type == "max"],
  y = curve00_nosus$pvs[curve00_nosus$type == "last"],
  sigma.x = sd(curve00_nosus$pvs[curve00_nosus$type == "max"], na.rm = TRUE),
  sigma.y = sd(curve00_nosus$pvs[curve00_nosus$type == "last"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ type, data = curve00_nosus)

shapiro_max <- shapiro.test(curve00_nosus$pvs[curve00_nosus$type == "max"])
shapiro_last <- shapiro.test(curve00_nosus$pvs[curve00_nosus$type == "last"])
wilcox <- wilcox.test(curve00_nosus$pvs[curve00_nosus$type == "max"], curve00_nosus$pvs[curve00_nosus$type == "last"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_nosus), pvs ~ type)
eff <- wilcox_effsize(as.data.frame(curve00_nosus), pvs ~ type)

# Last comparison 
curve00_last <- curve00[curve00$type == "last", ]
describeBy(curve00_last$pvs, curve00_last$treatment)
z.test(
  x = curve00_last$pvs[curve00_last$treatment == "Suspected Sepsis"],
  y = curve00_last$pvs[curve00_last$treatment == "No Suspected Sepsis"],
  sigma.x = sd(curve00_last$pvs[curve00_last$treatment == "Suspected Sepsis"], na.rm = TRUE),
  sigma.y = sd(curve00_last$pvs[curve00_last$treatment == "No Suspected Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)

cohensD(pvs ~ treatment, data = curve00_last)

shapiro_with <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "Suspected Sepsis"])
shapiro_without <- shapiro.test(curve00_last$pvs[curve00_last$treatment == "No Suspected Sepsis"])
wilcox <- wilcox.test(curve00_last$pvs[curve00_last$treatment == "Suspected Sepsis"], curve00_last$pvs[curve00_last$treatment == "No Suspected Sepsis"], exact = FALSE)

wilcox <- wilcox_test(as.data.frame(curve00_last), pvs ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve00_last), pvs ~ treatment)

# Panel Plot ###
maxlast <- pm1 +  pm3 + pm2
maxlast <- maxlast + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "top", plot.title = element_text(face = "bold"))
maxlast

ggsave(
  "figs/figure_2.tiff",
  plot = maxlast,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)


# 4.6 Multidimensional PVS-related septic states were different from non-SIRS patients --------



#Figure No SIRS and sepsis groups
curve <- final %>% filter(treatment == "NoSIRS") %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))

curve_sirs <- final %>% filter(treatment == "SIRS" & !(id %in% sus_ids) & !(id %in% proven_sepsis_ids)) %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))
curve_sirs$treatment <- as.factor("SIRS")  

curve_sus <- final %>% filter(id %in% sus_ids & !(id %in% proven_sepsis_ids)) %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))
curve_sus$treatment <- as.factor("Suspected Sepsis")

curve_prov <- final %>% filter(id %in% proven_sepsis_ids) %>%
  group_by(id, treatment) %>%
  summarise(max = max(as.numeric(unlist(rms)), na.rm = TRUE))
curve_prov$treatment <- as.factor("Proven Sepsis")

curve <- rbind(curve, curve_sirs, curve_sus, curve_prov)
curve <- curve %>% mutate(treatment = factor(treatment,
                                             levels = c("NoSIRS", "SIRS", "Suspected Sepsis", "Proven Sepsis"),
                                             labels = c("No SIRS", "SIRS", "Suspected Sepsis", "Proven Sepsis")))

p4 <- curve %>%
  ggplot(aes(x = treatment, y = max, color = treatment)) +
  #ggtitle("C") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 1),
             size = 3,
             alpha = 0.75) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(
    y_position = c(1.05, 1.15, 1.25), xmin = c(1,1,1), xmax = c(2,3,4),
    annotation = c("*", "***", "***"), tip_length = 0.03, color = "black", textsize = 6) +
  labs (x = "", y = expression('PVS'[max])) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 18)
  ) +
  scale_color_manual(values = c( "#001F5C", "#B30000", "#B30000", "#B30000"))
p4

ggsave(
  "figs/figure_3.tiff",
  plot = p4,
  dpi = 300,
  units = "in",
  compression = "lzw"
)

# Statistics NoSIRS SIRS
describeBy(curve$max, curve$treatment)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No SIRS"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "SIRS"])

curve0 <- curve %>% filter (treatment %in% c("No SIRS", "SIRS")) %>%
  mutate( treatment = factor(treatment, levels = c("No SIRS", "SIRS")))
wilcox <- wilcox_test(as.data.frame(curve0), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve0), max ~ treatment)

z.test(
  x = curve$max[curve$treatment == "No SIRS"],
  y = curve$max[curve$treatment == "SIRS"],
  sigma.x = sd(curve$max[curve$treatment == "No SIRS"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "SIRS"], na.rm = TRUE),
  alternative = "two.sided"
)
cohensD(max ~ treatment, data = curve0)

# Statistics NoSIRS Suspected Sepsis
describeBy(curve$max, curve$treatment)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No SIRS"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "Suspected Sepsis"])

curve0 <- curve %>% filter (treatment %in% c("No SIRS", "Suspected Sepsis")) %>%
  mutate( treatment = factor(treatment, levels = c("No SIRS", "Suspected Sepsis")))
wilcox <- wilcox_test(as.data.frame(curve0), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve0), max ~ treatment)

z.test(
  x = curve$max[curve$treatment == "No SIRS"],
  y = curve$max[curve$treatment == "Suspected Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No SIRS"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Suspected Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)
cohensD(max ~ treatment, data = curve0)

# Statistics NoSIRS Proven Sepsis
describeBy(curve$max, curve$treatment)

shapiro_without <- shapiro.test(curve$max[curve$treatment == "No SIRS"])
shapiro_with <- shapiro.test(curve$max[curve$treatment == "Proven Sepsis"])

curve0 <- curve %>% filter (treatment %in% c("No SIRS", "Proven Sepsis")) %>%
  mutate( treatment = factor(treatment, levels = c("No SIRS", "Proven Sepsis")))
wilcox <- wilcox_test(as.data.frame(curve0), max ~ treatment)
eff <- wilcox_effsize(as.data.frame(curve0), max ~ treatment)

z.test(
  x = curve$max[curve$treatment == "No SIRS"],
  y = curve$max[curve$treatment == "Proven Sepsis"],
  sigma.x = sd(curve$max[curve$treatment == "No SIRS"], na.rm = TRUE),
  sigma.y = sd(curve$max[curve$treatment == "Proven Sepsis"], na.rm = TRUE),
  alternative = "two.sided"
)
cohensD(max ~ treatment, data = curve0)

## Additional analyses






# Use Case Generator ------------------------------------------------------

#Settings
which_id <- "100750"
#which_id <- "101137"
#plot_names <- c("G", "H", "I", "J", "K", "L")
plot_names <- c("A", "B", "C", "D", "E", "F")
in_sus <- which_id %in% sus_ids
in_proven_sepsis <- which_id %in% proven_sepsis_ids
in_sirs <- which_id %in% sirs_ids
variable_name <- ifelse(in_sus, "sus_sepsis", ifelse(in_proven_sepsis, "proven_sepsis", ifelse(in_sirs, "sirs", NA)))
#variable_name <- "proven_sepsis"

line_size = 1.25
vline_size = 0.5
hline_size = 0.5

# Plot
curve0 <- final %>%
  filter(id == which_id) %>%
  mutate(count = row_number())

curve00 <- dat %>%
  filter(id == which_id) %>%
  mutate(count = row_number())

p01 <- ggplot() +
  ggtitle(plot_names[1]) +
  geom_line(
    data = curve0,
    aes(x = count, y = rms$rms, color = "PVS"),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "PVS", colour = "") +
  ylim(0, max(final$rms$rms, na.rm = TRUE)) +
  geom_hline(yintercept = 1, linetype = 2, size = hline_size) +
  scale_color_manual(values = c("#B30000")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
        )

p02 <- ggplot() +
  ggtitle(plot_names[2]) +
  geom_line(
    data = curve00,
    aes(x = count, y = hr),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "Heart Rate (bpm)") +
  #labs (x = "Consecutive Entries", y = expression("Heart Rate (min"^-1*")")) +
  ylim(min(dat$hr, na.rm = TRUE), max(dat$hr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 135, 145), linetype = 2, size = hline_size) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

## For ticks
p02 <- p02 + scale_x_continuous(breaks = seq(0, 649260, 100), expand = c(0, 0.2))

p03 <- ggplot() +
  ggtitle(plot_names[3]) +
  geom_line(
    data = curve00,
    aes(x = count, y = sao2),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "Oxygen Saturation (%)") +
  ylim(min(dat$sao2, na.rm = TRUE), max(dat$sao2, na.rm = TRUE)) +
  geom_hline(yintercept = 100, linetype = 2, size = hline_size) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

p04 <- ggplot() +
  ggtitle(plot_names[4]) +
  geom_line(
    data = curve00,
    aes(x = count, y = rr),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "Respiratory Rate (rpm)") +
  #labs (x = "Consecutive Entries", y = expression("Respiratory Rate (min"^-1*")")) +
  ylim(min(dat$rr, na.rm = TRUE), max(dat$rr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 35, 42.5), linetype = 2, size = hline_size) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

p05 <- ggplot() +
  ggtitle(plot_names[5]) +
  geom_line(
    data = curve00,
    aes(x = count, y = map),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "Mean Arterial Pressure (mmHg)") +
  ylim(min(dat$map, na.rm = TRUE), max(dat$map, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 70, 55), linetype = 2, size = hline_size) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

p06 <- ggplot() +
  ggtitle(plot_names[6]) +
  geom_line(
    data = curve00,
    aes(x = count, y = temperature),
    size = line_size
  ) +
  geom_vline(
    data = curve0 %>% filter(
      eval(parse(text = variable_name)) != lag(eval(parse(text = variable_name)))
    ),
    aes(xintercept = count),
    linetype = "dashed",
    size = vline_size
  ) +
  labs (x = "Consecutive Entries", y = "Temperature (°C)") +
  ylim(min(dat$temperature, na.rm = TRUE), max(dat$temperature, na.rm = TRUE)) +
  geom_hline(yintercept = 37.3, linetype = 2, size = hline_size) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank()
  )

usecase <- p01 +  p02 + p03 + p04 + p05 + p06
usecase <- usecase + plot_layout(ncol = 3) & theme(plot.title = element_text(face = "bold"),
                                                   text = element_text(size = 12),
                                                   plot.margin = margin(0.35, 0.35, 0.35, 0.35, "cm"))
usecase

ggsave(
  "figs/XXX.tiff",
  plot = usecase,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)
# Use Case Generator (real time) ------------------------------------------------------

#Settings
which_id <- "101130"
plot_names <- c("A", "B", "C", "D", "E", "F")

line_size = 1.5

# Plot
curve0 <- final %>%
  filter(id == which_id) %>%
  mutate(days = time/60/60/24)

curve00 <- dat %>%
  filter(id == which_id) %>%
  mutate(days = time/60/60/24)

line_size = 1.5

p01 <- ggplot() +
  ggtitle(plot_names[1]) +
  geom_line(
    data = curve0,
    aes(x = days, y = rms$rms, color = "PVS"),
    size = line_size
  ) +
  labs (x = "days", y = "PVS", colour = "") +
  ylim(0, 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_manual(values = c("#FF5251")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p02 <- ggplot() +
  ggtitle(plot_names[2]) +
  geom_line(
    data = curve00,
    aes(x = days, y = hr),
    size = line_size
  ) +
  labs (x = "days", y = "heart rate") +
  ylim(min(dat$hr, na.rm = TRUE), max(dat$hr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 135, 145), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p03 <- ggplot() +
  ggtitle(plot_names[3]) +
  geom_line(
    data = curve00,
    aes(x = days, y = sao2),
    size = line_size
  ) +
  labs (x = "days", y = "oxygen saturation") +
  ylim(min(dat$sao2, na.rm = TRUE), max(dat$sao2, na.rm = TRUE)) +
  geom_hline(yintercept = 100, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p04 <- ggplot() +
  ggtitle(plot_names[4]) +
  geom_line(
    data = curve00,
    aes(x = days, y = rr),
    size = line_size
  ) +
  labs (x = "days", y = "respiratory rate") +
  ylim(min(dat$rr, na.rm = TRUE), max(dat$rr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 35, 42.5), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p05 <- ggplot() +
  ggtitle(plot_names[5]) +
  geom_line(
    data = curve00,
    aes(x = days, y = map),
    size = line_size
  ) +
  labs (x = "days", y = "mean arterial pressure") +
  ylim(min(dat$map, na.rm = TRUE), max(dat$map, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 70, 55), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p06 <- ggplot() +
  ggtitle(plot_names[6]) +
  geom_line(
    data = curve00,
    aes(x = days, y = temperature),
    size = line_size
  ) +
  labs (x = "days", y = "temperature") +
  ylim(min(dat$temperature, na.rm = TRUE), max(dat$temperature, na.rm = TRUE)) +
  geom_hline(yintercept = 37.3, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

usecase <- p01 +  p02 + p03 + p04 + p05 + p06
usecase <- usecase + plot_layout(ncol = 3)
usecase

ggsave(
  "figs/suppmaterial5.tiff",
  plot = usecase,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)
# Use Case Generator (Fitting) ------------------------------------------------------

#Settings
which_id <- "101137"
plot_names <- c("A", "B", "C", "D", "E")

line_size = 1.25

# Plot
curve0 <- dat_before_rolling_means %>%
  filter(id == which_id) %>%
  mutate(count = row_number())

curve00 <- dat %>%
  filter(id == which_id) %>%
  mutate(count = curve0$count[match(paste(id, time), paste(curve0$id, curve0$time))])

p01 <- ggplot() +
  ggtitle(plot_names[1]) +
  geom_line(
    data = curve0,
    aes(x = count, y = hr),
    size = line_size,
    color = "black"
  ) +
  geom_line(
    data = curve00,
    aes(x = count, y = hr),
    size = line_size,
    color = "red"
  ) +
  labs (x = "Consecutive Entries", y = "Heart Rate (bpm)") +
  ylim(min(dat_before_rolling_means$hr, na.rm = TRUE), max(dat_before_rolling_means$hr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 135, 145), linetype = 2) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p02 <- ggplot() +
  ggtitle(plot_names[2]) +
  geom_line(
    data = curve0,
    aes(x = count, y = sao2),
    size = line_size,
    color = "black"
  ) +
  geom_line(
    data = curve00,
    aes(x = count, y = sao2),
    size = line_size,
    color = "red"
  ) +
  labs (x = "Consecutive Entries", y = "Oxygen Saturation (%)") +
  ylim(min(dat_before_rolling_means$sao2, na.rm = TRUE), max(dat_before_rolling_means$sao2, na.rm = TRUE)) +
  geom_hline(yintercept = 100, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p03 <- ggplot() +
  ggtitle(plot_names[3]) +
  geom_line(
    data = curve0,
    aes(x = count, y = rr),
    size = line_size,
    color = "black"
  ) +
  geom_line(
    data = curve00,
    aes(x = count, y = rr),
    size = line_size,
    color = "red"
  ) +
  labs (x = "Consecutive Entries", y = "Respiratory Rate (rpm)") +
  ylim(min(dat_before_rolling_means$rr, na.rm = TRUE), max(dat_before_rolling_means$rr, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 35, 42.5), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p04 <- ggplot() +
  ggtitle(plot_names[4]) +
  geom_line(
    data = curve0,
    aes(x = count, y = map),
    size = line_size,
    color = "black"
  ) +
  geom_line(
    data = curve00,
    aes(x = count, y = map),
    size = line_size,
    color = "red"
  ) +
  labs (x = "Consecutive Entries", y = "Mean Arterial Pressure (mmHg)") +
  ylim(min(dat_before_rolling_means$map, na.rm = TRUE), max(dat_before_rolling_means$map, na.rm = TRUE)) +
  geom_hline(yintercept = ifelse(mean(age_pre$age[age_pre$id == which_id]) > 0.547, 70, 55), linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

p05 <- ggplot() +
  ggtitle(plot_names[5]) +
  geom_line(
    data = curve0,
    aes(x = count, y = temperature),
    size = line_size,
    color = "black"
  ) +
  geom_line(
    data = curve00,
    aes(x = count, y = temperature),
    size = line_size,
    color = "red"
  ) +
  labs (x = "Consecutive Entries", y = "Temperature (°C)") +
  ylim(min(dat_before_rolling_means$temperature, na.rm = TRUE), max(dat_before_rolling_means$temperature, na.rm = TRUE)) +
  geom_hline(yintercept = 37.3, linetype = 2) +
  scale_color_manual(values = c("black")) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

usecase <- p01 + p02 + p03 + p04 + p05
usecase <- usecase + plot_layout(ncol = 3) & theme(plot.title = element_text(face = "bold"),
                                                   text = element_text(size = 12))
usecase

ggsave(
  "figs/supp_1.tiff",
  plot = usecase,
  dpi = 300,
  width = 14,
  height = 7,
  units = "in",
  compression = "lzw"
)


