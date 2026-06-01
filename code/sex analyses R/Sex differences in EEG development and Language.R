# EEG SEX DIFFERENCES IN INFANT DEVELOPMENT
#
## Author: Gabriel Blanco-Gomez

# Research Questions:

#   EEG SEX DIFFERENCES IN INFANT DEVELOPMENT
#   Q1,1: Do males and females differ in EEG development across various EEG measures (6 to 12 months)?
#   Q1.2: Does sex interact with ASD likelihood (TLA vs ELA)?
#   Q1.3 (Exploratory): Does sex interact with ASD diagnosis? (ELA only)

#   RELATIONSHIP BETWEEN EEG AND LANGUAGE

#   Q2.1: Does expressive language growth moderated by differences in speech connectivity?
#   Q2.2: Does expressive language growth moderated by differences in speech connectivity?


#   MULTIVARIATE EEG AND LANGUAGE

#   Q3.1: Does expressive language growth differ by sex and HC neurosubtype?
#   Q3.2: Does receptive language growth differ by sex and HC neurosubtype?


# SETUP

library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(effectsize)
library(flextable)
library(officer)
library(lmtest)
library(car)
library(sandwich)
library(interactions)
library(patchwork)         

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dir.create("../../figures/main",          recursive = TRUE, showWarnings = FALSE)
dir.create("../../figures/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/main",                 recursive = TRUE, showWarnings = FALSE)
dir.create("tables/supplementary",        recursive = TRUE, showWarnings = FALSE)


# DATA PREPARATION

raw <- read.csv(
  "../../datasets/sex_diff_main.csv",
  header = TRUE, sep = ","
)

# Check structure
glimpse(raw)
cat("N subjects:", nrow(raw), "\n")
cat("Columns:",    ncol(raw), "\n")

# Rename speech columns to drop the _left suffix so the pivot is clean
raw <- raw %>%
  rename(
    speech_con_6  = speech_con_6_left,
    speech_con_12 = speech_con_12_left
  )

# Wide -> Long: one row per subject x timepoint, 4 EEG columns
eeg_long <- raw %>%
  select(subject, sex, site, group_type, outcome,nonverbal_iq_6,
         front_gamma_6, front_gamma_12,
         auditory_con_6,  auditory_con_12,
         speech_con_6,    speech_con_12,
         gamma_lat_6,     gamma_lat_12) %>%
  pivot_longer(
    cols            = c(front_gamma_6,  front_gamma_12,
                        auditory_con_6,   auditory_con_12,
                        speech_con_6,     speech_con_12,
                        gamma_lat_6,      gamma_lat_12),
    names_to        = c(".value", "timepoint"),
    names_pattern   = "(.+)_(6|12)$"
  ) %>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6,          # centre: 0 = 6 mo, 6 = 12 mo
    sex         = factor(sex,        levels = c("M", "F")),
    site        = factor(site),
    iq          = as.numeric(nonverbal_iq_6), # for supplementary checks)
    group_type  = factor(group_type, levels = c("TLA", "ELA")),
    subject     = factor(subject)
  )

# Sample sizes per group
eeg_long %>%
  distinct(subject, sex, group_type) %>%
  count(sex, group_type)

# Missing data per EEG variable at each timepoint
eeg_long %>%
  group_by(timepoint) %>%
  summarise(
    n_front_gamma = sum(!is.na(front_gamma)),
    n_auditory_con  = sum(!is.na(auditory_con)),
    n_speech_con    = sum(!is.na(speech_con)),
    n_gamma_lat     = sum(!is.na(gamma_lat))
  )


# Q1: SEX x TIME: FRONTAL GAMMA POWER

q1_fg <- lmer(
  front_gamma ~ timepoint_c * sex + site + iq + (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1_fg)
performance::r2(q1_fg)
performance::check_model(q1_fg)

# Preliminary answer = no difference in frontal gamma power development


# Q1: SEX x TIME: AUDITORY NETWORK CONNECTIVITY

q1_ac <- lmer(
  auditory_con ~ timepoint_c * sex + site + iq + (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1_ac)
performance::r2(q1_ac)
performance::check_model(q1_ac)

# Preliminary results = no difference in auditory connectivity

# Check Homogenity of variance 
bptest(lm(auditory_con ~ timepoint_c * sex + site, data = eeg_long))

# HC3 robustness check: compare with original
q1_ac_lm <- lm(auditory_con ~ timepoint_c * sex + site, data = eeg_long)
coeftest(q1_ac_lm, vcov = vcovHC(q1_ac_lm, type = "HC3"))

# Q1: SEX x TIME: SPEECH NETWORK CONNECTIVITY

q1_sc <- lmer(
  speech_con ~ timepoint_c * sex + site + iq + (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1_sc)
performance::r2(q1_sc)
performance::check_model(q1_sc)

# Issue with singularity (random intercept variance estimated at 0)
# Rerun with singular model

q1_sc_lm <- lm(
  speech_con ~ timepoint_c * sex + site +iq ,
  data = eeg_long
)

summary(q1_sc_lm)
performance::r2(q1_sc_lm)
performance::check_model(q1_sc_lm)


#Preliminary results = sig difference in speech connectivity

# Q1: SEX x TIME: POWER LATERALIZATION

q1_gl <- lmer(
  gamma_lat ~ timepoint_c * sex + site + iq +(1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1_gl)
performance::r2(q1_gl)
performance::check_model(q1_gl)

# Singularity issue again, rerun without random effect
q1_gl_lm <- lm(
  gamma_lat ~ timepoint_c * sex + site + iq,
  data = eeg_long
)
summary(q1_gl_lm)
performance::r2(q1_gl_lm)
performance::check_model(q1_gl_lm)


#Q1 compare all four models  
performance::compare_performance(q1.3_fg, q1.3_ac, q1.3_sc, q1.3_gl)

# Q1.2: SEX x LIKELIHOOD x TIME: FRONTAL GAMMA POWER

q1.2_fg <- lmer(
  front_gamma ~ timepoint_c * sex * group_type + site + iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.2_fg)
performance::r2(q1.2_fg)
performance::check_model(q1.2_fg)


# Q1.2: SEX x LIKELIHOOD x TIME: AUDITORY NETWORK CONNECTIVITY

q1.2_ac <- lmer(
  auditory_con ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.2_ac)
performance::r2(q1.2_ac)
performance::check_model(q1.2_ac)

# Q1.2: SEX x LIKELIHOOD x TIME: SPEECH NETWORK CONNECTIVITY

q1.2_sc <- lmer(
  speech_con ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.2_sc)
performance::r2(q1.2_sc)
performance::check_model(q1.2_sc)

# Singularity issue again, rerun without random effect
q1.2_sc_lm <- lm(
  speech_con ~ timepoint_c * sex * group_type + site + iq,
  data = eeg_long
)

summary(q1.2_sc_lm)
performance::r2(q1.2_sc_lm)
performance::check_model(q1.2_sc_lm)
emmeans(q1.2_sc_lm, pairwise ~ sex | group_type)

# Q1.2: SEX x LIKELIHOOD x TIME: POWER LATERALIZATION

q1.2_gl <- lmer(
  gamma_lat ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.2_gl)
performance::r2(q1.2_gl)
performance::check_model(q1.2_gl)

# Singularity issue again, rerun without random effect
q1.2_gl_lm <- lm(
  gamma_lat ~ timepoint_c * sex * group_type + site +iq,
  data = eeg_long
)


summary(q1.2_gl_lm)
performance::r2(q1.2_gl_lm)
performance::check_model(q1.2_gl_lm)


# Q1.2: compare all four models
performance::compare_performance(q1.3_fg, q1.3_ac, q1.3_sc, q1.3_gl)


# Q1.3 PREP: ELA SUBSAMPLE + ASD OUTCOME CODING

eeg_ela <- eeg_long %>%
  filter(group_type == "ELA") %>%
  mutate(
    outcome_bin = factor(
      case_when(
        str_detect(tolower(outcome), "^asd")    ~ "ASD",
        str_detect(tolower(outcome), "no.?asd") ~ "no-ASD",
        TRUE                                     ~ NA_character_
      ),
      levels = c("no-ASD", "ASD")
    )
  ) %>%
  filter(!is.na(outcome_bin))

# N per cell: check before modelling
eeg_ela %>%
  distinct(subject, sex, outcome_bin) %>%
  count(sex, outcome_bin)


# Q1.3 EXPLORATORY: SEX x DIAGNOSIS x TIME: FRONTAL GAMMA POWER

q1.3_fg <- lmer(
  front_gamma ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.3_fg)
performance::r2(q1.3_fg)
performance::check_model(q1.3_fg)


# Q1.3 EXPLORATORY: SEX x DIAGNOSIS x TIME: AUDITORY NETWORK


q1.3_ac <- lmer(
  auditory_con ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.3_ac)
performance::r2(q1.3_ac)
performance::check_model(q1.3_ac)


 
# Q1.3: SEX x DIAGNOSIS x TIME: SPEECH NETWORK
 
q1.3_sc <- lmer(
  speech_con ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.3_sc)
performance::r2(q1.3_sc)
performance::check_model(q1.3_sc)

# Q1.3 : SEX x DIAGNOSIS x TIME: POWER LATERALIZATION

q1.3_gl <- lmer(
  gamma_lat ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1.3_gl)
performance::r2(q1.3_gl)
performance::check_model(q1.3_gl)

# Singularity issue again, rerun without random effect
q1.3_gl_lm <- lm(
  gamma_lat ~ timepoint_c * sex * outcome_bin + site + iq,
  data = eeg_ela
)
summary(q1.3_gl_lm)
performance::r2(q1.3_gl_lm)
performance::check_model(q1.3_gl_lm)

# Compare all four models
performance::compare_performance(q3_fg, q3_ac, q3_sc, q3_gl)


extract_p <- function(model, term) {
  coefs     <- as.data.frame(coef(summary(model)))
  coefs$term <- rownames(coefs)
  row       <- coefs[grepl(term, coefs$term, fixed = TRUE), ]
  if (nrow(row) == 0) return(NA_real_)
  row$`Pr(>|t|)`[1]
}

# Extract p for sexF (sex difference at 6 months = intercept)
p_sex_q1 <- c(
  front_gamma  = extract_p(q1_fg, "sexF"),
  auditory_con = extract_p(q1_ac, "sexF"),
  speech_con   = extract_p(q1_sc_lm, "sexF"),
  gamma_lat    = extract_p(q1_gl_lm, "sexF")
)

p_sex_q2 <- c(
  front_gamma  = extract_p(q1.2_fg, "sexF"),
  auditory_con = extract_p(q1.2_ac, "sexF"),
  speech_con   = extract_p(q1.2_sc_lm, "sexF"),
  gamma_lat    = extract_p(q1.2_gl_lm, "sexF")
)

p_sex_q3 <- c(
  front_gamma  = extract_p(q1.3_fg, "sexF"),
  auditory_con = extract_p(q1.3_ac, "sexF"),
  speech_con   = extract_p(q1.3_sc, "sexF"),
  gamma_lat    = extract_p(q1.3_gl, "sexF")
)
 
p_traj_q1 <- c(
  front_gamma  = extract_p(q1_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q1_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q1_sc_lm, "timepoint_c:sexF"),
  gamma_lat    = extract_p(q1_gl_lm, "timepoint_c:sexF")
)

# Q1: 8 tests, 4 intercept (sex at 6 mo) + 4 trajectory (sex x time)
p_q1_full      <- c(p_sex_q1, p_traj_q1)
p_q1_full_holm <- p.adjust(p_q1_full, method = "holm")
cat("Q1 full family Holm (intercept + trajectory):\n")
print(round(p_q1_full_holm, 4))

p_traj_q2 <- c(
  front_gamma  = extract_p(q1.2_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q1.2_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q1.2_sc_lm, "timepoint_c:sexF"),
  gamma_lat    = extract_p(q1.2_gl_lm, "timepoint_c:sexF")
)

p_q2_full      <- c(p_sex_q2, p_traj_q2)
p_q2_full_holm <- p.adjust(p_q2_full, method = "holm")
cat("Q2 full family Holm (intercept + trajectory):\n")
print(round(p_q2_full_holm, 4))

# Q1.3 : 8 tests, 4 intercept + 4 trajectory
p_traj_q3 <- c(
  front_gamma  = extract_p(q1.3_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q1.3_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q1.3_sc,    "timepoint_c:sexF"),
  gamma_lat    = extract_p(q1.3_gl_lm, "timepoint_c:sexF")
)

p_q3_full      <- c(p_sex_q3, p_traj_q3)
p_q3_full_holm <- p.adjust(p_q3_full, method = "holm")
p_q3_full_bh   <- p.adjust(p_q3_full, method = "BH")
cat("Q3 full family Holm (intercept + trajectory, exploratory):\n")
print(round(p_q3_full_holm, 4))
cat("Q3 full family BH:\n")
print(round(p_q3_full_bh, 4))

p_int_q2 <- c(
  front_gamma  = extract_p(q1.2_fg,    "sexF:group_typeELA"),
  auditory_con = extract_p(q1.2_ac,    "sexF:group_typeELA"),
  speech_con   = extract_p(q1.2_sc_lm, "sexF:group_typeELA"),
  gamma_lat    = extract_p(q1.2_gl_lm, "sexF:group_typeELA")
)
p_int_traj_q2 <- c(
  front_gamma  = extract_p(q1.2_fg,    "timepoint_c:sexF:group_typeELA"),
  auditory_con = extract_p(q1.2_ac,    "timepoint_c:sexF:group_typeELA"),
  speech_con   = extract_p(q1.2_sc_lm, "timepoint_c:sexF:group_typeELA"),
  gamma_lat    = extract_p(q1.2_gl_lm, "timepoint_c:sexF:group_typeELA")
)
p_q2_full      <- c(p_int_q2, p_int_traj_q2)
p_q2_full_holm <- p.adjust(p_q2_full, method = "holm")
cat("Q1.2 family Holm (sex × ASD-likelihood interaction):\n")
print(round(p_q2_full_holm, 4))

# Q1.3: Sex × ASD diagnosis interaction terms
p_int_q3 <- c(
  front_gamma  = extract_p(q1.3_fg,    "sexF:outcome_binASD"),
  auditory_con = extract_p(q1.3_ac,    "sexF:outcome_binASD"),
  speech_con   = extract_p(q1.3_sc,    "sexF:outcome_binASD"),
  gamma_lat    = extract_p(q1.3_gl_lm, "sexF:outcome_binASD")
)
p_int_traj_q3 <- c(
  front_gamma  = extract_p(q1.3_fg,    "timepoint_c:sexF:outcome_binASD"),
  auditory_con = extract_p(q1.3_ac,    "timepoint_c:sexF:outcome_binASD"),
  speech_con   = extract_p(q1.3_sc,    "timepoint_c:sexF:outcome_binASD"),
  gamma_lat    = extract_p(q1.3_gl_lm, "timepoint_c:sexF:outcome_binASD")
)
p_q3_full      <- c(p_int_q3, p_int_traj_q3)
p_q3_full_holm <- p.adjust(p_q3_full, method = "holm")
p_q3_full_bh   <- p.adjust(p_q3_full, method = "BH")


# EMMEANS MODEL
# Corrected pairwise contrasts for speech connectivity sex differences

# 1. Speech connectivity
emm_q1_sc <- emmeans(q1_sc_lm, ~ sex | timepoint_c,
                     at = list(timepoint_c = c(0, 6)))

# Holm adjustment within this contrast family
pairs(emm_q1_sc, simple = "sex", adjust = "holm")


# Auditory connectivity: sex differences by group (TLA vs ELA)
emm_q1.2_ac_6mo <- emmeans(q1.2_ac, ~ sex | group_type,
                           at = list(timepoint_c = 0))
pairs(emm_q1.2_ac_6mo, simple = "sex", adjust = "holm")



# Marginal R2 for each model
# Marginal R2: use lm $r.squared for singular-replaced models
r2_q1 <- c(
  performance::r2(q1_fg)$R2_marginal,
  performance::r2(q1_ac)$R2_marginal,
  summary(q1_sc_lm)$r.squared,        # was singular lmer
  summary(q1_gl_lm)$r.squared         # was singular lmer
)

r2_q2 <- c(
  performance::r2(q1.2_fg)$R2_marginal,
  performance::r2(q1.2_ac)$R2_marginal,
  summary(q1.2_sc_lm)$r.squared,        # was singular lmer
  summary(q1.2_gl_lm)$r.squared         # was singular lmer
)

r2_q3 <- c(
  performance::r2(q1.3_fg)$R2_marginal,
  performance::r2(q1.3_ac)$R2_marginal,
  performance::r2(q1.3_sc)$R2_marginal,
  summary(q1.3_gl_lm)$r.squared         # was singular lmer
)

# Assemble summary table

eeg_labels <- c("Frontal Gamma Power", "Auditory Network",
                "Speech Network",      "Power Lateralization")
# Model type label: flag which outcomes used lm() vs lmer()
model_type_q1 <- c("lmer", "lmer", "lm*", "lm*")
model_type_q2 <- c("lmer", "lmer", "lm*", "lm*")
model_type_q3 <- c("lmer", "lmer", "lmer", "lm*")

summary_table <- data.frame(
  EEG_Variable = rep(rep(eeg_labels, each = 1), 6),
  Question = c(
    rep("Sex at 6 mo (Q1)",          4),
    rep("Sex × Time (Q1)",           4),
    rep("Sex × ASD-likelihood at 6 mo (Q1.2)",  4),
    rep("Sex × ASD-likelihood × Time (Q1.2)",   4),
    rep("Sex × ASD-diagnosis at 6 mo (Q1.3)",   4),
    rep("Sex × ASD-diagnosis × Time (Q1.3)",    4)
  ),
  Test = c(
    rep("Intercept",    4), rep("Trajectory", 4),   # Q1
    rep("Interaction",  4), rep("3-way",      4),   # Q1.2
    rep("Interaction",  4), rep("3-way",      4)    # Q1.3
  ),
  Model_type = c(model_type_q1, model_type_q1,
                 model_type_q2, model_type_q2,
                 model_type_q3, model_type_q3),
  p_raw = round(c(
    p_q1_full,
    p_q2_full,
    p_q3_full
  ), 4),
  p_holm = round(c(
    p_q1_full_holm,
    p_q2_full_holm,
    p_q3_full_holm
  ), 4),
  p_BH = c(
    rep(NA, 16),
    round(p_q3_full_bh, 4)
  ),
  marginal_R2 = round(c(
    r2_q1, r2_q1,
    r2_q2, r2_q2,
    r2_q3, r2_q3
  ), 4)
)

# Add significance stars based on Holm-corrected p
summary_table$sig_holm <- case_when(
  summary_table$p_holm < 0.001 ~ "***",
  summary_table$p_holm < 0.01  ~ "**",
  summary_table$p_holm < 0.05  ~ "*",
  summary_table$p_holm < 0.10  ~ ".",
  TRUE                          ~ "ns"
)

print(summary_table)
# Save as Word via flextable

ft_summary <- flextable(summary_table) %>%
  set_header_labels(
    EEG_Variable = "EEG Variable",
    Question     = "Comparison",
    Test         = "Test type",
    Model_type   = "Model",
    p_raw        = "p (raw)",
    p_holm       = "p (Holm)",
    p_BH         = "p (BH)",
    marginal_R2  = "R²",
    sig_holm     = "Sig."
  ) %>%
  bold(i = ~ p_holm < 0.05, j = ~ p_holm + sig_holm) %>%
  color(i = ~ p_holm < 0.05, j = ~ p_holm, color = "#A85848") %>%
  italic(i = ~ Test == "Intercept", j = ~ Question) %>%
  bg(i = ~ Test == "Intercept", bg = "#F5F5F5") %>%
  autofit() %>%
  add_header_lines("Supplementary Table X. Sex differences in EEG metrics — full corrected results") %>%
  add_footer_lines(
    "Holm correction applied within each question family (8 tests per family). Q1: 4 sex-at-6mo + 4 sex×time. Q1.2: 4 sex×ASD-likelihood + 4 sex×ASD-likelihood×time. Q1.3: 4 sex×ASD-diagnosis + 4 sex×ASD-diagnosis×time."
  ) %>%
  add_footer_lines(
    "p (BH) shown for Q1.3 only (exploratory). lm* = singular lmer replaced with ordinary linear model (random intercept variance = 0)."
  ) %>%
  add_footer_lines(
    "R² = marginal R² for lmer models; ordinary R² for lm* models. Q1.2 Interaction = sexF:group_typeELA; Q1.2 3-way = timepoint_c:sexF:group_typeELA. Q1.3 Interaction = sexF:outcome_binASD; Q1.3 3-way = timepoint_c:sexF:outcome_binASD."
  )
ft_summary

save_as_docx(ft_summary, path = "tables/EEG_sex_diff_corrected_results MAY.docx")






### SECTION TWO: EEG AND LANGUAGE RELATIONSHIP

raw <- read.csv(
  "../../datasets/sex_diff_main.csv",
  header = TRUE, sep = ","
)


# Factor key variables
raw <- raw %>%
  mutate(
    
    sex        = factor(sex,        levels = c("F", "M")),
    site       = factor(site),
    group_type = factor(group_type, levels = c("TLA", "ELA")),
    subject    = factor(subject),
    outcome    = factor(tolower(outcome), levels = c("no-asd", "asd"))
  )

# Remove those without speech connectivity data 
raw <- raw %>% filter(!is.na(speech_con_6_left))


# Select data for speech connectivity and language outcomes
eeg_long <- raw %>%
  select(subject, sex, site, group_type, outcome,nonverbal_iq_6,
         speech_con_6_left, speech_con_12_left) %>%
  pivot_longer(
    cols            = c(speech_con_6_left,speech_con_12_left),
    names_to = "timepoint",
    values_to = "speech_con_left",
    names_pattern = "speech_con_(6|12)_left"
  )%>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6)       # centre: 0 = 6 mo, 6 = 12 mo)


# LONG FORMAT: EXPRESSIVE LANGUAGE

exp_long <- raw %>%
  select(subject, sex, site, group_type, speech_con_6_left, outcome,nonverbal_iq_6,
         expressive_6, expressive_12, expressive_18, expressive_24, expressive_36) %>%
  pivot_longer(
    cols         = starts_with("expressive_"),
    names_to     = "timepoint",
    names_prefix = "expressive_",
    values_to    = "expressive"
  ) %>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6          # 0 = 6 mo, 6 = 12 mo, 12 = 18 mo, 18 = 24 mo, 30 = 36 mo
  )

cat("\nExpressive language: observations per timepoint\n")
exp_long %>%
  group_by(timepoint) %>%
  summarise(n = sum(!is.na(expressive))) %>%
  print()



# LONG FORMAT: RECEPTIVE LANGUAGE
rec_long <- raw %>%
  select(subject, sex, site, group_type, speech_con_6_left, outcome,nonverbal_iq_6,
         receptive_6, receptive_12, receptive_18, receptive_24, receptive_36) %>%
  pivot_longer(
    cols         = starts_with("receptive_"),
    names_to     = "timepoint",
    names_prefix = "receptive_",
    values_to    = "receptive"
  ) %>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6
  )

cat("\nReceptive language: observations per timepoint\n")
rec_long %>%
  group_by(timepoint) %>%
  summarise(n = sum(!is.na(receptive))) %>%
  print()



# Sample sizes per group
eeg_long %>%
  distinct(subject, sex, group_type) %>%
  count(sex, group_type)

# Missing data per EEG variable at each timepoint
eeg_long %>%
  group_by(timepoint) %>%
  summarise(
    n_speech_con    = sum(!is.na(speech_con_left)))


# SPEECH CONNECTIVITY DESCRIPTIVE STATISTICS

sc_mean <- mean(exp_long$speech_con_6_left, na.rm = TRUE)
sc_sd   <- sd(exp_long$speech_con_6_left,   na.rm = TRUE)
cat(sprintf("Speech connectivity: M = %.3f, SD = %.3f\n", sc_mean, sc_sd))

# CENTRE SPEECH CONNECTIVITY
exp_long <- exp_long %>% mutate(speech_con_c = speech_con_6_left - sc_mean)
rec_long <- rec_long %>% mutate(speech_con_c = speech_con_6_left - sc_mean)

# Q2.1: EXPRESSIVE: Speech Con × Sex × Time
q2.1_exp <- lmer(
  expressive ~ timepoint_c * speech_con_c * sex + nonverbal_iq_6 + site +(1 | subject),
  data    = exp_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)


summary(q2.1_exp)
performance::r2(q2.1_exp)
performance::check_model(q2.1_exp)


# Q2.1: RECEPTIVE: Speech Con × Sex × Time

q2.1_rec <- lmer(
  receptive ~ timepoint_c * sex * speech_con_c + site +nonverbal_iq_6+ (1 | subject),
  data    = rec_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q2.1_rec)
performance::r2(q2.1_rec)
performance::check_model(q2.1_rec)

# 2.1.1 Does speech connectivity at 6 months predict receptive language at 36 months?

speech_outcomes_df <- raw %>%
  select(subject, sex, site, group_type, outcome,nonverbal_iq_6,
         expressive_36, receptive_36,expressive_6, expressive_12,
         receptive_6, receptive_12) %>%
  left_join(
    eeg_long %>% filter(timepoint == 6) %>% select(subject, speech_con_left),
    by = "subject"
  ) %>% filter(!is.na(receptive_36))

glimpse(speech_outcomes_df)


receptive_speech_36_outcomes <- lm(
  receptive_36 ~ speech_con_left * sex  + nonverbal_iq_6 ,
  data = speech_outcomes_df
)

summary(receptive_speech_36_outcomes)

# Low / mean / high values for marginal predictions later
sc_vals        <- c(sc_mean - sc_sd, sc_mean, sc_mean + sc_sd)
sc_labels      <- c("Low (-1 SD)", "Mean", "High (+1 SD)")
sc_vals_labels <- setNames(sc_vals, sc_labels)


# Multiple Comparisons: Holm correction across 4 tests

p_exp_growth <- extract_p(q2.1_exp, "timepoint_c:speech_con_c:sexM")
p_rec_growth <- extract_p(q2.1_rec, "timepoint_c:sexM:speech_con_c")
# Sex x timepoint interaction p-values from original models (uncentered)
p_exp_sex_time <- summary(q2.1_exp)$coefficients["timepoint_c:sexM", "Pr(>|t|)"]
p_rec_time <- summary(q2.1_rec)$coefficients["timepoint_c:sexM", "Pr(>|t|)"]

# 36-month outcome models (from original code, uncentered)
 p_rec36 <- summary(receptive_speech_36_outcomes)$coefficients["speech_con_left:sexM", "Pr(>|t|)"]

p_family <- c(
  "Expressive growth (3-way)"  = p_exp_growth,
  "Receptive growth (3-way)"   = p_rec_growth,
   "Receptive 36mo (interact)"  = p_rec36, 
  "Expressive sex x time" = p_exp_sex_time,
  "Receptive sex x time"  = p_rec_time
)

correction_tbl <- tibble(
  Test          = names(p_family),
  p_raw         = unname(p_family),
  p_holm        = p.adjust(p_family, method = "holm"),
  p_BH          = p.adjust(p_family, method = "BH")
)
print(correction_tbl, n = Inf)

correction_tbl

# EFFECT SIZES
# Standardised effect for the 3-way interaction (receptive)
cat("Effect sizes for key interactions")
standardize_parameters(q2.1_rec, method = "refit") %>%
  filter(grepl("timepoint_c", Parameter)) %>%
  print()

# Separate models 

rec_female <- rec_long %>% filter(sex == "F")
rec_male   <- rec_long %>% filter(sex == "M")

lme_f <- lmer(
  receptive ~ timepoint_c * speech_con_c + nonverbal_iq_6 + site + (1 | subject),
  data = rec_female, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

lme_m <- lmer(
  receptive ~ timepoint_c * speech_con_c + nonverbal_iq_6 + site + (1 | subject),
  data = rec_male, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(lme_f)
summary(lme_m)

# Effect sizes per stratum
standardize_parameters(lme_f, method = "refit")
standardize_parameters(lme_m, method = "refit")

# Simple slopes within each sex model
# Standard simple slopes at ±1 SD and mean of speech_con_c

sim_slopes(lme_f, pred = timepoint_c, modx = speech_con_c,
           modx.values = "plus-minus", jnplot = TRUE)
sim_slopes(lme_m, pred = timepoint_c, modx = speech_con_c,
           modx.values = "plus-minus", jnplot = TRUE)



# PLOTS 

# Shared style constants (matched to Python notebook palette)
F_DARK    <- "#A85848"
F_LIGHT   <- "#D4A898"      # light salmon
F_MID     <- "#BE8070"      # midpoint between F_LIGHT and F_DARK
M_DARK    <- "#5A6848"
M_LIGHT   <- "#B0BCA0"      # light sage
M_MID     <- "#859274"      # midpoint between M_LIGHT and M_DARK

pal_sex   <- c("F" = F_DARK,   "M" = M_DARK)
shape_sex <- c("F" = 15,       "M" = 19)
lty_sex   <- c("F" = "dashed", "M" = "solid")

theme_pub <- theme_classic(base_size = 13, base_family = "sans") +
  theme(
    plot.title       = element_text(face = "bold",   size = 14),
    plot.subtitle    = element_text(colour = "#999999", size = 11),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 12),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 12),
    legend.position  = "bottom",
    legend.key.width = unit(1.8, "lines"),
    strip.background = element_blank(),                # ← removes grey boxes
    strip.text       = element_text(face = "bold", size = 12, colour = "#333333"),
    panel.spacing    = unit(1.2, "lines"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white",  colour = NA)
  )


# Q2.1 Predicted receptive language trajectories: Sex × SpeechCon
pred_grid_rec <- expand.grid(
  timepoint_c      = c(0, 6, 12, 18, 30),
  sex              = factor(c("F", "M"), levels = c("F", "M")),
  speech_con_6_left = sc_vals,
  site             = factor("london"),
  nonverbal_iq_6   = mean(rec_long$nonverbal_iq_6, na.rm = TRUE)
) %>%
  mutate(
    speech_con_c  = speech_con_6_left - sc_mean,
    timepoint     = timepoint_c + 6,
    sc_group      = factor(speech_con_6_left,
                           levels = sc_vals,
                           labels = sc_labels)
  )

pred_grid_rec$predicted <- predict(q2.1_rec,
                                   newdata  = pred_grid_rec,
                                   re.form  = NA)

p_rec_traj <- ggplot(pred_grid_rec,
                     aes(x        = timepoint,
                         y        = predicted,
                         colour   = sex,
                         linetype = sex,
                         shape    = sex,
                         group    = sex)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.0) +
  facet_wrap(~ sc_group, nrow = 1,
             labeller = labeller(sc_group = function(x) paste("Speech connectivity:", x))) +
  scale_colour_manual(values = pal_sex,    labels = c("Female", "Male")) +
  scale_linetype_manual(values = lty_sex,  labels = c("Female", "Male")) +
  scale_shape_manual(values = shape_sex,   labels = c("Female", "Male")) +
  scale_x_continuous(breaks = c(6, 12, 18, 24, 36)) +
  labs(
    x        = "Age (months)",
    y        = "Predicted Receptive Language Score",
    colour   = "Sex",
    linetype = "Sex",
    shape    = "Sex",
    title    = "Sex \u00d7 Speech Connectivity \u00d7 Time on Receptive Language",
  ) +
  theme_pub

print(p_rec_traj)

ggsave("../../figures/main/Fig2_receptive_sex_speech_trajectories.tiff",
       p_rec_traj, width = 10, height = 5, bg = "white",
       dpi = 300, compression = "lzw")


# Q2.1 Predicted expressive language trajectories: Sex × SpeechCon (supplemental)
pred_grid_exp <- expand.grid(
  timepoint_c       = c(0, 6, 12, 18, 30),
  sex               = factor(c("F", "M"), levels = c("F", "M")),
  speech_con_6_left = sc_vals,
  site              = factor("london"),
  nonverbal_iq_6    = mean(exp_long$nonverbal_iq_6, na.rm = TRUE)
) %>%
  mutate(
    speech_con_c = speech_con_6_left - sc_mean,
    timepoint    = timepoint_c + 6,
    sc_group     = factor(speech_con_6_left, levels = sc_vals, labels = sc_labels)
  )

pred_grid_exp$predicted <- predict(q2.1_exp,
                                   newdata = pred_grid_exp,
                                   re.form = NA)

p_exp_traj <- ggplot(pred_grid_exp,
                     aes(x        = timepoint,
                         y        = predicted,
                         colour   = sex,
                         linetype = sex,
                         shape    = sex,
                         group    = sex)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.0) +
  facet_wrap(~ sc_group, nrow = 1) +
  scale_colour_manual(values = pal_sex,   labels = c("Female", "Male")) +
  scale_linetype_manual(values = lty_sex, labels = c("Female", "Male")) +
  scale_shape_manual(values = shape_sex,  labels = c("Female", "Male")) +
  scale_x_continuous(breaks = c(6, 12, 18, 24, 36)) +
  labs(
    x        = "Age (months)",
    y        = "Predicted Expressive Language Score",
    colour   = "Sex",
    linetype = "Sex",
    shape    = "Sex",
    title    = "Sex \u00d7 Speech Connectivity \u00d7 Time on Expressive Language",
    subtitle = "3-way interaction not significant after correction (p = 0.057, p-adj > 0.05)"
  ) +
  theme_pub

print(p_exp_traj)
ggsave("../../figures/supplementary/S7_expressive_sex_speech_trajectories.pdf",
       p_exp_traj, width = 10, height = 5, bg = "white")


# Scatter plot: Speech connectivity vs 36-month receptive, by sex (EXPLORATORY)
p_scatter_36 <- ggplot(speech_outcomes_df,
                       aes(x      = speech_con_left,
                           y      = receptive_36,
                           colour = sex,
                           fill   = sex,
                           shape  = sex)) +
  geom_point(alpha = 0.80, size = 3.0) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 1.4) +
  scale_colour_manual(values = pal_sex,   labels = c("Female", "Male")) +
  scale_fill_manual(values   = pal_sex,   labels = c("Female", "Male"),
                    guide    = "none") +
  scale_shape_manual(values  = shape_sex, labels = c("Female", "Male")) +
  labs(
    x        = "Left-hemisphere speech connectivity at 6 months",
    y        = "Receptive language score at 36 months",
    colour   = "Sex",
    shape    = "Sex",
    title    = "Speech Connectivity Predicts 36-Month Receptive Language",
  ) +
  theme_pub

print(p_scatter_36)
ggsave("../../figures/supplementary/S8_speech_con_receptive_36mo_scatter.tiff",
       p_scatter_36, width = 7, height = 6, bg = "white")







# SECTION 3 CLUSTERING AND LANGUAGE ANALYSES

raw <- read.csv(
  "../../datasets/sex_diff_main.csv",
  header = TRUE, sep = ","
)

# Factor key variables
raw <- raw %>%
  mutate(
    sex        = factor(sex,        levels = c("M", "F")),
    site       = factor(site),
    group_type = factor(group_type, levels = c("TLA", "ELA")),
    hc_class   = factor(as.integer(hc_class), levels = c(0L, 1L, 2L)),
    subject    = factor(subject),
    outcome    = factor(tolower(outcome), levels = c("no-asd", "asd"))
  )

# Remove those without HC 
raw <- raw %>% filter(!is.na(hc_class))

# DEMOGRAPHICS: SEX × HC CLASS DISTRIBUTION


cat("\nN per HC class × sex:\n")
raw %>%
  count(hc_class, sex) %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  print()

cat("\nN per HC class × ASD likelihood:\n")
raw %>%
  count(hc_class, group_type) %>%
  pivot_wider(names_from = group_type, values_from = n, values_fill = 0) %>%
  print()

cat("\nN per HC class × sex × ASD likelihood:\n")
raw %>%
  count(hc_class, sex, group_type) %>%
  pivot_wider(names_from = c(sex, group_type), values_from = n, values_fill = 0) %>%
  print()

# Chi-square: sex composition across HC classes
sex_class_table <- table(raw$hc_class, raw$sex)
cat("\nChi-square: sex ~ HC class\n")
print(chisq.test(sex_class_table))

# Check average expressive language scores for each subgroup at timepoint 36
raw %>%
  group_by(hc_class,sex) %>%
  summarise(mean_expressive_36 = mean(expressive_36, na.rm = TRUE),
            sd_expressive_36 = sd(expressive_36, na.rm = TRUE),
            n = sum(!is.na(expressive_36)))

# Check average receotive language scores for each subgroup at timepoint 36
raw %>%
  group_by(hc_class) %>%
  summarise(mean_receptive_36 = mean(receptive_36, na.rm = TRUE),
            sd_receptive_36 = sd(receptive_36, na.rm = TRUE),
            n = sum(!is.na(receptive_36)))


# LONG FORMAT: EXPRESSIVE LANGUAGE

exp_long <- raw %>%
  select(subject, sex, site, group_type, hc_class, outcome,nonverbal_iq_6,
         expressive_6, expressive_12, expressive_18, expressive_24, expressive_36) %>%
  pivot_longer(
    cols         = starts_with("expressive_"),
    names_to     = "timepoint",
    names_prefix = "expressive_",
    values_to    = "expressive"
  ) %>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6          # 0 = 6 mo, 6 = 12 mo, 12 = 18 mo, 18 = 24 mo, 30 = 36 mo
  )

# LONG FORMAT: RECEPTIVE LANGUAGE

rec_long <- raw %>%
  select(subject, sex, site, group_type, hc_class, outcome,nonverbal_iq_6,
         receptive_6, receptive_12, receptive_18, receptive_24, receptive_36) %>%
  pivot_longer(
    cols         = starts_with("receptive_"),
    names_to     = "timepoint",
    names_prefix = "receptive_",
    values_to    = "receptive"
  ) %>%
  mutate(
    timepoint   = as.numeric(timepoint),
    timepoint_c = timepoint - 6
  )

cat("\nReceptive language: observations per timepoint\n")
rec_long %>%
  group_by(timepoint) %>%
  summarise(n = sum(!is.na(receptive))) %>%
  print()


# Q3.1 EXPRESSIVE: HC class × Sex × Time

# (Do different HC classes show different sex gaps in language development?)

q3_exp <- lmer(
  expressive ~ timepoint_c * sex * hc_class + site + nonverbal_iq_6 + (1 | subject),
  data    = exp_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)


summary(q3_exp)
anova(q3_exp)                  # Type III F-tests (Satterthwaite df)
performance::r2(q3_exp)performance::check_model(q3_exp)



# Q3.1: RECEPTIVE: HC class × Sex × Time

q3_rec <- lmer(
  receptive ~ timepoint_c * sex * hc_class + site + nonverbal_iq_6 + (1 | subject),
  data    = rec_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_rec)
anova(q3_rec)
performance::r2(q3_rec)
performance::check_model(q3_rec)



# 8. MULTIPLE TESTING CORRECTIONS & SUMMARY TABLE


# Helper: extract p from anova() row matching a term
extract_anova_p <- function(model, term_pattern) {
  aov_tbl      <- as.data.frame(anova(model))
  aov_tbl$term <- rownames(aov_tbl)
  row          <- aov_tbl[aov_tbl$term == term_pattern, ]
  if (nrow(row) == 0) return(NA_real_)
  row$`Pr(>F)`[1]
}
# Raw p-values
# Family 1: all sex-related omnibus terms from Q1 and Q3
p_fam_primary <- c(
  exp_sex_t        = extract_anova_p(q3_exp, "timepoint_c:sex"),
  exp_sex_class    = extract_anova_p(q3_exp, "sex:hc_class"),
  exp_sex_class_t  = extract_anova_p(q3_exp, "timepoint_c:sex:hc_class"),
  rec_sex_t        = extract_anova_p(q3_rec, "timepoint_c:sex"),
  rec_sex_class    = extract_anova_p(q3_rec, "sex:hc_class"),
  rec_sex_class_t  = extract_anova_p(q3_rec, "timepoint_c:sex:hc_class")
)

p_fam_primary_holm <- p.adjust(p_fam_primary, method = "holm")
p_fam_primary_holm


# Post-hoc pairwise comparisons for primary three-way models.
emm_exp_pairwise <- emmeans(q3_exp, ~ sex * hc_class | timepoint_c,
                            at = list(timepoint_c = c(0, 30)))
cat("\nExpressive language: sex contrasts by cluster at 6 and 36 months\n")
pairs(emm_exp_pairwise, simple = "sex", adjust = "holm")

emm_rec_pairwise <- emmeans(q3_rec, ~ sex * hc_class | timepoint_c,
                            at = list(timepoint_c = c(0, 30)))
cat("\nReceptive language: sex contrasts by cluster at 6 and 36 months\n")
pairs(emm_rec_pairwise, simple = "sex", adjust = "holm")

# Summary table

summary_table <- data.frame(
  Language = c("Expressive", "Expressive", "Expressive",
               "Receptive",  "Receptive",  "Receptive"),
  Term = c("Sex × Time",
           "Sex × Cluster",
           "Sex × Cluster × Time",
           "Sex × Time",
           "Sex × Cluster",
           "Sex × Cluster × Time"),
  p_raw  = round(p_fam_primary, 4),
  p_holm = round(p_fam_primary_holm, 4),
  marginal_R2 = round(c(
    rep(performance::r2(q3_exp)$R2_marginal, 3),
    rep(performance::r2(q3_rec)$R2_marginal, 3)
  ), 4),
  sig = case_when(
    p_fam_primary_holm < 0.001 ~ "***",
    p_fam_primary_holm < 0.01  ~ "**",
    p_fam_primary_holm < 0.05  ~ "*",
    p_fam_primary_holm < 0.10  ~ ".",
    TRUE                       ~ "ns"
  )
)

print(summary_table)

# Save as Word via flextable
ft_summary <- flextable(summary_table) %>%
  set_header_labels(
    Language    = "Language",
    Term        = "Interaction Term",
    p_raw       = "p (raw)",
    p_holm      = "p (Holm)",
    marginal_R2 = "Marginal R2",
    sig         = "Sig. (Holm)"
  ) %>%
  bold(~ p_holm < 0.05, ~ p_holm + sig) %>%
  color(~ p_holm < 0.05, ~ p_holm, color = "#A85848") %>%
  theme_zebra() %>%
  autofit() %>%
  add_header_lines("Language Growth: HC Class × Sex Interactions") %>%
  add_footer_lines(
    "Holm FWER correction applied within each family. Family 1: Q1+Q3 (sex×class×time)"
  )

ft_summary
save_as_docx(ft_summary, path = "tables/main/Language_HCclass_sex_corrected_results.docx")


# Supplemmentary analyses

# Check if were are differences in non-verbal IQ between males and females at 6mo

non_verbal_df <- raw %>%
  select(subject, sex, group_type, outcome, site,
         nonverbal_iq_6) %>%
  filter(!is.na(nonverbal_iq_6))

# Model: non-verbal IQ ~ sex
nonverbal_iq_model <- lm(
  nonverbal_iq_6 ~ sex + site,
  data = non_verbal_df
)
summary(nonverbal_iq_model)

