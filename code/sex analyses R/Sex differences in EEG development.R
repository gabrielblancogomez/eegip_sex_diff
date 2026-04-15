# EEG SEX DIFFERENCES IN INFANT DEVELOPMENT 

# Author Gabriel Blanco-Gomez

# Questions:
#   Main Q1: Do males and females differ in EEG development across various EEG measures (6->12 months)?
#   Supplemetary Q3: Does sex interact with ASD likelihood (TLA vs ELA)?
#   Supplemetary Q3: Does sex interact with ASD diagnosis? (ELA only)
#


# 1. SETUP

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



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dir.create("figures/main",          recursive = TRUE, showWarnings = FALSE)
dir.create("figures/supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/main",           recursive = TRUE, showWarnings = FALSE)
dir.create("tables/supplementary",  recursive = TRUE, showWarnings = FALSE)


# 2. DATA PREPARATION

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


# 3a. Q1 — SEX x TIME: FRONTAL GAMMA POWER

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



# 3b. Q1 — SEX x TIME: AUDITORY NETWORK CONNECTIVITY

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
leveneTest(residuals(q1_ac) ~ eeg_long$sex[!is.na(eeg_long$auditory_con)])

# HC3 robustness check — compare with original
q1_ac_lm <- lm(auditory_con ~ timepoint_c * sex + site, data = eeg_long)
coeftest(q1_ac_lm, vcov = vcovHC(q1_ac_lm, type = "HC3"))

# 3c. Q1 — SEX x TIME: SPEECH NETWORK CONNECTIVITY

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
# Rerun with singurla model

q1_sc_lm <- lm(
  speech_con ~ timepoint_c * sex + site +iq ,
  data = eeg_long
)

summary(q1_sc_lm)
performance::r2(q1_sc_lm)
performance::check_model(q1_sc_lm)


#Preliminary results = sig difference in speech connectivity

# 3d. Q1 — SEX x TIME: POWER LATERALIZATION

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
performance::compare_performance(q1_fg, q1_ac, q1_sc, q1_gl)



# 4a. Q2 — SEX x LIKELIHOOD x TIME: FRONTAL GAMMA POWER

q2_fg <- lmer(
  front_gamma ~ timepoint_c * sex * group_type + site + iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q2_fg)
performance::r2(q2_fg)
performance::check_model(q2_fg)


# 4b. Q2 — SEX x LIKELIHOOD x TIME: AUDITORY NETWORK CONNECTIVITY

q2_ac <- lmer(
  auditory_con ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q2_ac)
performance::r2(q2_ac)
performance::check_model(q2_ac)


# 4c. Q2 — SEX x LIKELIHOOD x TIME: SPEECH NETWORK CONNECTIVITY

q2_sc <- lmer(
  speech_con ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q2_sc)
performance::r2(q2_sc)
performance::check_model(q2_sc)

# Singularity issue again, rerun without random effect
q2_sc_lm <- lm(
  speech_con ~ timepoint_c * sex * group_type + site + iq,
  data = eeg_long
)

summary(q2_sc_lm)
performance::r2(q2_sc_lm)
performance::check_model(q2_sc_lm)
emmeans(q2_sc_lm, pairwise ~ sex | group_type)
# 4d. Q2 — SEX x LIKELIHOOD x TIME: POWER LATERALIZATION

q2_gl <- lmer(
  gamma_lat ~ timepoint_c * sex * group_type + site +iq +
    (1 | subject),
  data    = eeg_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q2_gl)
performance::r2(q2_gl)
performance::check_model(q2_gl)

# Singularity issue again, rerun without random effect
q2_gl_lm <- lm(
  gamma_lat ~ timepoint_c * sex * group_type + site +iq,
  data = eeg_long
)


summary(q2_gl_lm)
performance::r2(q2_gl_lm)
performance::check_model(q2_gl_lm)


# -- Q2 compare all four models -----------------------------------------------
performance::compare_performance(q2_fg, q2_ac, q2_sc, q2_gl)


# 5. Q3 PREP — ELA SUBSAMPLE + ASD OUTCOME CODING  [EXPLORATORY]

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

# N per cell — check before modelling
eeg_ela %>%
  distinct(subject, sex, outcome_bin) %>%
  count(sex, outcome_bin)


# 5a. Q3 [EXPLORATORY] — SEX x DIAGNOSIS x TIME: FRONTAL GAMMA POWER

q3_fg <- lmer(
  front_gamma ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_fg)
performance::r2(q3_fg)
performance::check_model(q3_fg)


# 5b. Q3 [EXPLORATORY] — SEX x DIAGNOSIS x TIME: AUDITORY NETWORK


q3_ac <- lmer(
  auditory_con ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_ac)
performance::r2(q3_ac)
performance::check_model(q3_ac)


 
# 5c. Q3 [EXPLORATORY] — SEX x DIAGNOSIS x TIME: SPEECH NETWORK
 
q3_sc <- lmer(
  speech_con ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_sc)
performance::r2(q3_sc)
performance::check_model(q3_sc)

# 5d. Q3 [EXPLORATORY] — SEX x DIAGNOSIS x TIME: POWER LATERALIZATION

q3_gl <- lmer(
  gamma_lat ~ timepoint_c * sex * outcome_bin + site +iq +
    (1 | subject),
  data    = eeg_ela,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_gl)
performance::r2(q3_gl)
performance::check_model(q3_gl)

# Singularity issue again, rerun without random effect
q3_gl_lm <- lm(
  gamma_lat ~ timepoint_c * sex * outcome_bin + site + iq,
  data = eeg_ela
)
summary(q3_gl_lm)
performance::r2(q3_gl_lm)
performance::check_model(q3_gl_lm)

#  Q3 compare all four models 
performance::compare_performance(q3_fg, q3_ac, q3_sc, q3_gl)


extract_p <- function(model, term) {
  coefs     <- as.data.frame(coef(summary(model)))
  coefs$term <- rownames(coefs)
  row       <- coefs[grepl(term, coefs$term, fixed = TRUE), ]
  if (nrow(row) == 0) return(NA_real_)
  row$`Pr(>|t|)`[1]
}

#  Extract p for sexF (sex difference at 6 months = intercept)
p_sex_q1 <- c(
  front_gamma  = extract_p(q1_fg, "sexF"),
  auditory_con = extract_p(q1_ac, "sexF"),
  speech_con   = extract_p(q1_sc, "sexF"),
  gamma_lat    = extract_p(q1_gl, "sexF")
)

p_sex_q2 <- c(
  front_gamma  = extract_p(q2_fg, "sexF"),
  auditory_con = extract_p(q2_ac, "sexF"),
  speech_con   = extract_p(q2_sc, "sexF"),
  gamma_lat    = extract_p(q2_gl, "sexF")
)

p_sex_q3 <- c(
  front_gamma  = extract_p(q3_fg, "sexF"),
  auditory_con = extract_p(q3_ac, "sexF"),
  speech_con   = extract_p(q3_sc, "sexF"),
  gamma_lat    = extract_p(q3_gl, "sexF")
)
 
p_traj_q1 <- c(
  front_gamma  = extract_p(q1_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q1_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q1_sc_lm, "timepoint_c:sexF"),
  gamma_lat    = extract_p(q1_gl_lm, "timepoint_c:sexF")
)

# Q1: 8 tests — 4 intercept (sex at 6 mo) + 4 trajectory (sex x time)
p_q1_full      <- c(p_sex_q1, p_traj_q1)
p_q1_full_holm <- p.adjust(p_q1_full, method = "holm")
cat("Q1 full family Holm (intercept + trajectory):\n")
print(round(p_q1_full_holm, 4))

p_traj_q2 <- c(
  front_gamma  = extract_p(q2_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q2_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q2_sc_lm, "timepoint_c:sexF"),
  gamma_lat    = extract_p(q2_gl_lm, "timepoint_c:sexF")
)

p_q2_full      <- c(p_sex_q2, p_traj_q2)
p_q2_full_holm <- p.adjust(p_q2_full, method = "holm")
cat("Q2 full family Holm (intercept + trajectory):\n")
print(round(p_q2_full_holm, 4))

# Q3 [EXPLORATORY]: 8 tests — 4 intercept + 4 trajectory
p_traj_q3 <- c(
  front_gamma  = extract_p(q3_fg,    "timepoint_c:sexF"),
  auditory_con = extract_p(q3_ac,    "timepoint_c:sexF"),
  speech_con   = extract_p(q3_sc,    "timepoint_c:sexF"),
  gamma_lat    = extract_p(q3_gl_lm, "timepoint_c:sexF")
)

p_q3_full      <- c(p_sex_q3, p_traj_q3)
p_q3_full_holm <- p.adjust(p_q3_full, method = "holm")
p_q3_full_bh   <- p.adjust(p_q3_full, method = "BH")
cat("Q3 full family Holm (intercept + trajectory, exploratory):\n")
print(round(p_q3_full_holm, 4))
cat("Q3 full family BH:\n")
print(round(p_q3_full_bh, 4))

# **** EMMEANS MODEL ****
# Corrected pairwise contrasts for speech connectivity sex differences

# 1. Speech connectivity
emm_q1_sc <- emmeans(q1_sc_lm, ~ sex | timepoint_c,
                     at = list(timepoint_c = c(0, 6)))

# Holm adjustment within this contrast family
pairs(emm_q1_sc, simple = "sex", adjust = "holm")

# 2. Frontal Gammma in ELA
# Corrected pairwise contrasts for frontal gamma in ELA
emm_q3_fg <- emmeans(q3_fg, ~ sex * outcome_bin | timepoint_c,
                     at = list(timepoint_c = c(0, 6)))

pairs(emm_q3_fg, simple = "sex", adjust = "holm")



#  Marginal R2 for each model 
#  Marginal R2: use lm $r.squared for singular-replaced models -
r2_q1 <- c(
  performance::r2(q1_fg)$R2_marginal,
  performance::r2(q1_ac)$R2_marginal,
  summary(q1_sc_lm)$r.squared,        # was singular lmer
  summary(q1_gl_lm)$r.squared         # was singular lmer
)

r2_q2 <- c(
  performance::r2(q2_fg)$R2_marginal,
  performance::r2(q2_ac)$R2_marginal,
  summary(q2_sc_lm)$r.squared,        # was singular lmer
  summary(q2_gl_lm)$r.squared         # was singular lmer
)

r2_q3 <- c(
  performance::r2(q3_fg)$R2_marginal,
  performance::r2(q3_ac)$R2_marginal,
  performance::r2(q3_sc)$R2_marginal,
  summary(q3_gl_lm)$r.squared         # was singular lmer
)

#  Assemble summary table 

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
    rep("Sex at 6 mo (Q2, TLA ref)", 4),
    rep("Sex × Time (Q2)",           4),
    rep("Sex at 6 mo (Q3, no-ASD ref)", 4),
    rep("Sex × Time (Q3)",           4)
  ),
  Test = c(
    rep("Intercept",   4), rep("Trajectory", 4),
    rep("Intercept",   4), rep("Trajectory", 4),
    rep("Intercept",   4), rep("Trajectory", 4)
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
#- Save as Word via flextable 

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
    "Holm correction applied within each question family (8 tests per family: 4 intercept + 4 trajectory)."
  ) %>%
  add_footer_lines(
    "p (BH) shown for Q3 only (exploratory). lm* = singular lmer replaced with ordinary linear model (random intercept variance = 0)."
  ) %>%
  add_footer_lines(
    "R² = marginal R² for lmer models; ordinary R² for lm* models. Q2 and Q3 intercept tests reflect sex difference in reference group (TLA and no-ASD respectively)."
  )

ft_summary

save_as_docx(ft_summary, path = "tables/EEG_sex_diff_corrected_results.docx")


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

