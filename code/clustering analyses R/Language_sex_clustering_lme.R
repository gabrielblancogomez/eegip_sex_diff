# LANGUAGE GROWTH: HC CLASS × SEX INTERACTIONS
#
# Questions:
#   Q1: Does HC class interact with sex to predict expressive language growth?
#   Q2: Same as Q1 for receptive language.
#   Q3: Exploratory Does HC class × sex × ASD likelihood interact to predict expressive language growth?
#   Q4: Same as Q3 for receptive language.
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

glimpse(raw)
cat("N subjects:", nrow(raw), "\n")
cat("Columns:",    ncol(raw), "\n")

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


# 3. DEMOGRAPHICS — SEX × HC CLASS DISTRIBUTION


cat("\n-- N per HC class × sex ------------------------------------------\n")
raw %>%
  count(hc_class, sex) %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  print()

cat("\n-- N per HC class × ASD likelihood ----------------------------------\n")
raw %>%
  count(hc_class, group_type) %>%
  pivot_wider(names_from = group_type, values_from = n, values_fill = 0) %>%
  print()

cat("\n-- N per HC class × sex × ASD likelihood ----------------------------\n")
raw %>%
  count(hc_class, sex, group_type) %>%
  pivot_wider(names_from = c(sex, group_type), values_from = n, values_fill = 0) %>%
  print()

# Chi-square: sex composition across HC classes
sex_class_table <- table(raw$hc_class, raw$sex)
cat("\n-- Chi-square: sex ~ HC class ----------------------------------------\n")
print(chisq.test(sex_class_table))

# Check average expressive language scores for each subgroup at timepoint 36
raw %>%
  group_by(hc_class) %>%
  summarise(mean_expressive_36 = mean(expressive_36, na.rm = TRUE),
            sd_expressive_36 = sd(expressive_36, na.rm = TRUE),
            n = sum(!is.na(expressive_36)))

# Check average receotive language scores for each subgroup at timepoint 36
raw %>%
  group_by(hc_class) %>%
  summarise(mean_receptive_36 = mean(receptive_36, na.rm = TRUE),
            sd_receptive_36 = sd(receptive_36, na.rm = TRUE),
            n = sum(!is.na(receptive_36)))


# 4. LONG FORMAT — EXPRESSIVE LANGUAGE

exp_long <- raw %>%
  select(subject, sex, site, group_type, hc_class, outcome,
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

cat("\n-- Expressive language: observations per timepoint -------------------\n")
exp_long %>%
  group_by(timepoint) %>%
  summarise(n = sum(!is.na(expressive))) %>%
  print()



# 5. LONG FORMAT — RECEPTIVE LANGUAGE


rec_long <- raw %>%
  select(subject, sex, site, group_type, hc_class, outcome,
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

cat("\n-- Receptive language: observations per timepoint --------------------\n")
rec_long %>%
  group_by(timepoint) %>%
  summarise(n = sum(!is.na(receptive))) %>%
  print()


# 6a. Q1 — EXPRESSIVE: HC class × Sex × Time

# (Do different HC classes show different sex gaps in language development?)

q1_exp <- lmer(
  expressive ~ timepoint_c * sex * hc_class + site + (1 | subject),
  data    = exp_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q1_exp)
anova(q1_exp)                  # Type III F-tests (Satterthwaite df)
performance::r2(q1_exp)
performance::check_model(q1_exp)

 

# So far, it looks like there a strong sex difference but weak evidence for a 
# sex x class difference. Pairs emmeans test suggest sex differences in languae at
# 36 months within HC class B but not at earlier timepoints

# 6b. Q2 — EXPRESSIVE: HC class × Sex × ASD Likelihood × Time
# Key term: timepoint_c : sex : hc_class : group_type
# (Is the sex × class pattern in language moderated by ASD likelihood?)

q2_exp <- lmer(
  expressive ~ timepoint_c * sex * hc_class * group_type + site + (1 | subject),
  data    = exp_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
# Check multi collinearity


summary(q2_exp)
anova(q2_exp)
performance::r2(q2_exp)
performance::check_collinearity(q2_exp)
performance::check_model(q2_exp)

emm_q2 <- emmeans(q2_exp, ~ sex * hc_class * group_type | timepoint_c,
                  at = list(timepoint_c = c(0, 30)))
pairs(emm_q2, simple = "sex")


# -- Compare Q1 and Q2 expressive models --------------------------------------
performance::compare_performance(q1_exp, q2_exp)


# 7a. Q3 — RECEPTIVE: HC class × Sex × Time

q3_rec <- lmer(
  receptive ~ timepoint_c * sex * hc_class + site + (1 | subject),
  data    = rec_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q3_rec)
anova(q3_rec)
performance::r2(q3_rec)
performance::check_model(q3_rec)
 


# 7b. Q4 — RECEPTIVE: HC class × Sex × ASD Likelihood × Time


q4_rec <- lmer(
  receptive ~ timepoint_c * sex * hc_class * group_type + site + (1 | subject),
  data    = rec_long,
  REML    = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(q4_rec)
anova(q4_rec)
performance::r2(q4_rec)
performance::check_model(q4_rec)


# -- Compare Q3 and Q4 receptive models ---------------------------------------
performance::compare_performance(q3_rec, q4_rec)


# 8. MULTIPLE TESTING CORRECTIONS & SUMMARY TABLE


# -- Helper: extract p from anova() row matching a term ----------------------
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
  exp_sex_time     = extract_anova_p(q1_exp, "timepoint_c:sex"),
  exp_sex_class    = extract_anova_p(q1_exp, "sex:hc_class"),
  exp_sex_class_t  = extract_anova_p(q1_exp, "timepoint_c:sex:hc_class"),
  rec_sex_time     = extract_anova_p(q3_rec, "timepoint_c:sex"),
  rec_sex_class    = extract_anova_p(q3_rec, "sex:hc_class"),
  rec_sex_class_t  = extract_anova_p(q3_rec, "timepoint_c:sex:hc_class")
)

p_fam_primary_holm <- p.adjust(p_fam_primary, method = "holm")

# SUPPLEMENTARY/EXPLORATORY: four-way model corrections reported separately
# and clearly labelled as not primary due to VIF and cell-size concerns.
p_fam_exploratory <- c(
  exp_four_way = extract_anova_p(q2_exp, "timepoint_c:sex:hc_class:group_type"),
  rec_four_way = extract_anova_p(q4_rec, "timepoint_c:sex:hc_class:group_type")
)
p_fam_exploratory_holm <- p.adjust(p_fam_exploratory, method = "holm")

cat("\n-- PRIMARY Holm-corrected p-values (Q1 + Q3) ----------------------------\n")
print(round(p_fam_primary_holm, 4))

cat("\n-- EXPLORATORY four-way p-values (Q2 + Q4, interpret with caution) ------\n")
cat("NOTE: VIF > 300 for key terms; minimum cell n = 4. For reference only.\n")
print(round(p_fam_exploratory_holm, 4))

cat("\n-- Raw p-values (primary family) ----------------------------------------\n")
print(round(p_fam_primary, 4))
# Marginal R2 for each model
safe_r2_marginal <- function(model) {
  r2_result <- tryCatch(performance::r2(model), error = function(e) NULL)
  if (is.null(r2_result) || is.null(r2_result$R2_marginal)) return(NA_real_)
  r2_result$R2_marginal
}

r2_q1 <- safe_r2_marginal(q1_exp)
r2_q2 <- safe_r2_marginal(q2_exp)   # kept for reference
r2_q3 <- safe_r2_marginal(q3_rec)
r2_q4 <- safe_r2_marginal(q4_rec)   # kept for reference


# Pairwise comparisons
# Four-way models are retained for reference onl because: 
#   (1) severe multicollinearity (VIF > 300 for key interaction terms), and
#   (2) insufficient cell sizes in some hc_class × sex × group_type cells
#       (minimum n = 4; see demographics table above).
# Primary pairwise comparisons are drawn from Q1 (expressive) and Q3 (receptive).

# Post-hoc pairwise comparisons for primary three-way models.
emm_exp_pairwise <- emmeans(q1_exp, ~ sex * hc_class | timepoint_c,
                            at = list(timepoint_c = c(0, 30)))
cat("\n-- Expressive language: sex contrasts by cluster at 6 and 36 months ----\n")
pairs(emm_exp_pairwise, simple = "sex", adjust = "holm")

emm_rec_pairwise <- emmeans(q3_rec, ~ sex * hc_class | timepoint_c,
                            at = list(timepoint_c = c(0, 30)))
cat("\n-- Receptive language: sex contrasts by cluster at 6 and 36 months -----\n")
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
  marginal_R2 = round(c(rep(r2_q1, 3), rep(r2_q3, 3)), 4),
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

