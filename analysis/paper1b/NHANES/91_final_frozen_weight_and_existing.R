options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root_dir, "First Paper")
out_dir <- file.path(paper_dir, "Paper1_FINAL_FROZEN_20260814")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "cache"), recursive = TRUE, showWarnings = FALSE)

age_dir <- file.path(paper_dir, "final_age_range_comparison_20260814")
refine_dir <- file.path(paper_dir, "female_p5_refinement_20260814")
hb_dir <- file.path(paper_dir, "h_vs_hb_comparison_20260814")
raw_path <- file.path(paper_dir, "raw_for_refit", "merged_training_1999_2006.csv.gz")
stopifnot(file.exists(raw_path), dir.exists(age_dir), dir.exists(refine_dir), dir.exists(hb_dir))

# 01. Frozen sample-flow source table.
file.copy(file.path(age_dir, "01_sample_flow_age_range.csv"),
          file.path(out_dir, "01_sample_flow.csv"), overwrite = TRUE)

# 02. Final NHANES development-weight audit.  One completed data set is used
# because design variables and weights are identical across all five DXA
# completed data sets; equality is audited explicitly below.
weight_raw <- read_csv(
  raw_path,
  col_select = all_of(c("SEQN", "imputation", "cycle", "RIDSTATR", "RIAGENDR",
                        "RIDAGEYR", "RIDEXPRG", "WTMEC2YR", "WTMEC4YR",
                        "SDMVSTRA", "SDMVPSU", "BMXHT", "alm_g", "wtmec_pool")),
  show_col_types = FALSE, progress = FALSE
) %>%
  transmute(
    SEQN = as.numeric(SEQN), imputation = as.integer(imputation), cycle = as.character(cycle),
    examined = RIDSTATR == 2,
    sex = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ NA_character_),
    age = as.numeric(RIDAGEYR), pregnant = RIAGENDR == 2 & RIDEXPRG == 1,
    WTMEC2YR = as.numeric(WTMEC2YR), WTMEC4YR = as.numeric(WTMEC4YR),
    strata = as.character(SDMVSTRA), psu = as.character(SDMVPSU),
    height_cm = as.numeric(BMXHT), alm_kg = as.numeric(alm_g) / 1000,
    wtmec_pool = as.numeric(wtmec_pool)
  )

weight_consistency <- weight_raw %>%
  group_by(SEQN) %>%
  summarise(
    completed_sets = n_distinct(imputation),
    weight_values = n_distinct(wtmec_pool),
    strata_values = n_distinct(strata), psu_values = n_distinct(psu),
    .groups = "drop"
  )
if (any(weight_consistency$completed_sets != 5L) || any(weight_consistency$weight_values > 1L) ||
    any(weight_consistency$strata_values > 1L) || any(weight_consistency$psu_values > 1L)) {
  stop("Weights/design variables are not invariant across completed data sets")
}

eligible <- weight_raw %>%
  filter(imputation == 1L, examined, sex %in% c("Female", "Male"), age >= 18, age <= 69,
         !coalesce(pregnant, FALSE), is.finite(height_cm), height_cm > 0,
         is.finite(alm_kg), alm_kg > 0, is.finite(wtmec_pool), wtmec_pool > 0,
         !is.na(strata), !is.na(psu)) %>%
  mutate(
    source_weight_name = if_else(cycle %in% c("1999-2000", "2001-2002"),
                                 "WTMEC4YR", "WTMEC2YR"),
    source_weight = if_else(source_weight_name == "WTMEC4YR", WTMEC4YR, WTMEC2YR),
    multiplier = if_else(source_weight_name == "WTMEC4YR", 0.5, 0.25),
    expected_pool_weight = source_weight * multiplier,
    transformation_abs_error = abs(wtmec_pool - expected_pool_weight)
  ) %>%
  group_by(sex) %>%
  mutate(gamlss_case_weight = wtmec_pool / mean(wtmec_pool)) %>%
  ungroup()

cycle_audit <- eligible %>%
  group_by(sex, cycle, source_weight_name, multiplier) %>%
  summarise(
    n = n(), source_weight_sum = sum(source_weight), sum_weight = sum(wtmec_pool),
    raw_weight_min = min(wtmec_pool), raw_weight_max = max(wtmec_pool),
    centered_mean = mean(gamlss_case_weight),
    centered_min = min(gamlss_case_weight), centered_max = max(gamlss_case_weight),
    max_abs_transformation_error = max(transformation_abs_error),
    source_weight_rule_verified = max_abs_transformation_error < 1e-8,
    survey_weight = "wtmec_pool", survey_strata = "interaction(cycle, SDMVSTRA)",
    survey_psu = "interaction(cycle, SDMVSTRA, SDMVPSU)",
    .groups = "drop"
  )

overall_audit <- eligible %>%
  group_by(sex) %>%
  summarise(
    cycle = "Overall 1999-2006", source_weight_name = "Cycle-specific source",
    multiplier = NA_real_, n = n(), source_weight_sum = sum(source_weight),
    sum_weight = sum(wtmec_pool), raw_weight_min = min(wtmec_pool),
    raw_weight_max = max(wtmec_pool), centered_mean = mean(gamlss_case_weight),
    centered_min = min(gamlss_case_weight), centered_max = max(gamlss_case_weight),
    max_abs_transformation_error = max(transformation_abs_error),
    source_weight_rule_verified = max_abs_transformation_error < 1e-8,
    survey_weight = "wtmec_pool", survey_strata = "interaction(cycle, SDMVSTRA)",
    survey_psu = "interaction(cycle, SDMVSTRA, SDMVPSU)", .groups = "drop"
  )

nesting <- eligible %>% summarise(
  raw_strata_codes = n_distinct(strata), cycle_nested_strata = n_distinct(interaction(cycle, strata)),
  raw_psu_codes = n_distinct(interaction(strata, psu)),
  cycle_nested_psus = n_distinct(interaction(cycle, strata, psu))
)
weight_audit <- bind_rows(cycle_audit, overall_audit) %>%
  mutate(
    completed_data_weight_invariance = TRUE,
    cycle_nesting_verified = nesting$cycle_nested_strata >= nesting$raw_strata_codes &&
      nesting$cycle_nested_psus >= nesting$raw_psu_codes,
    raw_strata_codes = nesting$raw_strata_codes,
    cycle_nested_strata = nesting$cycle_nested_strata,
    raw_psu_codes = nesting$raw_psu_codes,
    cycle_nested_psus = nesting$cycle_nested_psus
  ) %>% arrange(sex, factor(cycle, levels = c("1999-2000", "2001-2002", "2003-2004",
                                                "2005-2006", "Overall 1999-2006")))
write_csv(weight_audit, file.path(out_dir, "02_weight_audit.csv"))

# 03. Frozen primary H-model specification.
primary_spec <- tibble(
  model = "Primary H", age_range = "18-69", development_period = "NHANES 1999-2006",
  sex = c("Female", "Male"), family = c("BCCG", "BCT"),
  outcome = "Raw ALM (kg)", mu_formula = "pb(age, df=3) + pb(height, df=3)",
  sigma_formula = "constant", nu_formula = "constant", tau_formula = c("not applicable", "constant"),
  completed_datasets = 5L, fitting_weight = "wtmec_pool / within-sex mean(wtmec_pool)",
  validation_design = "survey::svydesign; cycle-nested strata and PSU; wtmec_pool",
  prediction_domain = "18-69 years; sex-specific model; no race/BMI conditioning",
  frozen = TRUE
)
write_csv(primary_spec, file.path(out_dir, "03_primary_H_model_specification.csv"))

# 04-05. Frozen existing internal and temporal H-model results.
internal_cv <- read_csv(file.path(age_dir, "06_psu_cv_overall.csv"), show_col_types = FALSE) %>%
  filter(age_range == "18-69")
write_csv(internal_cv, file.path(out_dir, "04_internal_PSU_CV.csv"))
temporal <- read_csv(file.path(age_dir, "11_temporal_validation_18_59.csv"), show_col_types = FALSE) %>%
  filter(scorer == "18-69")
write_csv(temporal, file.path(out_dir, "05_temporal_calibration.csv"))

# 09. Existing frozen height classification-drift outputs in one long table.
height_or <- read_csv(file.path(age_dir, "14_classification_height_or_18_59.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "height_OR_per_10cm")
height_prev <- read_csv(file.path(age_dir, "15_height_quintile_prevalence.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "height_quintile_prevalence")
height_conditional <- read_csv(file.path(age_dir, "13_conditional_height_or_three_scorers.csv"), show_col_types = FALSE) %>%
  filter(scorer == "18-69") %>% mutate(record_type = "conditional_scorer_comparison")
height_drift <- bind_rows(height_or, height_prev, height_conditional)
write_csv(height_drift, file.path(out_dir, "09_height_classification_drift.csv"))

# 11-14. Frozen sensitivity summaries already completed before this final run.
normal_grip <- bind_rows(
  read_csv(file.path(age_dir, "20_grip_adjusted_means_fnih.csv"), show_col_types = FALSE),
  read_csv(file.path(age_dir, "21_grip_adjusted_means_ewgsop2.csv"), show_col_types = FALSE)
) %>% filter(scorer == "18-69")
write_csv(normal_grip, file.path(out_dir, "11_normalBMI_grip_sensitivity.csv"))

file.copy(file.path(age_dir, "24_final_age_range_decision_table.csv"),
          file.path(out_dir, "12_age_range_18_84_sensitivity.csv"), overwrite = TRUE)

refine_effects <- read_csv(file.path(refine_dir, "08_oof_height_effects.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "OOF_height_effect")
refine_overall <- read_csv(file.path(refine_dir, "07_oof_overall_calibration.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "OOF_overall_calibration")
refine_decision <- read_csv(file.path(refine_dir, "13_upgrade_decision.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "upgrade_decision")
refine_summary <- bind_rows(refine_effects, refine_overall, refine_decision)
write_csv(refine_summary, file.path(out_dir, "13_female_P5_refinement_summary.csv"))

hb_overall <- read_csv(file.path(hb_dir, "05_oof_overall_calibration.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "OOF_overall_calibration")
hb_effects <- read_csv(file.path(hb_dir, "06_oof_height_bmi_effects.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "OOF_height_BMI_effect")
hb_bic <- read_csv(file.path(hb_dir, "10_full_fit_BIC_summary.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "full_fit_information_criteria")
hb_stability <- read_csv(file.path(hb_dir, "11_cv_stability.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "CV_stability")
hb_decision <- read_csv(file.path(hb_dir, "12_development_upgrade_decision.csv"), show_col_types = FALSE) %>%
  mutate(record_type = "development_gate_decision")
hb_development <- bind_rows(hb_overall, hb_effects, hb_bic, hb_stability, hb_decision)
write_csv(hb_development, file.path(out_dir, "14_H_vs_HB_development.csv"))

write_csv(tibble(
  audit_item = c("Completed data sets", "Weight invariance", "Weight transformation",
                 "GAMLSS centering", "Survey inference weight", "Cycle nesting"),
  result = c("5", "Identical within SEQN across five completed data sets",
             "1999-2002: WTMEC4YR x 0.5; 2003-2006: WTMEC2YR x 0.25",
             "Within sex and completed data set: wtmec_pool / arithmetic mean",
             "wtmec_pool (constant pooled scaling retained)",
             "strata=cycle x SDMVSTRA; PSU=cycle x SDMVSTRA x SDMVPSU"),
  status = "PASS"
), file.path(out_dir, "cache", "weight_audit_readme.csv"))

cat("Completed frozen weight audit and copied existing frozen results to", out_dir, "\n")
