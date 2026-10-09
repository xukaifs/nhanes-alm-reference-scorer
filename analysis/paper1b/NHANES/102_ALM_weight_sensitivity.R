# =============================================================================
# 102_ALM_weight_sensitivity.R
#
# Purpose:
#   Add ALM/weight definitions from Reiter et al. (Clinical Nutrition, 2025)
#   to the SAME 2011-2018 NHANES validation framework already used in Paper 1.
#
# Main comparison (same common analysis sample):
#   1) Corrected conditional H P5 (18-69 reference)
#   2) EWGSOP2 low ALMI
#   3) FNIH low ALM/BMI
#   4) ALM/weight Definition 1
#   5) ALM/weight Definition 2
#
# Outputs:
#   - overall survey-weighted prevalence
#   - adjusted height OR per +10 cm
#   - sex-specific height-quintile prevalence
#   - BMI phenotype (low vs non-low)
#   - one combined summary table
#
# IMPORTANT:
#   The current frozen validation cache does not carry measured body weight.
#   Therefore body weight is reconstructed as BMI * height^2. This is expected
#   to be extremely close to the measured weight from which NHANES BMI was
#   calculated, but minor differences can arise from BMI rounding.
#
# This script DOES NOT refit any GAMLSS model and DOES NOT overwrite primary
# Paper 1 outputs.
# =============================================================================

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(splines)
})

set.seed(20260817)

# -----------------------------------------------------------------------------
# 0. Locate project and input
# -----------------------------------------------------------------------------
find_project_root <- function(start_dir = getwd()) {
  d <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 0:10) {
    candidate <- file.path(
      d, "First Paper", "raw_for_refit", "merged_training_1999_2006.csv.gz"
    )
    if (file.exists(candidate)) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop(
    "Could not locate project root. Expected to find:\n",
    "First Paper/raw_for_refit/merged_training_1999_2006.csv.gz\n",
    "above current working directory:\n", getwd(),
    call. = FALSE
  )
}

root_dir <- find_project_root()
paper_dir <- file.path(root_dir, "First Paper")
frozen_dir <- file.path(paper_dir, "Paper1_FINAL_FROZEN_20260814")
cache_path <- file.path(frozen_dir, "cache", "validation_corrected_H_18_59.rds")

out_dir <- file.path(paper_dir, "ALM_weight_sensitivity_20260817")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_path)) {
  stop(
    "Final corrected validation cache was not found:\n", cache_path,
    "\nPlease confirm that 93_corrected_H_validation_analyses.R has been run.",
    call. = FALSE
  )
}

cat("Project root:\n", root_dir, "\n")
cat("Input:\n", cache_path, "\n")
cat("Output:\n", out_dir, "\n\n")

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------
make_design <- function(dat) {
  svydesign(
    ids = ~psu_pool,
    strata = ~strata_pool,
    weights = ~wtmec_pool,
    nest = TRUE,
    data = dat
  )
}

weighted_quantile <- function(x, w, probabilities) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative <- cumsum(w) / sum(w)
  vapply(
    probabilities,
    function(p) x[which(cumulative >= p)[1]],
    numeric(1)
  )
}

pool_scalar <- function(estimates, standard_errors, null = 0) {
  ok <- is.finite(estimates) & is.finite(standard_errors) & standard_errors >= 0
  estimates <- estimates[ok]
  standard_errors <- standard_errors[ok]
  m <- length(estimates)

  if (!m) {
    return(tibble(
      m = 0, estimate = NA_real_, standard_error = NA_real_,
      df = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_
    ))
  }

  qbar <- mean(estimates)
  within <- mean(standard_errors^2)
  between <- if (m > 1) var(estimates) else 0
  total <- within + (1 + 1 / m) * between
  standard_error <- sqrt(total)

  df <- if (m > 1 && between > 0) {
    (m - 1) * (1 + within / ((1 + 1 / m) * between))^2
  } else {
    Inf
  }

  critical <- qt(0.975, df)
  statistic <- if (standard_error > 0) (qbar - null) / standard_error else NA_real_

  tibble(
    m = m,
    estimate = qbar,
    standard_error = standard_error,
    df = df,
    conf_low = qbar - critical * standard_error,
    conf_high = qbar + critical * standard_error,
    p_value = if (is.finite(statistic)) {
      2 * pt(abs(statistic), df = df, lower.tail = FALSE)
    } else {
      NA_real_
    }
  )
}

pool_table <- function(dat, groups, estimate_col = "estimate", se_col = "standard_error") {
  dat %>%
    group_by(across(all_of(groups))) %>%
    group_modify(~ pool_scalar(.x[[estimate_col]], .x[[se_col]])) %>%
    ungroup()
}

safe_mean <- function(design, variable) {
  fit <- try(
    svymean(as.formula(paste0("~", variable)), design, na.rm = TRUE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(c(estimate = NA_real_, standard_error = NA_real_))
  }
  c(
    estimate = as.numeric(coef(fit)[1]),
    standard_error = as.numeric(SE(fit)[1])
  )
}

fit_height_effect <- function(dat, outcome) {
  formula <- as.formula(
    paste(
      outcome,
      "~ I(height_cm/10) + ns(age, df=3) + race + cycle"
    )
  )
  try(
    svyglm(
      formula,
      design = make_design(dat),
      family = quasibinomial()
    ),
    silent = TRUE
  )
}

# -----------------------------------------------------------------------------
# 2. Read final corrected validation sample and derive ALM/weight
# -----------------------------------------------------------------------------
validation <- readRDS(cache_path)

# Prefer measured body weight if it is available either in the frozen cache or
# in the five processed temporal-validation source files. Otherwise fall back
# to BMI * height^2.
weight_source <- "BMI x height_m^2 reconstructed weight"
validation$weight_kg_for_analysis <- NA_real_

cache_weight_candidates <- c("weight_kg", "BMXWT", "body_weight_kg", "weight")
cache_weight_hit <- intersect(cache_weight_candidates, names(validation))

if (length(cache_weight_hit) > 0) {
  chosen_weight <- cache_weight_hit[[1]]
  validation$weight_kg_for_analysis <- as.numeric(validation[[chosen_weight]])
  weight_source <- paste0("Measured weight from frozen cache: ", chosen_weight)
} else {
  raw_validation_paths <- file.path(
    root_dir, "新建文件夹",
    paste0("07_validation_reference_imputation", 1:5, ".csv.gz")
  )

  if (all(file.exists(raw_validation_paths))) {
    first_names <- names(
      suppressMessages(
        read_csv(raw_validation_paths[[1]], n_max = 0, show_col_types = FALSE)
      )
    )
    raw_weight_hit <- intersect(
      c("weight_kg", "BMXWT", "body_weight_kg", "weight"),
      first_names
    )

    if (length(raw_weight_hit) > 0) {
      chosen_weight <- raw_weight_hit[[1]]

      measured_weight <- map2_dfr(
        raw_validation_paths, 1:5,
        function(path, imp) {
          read_csv(
            path,
            col_select = all_of(c("SEQN", chosen_weight)),
            show_col_types = FALSE,
            progress = FALSE
          ) %>%
            transmute(
              SEQN = as.numeric(SEQN),
              imputation = imp,
              weight_kg_measured = as.numeric(.data[[chosen_weight]])
            ) %>%
            distinct(SEQN, .keep_all = TRUE)
        }
      )

      validation <- validation %>%
        left_join(measured_weight, by = c("SEQN", "imputation")) %>%
        mutate(weight_kg_for_analysis = weight_kg_measured)

      weight_source <- paste0(
        "Measured weight from temporal-validation source files: ",
        chosen_weight
      )
    }
  }
}

validation <- validation %>%
  mutate(
    sex = factor(as.character(sex), levels = c("Female", "Male")),
    race = factor(as.character(race)),
    cycle = factor(as.character(cycle)),
    strata_pool = factor(strata_pool),
    psu_pool = factor(psu_pool),
    age = as.numeric(age),
    height_m = as.numeric(height_m),
    height_cm = as.numeric(height_cm),
    bmi = as.numeric(bmi),
    alm_kg = as.numeric(alm_kg),
    almi = as.numeric(almi),
    wtmec_pool = as.numeric(wtmec_pool),

    weight_kg_reconstructed = bmi * height_m^2,
    weight_kg_for_analysis = if_else(
      is.finite(weight_kg_for_analysis) & weight_kg_for_analysis > 0,
      weight_kg_for_analysis,
      weight_kg_reconstructed
    ),
    alm_weight_pct = 100 * alm_kg / weight_kg_for_analysis,

    # Recalculate all comparator flags from their definitions to keep this
    # sensitivity analysis self-contained.
    low_ewgsop2_recalc = case_when(
      sex == "Male"   ~ as.numeric(almi < 7.0),
      sex == "Female" ~ as.numeric(almi < 5.5),
      TRUE ~ NA_real_
    ),

    alm_bmi_ratio_recalc = alm_kg / bmi,
    low_fnih_recalc = case_when(
      sex == "Male"   ~ as.numeric(alm_bmi_ratio_recalc < 0.789),
      sex == "Female" ~ as.numeric(alm_bmi_ratio_recalc < 0.512),
      TRUE ~ NA_real_
    ),

    # Reiter et al. 2025 / ALST(weight) definitions.
    low_weight_def1 = case_when(
      sex == "Male"   ~ as.numeric(alm_weight_pct < 25.7),
      sex == "Female" ~ as.numeric(alm_weight_pct < 19.4),
      TRUE ~ NA_real_
    ),
    low_weight_def2 = case_when(
      sex == "Male"   ~ as.numeric(alm_weight_pct < 28.27),
      sex == "Female" ~ as.numeric(alm_weight_pct < 23.47),
      TRUE ~ NA_real_
    ),

    low_conditional_p5 = as.numeric(low_p5_H_corrected)
  ) %>%
  filter(
    age >= 18, age <= 59,
    sex %in% c("Female", "Male"),
    is.finite(height_m), height_m > 0,
    is.finite(height_cm), height_cm > 0,
    is.finite(bmi), bmi > 0,
    is.finite(alm_kg), alm_kg > 0,
    is.finite(almi),
    is.finite(weight_kg_for_analysis), weight_kg_for_analysis > 0,
    is.finite(alm_weight_pct),
    is.finite(wtmec_pool), wtmec_pool > 0,
    is.finite(low_conditional_p5),
    !is.na(race), !is.na(cycle),
    !is.na(strata_pool), !is.na(psu_pool)
  )

# Definition metadata
definitions <- tribble(
  ~definition_order, ~definition, ~outcome, ~normalization, ~cutoff_text,
  1L, "Conditional H P5", "low_conditional_p5",
      "ALM | age + height", "Percentile < 5",
  2L, "EWGSOP2 low ALMI", "low_ewgsop2_recalc",
      "ALM / height^2", "Men <7.0 kg/m^2; women <5.5 kg/m^2",
  3L, "FNIH low ALM/BMI", "low_fnih_recalc",
      "ALM / BMI", "Men <0.789; women <0.512",
  4L, "ALM/weight Definition 1", "low_weight_def1",
      "100 x ALM / weight", "Men <25.7%; women <19.4%",
  5L, "ALM/weight Definition 2", "low_weight_def2",
      "100 x ALM / weight", "Men <28.27%; women <23.47%"
)

write_csv(
  definitions %>%
    mutate(
      analysis_domain = "NHANES 2011-2018, age 18-59",
      survey_design = "WTMEC pooled weight + cycle-nested strata/PSU",
      weight_source_for_ALM_weight = weight_source
    ),
  file.path(out_dir, "00_definition_table.csv")
)

# QA: compare stored vs recalculated existing comparator flags
qa_existing <- validation %>%
  summarise(
    rows = n(),
    stored_ewgsop2_available = sum(is.finite(low_ewgsop2)),
    ewgsop2_discordant = sum(
      is.finite(low_ewgsop2) &
      low_ewgsop2 != low_ewgsop2_recalc
    ),
    stored_fnih_available = sum(is.finite(low_fnih)),
    fnih_discordant = sum(
      is.finite(low_fnih) &
      low_fnih != low_fnih_recalc
    )
  )

sample_qa <- validation %>%
  filter(imputation == 1) %>%
  group_by(sex) %>%
  summarise(
    unweighted_n = n(),
    age_min = min(age),
    age_max = max(age),
    mean_weight_kg_for_analysis = weighted.mean(
      weight_kg_for_analysis, wtmec_pool
    ),
    mean_reconstructed_weight_kg = weighted.mean(
      weight_kg_reconstructed, wtmec_pool
    ),
    mean_alm_weight_pct = weighted.mean(alm_weight_pct, wtmec_pool),
    .groups = "drop"
  )

write_csv(qa_existing, file.path(out_dir, "01a_existing_definition_QA.csv"))
write_csv(sample_qa, file.path(out_dir, "01b_analysis_sample_QA.csv"))

# -----------------------------------------------------------------------------
# 3. Fixed sex-specific height quintiles from imputation 1
# -----------------------------------------------------------------------------
height_cutpoints <- map_dfr(c("Female", "Male"), function(sex_value) {
  dat <- validation %>%
    filter(imputation == 1, as.character(sex) == sex_value)

  cuts <- weighted_quantile(
    dat$height_m,
    dat$wtmec_pool,
    seq(0, 1, 0.2)
  )

  tibble(
    sex = sex_value,
    boundary = 0:5,
    probability = seq(0, 1, 0.2),
    height_m = cuts,
    height_cm = 100 * cuts
  )
})

write_csv(
  height_cutpoints,
  file.path(out_dir, "04_height_quintile_cutpoints.csv")
)

assign_height_quintile <- function(dat) {
  map_dfr(c("Female", "Male"), function(sex_value) {
    tmp <- dat %>% filter(as.character(sex) == sex_value)

    cuts <- height_cutpoints %>%
      filter(sex == sex_value) %>%
      arrange(boundary) %>%
      pull(height_m)

    cuts[c(1, length(cuts))] <- c(-Inf, Inf)

    tmp %>%
      mutate(
        height_quintile = cut(
          height_m,
          breaks = cuts,
          labels = paste0("Q", 1:5),
          include.lowest = TRUE
        )
      )
  })
}

validation <- assign_height_quintile(validation)

# -----------------------------------------------------------------------------
# 4. Overall prevalence: same common sample, all five definitions
# -----------------------------------------------------------------------------
prevalence_single <- list()
counter <- 1L

for (imp in 1:5) {
  for (sex_value in c("Female", "Male")) {
    dat <- validation %>%
      filter(imputation == imp, as.character(sex) == sex_value)

    design <- make_design(dat)

    for (i in seq_len(nrow(definitions))) {
      spec <- definitions[i, ]
      result <- safe_mean(design, spec$outcome)

      prevalence_single[[counter]] <- tibble(
        imputation = imp,
        sex = sex_value,
        definition_order = spec$definition_order,
        definition = spec$definition,
        analysis_n = nrow(dat),
        estimate = result["estimate"],
        standard_error = result["standard_error"]
      )
      counter <- counter + 1L
    }
  }
}

prevalence_single <- bind_rows(prevalence_single)

prevalence <- pool_table(
  prevalence_single,
  c("sex", "definition_order", "definition")
) %>%
  left_join(
    prevalence_single %>%
      group_by(sex, definition_order, definition) %>%
      summarise(analysis_n = mean(analysis_n), .groups = "drop"),
    by = c("sex", "definition_order", "definition")
  ) %>%
  mutate(
    prevalence_percent = 100 * estimate,
    prevalence_low_percent = 100 * pmax(0, conf_low),
    prevalence_high_percent = 100 * pmin(1, conf_high)
  ) %>%
  arrange(sex, definition_order)

write_csv(
  prevalence,
  file.path(out_dir, "02_overall_prevalence_by_definition.csv")
)

# -----------------------------------------------------------------------------
# 5. Adjusted classification OR per +10 cm height
#    Same model as the existing Paper 1 classification analysis:
#    low status ~ height/10 + ns(age,3) + race + cycle
# -----------------------------------------------------------------------------
height_single <- list()
counter <- 1L

for (imp in 1:5) {
  for (sex_value in c("Female", "Male")) {
    dat <- validation %>%
      filter(imputation == imp, as.character(sex) == sex_value)

    for (i in seq_len(nrow(definitions))) {
      spec <- definitions[i, ]
      fit <- fit_height_effect(dat, spec$outcome)
      term <- "I(height_cm/10)"

      height_single[[counter]] <- tibble(
        imputation = imp,
        sex = sex_value,
        definition_order = spec$definition_order,
        definition = spec$definition,
        analysis_n = nrow(dat),
        estimate = if (inherits(fit, "try-error")) {
          NA_real_
        } else {
          unname(coef(fit)[term])
        },
        standard_error = if (inherits(fit, "try-error")) {
          NA_real_
        } else {
          unname(SE(fit)[term])
        }
      )
      counter <- counter + 1L
    }
  }
}

height_single <- bind_rows(height_single)

height_or <- pool_table(
  height_single,
  c("sex", "definition_order", "definition")
) %>%
  left_join(
    height_single %>%
      group_by(sex, definition_order, definition) %>%
      summarise(analysis_n = mean(analysis_n), .groups = "drop"),
    by = c("sex", "definition_order", "definition")
  ) %>%
  mutate(
    odds_ratio_per_10cm = exp(estimate),
    or_conf_low = exp(conf_low),
    or_conf_high = exp(conf_high),
    percent_change_in_odds = 100 * (odds_ratio_per_10cm - 1)
  ) %>%
  arrange(sex, definition_order)

write_csv(
  height_or,
  file.path(out_dir, "03_height_OR_per10cm.csv")
)

# -----------------------------------------------------------------------------
# 6. Height-quintile prevalence
# -----------------------------------------------------------------------------
quintile_single <- list()
counter <- 1L

for (imp in 1:5) {
  for (sex_value in c("Female", "Male")) {
    dat <- validation %>%
      filter(imputation == imp, as.character(sex) == sex_value)

    design <- make_design(dat)

    for (i in seq_len(nrow(definitions))) {
      spec <- definitions[i, ]

      for (q in paste0("Q", 1:5)) {
        keep <- as.character(dat$height_quintile) == q
        subdesign <- design[keep, ]
        result <- safe_mean(subdesign, spec$outcome)

        quintile_single[[counter]] <- tibble(
          imputation = imp,
          sex = sex_value,
          definition_order = spec$definition_order,
          definition = spec$definition,
          height_quintile = q,
          unweighted_n = sum(keep),
          height_low_cm = min(dat$height_cm[keep], na.rm = TRUE),
          height_high_cm = max(dat$height_cm[keep], na.rm = TRUE),
          estimate = result["estimate"],
          standard_error = result["standard_error"]
        )
        counter <- counter + 1L
      }
    }
  }
}

quintile_single <- bind_rows(quintile_single)

quintile_prev <- pool_table(
  quintile_single,
  c("sex", "definition_order", "definition", "height_quintile")
) %>%
  left_join(
    quintile_single %>%
      group_by(sex, definition_order, definition, height_quintile) %>%
      summarise(
        unweighted_n = mean(unweighted_n),
        height_low_cm = min(height_low_cm),
        height_high_cm = max(height_high_cm),
        .groups = "drop"
      ),
    by = c("sex", "definition_order", "definition", "height_quintile")
  ) %>%
  mutate(
    prevalence_percent = 100 * estimate,
    prevalence_low_percent = 100 * pmax(0, conf_low),
    prevalence_high_percent = 100 * pmin(1, conf_high)
  ) %>%
  arrange(sex, definition_order, height_quintile)

write_csv(
  quintile_prev,
  file.path(out_dir, "05_height_quintile_prevalence.csv")
)

# -----------------------------------------------------------------------------
# 7. BMI phenotype of low vs non-low groups
# -----------------------------------------------------------------------------
bmi_mean_single <- list()
bmi_diff_single <- list()
mean_counter <- 1L
diff_counter <- 1L

for (imp in 1:5) {
  for (sex_value in c("Female", "Male")) {
    dat <- validation %>%
      filter(imputation == imp, as.character(sex) == sex_value)

    for (i in seq_len(nrow(definitions))) {
      spec <- definitions[i, ]
      outcome <- spec$outcome

      # Means by low/non-low group
      for (status in c(0, 1)) {
        subdat <- dat %>% filter(.data[[outcome]] == status)

        if (nrow(subdat) > 0) {
          subdesign <- make_design(subdat)
          result <- safe_mean(subdesign, "bmi")
        } else {
          result <- c(estimate = NA_real_, standard_error = NA_real_)
        }

        bmi_mean_single[[mean_counter]] <- tibble(
          imputation = imp,
          sex = sex_value,
          definition_order = spec$definition_order,
          definition = spec$definition,
          low_status = status,
          status_label = if_else(status == 1, "Low", "Not low"),
          unweighted_n = nrow(subdat),
          estimate = result["estimate"],
          standard_error = result["standard_error"]
        )
        mean_counter <- mean_counter + 1L
      }

      # Difference: Low minus not-low
      formula <- as.formula(paste("bmi ~", outcome))
      fit <- try(
        svyglm(formula, design = make_design(dat)),
        silent = TRUE
      )

      term <- outcome

      bmi_diff_single[[diff_counter]] <- tibble(
        imputation = imp,
        sex = sex_value,
        definition_order = spec$definition_order,
        definition = spec$definition,
        estimate = if (inherits(fit, "try-error")) {
          NA_real_
        } else {
          unname(coef(fit)[term])
        },
        standard_error = if (inherits(fit, "try-error")) {
          NA_real_
        } else {
          unname(SE(fit)[term])
        }
      )
      diff_counter <- diff_counter + 1L
    }
  }
}

bmi_mean_single <- bind_rows(bmi_mean_single)
bmi_diff_single <- bind_rows(bmi_diff_single)

bmi_means <- pool_table(
  bmi_mean_single,
  c("sex", "definition_order", "definition", "low_status", "status_label")
) %>%
  left_join(
    bmi_mean_single %>%
      group_by(sex, definition_order, definition, low_status, status_label) %>%
      summarise(unweighted_n = mean(unweighted_n), .groups = "drop"),
    by = c("sex", "definition_order", "definition", "low_status", "status_label")
  ) %>%
  rename(
    mean_bmi = estimate,
    mean_bmi_se = standard_error,
    mean_bmi_low95 = conf_low,
    mean_bmi_high95 = conf_high
  )

bmi_diff <- pool_table(
  bmi_diff_single,
  c("sex", "definition_order", "definition")
) %>%
  rename(
    bmi_difference_low_minus_not_low = estimate,
    difference_se = standard_error,
    difference_low95 = conf_low,
    difference_high95 = conf_high,
    difference_p_value = p_value
  )

bmi_phenotype <- bmi_means %>%
  select(
    sex, definition_order, definition, low_status, status_label,
    unweighted_n, mean_bmi, mean_bmi_se, mean_bmi_low95, mean_bmi_high95
  ) %>%
  left_join(
    bmi_diff %>%
      select(
        sex, definition_order, definition,
        bmi_difference_low_minus_not_low,
        difference_se, difference_low95, difference_high95,
        difference_p_value
      ),
    by = c("sex", "definition_order", "definition")
  ) %>%
  arrange(sex, definition_order, low_status)

write_csv(
  bmi_phenotype,
  file.path(out_dir, "06_BMI_phenotype_low_vs_not_low.csv")
)

# -----------------------------------------------------------------------------
# 8. Main compact summary table
# -----------------------------------------------------------------------------
bmi_wide <- bmi_means %>%
  select(sex, definition_order, definition, status_label, mean_bmi) %>%
  pivot_wider(
    names_from = status_label,
    values_from = mean_bmi,
    names_prefix = "mean_BMI_"
  )

summary_table <- prevalence %>%
  select(
    sex, definition_order, definition, analysis_n,
    prevalence_percent, prevalence_low_percent, prevalence_high_percent
  ) %>%
  left_join(
    height_or %>%
      select(
        sex, definition_order, definition,
        odds_ratio_per_10cm, or_conf_low, or_conf_high
      ),
    by = c("sex", "definition_order", "definition")
  ) %>%
  left_join(
    bmi_wide,
    by = c("sex", "definition_order", "definition")
  ) %>%
  left_join(
    bmi_diff %>%
      select(
        sex, definition_order, definition,
        bmi_difference_low_minus_not_low,
        difference_low95, difference_high95
      ),
    by = c("sex", "definition_order", "definition")
  ) %>%
  arrange(sex, definition_order)

write_csv(
  summary_table,
  file.path(out_dir, "07_MAIN_summary_all_definitions.csv")
)

# -----------------------------------------------------------------------------
# 9. Save analysis-ready participant data used for this sensitivity only
# -----------------------------------------------------------------------------
write_csv(
  validation %>%
    select(
      SEQN, imputation, cycle, sex, race, age, height_cm, height_m,
      bmi, alm_kg, almi, weight_kg_for_analysis, weight_kg_reconstructed, alm_weight_pct,
      low_conditional_p5, low_ewgsop2_recalc, low_fnih_recalc,
      low_weight_def1, low_weight_def2,
      height_quintile, wtmec_pool, strata, psu, strata_pool, psu_pool
    ),
  file.path(out_dir, "08_analysis_ready_validation_data.csv.gz")
)

cat("\n============================================================\n")
cat("ALM/weight sensitivity analysis completed successfully.\n")
cat("============================================================\n\n")

cat("Main summary:\n")
print(summary_table)

cat("\nResults folder:\n", out_dir, "\n\n")
cat("Please upload these files back to ChatGPT:\n")
cat("  00_definition_table.csv\n")
cat("  02_overall_prevalence_by_definition.csv\n")
cat("  03_height_OR_per10cm.csv\n")
cat("  05_height_quintile_prevalence.csv\n")
cat("  06_BMI_phenotype_low_vs_not_low.csv\n")
cat("  07_MAIN_summary_all_definitions.csv\n")
