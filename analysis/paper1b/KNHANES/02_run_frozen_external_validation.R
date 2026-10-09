# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(survey)
  library(gamlss)
  library(gamlss.dist)
})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
raw_all <- file.path(root, "raw_all_sav")
raw_dxa <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/Knhanes\u6570\u636e\u5e93"
paper_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/First Paper/Final Paper/\u65b0\u5efa\u6587\u4ef6\u5939"
scorer_root <- file.path(paper_root, "nhanes-alm-reference-scorer-v2.1.0")
age_only_root <- file.path(paper_root, "age_only_comparator_20260820", "models")

dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "logs"), showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(root, "logs", "02_run_log.txt"), open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("KNHANES frozen external validation started:", format(Sys.time()), "\n")

years <- 2008:2011
all_files <- setNames(
  file.path(raw_all, c("HN08_ALL.sav", "HN09_all.sav", "HN10_all.sav", "HN11_all.sav")),
  years
)
dxa_files <- setNames(file.path(raw_dxa, sprintf("hn%02d_dxa.sav", 8:11)), years)
weight_map <- c(`2008` = "wt_ex1", `2009` = "wt_dw", `2010` = "wt_itvex", `2011` = "wt_ex1")
pool_numerator <- c(`2008` = 108, `2009` = 199, `2010` = 192, `2011` = 80)

all_wanted <- c("id", "psu", "sex", "age", "kstrata", "he_ht", "he_wt", "he_bmi",
                "wt_ex1", "wt_dw", "wt_itvex")
dxa_wanted <- c("id", "dw_lrm_ln", "dw_rrm_ln", "dw_llg_ln", "dw_rlg_ln",
                "dw_lrm_bmc", "dw_rrm_bmc", "dw_llg_bmc", "dw_rlg_bmc")

read_selected <- function(path, wanted_lower) {
  header <- read_sav(path, n_max = 1)
  selected <- names(header)[tolower(names(header)) %in% wanted_lower]
  missing <- setdiff(wanted_lower, tolower(selected))
  required_missing <- setdiff(missing, c("wt_ex1", "wt_dw", "wt_itvex"))
  if (length(required_missing)) {
    stop("Missing expected variables in ", basename(path), ": ", paste(required_missing, collapse = ", "))
  }
  x <- read_sav(path, col_select = all_of(selected))
  names(x) <- tolower(names(x))
  as.data.frame(x)
}

as_clean_numeric <- function(x) {
  suppressWarnings(as.numeric(zap_labels(zap_missing(x))))
}

not_special <- function(x) {
  is.finite(x) &
    abs(x - 9999) > 0.02 & abs(x - 9999.9) > 0.02 &
    abs(x - 9999.99) > 0.02 & abs(x - 99999) > 0.02
}

year_data <- list()
flow_rows <- list()
for (year in years) {
  cat("Reading and harmonizing", year, "\n")
  all_dat <- read_selected(all_files[[as.character(year)]], all_wanted)
  dxa_dat <- read_selected(dxa_files[[as.character(year)]], dxa_wanted)
  all_dat$id <- trimws(as.character(all_dat$id))
  dxa_dat$id <- trimws(as.character(dxa_dat$id))

  if (anyDuplicated(all_dat$id)) stop("Duplicate ALL IDs in ", year)
  if (anyDuplicated(dxa_dat$id)) stop("Duplicate DXA IDs in ", year)

  matched <- dxa_dat$id %in% all_dat$id
  dat <- merge(dxa_dat, all_dat, by = "id", all.x = TRUE, sort = FALSE)
  dat$year <- year
  numeric_vars <- setdiff(names(dat), c("id", "psu"))
  for (v in numeric_vars) dat[[v]] <- as_clean_numeric(dat[[v]])

  ln_vars <- c("dw_lrm_ln", "dw_rrm_ln", "dw_llg_ln", "dw_rlg_ln")
  bmc_vars <- c("dw_lrm_bmc", "dw_rrm_bmc", "dw_llg_bmc", "dw_rlg_bmc")
  valid_ln_matrix <- vapply(dat[ln_vars], function(x) not_special(x) & x > 0 & x < 20000, logical(nrow(dat)))
  valid_bmc_matrix <- vapply(dat[bmc_vars], function(x) not_special(x) & x > 0 & x < 900, logical(nrow(dat)))
  dat$valid_four_ln <- rowSums(valid_ln_matrix) == 4L
  dat$valid_four_bmc <- rowSums(valid_bmc_matrix) == 4L

  dat$alm_raw_kg <- NA_real_
  # KNHANES LN is the limb lean/muscle compartment itself: MS = LN + FT.
  # BMC is reported separately and must not be subtracted from LN.  The 0.946
  # cross-calibration is applied to the sum of the four limb LN values.
  component_ok <- dat$valid_four_ln
  dat$alm_raw_kg[component_ok] <-
    rowSums(dat[component_ok, ln_vars, drop = FALSE]) / 1000
  dat$valid_alm_raw <- is.finite(dat$alm_raw_kg) & dat$alm_raw_kg > 5 & dat$alm_raw_kg < 70
  dat$alm_kg <- ifelse(dat$valid_alm_raw, dat$alm_raw_kg * 0.946, NA_real_)

  dat$height_cm <- ifelse(is.finite(dat$he_ht) & dat$he_ht >= 120 & dat$he_ht <= 220, dat$he_ht, NA_real_)
  dat$height_m <- dat$height_cm / 100
  dat$weight_kg <- ifelse(is.finite(dat$he_wt) & dat$he_wt >= 20 & dat$he_wt <= 300, dat$he_wt, NA_real_)
  calculated_bmi <- dat$weight_kg / dat$height_m^2
  dat$bmi <- ifelse(is.finite(dat$he_bmi) & dat$he_bmi >= 10 & dat$he_bmi <= 80,
                    dat$he_bmi, calculated_bmi)
  dat$bmi[!is.finite(dat$bmi) | dat$bmi < 10 | dat$bmi > 80] <- NA_real_
  dat$sex_label <- ifelse(dat$sex == 1, "Male", ifelse(dat$sex == 2, "Female", NA_character_))
  dat$wt_year <- dat[[weight_map[[as.character(year)]]]]
  dat$wt_year[!is.finite(dat$wt_year) | dat$wt_year <= 0] <- NA_real_

  design_complete <- is.finite(dat$wt_year) & !is.na(dat$psu) & nzchar(dat$psu) & is.finite(dat$kstrata)
  base_target <- dat$age >= 19 & dat$age <= 69 & !is.na(dat$sex_label) &
    dat$valid_alm_raw & is.finite(dat$height_m)
  core_target <- base_target & is.finite(dat$bmi) & design_complete

  flow_rows[[as.character(year)]] <- data.frame(
    year = year,
    dxa_n = nrow(dxa_dat),
    matched_all_n = sum(matched),
    valid_four_limb_ln_n = sum(dat$valid_four_ln),
    valid_four_limb_bmc_n = sum(dat$valid_four_bmc),
    valid_four_limb_ln_alm_n = sum(dat$valid_alm_raw),
    valid_height_n = sum(is.finite(dat$height_m)),
    valid_bmi_n = sum(is.finite(dat$bmi)),
    age_19_69_alm_height_n = sum(base_target, na.rm = TRUE),
    age_19_69_core_design_n = sum(core_target, na.rm = TRUE),
    male_core_n = sum(core_target & dat$sex_label == "Male", na.rm = TRUE),
    female_core_n = sum(core_target & dat$sex_label == "Female", na.rm = TRUE),
    age_19_59_core_design_n = sum(core_target & dat$age <= 59, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  year_data[[as.character(year)]] <- dat
  rm(all_dat, dxa_dat, dat)
  gc()
}

flow <- do.call(rbind, flow_rows)
write.csv(flow, file.path(root, "tables", "01_sample_flow.csv"), row.names = FALSE)
print(flow, row.names = FALSE)

dat <- bind_rows(year_data)
rm(year_data)
dat$pool_coefficient <- unname(pool_numerator[as.character(dat$year)]) / 579
dat$wt_pool <- dat$wt_year * dat$pool_coefficient
dat$strata_pool <- interaction(dat$year, dat$kstrata, drop = TRUE)
dat$psu_pool <- interaction(dat$year, dat$psu, drop = TRUE)
dat$row_key <- seq_len(nrow(dat))

score_candidates <- is.finite(dat$age) & dat$age >= 19 & dat$age <= 69 &
  !is.na(dat$sex_label) & dat$valid_alm_raw & is.finite(dat$height_m)

cat("Loading frozen age+height scorer and models\n")
source(file.path(scorer_root, "R", "score_conditional_alm.R"), local = .GlobalEnv)
h_model_files <- file.path(
  scorer_root, "models",
  c("corrected_H_18_69_Female_bundles.rds", "corrected_H_18_69_Male_bundles.rds")
)
h_bundles <- load_alm_reference(h_model_files)

score_input <- data.frame(
  row_key = dat$row_key[score_candidates],
  sex = dat$sex_label[score_candidates],
  age = dat$age[score_candidates],
  height_m = dat$height_m[score_candidates],
  alm_kg = dat$alm_kg[score_candidates],
  stringsAsFactors = FALSE
)
h_scored <- score_conditional_alm(score_input, bundles = h_bundles, warn = TRUE)$pooled
dat$z_conditional <- NA_real_
dat$z_conditional[match(h_scored$row_key, dat$row_key)] <- h_scored$model_averaged_z
dat$scorer_valid <- FALSE
dat$scorer_valid[match(h_scored$row_key, dat$row_key)] <- h_scored$score_valid
dat$scorer_caution <- FALSE
dat$scorer_caution[match(h_scored$row_key, dat$row_key)] <- h_scored$use_caution
cat("Conditional score complete; finite z:", sum(is.finite(dat$z_conditional)), "\n")

cat("Loading and applying frozen age-only comparator\n")
age_only_bundles <- c(
  readRDS(file.path(age_only_root, "age_only_18_69_Female_bundles.rds")),
  readRDS(file.path(age_only_root, "age_only_18_69_Male_bundles.rds"))
)
age_scored <- score_conditional_alm(score_input, bundles = age_only_bundles, warn = FALSE)$pooled
dat$z_age_only <- NA_real_
dat$z_age_only[match(age_scored$row_key, dat$row_key)] <- age_scored$model_averaged_z
cat("Age-only score complete; finite z:", sum(is.finite(dat$z_age_only)), "\n")

dat$almi <- dat$alm_kg / dat$height_m^2
dat$alm_bmi <- dat$alm_kg / dat$bmi
dat$low_p5 <- is.finite(dat$z_conditional) & dat$z_conditional < qnorm(0.05)
dat$low_p10 <- is.finite(dat$z_conditional) & dat$z_conditional < qnorm(0.10)
dat$low_ewgsop2 <- ifelse(dat$sex_label == "Male", dat$almi < 7.0,
                          ifelse(dat$sex_label == "Female", dat$almi < 5.5, NA))
dat$low_fnih <- ifelse(dat$sex_label == "Male", dat$alm_bmi < 0.789,
                       ifelse(dat$sex_label == "Female", dat$alm_bmi < 0.512, NA))
dat$height10 <- dat$height_cm / 10

analytic <- dat$age >= 19 & dat$age <= 69 & !is.na(dat$sex_label) &
  is.finite(dat$z_conditional) & is.finite(dat$z_age_only) &
  is.finite(dat$bmi) & is.finite(dat$wt_pool) & dat$wt_pool > 0 &
  !is.na(dat$psu_pool) & !is.na(dat$strata_pool)
analysis_dat <- dat[analytic, ]

flow$scorer_domain_valid_n <- vapply(years, function(y) {
  sum(dat$year == y & score_candidates & dat$scorer_valid, na.rm = TRUE)
}, numeric(1))
flow$scorer_caution_n <- vapply(years, function(y) {
  sum(dat$year == y & score_candidates & dat$scorer_caution, na.rm = TRUE)
}, numeric(1))
flow$final_common_analytic_19_69_n <- vapply(years, function(y) sum(analysis_dat$year == y), numeric(1))
write.csv(flow, file.path(root, "tables", "01_sample_flow.csv"), row.names = FALSE)

saveRDS(
  analysis_dat[, c("id", "year", "sex_label", "age", "height_cm", "height_m", "weight_kg", "bmi",
                   "alm_raw_kg", "alm_kg", "almi", "alm_bmi", "z_conditional", "z_age_only",
                   "low_p5", "low_p10", "low_ewgsop2", "low_fnih", "wt_year", "wt_pool",
                   "kstrata", "psu", "strata_pool", "psu_pool", "scorer_caution")],
  file.path(root, "data", "knhanes_2008_2011_analysis_19_69.rds"),
  compress = "xz"
)

write.csv(
  data.frame(
    item = c("DXA records", "19-69 ALM+height candidates", "frozen scorer valid",
             "frozen scorer outside observed domain", "recommended-range caution", "final common analytic sample"),
    n = c(nrow(dat), sum(score_candidates), sum(dat$scorer_valid & score_candidates),
          sum(score_candidates & !dat$scorer_valid), sum(dat$scorer_caution & score_candidates), nrow(analysis_dat))
  ),
  file.path(root, "tables", "01b_overall_flow.csv"), row.names = FALSE
)

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  x[pmin(length(x), findInterval(probs, cumsum(w) / sum(w)) + 1L)]
}

qc_summary_rows <- list()
for (year in years) for (sex_value in c("Female", "Male")) {
  z <- dat[dat$year == year & dat$sex_label == sex_value & dat$age >= 19 & dat$age <= 69 &
             is.finite(dat$alm_raw_kg) & is.finite(dat$wt_year) & dat$wt_year > 0, ]
  for (variable in c("alm_raw_kg", "alm_kg", "height_cm", "bmi")) {
    x <- z[[variable]]; w <- z$wt_year
    ok <- is.finite(x) & is.finite(w) & w > 0
    q <- weighted_quantile(x, w, c(.05, .10, .50, .90, .95))
    m <- weighted.mean(x[ok], w[ok])
    sd_w <- sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok]))
    qc_summary_rows[[length(qc_summary_rows) + 1L]] <- data.frame(
      year = year, sex = sex_value, variable = variable, n = sum(ok), mean = m, sd = sd_w,
      p05 = q[1], p10 = q[2], median = q[3], p90 = q[4], p95 = q[5]
    )
  }
}
qc_summary <- bind_rows(qc_summary_rows)
write.csv(qc_summary, file.path(root, "tables", "02_qc_weighted_distributions.csv"), row.names = FALSE)

make_design <- function(x) {
  svydesign(ids = ~psu_pool, strata = ~strata_pool, weights = ~wt_pool,
            nest = TRUE, data = x)
}

design_19_69 <- make_design(analysis_dat)

safe_prop <- function(design, variable) {
  values <- design$variables[[variable]]
  values <- values[!is.na(values)]
  if (!length(values)) return(c(estimate = NA_real_, lcl = NA_real_, ucl = NA_real_))
  if (all(!values)) return(c(estimate = 0, lcl = 0, ucl = 0))
  if (all(values)) return(c(estimate = 1, lcl = 1, ucl = 1))
  f <- as.formula(paste0("~", variable))
  ans <- try(svyciprop(f, design, method = "logit", level = 0.95, na.rm = TRUE), silent = TRUE)
  if (inherits(ans, "try-error")) {
    ans <- svymean(f, design, na.rm = TRUE)
  }
  c(estimate = as.numeric(coef(ans))[1], lcl = as.numeric(confint(ans))[1], ucl = as.numeric(confint(ans))[2])
}

calibration_one <- function(design, data_frame, scope, sex_value, year_value = NA_integer_) {
  d <- subset(design, sex_label == sex_value)
  idx <- data_frame$sex_label == sex_value
  if (is.finite(year_value)) {
    d <- subset(d, year == year_value)
    idx <- idx & data_frame$year == year_value
  }
  mz <- svymean(~z_conditional, d, na.rm = TRUE)
  vz <- svyvar(~z_conditional, d, na.rm = TRUE)
  p5 <- safe_prop(d, "low_p5")
  p10 <- safe_prop(d, "low_p10")
  data.frame(
    scope = scope, year = ifelse(is.finite(year_value), year_value, NA), sex = sex_value,
    n = sum(idx), weighted_mean_z = coef(mz)[1], mean_z_lcl = confint(mz)[1], mean_z_ucl = confint(mz)[2],
    weighted_sd_z = sqrt(as.numeric(vz)[1]),
    p5_percent = 100 * p5[1], p5_lcl = 100 * p5[2], p5_ucl = 100 * p5[3],
    p10_percent = 100 * p10[1], p10_lcl = 100 * p10[2], p10_ucl = 100 * p10[3],
    stringsAsFactors = FALSE
  )
}

calibration_rows <- list()
for (sex_value in c("Female", "Male")) {
  calibration_rows[[length(calibration_rows) + 1L]] <- calibration_one(
    design_19_69, analysis_dat, "Primary 19-69 pooled", sex_value
  )
  for (year in years) calibration_rows[[length(calibration_rows) + 1L]] <- calibration_one(
    design_19_69, analysis_dat, "Primary 19-69 by year", sex_value, year
  )
}
analysis_19_59 <- analysis_dat[analysis_dat$age <= 59, ]
design_19_59 <- make_design(analysis_19_59)
for (sex_value in c("Female", "Male")) calibration_rows[[length(calibration_rows) + 1L]] <- calibration_one(
  design_19_59, analysis_19_59, "Sensitivity 19-59 pooled", sex_value
)
analysis_60_69 <- analysis_dat[analysis_dat$age >= 60, ]
design_60_69 <- make_design(analysis_60_69)
for (sex_value in c("Female", "Male")) calibration_rows[[length(calibration_rows) + 1L]] <- calibration_one(
  design_60_69, analysis_60_69, "Age extension 60-69", sex_value
)
calibration <- bind_rows(calibration_rows)
write.csv(calibration, file.path(root, "tables", "03_frozen_calibration.csv"), row.names = FALSE)

extract_linear <- function(fit, term = "height10") {
  b <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  data.frame(estimate = b, se = se, lcl = b - qnorm(.975) * se, ucl = b + qnorm(.975) * se,
             p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}

height_slope_rows <- list()
run_height_slopes <- function(design, scope, year_value = NA_integer_) {
  out <- list()
  for (sex_value in c("Female", "Male")) for (score in c("z_age_only", "z_conditional")) {
    d <- subset(design, sex_label == sex_value)
    formula_text <- paste0(score, " ~ height10 + age")
    if (!is.finite(year_value)) formula_text <- paste0(formula_text, " + factor(year)")
    fit <- svyglm(as.formula(formula_text), design = d)
    e <- extract_linear(fit)
    out[[length(out) + 1L]] <- data.frame(scope = scope, year = ifelse(is.finite(year_value), year_value, NA),
                                          sex = sex_value, score = score, e)
  }
  bind_rows(out)
}
height_slope_rows[[1]] <- run_height_slopes(design_19_69, "Primary 19-69 pooled")
height_slope_rows[[2]] <- run_height_slopes(design_19_59, "Sensitivity 19-59 pooled")
height_slope_rows[[3]] <- run_height_slopes(design_60_69, "Age extension 60-69")
for (year_value in years) {
  dy <- subset(design_19_69, year == year_value)
  height_slope_rows[[length(height_slope_rows) + 1L]] <- run_height_slopes(dy, "Primary 19-69 by year", year_value)
}
height_slopes <- bind_rows(height_slope_rows)
write.csv(height_slopes, file.path(root, "tables", "04_z_height_slopes.csv"), row.names = FALSE)

extract_or <- function(fit, term = "height10") {
  x <- extract_linear(fit, term)
  data.frame(log_or = x$estimate, se = x$se, or_per_10cm = exp(x$estimate),
             or_lcl = exp(x$lcl), or_ucl = exp(x$ucl), p_value = x$p_value)
}

outcome_labels <- c(low_p5 = "Conditional P5", low_p10 = "Conditional P10",
                    low_ewgsop2 = "EWGSOP2", low_fnih = "FNIH")
run_or_models <- function(design, scope, year_value = NA_integer_) {
  out <- list()
  for (sex_value in c("Female", "Male")) for (outcome in names(outcome_labels)) {
    d <- subset(design, sex_label == sex_value)
    formula_text <- paste0(outcome, " ~ height10 + age")
    if (!is.finite(year_value)) formula_text <- paste0(formula_text, " + factor(year)")
    fit <- try(svyglm(as.formula(formula_text), design = d, family = quasibinomial()), silent = TRUE)
    e <- if (inherits(fit, "try-error") || !"height10" %in% names(coef(fit))) {
      data.frame(log_or = NA, se = NA, or_per_10cm = NA, or_lcl = NA, or_ucl = NA, p_value = NA)
    } else extract_or(fit)
    out[[length(out) + 1L]] <- data.frame(scope = scope, year = ifelse(is.finite(year_value), year_value, NA),
                                          sex = sex_value, definition = outcome_labels[[outcome]], e)
  }
  bind_rows(out)
}

or_rows <- list(
  run_or_models(design_19_69, "Primary 19-69 pooled"),
  run_or_models(design_19_59, "Sensitivity 19-59 pooled"),
  run_or_models(design_60_69, "Age extension 60-69")
)
for (year_value in years) {
  dy <- subset(design_19_69, year == year_value)
  or_rows[[length(or_rows) + 1L]] <- run_or_models(dy, "Primary 19-69 by year", year_value)
}
height_or <- bind_rows(or_rows)
write.csv(height_or, file.path(root, "tables", "05_low_alm_height_or.csv"), row.names = FALSE)

# Sex-specific survey-weighted height quintiles in the primary 19-69 sample.
analysis_dat$height_quintile <- NA_character_
quintile_cut_rows <- list()
for (sex_value in c("Female", "Male")) {
  idx <- analysis_dat$sex_label == sex_value
  cuts <- weighted_quantile(analysis_dat$height_cm[idx], analysis_dat$wt_pool[idx], c(.2, .4, .6, .8))
  analysis_dat$height_quintile[idx] <- as.character(cut(
    analysis_dat$height_cm[idx], breaks = c(-Inf, cuts, Inf), labels = paste0("Q", 1:5), include.lowest = TRUE
  ))
  quintile_cut_rows[[sex_value]] <- data.frame(sex = sex_value, boundary = c("Q20", "Q40", "Q60", "Q80"), height_cm = cuts)
}
analysis_dat$height_quintile <- factor(analysis_dat$height_quintile, levels = paste0("Q", 1:5))
write.csv(bind_rows(quintile_cut_rows), file.path(root, "tables", "06a_height_quintile_cutpoints.csv"), row.names = FALSE)
design_19_69_q <- make_design(analysis_dat)

prev_rows <- list()
for (sex_value in c("Female", "Male")) for (q in paste0("Q", 1:5)) for (outcome in names(outcome_labels)) {
  d <- subset(design_19_69_q, sex_label == sex_value & height_quintile == q)
  p <- safe_prop(d, outcome)
  n_cell <- sum(analysis_dat$sex_label == sex_value & analysis_dat$height_quintile == q)
  prev_rows[[length(prev_rows) + 1L]] <- data.frame(
    sex = sex_value, height_quintile = q, definition = outcome_labels[[outcome]], n = n_cell,
    prevalence_percent = 100 * p[1], lcl = 100 * p[2], ucl = 100 * p[3]
  )
}
height_prev <- bind_rows(prev_rows)
write.csv(height_prev, file.path(root, "tables", "06b_height_quintile_prevalence.csv"), row.names = FALSE)

discordance_rows <- list()
for (sex_value in c("Female", "Male")) for (comparator in c("low_fnih", "low_ewgsop2")) {
  comp_label <- outcome_labels[[comparator]]
  sex_design <- subset(design_19_69_q, sex_label == sex_value)
  sex_data <- analysis_dat[analysis_dat$sex_label == sex_value, ]
  conditional <- sex_data$low_p5
  comp <- sex_data[[comparator]]
  group <- ifelse(!conditional & !comp, "Neither",
                  ifelse(!conditional & comp, paste0(comp_label, " only"),
                         ifelse(conditional & !comp, "Conditional only", "Both")))
  analysis_dat$discordance_tmp <- NA_character_
  analysis_dat$discordance_tmp[analysis_dat$sex_label == sex_value] <- group
  dtmp <- make_design(analysis_dat)
  dsex <- subset(dtmp, sex_label == sex_value)
  for (g in c("Neither", paste0(comp_label, " only"), "Conditional only", "Both")) {
    n_cell <- sum(group == g)
    indicator_name <- "discordance_indicator"
    dsex$variables[[indicator_name]] <- dsex$variables$discordance_tmp == g
    p <- safe_prop(dsex, indicator_name)
    dgroup <- subset(dsex, discordance_tmp == g)
    mean_vars <- c("age", "height_cm", "bmi", "alm_kg", "almi")
    means <- if (n_cell > 0) coef(svymean(as.formula(paste("~", paste(mean_vars, collapse = "+"))), dgroup, na.rm = TRUE)) else rep(NA, 5)
    discordance_rows[[length(discordance_rows) + 1L]] <- data.frame(
      sex = sex_value, comparison = paste("Conditional P5 vs", comp_label), group = g,
      unweighted_n = n_cell, weighted_percent = 100 * p[1], percent_lcl = 100 * p[2], percent_ucl = 100 * p[3],
      mean_age = means["age"], mean_height_cm = means["height_cm"], mean_bmi = means["bmi"],
      mean_alm_kg = means["alm_kg"], mean_almi = means["almi"]
    )
  }
}
discordance <- bind_rows(discordance_rows)
write.csv(discordance, file.path(root, "tables", "07_classification_discordance.csv"), row.names = FALSE)

write.csv(
  data.frame(
    component = c("R", "haven", "survey", "gamlss", "gamlss.dist"),
    version = c(as.character(getRversion()), as.character(packageVersion("haven")),
                as.character(packageVersion("survey")), as.character(packageVersion("gamlss")),
                as.character(packageVersion("gamlss.dist")))
  ),
  file.path(root, "tables", "99_software_versions.csv"), row.names = FALSE
)

cat("KNHANES frozen external validation completed:", format(Sys.time()), "\n")
cat("Final primary analytic N:", nrow(analysis_dat), "\n")
cat("Final 19-59 sensitivity N:", nrow(analysis_19_59), "\n")
cat("Final 60-69 age-extension N:", nrow(analysis_60_69), "\n")
print(calibration[calibration$scope != "Primary 19-69 by year", ], row.names = FALSE)
print(height_slopes[height_slopes$scope != "Primary 19-69 by year", ], row.names = FALSE)
print(height_or[height_or$scope != "Primary 19-69 by year", ], row.names = FALSE)
