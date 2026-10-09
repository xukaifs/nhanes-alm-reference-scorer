# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages(library(survey))

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
out <- file.path(root, "tables_official_weight")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

dat <- readRDS(file.path(root, "data", "knhanes_2008_2011_analysis_19_69.rds"))
pool_numerator <- c(`2008` = 108, `2009` = 199, `2010` = 192, `2011` = 80)
dat$wt_pool_equal_year <- dat$wt_year / 4
dat$pool_coefficient <- unname(pool_numerator[as.character(dat$year)]) / 579
dat$wt_pool <- dat$wt_year * dat$pool_coefficient
dat$height10 <- dat$height_cm / 10
dat$height_recommended <- ifelse(
  dat$sex_label == "Female",
  dat$height_m >= 1.458 & dat$height_m <= 1.781,
  dat$height_m >= 1.570 & dat$height_m <= 1.943
)

write.csv(
  data.frame(year = as.integer(names(pool_numerator)), survey_areas = unname(pool_numerator),
             denominator = 579, coefficient = unname(pool_numerator) / 579),
  file.path(out, "00_official_pooling_coefficients.csv"), row.names = FALSE
)

saveRDS(dat, file.path(root, "data", "knhanes_2008_2011_analysis_19_69_official_weight.rds"), compress = "xz")

make_design <- function(x) svydesign(
  ids = ~psu_pool, strata = ~strata_pool, weights = ~wt_pool,
  nest = TRUE, data = x
)

safe_prop <- function(design, variable) {
  values <- design$variables[[variable]]
  values <- values[!is.na(values)]
  if (!length(values)) return(c(estimate = NA_real_, lcl = NA_real_, ucl = NA_real_))
  if (all(!values)) return(c(estimate = 0, lcl = 0, ucl = 0))
  if (all(values)) return(c(estimate = 1, lcl = 1, ucl = 1))
  f <- as.formula(paste0("~", variable))
  ans <- try(svyciprop(f, design, method = "logit", level = .95, na.rm = TRUE), silent = TRUE)
  if (inherits(ans, "try-error")) ans <- svymean(f, design, na.rm = TRUE)
  c(estimate = as.numeric(coef(ans))[1], lcl = as.numeric(confint(ans))[1], ucl = as.numeric(confint(ans))[2])
}

scope_data <- list(
  `Primary 19-69 pooled` = dat,
  `Sensitivity 19-59 pooled` = dat[dat$age <= 59, ],
  `Age extension 60-69` = dat[dat$age >= 60, ],
  `Recommended height 19-69` = dat[dat$height_recommended, ]
)

calibration_one <- function(x, scope, sex_value) {
  d <- subset(make_design(x), sex_label == sex_value)
  mz <- svymean(~z_conditional, d, na.rm = TRUE)
  p5 <- safe_prop(d, "low_p5")
  p10 <- safe_prop(d, "low_p10")
  data.frame(
    scope = scope, year = NA_integer_, sex = sex_value, n = sum(x$sex_label == sex_value),
    weighted_mean_z = coef(mz)[1], mean_z_lcl = confint(mz)[1], mean_z_ucl = confint(mz)[2],
    weighted_sd_z = sqrt(as.numeric(svyvar(~z_conditional, d, na.rm = TRUE))[1]),
    p5_percent = 100 * p5[1], p5_lcl = 100 * p5[2], p5_ucl = 100 * p5[3],
    p10_percent = 100 * p10[1], p10_lcl = 100 * p10[2], p10_ucl = 100 * p10[3]
  )
}

cal_rows <- list()
for (scope in names(scope_data)) for (sex_value in c("Female", "Male")) {
  cal_rows[[length(cal_rows) + 1L]] <- calibration_one(scope_data[[scope]], scope, sex_value)
}
for (year_value in 2008:2011) for (sex_value in c("Female", "Male")) {
  x <- dat[dat$year == year_value, ]
  row <- calibration_one(x, "Primary 19-69 by year", sex_value)
  row$year <- year_value
  cal_rows[[length(cal_rows) + 1L]] <- row
}
calibration <- do.call(rbind, cal_rows)
write.csv(calibration, file.path(out, "01_calibration.csv"), row.names = FALSE)

extract_linear <- function(fit, term = "height10") {
  b <- coef(fit)[term]
  se <- sqrt(vcov(fit)[term, term])
  data.frame(estimate = b, se = se, lcl = b - qnorm(.975) * se, ucl = b + qnorm(.975) * se,
             p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}

run_slopes <- function(x, scope, year_adjust = TRUE) {
  design <- make_design(x)
  rows <- list()
  for (sex_value in c("Female", "Male")) for (score in c("z_age_only", "z_conditional")) {
    d <- subset(design, sex_label == sex_value)
    formula_text <- paste0(score, " ~ height10 + age", if (year_adjust) " + factor(year)" else "")
    e <- extract_linear(svyglm(as.formula(formula_text), design = d))
    rows[[length(rows) + 1L]] <- data.frame(scope = scope, year = NA_integer_, sex = sex_value, score = score, e)
  }
  do.call(rbind, rows)
}

slope_rows <- lapply(names(scope_data), function(scope) run_slopes(scope_data[[scope]], scope, TRUE))
for (year_value in 2008:2011) {
  z <- run_slopes(dat[dat$year == year_value, ], "Primary 19-69 by year", FALSE)
  z$year <- year_value
  slope_rows[[length(slope_rows) + 1L]] <- z
}
height_slopes <- do.call(rbind, slope_rows)
write.csv(height_slopes, file.path(out, "02_z_height_slopes.csv"), row.names = FALSE)

outcome_labels <- c(low_p5 = "Conditional P5", low_p10 = "Conditional P10",
                    low_ewgsop2 = "EWGSOP2", low_fnih = "FNIH")
run_or <- function(x, scope, year_adjust = TRUE) {
  design <- make_design(x)
  rows <- list()
  for (sex_value in c("Female", "Male")) for (outcome in names(outcome_labels)) {
    d <- subset(design, sex_label == sex_value)
    formula_text <- paste0(outcome, " ~ height10 + age", if (year_adjust) " + factor(year)" else "")
    fit <- svyglm(as.formula(formula_text), design = d, family = quasibinomial())
    e <- extract_linear(fit)
    rows[[length(rows) + 1L]] <- data.frame(
      scope = scope, year = NA_integer_, sex = sex_value, definition = outcome_labels[[outcome]],
      log_or = e$estimate, se = e$se, or_per_10cm = exp(e$estimate),
      or_lcl = exp(e$lcl), or_ucl = exp(e$ucl), p_value = e$p_value
    )
  }
  do.call(rbind, rows)
}

or_rows <- lapply(names(scope_data), function(scope) run_or(scope_data[[scope]], scope, TRUE))
for (year_value in 2008:2011) {
  z <- run_or(dat[dat$year == year_value, ], "Primary 19-69 by year", FALSE)
  z$year <- year_value
  or_rows[[length(or_rows) + 1L]] <- z
}
height_or <- do.call(rbind, or_rows)
write.csv(height_or, file.path(out, "03_low_alm_height_or.csv"), row.names = FALSE)

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  ord <- order(x); x <- x[ord]; w <- w[ord]
  x[pmin(length(x), findInterval(probs, cumsum(w) / sum(w)) + 1L)]
}

dat$height_quintile <- NA_character_
cut_rows <- list()
for (sex_value in c("Female", "Male")) {
  idx <- dat$sex_label == sex_value
  cuts <- weighted_quantile(dat$height_cm[idx], dat$wt_pool[idx], c(.2, .4, .6, .8))
  dat$height_quintile[idx] <- as.character(cut(
    dat$height_cm[idx], c(-Inf, cuts, Inf), labels = paste0("Q", 1:5), include.lowest = TRUE
  ))
  cut_rows[[sex_value]] <- data.frame(sex = sex_value, boundary = c("Q20", "Q40", "Q60", "Q80"), height_cm = cuts)
}
dat$height_quintile <- factor(dat$height_quintile, levels = paste0("Q", 1:5))
write.csv(do.call(rbind, cut_rows), file.path(out, "04_height_quintile_cutpoints.csv"), row.names = FALSE)
design_q <- make_design(dat)

prev_rows <- list()
for (sex_value in c("Female", "Male")) for (quintile in paste0("Q", 1:5)) for (outcome in names(outcome_labels)) {
  d <- subset(design_q, sex_label == sex_value & height_quintile == quintile)
  p <- safe_prop(d, outcome)
  prev_rows[[length(prev_rows) + 1L]] <- data.frame(
    sex = sex_value, height_quintile = quintile, definition = outcome_labels[[outcome]],
    n = sum(dat$sex_label == sex_value & dat$height_quintile == quintile),
    prevalence_percent = 100 * p[1], lcl = 100 * p[2], ucl = 100 * p[3]
  )
}
height_prev <- do.call(rbind, prev_rows)
write.csv(height_prev, file.path(out, "05_height_quintile_prevalence.csv"), row.names = FALSE)

discordance_rows <- list()
for (sex_value in c("Female", "Male")) for (comparator in c("low_fnih", "low_ewgsop2")) {
  comp_label <- outcome_labels[[comparator]]
  x <- dat[dat$sex_label == sex_value, ]
  group <- ifelse(!x$low_p5 & !x[[comparator]], "Neither",
                  ifelse(!x$low_p5 & x[[comparator]], paste0(comp_label, " only"),
                         ifelse(x$low_p5 & !x[[comparator]], "Conditional only", "Both")))
  x$discordance_group <- group
  d <- make_design(x)
  for (g in c("Neither", paste0(comp_label, " only"), "Conditional only", "Both")) {
    d$variables$discordance_indicator <- d$variables$discordance_group == g
    p <- safe_prop(d, "discordance_indicator")
    n_cell <- sum(group == g)
    means <- if (n_cell) coef(svymean(~age + height_cm + bmi + alm_kg + almi,
                                      subset(d, discordance_group == g), na.rm = TRUE)) else rep(NA_real_, 5)
    discordance_rows[[length(discordance_rows) + 1L]] <- data.frame(
      sex = sex_value, comparison = paste("Conditional P5 vs", comp_label), group = g,
      unweighted_n = n_cell, weighted_percent = 100 * p[1], percent_lcl = 100 * p[2], percent_ucl = 100 * p[3],
      mean_age = means[1], mean_height_cm = means[2], mean_bmi = means[3], mean_alm_kg = means[4], mean_almi = means[5]
    )
  }
}
discordance <- do.call(rbind, discordance_rows)
write.csv(discordance, file.path(out, "06_classification_discordance.csv"), row.names = FALSE)

# Direct comparison of equal-year pooling against the official survey-area coefficients.
dat_equal <- dat
dat_equal$wt_pool <- dat_equal$wt_pool_equal_year
new_cal <- subset(calibration, scope == "Primary 19-69 pooled")
old_cal <- do.call(rbind, lapply(c("Female", "Male"), function(s) {
  calibration_one(dat_equal, "Primary 19-69 pooled", s)
}))
new_slope <- subset(height_slopes, scope == "Primary 19-69 pooled" & score == "z_conditional")
old_slope <- subset(run_slopes(dat_equal, "Primary 19-69 pooled", TRUE), score == "z_conditional")
new_or <- subset(height_or, scope == "Primary 19-69 pooled" & definition %in% c("Conditional P5", "FNIH"))
old_or <- subset(run_or(dat_equal, "Primary 19-69 pooled", TRUE), definition %in% c("Conditional P5", "FNIH"))
comparison <- rbind(
  data.frame(sex = new_cal$sex, metric = "mean_z", equal_year = old_cal$weighted_mean_z, official = new_cal$weighted_mean_z),
  data.frame(sex = new_cal$sex, metric = "p5_percent", equal_year = old_cal$p5_percent, official = new_cal$p5_percent),
  data.frame(sex = new_slope$sex, metric = "conditional_z_slope", equal_year = old_slope$estimate, official = new_slope$estimate),
  data.frame(sex = new_or$sex, metric = paste0(new_or$definition, "_OR"), equal_year = old_or$or_per_10cm, official = new_or$or_per_10cm)
)
comparison$absolute_change <- comparison$official - comparison$equal_year
write.csv(comparison, file.path(out, "07_equal_year_vs_official_comparison.csv"), row.names = FALSE)

cat("Official pooled-weight survey-layer rerun complete.\n")
print(new_cal, row.names = FALSE)
print(subset(height_slopes, scope %in% c("Primary 19-69 pooled", "Recommended height 19-69") & score == "z_conditional"), row.names = FALSE)
print(subset(height_or, scope %in% c("Primary 19-69 pooled", "Recommended height 19-69") & definition %in% c("Conditional P5", "FNIH")), row.names = FALSE)
print(comparison, row.names = FALSE)
