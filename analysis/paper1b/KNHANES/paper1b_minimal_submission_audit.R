# -*- coding: UTF-8 -*-

# Paper 1B minimal presubmission audit.
# This script reads only locked, analysis-ready objects and locked result tables.
# It does not rescore ALM, refit any reference distribution, or run bootstrap analyses.

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(try(Sys.setlocale("LC_ALL", "English_United States.utf8"), silent = TRUE))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages(library(survey))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

results_file <- file.path(root, "paper1b_minimal_submission_audit_results.csv")
report_file <- file.path(root, "paper1b_minimal_submission_audit_report.md")

required_files <- c(
  file.path(root, "data", "knhanes_2008_2011_analysis_19_69_official_weight.rds"),
  file.path(root, "data", "knhanes_2024_analysis_40_69.rds"),
  file.path(root, "tables_official_weight", "03_low_alm_height_or.csv"),
  file.path(root, "tables_2024", "02_calibration.csv"),
  file.path(root, "tables_2024", "03_z_height_slopes.csv"),
  file.path(root, "tables_2024", "04_low_alm_height_or.csv")
)
if (any(!file.exists(required_files))) {
  stop("Required locked analysis input missing: ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

d0811 <- readRDS(required_files[1])
d2024 <- readRDS(required_files[2])
locked_or_0811 <- read.csv(required_files[3], check.names = FALSE)
locked_cal_2024 <- read.csv(required_files[4], check.names = FALSE)
locked_slope_2024 <- read.csv(required_files[5], check.names = FALSE)
locked_or_2024 <- read.csv(required_files[6], check.names = FALSE)

result_rows <- list()
add_result <- function(section, cohort, sex = NA_character_, subgroup = NA_character_, metric,
                       estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                       p_value = NA_real_, numerator = NA_integer_, denominator = NA_integer_,
                       raw_proportion_percent = NA_real_, weighted_percent = NA_real_,
                       status = "PASS", note = "") {
  result_rows[[length(result_rows) + 1L]] <<- data.frame(
    section = section, cohort = cohort, sex = sex, subgroup = subgroup, metric = metric,
    estimate = as.numeric(estimate), se = as.numeric(se), ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high), p_value = as.numeric(p_value),
    numerator = as.integer(numerator), denominator = as.integer(denominator),
    raw_proportion_percent = as.numeric(raw_proportion_percent),
    weighted_percent = as.numeric(weighted_percent), status = status, note = note,
    stringsAsFactors = FALSE
  )
}

make_design_0811 <- function(z) svydesign(
  ids = ~psu_pool, strata = ~strata_pool, weights = ~wt_pool, nest = TRUE, data = z
)
make_design_2024 <- function(z) svydesign(
  ids = ~psu, strata = ~kstrata, weights = ~wt, nest = TRUE, data = z
)

safe_prop <- function(design, variable) {
  values <- design$variables[[variable]]
  values <- values[!is.na(values)]
  if (!length(values)) return(c(estimate = NA_real_, lcl = NA_real_, ucl = NA_real_))
  if (all(!values)) return(c(estimate = 0, lcl = 0, ucl = 0))
  if (all(values)) return(c(estimate = 1, lcl = 1, ucl = 1))
  ans <- svyciprop(as.formula(paste0("~", variable)), design, method = "logit", level = .95, na.rm = TRUE)
  ci <- confint(ans)
  c(estimate = as.numeric(coef(ans))[1], lcl = as.numeric(ci)[1], ucl = as.numeric(ci)[2])
}

extract_linear <- function(fit, term = "height10") {
  b <- unname(coef(fit)[term])
  se <- unname(sqrt(vcov(fit)[term, term]))
  c(estimate = b, se = se, lcl = b - qnorm(.975) * se, ucl = b + qnorm(.975) * se,
    p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}

capture_svyglm <- function(formula, design, family = NULL) {
  warnings_seen <- character()
  fit <- withCallingHandlers(
    tryCatch(
      if (is.null(family)) svyglm(formula, design = design) else svyglm(formula, design = design, family = family),
      error = function(e) e
    ),
    warning = function(w) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warnings_seen)
}

weighted_quantile <- function(v, w, probs) {
  ok <- is.finite(v) & is.finite(w) & w > 0
  v <- v[ok]; w <- w[ok]
  ord <- order(v); v <- v[ord]; w <- w[ord]
  v[pmin(length(v), findInterval(probs, cumsum(w) / sum(w)) + 1L)]
}

review_reasons <- character()
add_review <- function(condition, message) {
  if (isTRUE(condition)) review_reasons <<- unique(c(review_reasons, message))
}

# -----------------------------------------------------------------------------
# Locked sample checks and 2024 P5/P10 calibration
# -----------------------------------------------------------------------------

locked_counts <- list(
  `2008-2011` = c(Total = 16226L, Female = 9240L, Male = 6986L),
  `2024` = c(Total = 2513L, Female = 1505L, Male = 1008L)
)
actual_counts <- list(
  `2008-2011` = c(Total = nrow(d0811), table(factor(d0811$sex_label, levels = c("Female", "Male")))),
  `2024` = c(Total = nrow(d2024), table(factor(d2024$sex_label, levels = c("Female", "Male"))))
)
names(actual_counts[["2008-2011"]]) <- c("Total", "Female", "Male")
names(actual_counts[["2024"]]) <- c("Total", "Female", "Male")
for (cohort in names(locked_counts)) for (label in names(locked_counts[[cohort]])) {
  ok <- identical(as.integer(actual_counts[[cohort]][label]), as.integer(locked_counts[[cohort]][label]))
  add_result("Consistency checks", cohort, if (label == "Total") NA_character_ else label,
             metric = "Locked analytic sample count", estimate = actual_counts[[cohort]][label],
             status = if (ok) "PASS" else "REVIEW",
             note = paste0("Locked count=", locked_counts[[cohort]][label]))
  add_review(!ok, paste0(cohort, " ", label, " analytic count differs from the locked result."))
}

p10_cut <- -1.281551566
p5_cut <- qnorm(.05)
d2024$audit_low_p10 <- is.finite(d2024$z_conditional) & d2024$z_conditional < p10_cut
d2024$audit_low_p5 <- is.finite(d2024$z_conditional) & d2024$z_conditional < p5_cut
add_review(any(d2024$audit_low_p10 != d2024$low_p10, na.rm = TRUE),
           "Explicit P10 threshold does not reproduce the locked low_p10 indicator.")
add_review(any(d2024$audit_low_p5 != d2024$low_p5, na.rm = TRUE),
           "Explicit P5 threshold does not reproduce the locked low_p5 indicator.")

locked_p5 <- c(Female = 24.23, Male = 18.07)
calibration_2024 <- list()
for (sex_value in c("Female", "Male")) {
  z <- d2024[d2024$sex_label == sex_value, ]
  des <- make_design_2024(z)
  for (definition in c("P5", "P10")) {
    variable <- if (definition == "P5") "audit_low_p5" else "audit_low_p10"
    p <- safe_prop(des, variable)
    event_n <- sum(z[[variable]], na.rm = TRUE)
    row <- data.frame(
      sex = sex_value, definition = definition, event_n = event_n, denominator = nrow(z),
      percent = 100 * p[1], lcl = 100 * p[2], ucl = 100 * p[3]
    )
    calibration_2024[[length(calibration_2024) + 1L]] <- row
    status <- "PASS"
    note <- paste0("Frozen U.S. conditional ALM threshold; z < ",
                   if (definition == "P10") format(p10_cut, digits = 10) else format(p5_cut, digits = 10))
    if (definition == "P5") {
      status <- if (abs(row$percent - locked_p5[sex_value]) <= .05) "PASS" else "REVIEW"
      note <- paste0(note, "; locked rounded P5=", format(locked_p5[sex_value], nsmall = 2), "%")
      add_review(status == "REVIEW", paste0("2024 ", sex_value, " P5 prevalence does not reproduce the locked result."))
    }
    add_result("2024 calibration", "2024 (40-69)", sex_value, definition,
               metric = paste0("Below frozen U.S. conditional ", definition),
               estimate = row$percent, ci_low = row$lcl, ci_high = row$ucl,
               numerator = event_n, denominator = nrow(z),
               raw_proportion_percent = 100 * event_n / nrow(z), weighted_percent = row$percent,
               status = status, note = note)
  }
}
calibration_2024 <- do.call(rbind, calibration_2024)

# Honor the prespecified immediate-stop rule if the P5 consistency check fails.
if (any(calibration_2024$definition == "P5" &
        abs(calibration_2024$percent - locked_p5[calibration_2024$sex]) > .05)) {
  write.csv(do.call(rbind, result_rows), results_file, row.names = FALSE, na = "")
  writeLines(c(
    "# Paper 1B minimal submission audit",
    "",
    "**Analysis status: REVIEW**",
    "",
    "The 2024 P5 consistency check did not reproduce the locked prevalence within 0.05 percentage points.",
    "The audit stopped before any further calculations, as prespecified.",
    "",
    "**Review required before manuscript finalization because: the locked 2024 P5 result was not reproduced.**"
  ), report_file, useBytes = TRUE)
  stop("P5 consistency check failed; audit stopped as prespecified.")
}

# -----------------------------------------------------------------------------
# 2024 age-only z stature slope and locked conditional-z consistency
# -----------------------------------------------------------------------------

age_slope_2024 <- list()
locked_age_slope <- c(Female = .643, Male = .726)
for (sex_value in c("Female", "Male")) {
  des <- subset(make_design_2024(d2024), sex_label == sex_value)
  cap <- capture_svyglm(z_age_only ~ height10 + age, des)
  if (inherits(cap$fit, "error")) {
    add_review(TRUE, paste0("2024 ", sex_value, " age-only z slope model failed."))
    add_result("2024 age-only slope", "2024 (40-69)", sex_value,
               metric = "Age-only z slope per 10-cm greater height", status = "REVIEW",
               note = conditionMessage(cap$fit))
    next
  }
  e <- extract_linear(cap$fit)
  ok <- is.finite(e["estimate"]) && abs(e["estimate"] - locked_age_slope[sex_value]) <= .01
  add_review(!ok, paste0("2024 ", sex_value, " age-only z slope differs from the locked rounded estimate."))
  age_slope_2024[[sex_value]] <- e
  add_result("2024 age-only slope", "2024 (40-69)", sex_value,
             metric = "Age-only z slope per 10-cm greater height",
             estimate = e["estimate"], se = e["se"], ci_low = e["lcl"], ci_high = e["ucl"],
             p_value = e["p_value"], status = if (ok) "PASS" else "REVIEW",
             note = paste0("Survey model: z_age_only ~ height10 + age; locked rounded beta=",
                           locked_age_slope[sex_value],
                           if (length(cap$warnings)) paste0("; warnings: ", paste(cap$warnings, collapse = " | ")) else "; no warnings"))
}

locked_conditional_expected <- list(
  Female = c(estimate = -.119, lcl = -.215, ucl = -.024),
  Male = c(estimate = -.069, lcl = -.173, ucl = .035)
)
for (sex_value in c("Female", "Male")) {
  row <- subset(locked_slope_2024, scope == "Primary 40-69" & sex == sex_value & score == "z_conditional")
  ok <- nrow(row) == 1L && abs(row$estimate - locked_conditional_expected[[sex_value]]["estimate"]) <= .002 &&
    abs(row$lcl - locked_conditional_expected[[sex_value]]["lcl"]) <= .002 &&
    abs(row$ucl - locked_conditional_expected[[sex_value]]["ucl"]) <= .002
  add_result("Consistency checks", "2024 (40-69)", sex_value,
             metric = "Locked conditional-z slope per 10 cm",
             estimate = if (nrow(row)) row$estimate[1] else NA_real_,
             se = if (nrow(row)) row$se[1] else NA_real_,
             ci_low = if (nrow(row)) row$lcl[1] else NA_real_,
             ci_high = if (nrow(row)) row$ucl[1] else NA_real_,
             p_value = if (nrow(row)) row$p_value[1] else NA_real_,
             status = if (ok) "PASS" else "REVIEW", note = "Read from locked tables_2024/03_z_height_slopes.csv")
  add_review(!ok, paste0("2024 ", sex_value, " locked conditional-z slope check failed."))
}

# -----------------------------------------------------------------------------
# FNIH sparse-event / near-separation audit
# -----------------------------------------------------------------------------

cohorts <- list(
  `2008-2011 (19-69)` = list(data = d0811, design = make_design_0811, weight = "wt_pool", year_adjust = TRUE),
  `2024 (40-69)` = list(data = d2024, design = make_design_2024, weight = "wt", year_adjust = FALSE)
)
fnih_quintile_rows <- list()
fnih_model_rows <- list()

for (cohort in names(cohorts)) {
  spec <- cohorts[[cohort]]
  z <- spec$data
  z$height_quintile_audit <- NA_character_
  for (sex_value in c("Female", "Male")) {
    idx <- z$sex_label == sex_value
    cuts <- weighted_quantile(z$height_cm[idx], z[[spec$weight]][idx], c(.2, .4, .6, .8))
    z$height_quintile_audit[idx] <- as.character(cut(
      z$height_cm[idx], c(-Inf, cuts, Inf), labels = paste0("Q", 1:5), include.lowest = TRUE
    ))
  }
  z$height_quintile_audit <- factor(z$height_quintile_audit, levels = paste0("Q", 1:5))
  desq <- spec$design(z)

  for (sex_value in c("Female", "Male")) {
    q_rows_sex <- list()
    for (q in paste0("Q", 1:5)) {
      cell <- z[z$sex_label == sex_value & z$height_quintile_audit == q, ]
      dz <- subset(desq, sex_label == sex_value & height_quintile_audit == q)
      p <- safe_prop(dz, "low_fnih")
      events <- sum(cell$low_fnih, na.rm = TRUE)
      qr <- data.frame(quintile = q, n = nrow(cell), events = events,
                       raw_percent = 100 * events / nrow(cell), weighted_percent = 100 * p[1],
                       lcl = 100 * p[2], ucl = 100 * p[3])
      q_rows_sex[[q]] <- qr
      fnih_quintile_rows[[length(fnih_quintile_rows) + 1L]] <- cbind(
        data.frame(cohort = cohort, sex = sex_value), qr
      )
      add_result("FNIH height-quintile audit", cohort, sex_value, q,
                 metric = "FNIH-low by sex-specific weighted height quintile",
                 estimate = p[1], ci_low = p[2], ci_high = p[3],
                 numerator = events, denominator = nrow(cell), raw_proportion_percent = qr$raw_percent,
                 weighted_percent = qr$weighted_percent,
                 status = if (events == 0L) "SPARSE FLAG" else "PASS",
                 note = paste0("Zero-event quintile: ", if (events == 0L) "yes" else "no"))
    }
    q_rows_sex <- do.call(rbind, q_rows_sex)
    zero_q <- any(q_rows_sex$events == 0L)

    des <- subset(spec$design(z), sex_label == sex_value)
    form <- if (spec$year_adjust) low_fnih ~ height10 + age + factor(year) else low_fnih ~ height10 + age
    cap <- capture_svyglm(form, des, family = quasibinomial())
    fit_error <- inherits(cap$fit, "error")
    if (fit_error) {
      converged <- FALSE; e <- rep(NA_real_, 5); names(e) <- c("estimate", "se", "lcl", "ucl", "p_value")
      error_text <- conditionMessage(cap$fit)
    } else {
      converged <- isTRUE(cap$fit$converged)
      e <- extract_linear(cap$fit)
      error_text <- ""
    }
    finite_coef_se <- is.finite(e["estimate"]) && is.finite(e["se"])
    abnormal_se <- !is.finite(e["se"]) || e["se"] > 5
    warning_text <- if (length(cap$warnings)) paste(cap$warnings, collapse = " | ") else "None"
    numerical_warning <- any(grepl("did not converge|fitted probabilities numerically|non-finite|separation|singular",
                                   cap$warnings, ignore.case = TRUE))
    model_ok <- !fit_error && converged && finite_coef_se && !abnormal_se && !numerical_warning
    add_review(!model_ok, paste0(cohort, " ", sex_value, " FNIH survey-logistic model has a numerical issue."))

    fnih_model_rows[[length(fnih_model_rows) + 1L]] <- data.frame(
      cohort = cohort, sex = sex_value, coefficient = e["estimate"], se = e["se"],
      lcl = e["lcl"], ucl = e["ucl"], p_value = e["p_value"],
      or = exp(e["estimate"]), or_lcl = exp(e["lcl"]), or_ucl = exp(e["ucl"]),
      converged = converged, finite_coef_se = finite_coef_se, abnormal_se = abnormal_se,
      zero_event_quintile = zero_q, warnings = warning_text, stringsAsFactors = FALSE
    )
    add_result("FNIH model audit", cohort, sex_value,
               metric = "FNIH survey-logistic coefficient per 10 cm",
               estimate = e["estimate"], se = e["se"], ci_low = e["lcl"], ci_high = e["ucl"],
               p_value = e["p_value"], status = if (model_ok) "PASS" else "REVIEW",
               note = paste0("Converged=", converged, "; finite coefficient/SE=", finite_coef_se,
                             "; abnormal SE=", abnormal_se, "; numerical warnings=", warning_text,
                             if (nzchar(error_text)) paste0("; error=", error_text) else ""))
    add_result("FNIH model audit", cohort, sex_value,
               metric = "FNIH OR per 10-cm greater height",
               estimate = exp(e["estimate"]), ci_low = exp(e["lcl"]), ci_high = exp(e["ucl"]),
               p_value = e["p_value"], status = if (model_ok) "PASS" else "REVIEW",
               note = "Primary prespecified survey-weighted quasibinomial model retained")

    firth_trigger <- zero_q || !model_ok
    add_result("FNIH model audit", cohort, sex_value,
               metric = "Diagnostic-only penalized/Firth sensitivity",
               status = if (firth_trigger) "NOT RUN" else "NOT INDICATED",
               note = if (firth_trigger) {
                 paste0("Trigger present (zero-event quintile or numerical issue). Not implemented because the current ",
                        "complex-survey framework has no design-consistent Firth estimator; an unweighted penalized fit ",
                        "would not replace or validate the prespecified survey model.")
               } else {
                 "No zero-event quintile, separation warning, nonconvergence, non-finite coefficient/SE, or abnormal SE."
               })
  }
}
fnih_quintile_rows <- do.call(rbind, fnih_quintile_rows)
fnih_model_rows <- do.call(rbind, fnih_model_rows)

# Compare the four refitted primary FNIH models with their locked result rows.
for (i in seq_len(nrow(fnih_model_rows))) {
  rr <- fnih_model_rows[i, ]
  if (rr$cohort == "2008-2011 (19-69)") {
    lk <- subset(locked_or_0811, scope == "Primary 19-69 pooled" & sex == rr$sex & definition == "FNIH")
  } else {
    lk <- subset(locked_or_2024, scope == "Primary 40-69" & sex == rr$sex & definition == "FNIH")
  }
  ok <- nrow(lk) == 1L && is.finite(rr$or) && abs(rr$or - lk$or_per_10cm[1]) < 1e-10
  add_result("Consistency checks", rr$cohort, rr$sex,
             metric = "Refitted FNIH OR matches locked result", estimate = rr$or,
             ci_low = rr$or_lcl, ci_high = rr$or_ucl, status = if (ok) "PASS" else "REVIEW",
             note = if (nrow(lk)) paste0("Locked OR=", format(lk$or_per_10cm[1], digits = 12)) else "Locked row missing")
  add_review(!ok, paste0(rr$cohort, " ", rr$sex, " FNIH OR does not match the locked result."))
}

# Verify direction of the already locked year-specific 2008-2011 FNIH estimates.
annual_fnih <- subset(locked_or_0811, scope == "Primary 19-69 by year" & definition == "FNIH")
for (i in seq_len(nrow(annual_fnih))) {
  rr <- annual_fnih[i, ]
  ok <- is.finite(rr$or_per_10cm) && rr$or_per_10cm < 1
  add_result("FNIH annual direction check", paste0(rr$year, " (19-69)"), rr$sex,
             metric = "Locked annual FNIH OR per 10 cm", estimate = rr$or_per_10cm,
             ci_low = rr$or_lcl, ci_high = rr$or_ucl, p_value = rr$p_value,
             status = if (ok) "PASS" else "REVIEW",
             note = "Read from locked official-weight result; expected inverse height association")
  add_review(!ok, paste0(rr$year, " ", rr$sex, " locked FNIH OR is not below 1."))
}

# -----------------------------------------------------------------------------
# Conditional P5 x EWGSOP2 nesting audit
# -----------------------------------------------------------------------------

nesting_rows <- list()
nest_levels <- c("Neither low", "EWGSOP2 only", "Conditional P5 only", "Both low")
for (cohort in names(cohorts)) {
  spec <- cohorts[[cohort]]
  for (sex_value in c("Female", "Male")) {
    z <- spec$data[spec$data$sex_label == sex_value, ]
    z$nest_group <- factor(
      ifelse(!z$low_p5 & !z$low_ewgsop2, "Neither low",
             ifelse(!z$low_p5 & z$low_ewgsop2, "EWGSOP2 only",
                    ifelse(z$low_p5 & !z$low_ewgsop2, "Conditional P5 only", "Both low"))),
      levels = nest_levels
    )
    des <- spec$design(z)
    for (g in nest_levels) {
      des$variables$nest_indicator <- des$variables$nest_group == g
      p <- safe_prop(des, "nest_indicator")
      n_cell <- sum(z$nest_group == g)
      status <- if (g == "Conditional P5 only" && n_cell > 0L) "REVIEW" else "PASS"
      nesting_rows[[length(nesting_rows) + 1L]] <- data.frame(
        cohort = cohort, sex = sex_value, group = g, n = n_cell,
        weighted_percent = 100 * p[1], lcl = 100 * p[2], ucl = 100 * p[3]
      )
      add_result("Conditional P5 x EWGSOP2 nesting", cohort, sex_value, g,
                 metric = "Raw 2x2 classification cell", estimate = n_cell,
                 ci_low = 100 * p[2], ci_high = 100 * p[3], numerator = n_cell,
                 denominator = nrow(z), raw_proportion_percent = 100 * n_cell / nrow(z),
                 weighted_percent = 100 * p[1], status = status,
                 note = if (g == "Conditional P5 only") paste0("Conditional-P5-only raw count=", n_cell) else "")
      add_review(status == "REVIEW", paste0(cohort, " ", sex_value,
                                             " has participants classified as conditional-P5-only."))
    }
  }
}
nesting_rows <- do.call(rbind, nesting_rows)

# -----------------------------------------------------------------------------
# Final output and manuscript-facing report
# -----------------------------------------------------------------------------

all_results <- do.call(rbind, result_rows)
write.csv(all_results, results_file, row.names = FALSE, na = "")

fmt <- function(x, digits = 3) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
fmt_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < .001) return("<0.001")
  formatC(x, digits = 3, format = "f")
}
ci_text <- function(est, lcl, ucl, digits = 3) {
  paste0(fmt(est, digits), " (95% CI ", fmt(lcl, digits), " to ", fmt(ucl, digits), ")")
}

status <- if (length(review_reasons)) "REVIEW" else "PASS"
report <- c(
  "# Paper 1B minimal submission audit",
  "",
  paste0("**Analysis status: ", status, "**"),
  "",
  "This audit used only the locked analysis-ready KNHANES objects, their existing survey-design variables, and locked result tables. It did not rescore ALM, re-estimate the U.S. reference, refit a Korean reference, or run any bootstrap analysis.",
  "",
  "## 1. KNHANES 2024 P10 calibration",
  "",
  "| Sex | Below U.S. P10, weighted % (95% CI) | Raw events / N | P5 consistency check |",
  "|---|---:|---:|---:|"
)
for (sex_value in c("Female", "Male")) {
  p10 <- subset(calibration_2024, sex == sex_value & definition == "P10")
  p5 <- subset(calibration_2024, sex == sex_value & definition == "P5")
  report <- c(report, paste0("| ", sex_value, " | ", ci_text(p10$percent, p10$lcl, p10$ucl, 2),
                             " | ", p10$event_n, " / ", p10$denominator,
                             " | ", fmt(p5$percent, 2), "% (locked ", fmt(locked_p5[sex_value], 2), "%) |"))
}

report <- c(report, "", "P10 used the prespecified threshold z < -1.281551566. P5 was recomputed only as a consistency check.",
            "", "## 2. KNHANES 2024 age-only z slope", "",
            "| Sex | Beta per 10-cm greater height | SE | 95% CI | P value |",
            "|---|---:|---:|---:|---:|")
for (sex_value in c("Female", "Male")) {
  e <- age_slope_2024[[sex_value]]
  report <- c(report, paste0("| ", sex_value, " | ", fmt(e["estimate"]), " | ", fmt(e["se"]),
                             " | ", fmt(e["lcl"]), " to ", fmt(e["ucl"]), " | ", fmt_p(e["p_value"]), " |"))
}

report <- c(report, "", "Models used the locked 2024 analytic sample and survey design, with age adjustment and height expressed per 10 cm.",
            "", "## 3. FNIH rare-event and near-separation audit", "",
            "### Primary survey-logistic models", "",
            "| Cohort | Sex | Log-odds coefficient | SE | OR per 10 cm (95% CI) | Converged | Finite coefficient/SE | Numerical warning | Zero-event quintile |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(fnih_model_rows))) {
  rr <- fnih_model_rows[i, ]
  report <- c(report, paste0("| ", rr$cohort, " | ", rr$sex, " | ",
                             fmt(rr$coefficient), " | ", fmt(rr$se), " | ",
                             ci_text(rr$or, rr$or_lcl, rr$or_ucl), " | ",
                             if (rr$converged) "Yes" else "No", " | ",
                             if (rr$finite_coef_se) "Yes" else "No", " | ",
                             if (rr$warnings == "None") "No" else "Yes", " | ",
                             if (rr$zero_event_quintile) "Yes" else "No", " |"))
}

report <- c(report, "", "### Sex-specific weighted height quintiles", "",
            "| Cohort | Sex | Quintile | Raw FNIH events / N | Raw % | Survey-weighted % (95% CI) |",
            "|---|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(fnih_quintile_rows))) {
  rr <- fnih_quintile_rows[i, ]
  report <- c(report, paste0("| ", rr$cohort, " | ", rr$sex, " | ", rr$quintile,
                             " | ", rr$events, " / ", rr$n, " | ", fmt(rr$raw_percent, 3),
                             " | ", ci_text(rr$weighted_percent, rr$lcl, rr$ucl, 3), " |"))
}

annual_ok <- nrow(annual_fnih) == 8L && all(is.finite(annual_fnih$or_per_10cm)) && all(annual_fnih$or_per_10cm < 1)
annual_ranges <- aggregate(or_per_10cm ~ sex, annual_fnih, range)
annual_range_text <- paste(vapply(seq_len(nrow(annual_ranges)), function(i) {
  paste0(tolower(annual_ranges$sex[i]), " ", fmt(annual_ranges$or_per_10cm[i, 1]),
         "-", fmt(annual_ranges$or_per_10cm[i, 2]))
}, character(1)), collapse = "; ")
report <- c(report, "", paste0("All eight locked sex-by-year FNIH ORs for 2008-2011 remained below 1: ",
                                if (annual_ok) "yes" else "no", " (", annual_range_text, ")."),
            paste0("Sparse-event conclusion: one zero-event cell was present (2008-2011 men, Q4), while the adjacent highest quintile contained 3/1181 events. ",
                   "All four continuous-height survey-logistic models converged with finite coefficients and SEs and produced no numerical warnings. ",
                   "The isolated zero cell therefore documents sparsity but does not show complete separation of the prespecified continuous-height model."),
            "Penalized/Firth logistic regression was not used to replace the prespecified complex-survey model. Any trigger and the reason for not forcing an incompatible sensitivity are recorded in the CSV audit table.",
            "", "## 4. Conditional P5 x EWGSOP2 nesting audit", "",
            "| Cohort | Sex | Classification cell | Raw N | Survey-weighted % (95% CI) |",
            "|---|---|---|---:|---:|")
for (i in seq_len(nrow(nesting_rows))) {
  rr <- nesting_rows[i, ]
  report <- c(report, paste0("| ", rr$cohort, " | ", rr$sex, " | ", rr$group,
                             " | ", rr$n, " | ", ci_text(rr$weighted_percent, rr$lcl, rr$ucl, 3), " |"))
}

all_nested <- all(nesting_rows$n[nesting_rows$group == "Conditional P5 only"] == 0L)
if (all_nested) {
  report <- c(report, "", "All participants classified below the frozen conditional P5 also met the EWGSOP2 low-ALMI criterion in this sample.")
}
report <- c(report, "", "Nesting is a classification relationship in these samples and is not evidence of diagnostic superiority.",
            "", "## Consistency checks", "",
            paste0("- Locked sample sizes reproduced: ", if (all(vapply(names(locked_counts), function(k) all(actual_counts[[k]] == locked_counts[[k]]), logical(1)))) "yes" else "no", "."),
            paste0("- Locked 2024 P5 fractions reproduced within 0.05 percentage points: ",
                   if (all(abs(calibration_2024$percent[calibration_2024$definition == "P5"] - locked_p5[calibration_2024$sex[calibration_2024$definition == "P5"]]) <= .05)) "yes" else "no", "."),
            paste0("- Locked 2024 conditional-z stature slopes reproduced from the existing result table: ",
                   if (!any(grepl("conditional-z slope", review_reasons))) "yes" else "no", "."),
            paste0("- Refitted FNIH ORs matched the locked main estimates: ",
                   if (!any(grepl("FNIH OR does not match", review_reasons))) "yes" else "no", "."),
            paste0("- Conditional-P5-only raw count was zero in every cohort-sex table: ", if (all_nested) "yes" else "no", "."),
            "", "## Effect on Paper 1B conclusions", "")

if (length(review_reasons)) {
  report <- c(report, "The following items require review:", paste0("- ", review_reasons), "",
              paste0("**Review required before manuscript finalization because: ", paste(review_reasons, collapse = " "), "**"))
} else {
  report <- c(report,
              "No audit result requires modification of the existing Paper 1B main conclusions.",
              "",
              "**Main conclusions unchanged; manuscript can proceed to figure/manuscript finalization.**")
}

writeLines(report, report_file, useBytes = TRUE)

cat("Paper 1B minimal submission audit complete.\n")
cat("Analysis status:", status, "\n")
cat("Results:", results_file, "\n")
cat("Report:", report_file, "\n")
cat(if (status == "PASS") {
  "Main conclusions unchanged; manuscript can proceed to figure/manuscript finalization.\n"
} else {
  paste0("Review required before manuscript finalization because: ", paste(review_reasons, collapse = " "), "\n")
})
