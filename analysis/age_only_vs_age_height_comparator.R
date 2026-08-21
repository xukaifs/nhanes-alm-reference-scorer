# Age-only versus age+height ALM reference sensitivity analysis
# Repository release: v2.1.0
#
# Purpose
# -------
# Reproduce the controlled sensitivity analysis used in the manuscript to isolate
# the incremental contribution of explicit stature conditioning. The age-only
# comparator keeps the sex-specific distribution family, age smooth, scale and
# shape structures, corrected development weights, completed DXA datasets, and
# deterministic PSU folds identical to the frozen age+height model. The only
# model change is removal of height from the GAMLSS location predictor.
#
# This script does NOT alter the public primary scorer or frozen age+height model
# objects. It is a manuscript-analysis script.
#
# Inputs
# ------
# 1) The v2.0.0/v2.1.0 frozen model bundles in ../models/. Their embedded
#    train_data frames contain the final corrected-weight 1999-2006 development
#    analysis data used for each completed DXA dataset.
# 2) A user-supplied analysis-ready NHANES 2011-2018 validation file, either RDS
#    or CSV/CSV.GZ. Set environment variable ALM_VALIDATION_FILE or pass the file
#    as the first command-line argument. Required variables are documented in
#    analysis/README.md. Participant-level validation data are not redistributed
#    in this repository.
#
# Outputs
# -------
# analysis/output/age_only_full_fit_audit.csv
# analysis/output/age_only_oof_calibration.csv
# analysis/output/age_only_vs_age_height_calibration.csv
# analysis/output/age_only_vs_age_height_temporal.csv
# analysis/output/age_only_vs_age_height_spline_audit.csv
# analysis/output/age_only_vs_age_height_reclassification.csv
# analysis/output/age_only_vs_age_height_summary.csv
#
# The expected manuscript-level results are archived separately in
# analysis/expected/age_only_vs_age_height_expected.csv.

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(gamlss)
  library(gamlss.dist)
  library(survey)
  library(splines)
})

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "R", "score_conditional_alm.R"))) {
  stop("Run this script from the repository root.")
}

out_dir <- file.path(repo_root, "analysis", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
validation_path <- if (length(args) >= 1L && nzchar(args[[1]])) {
  args[[1]]
} else {
  Sys.getenv("ALM_VALIDATION_FILE", unset = "")
}
if (!nzchar(validation_path) || !file.exists(validation_path)) {
  stop(
    "Validation file not found. Supply an analysis-ready NHANES 2011-2018 file ",
    "as the first argument or set ALM_VALIDATION_FILE. See analysis/README.md."
  )
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
read_validation <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    read_csv(path, show_col_types = FALSE, progress = FALSE)
  }
}

rubin_pool <- function(q, se, null = 0) {
  ok <- is.finite(q) & is.finite(se) & se >= 0
  q <- q[ok]; se <- se[ok]; m <- length(q)
  if (!m) {
    return(tibble(
      m = 0, estimate = NA_real_, standard_error = NA_real_, df = NA_real_,
      conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_
    ))
  }
  qbar <- mean(q)
  u <- mean(se^2)
  b <- if (m > 1) var(q) else 0
  tvar <- u + (1 + 1 / m) * b
  s <- sqrt(tvar)
  df <- if (m > 1 && b > 0) (m - 1) * (1 + u / ((1 + 1 / m) * b))^2 else Inf
  crit <- qt(.975, df)
  tibble(
    m = m, estimate = qbar, standard_error = s, df = df,
    conf_low = qbar - crit * s, conf_high = qbar + crit * s,
    p_value = 2 * pt(abs((qbar - null) / s), df = df, lower.tail = FALSE)
  )
}

pool_table <- function(dat, groups) {
  dat %>%
    group_by(across(all_of(groups))) %>%
    group_modify(~ rubin_pool(.x$estimate, .x$standard_error)) %>%
    ungroup()
}

make_design <- function(dat) {
  svydesign(
    ids = ~psu_pool, strata = ~strata_pool, weights = ~wtmec_pool,
    nest = TRUE, data = dat
  )
}

safe_sd <- function(design, var) {
  ans <- try(svyvar(as.formula(paste0("~", var)), design, na.rm = TRUE), silent = TRUE)
  if (inherits(ans, "try-error")) return(NA_real_)
  sqrt(pmax(as.numeric(coef(ans)[1]), 0))
}

predict_parameters <- function(model, train_data, new_data) {
  nd <- as.data.frame(new_data[, intersect(names(new_data), names(train_data)), drop = FALSE])
  ans <- try(
    suppressWarnings(predictAll(model, newdata = nd, data = as.data.frame(train_data), type = "response")),
    silent = TRUE
  )
  if (inherits(ans, "try-error")) {
    ans <- suppressWarnings(predictAll(model, newdata = nd, type = "response"))
  }
  as.data.frame(ans)
}

score_fit <- function(model, family_name, train_data, new_data) {
  pars <- predict_parameters(model, train_data, new_data)
  pfun <- get(paste0("p", family_name), envir = asNamespace("gamlss.dist"))
  p_args <- list(q = new_data$alm_kg, mu = pars$mu, sigma = pars$sigma)
  if ("nu" %in% names(pars)) p_args$nu <- pars$nu
  if ("tau" %in% names(pars)) p_args$tau <- pars$tau
  p <- as.numeric(do.call(pfun, p_args))
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  tibble(percentile = p, z = qnorm(p))
}

fit_age_only <- function(train_data, sex_value, family_name,
                         mu_step = 1, algorithm = "RS", start_mode = "lm") {
  family_function <- get(family_name, envir = asNamespace("gamlss.dist"))
  formula_text <- if (sex_value == "Female") {
    "alm_kg ~ pb(age, df=3, inter=10)"
  } else {
    "alm_kg ~ pb(age, df=3)"
  }
  lm_start <- lm(alm_kg ~ age, data = train_data, weights = w_model)
  mu_start <- if (start_mode == "lm") {
    pmax(as.numeric(predict(lm_start, train_data)), .5)
  } else {
    NULL
  }
  method_call <- switch(
    algorithm,
    RS = quote(RS()),
    CG = quote(CG()),
    mixed = quote(mixed()),
    stop("Unknown algorithm: ", algorithm)
  )
  fit_call <- substitute(
    gamlss(
      as.formula(FORMULA_TEXT),
      sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
      family = family_function(), data = train_data, weights = w_model,
      mu.start = START, method = METHOD,
      control = gamlss.control(n.cyc = 500, trace = FALSE, mu.step = STEP)
    ),
    list(FORMULA_TEXT = formula_text, START = mu_start, METHOD = method_call, STEP = mu_step)
  )
  fit <- eval(fit_call)
  fit$call$family <- as.call(list(as.name(family_name)))
  fit$call$data <- quote(train_data)
  fit
}

fit_age_only_with_retries <- function(train_data, sex_value, family_name) {
  attempts <- tribble(
    ~mu_step, ~algorithm, ~start_mode,
    1.00, "RS", "lm",
    0.50, "RS", "lm",
    0.25, "RS", "lm",
    0.10, "RS", "lm",
    0.50, "mixed", "lm",
    0.25, "mixed", "lm",
    0.25, "CG", "lm",
    0.50, "RS", "default"
  )
  last <- NULL
  for (i in seq_len(nrow(attempts))) {
    fit <- try(
      fit_age_only(
        train_data, sex_value, family_name,
        attempts$mu_step[i], attempts$algorithm[i], attempts$start_mode[i]
      ),
      silent = TRUE
    )
    last <- fit
    if (!inherits(fit, "try-error") && isTRUE(fit$converged)) {
      return(list(fit = fit, attempt = i))
    }
  }
  list(fit = last, attempt = NA_integer_)
}

make_age_only_anchors <- function(dat, sex_value) {
  # Clone complete training rows so auxiliary columns remain non-missing in the
  # GAMLSS data frame; then overwrite only the fields needed for numerical
  # support. This mirrors the purpose of the primary H-model support anchors.
  lm_fit <- lm(alm_kg ~ age, data = dat, weights = wtmec_correct)
  ages <- c(18, 69, median(dat$age), median(dat$age))
  template <- dat[rep(1L, length(ages)), , drop = FALSE]
  template$SEQN <- -seq_along(ages)
  template$imputation <- unique(dat$imputation)[1]
  template$cycle <- "anchor"
  template$sex <- sex_value
  template$age <- ages
  template$height_m <- median(dat$height_m, na.rm = TRUE)
  template$height_cm <- 100 * template$height_m
  if ("bmi" %in% names(template)) template$bmi <- median(dat$bmi, na.rm = TRUE)
  template$alm_kg <- pmax(as.numeric(predict(lm_fit, newdata = data.frame(age = ages))), .5)
  template$wtmec_correct <- mean(dat$wtmec_correct) * 1e-4
  if ("wtmec_pool" %in% names(template)) template$wtmec_pool <- template$wtmec_correct
  if ("w_model" %in% names(template)) template$w_model <- NA_real_
  template$strata <- "anchor"
  template$psu <- paste0("anchor", seq_along(ages))
  template$strata_pool <- factor("anchor")
  template$psu_pool <- factor(paste0("anchor", seq_along(ages)))
  template
}

prepare_development_from_frozen_bundles <- function() {
  rows <- list(); k <- 1L
  for (sex_value in c("Female", "Male")) {
    bundle_path <- file.path(
      repo_root, "models", paste0("corrected_H_18_69_", sex_value, "_bundles.rds")
    )
    bundles <- readRDS(bundle_path)
    for (imp in 1:5) {
      key <- paste0("imp", imp, "_", sex_value)
      b <- bundles[[key]]
      dat <- as.data.frame(b$train_data)
      needed <- c("SEQN", "imputation", "cycle", "sex", "age", "height_m", "alm_kg",
                  "wtmec_correct", "strata", "psu", "strata_pool", "psu_pool")
      miss <- setdiff(needed, names(dat))
      if (length(miss)) stop("Frozen train_data missing: ", paste(miss, collapse = ", "))
      if (!"height_cm" %in% names(dat)) dat$height_cm <- 100 * dat$height_m
      if (!"bmi" %in% names(dat)) dat$bmi <- NA_real_
      dat$imputation <- imp
      dat$sex <- sex_value
      rows[[k]] <- dat
      k <- k + 1L
    }
  }
  bind_rows(rows) %>%
    mutate(
      sex = as.character(sex),
      strata_pool = factor(strata_pool),
      psu_pool = factor(psu_pool)
    )
}

prepare_validation <- function(dat) {
  dat <- as.data.frame(dat)
  if (!"height_m" %in% names(dat) && "height_cm" %in% names(dat)) {
    dat$height_m <- dat$height_cm / 100
  }
  if (!"height_cm" %in% names(dat) && "height_m" %in% names(dat)) {
    dat$height_cm <- dat$height_m * 100
  }
  if (!"imputation" %in% names(dat)) {
    dat <- tidyr::crossing(imputation = 1:5, dat)
  }
  required <- c(
    "imputation", "sex", "age", "height_m", "height_cm", "alm_kg",
    "wtmec_pool", "strata_pool", "psu_pool", "race", "cycle"
  )
  miss <- setdiff(required, names(dat))
  if (length(miss)) stop("Validation file missing required variables: ", paste(miss, collapse = ", "))
  dat %>%
    filter(age >= 18, age <= 59, sex %in% c("Female", "Male")) %>%
    mutate(
      imputation = as.integer(imputation),
      sex = factor(as.character(sex), levels = c("Female", "Male")),
      race = factor(as.character(race)), cycle = factor(as.character(cycle)),
      strata_pool = factor(strata_pool), psu_pool = factor(psu_pool)
    )
}

# -----------------------------------------------------------------------------
# Load final development data from the frozen H bundles and the user-supplied
# later-period validation data.
# -----------------------------------------------------------------------------
development <- prepare_development_from_frozen_bundles()
validation <- prepare_validation(read_validation(validation_path))

# Score the frozen H reference on the validation file if H scores are not already
# present. This makes the comparator independent of legacy cache column names.
if (!all(c("z_H", "p_H") %in% names(validation))) {
  scored_h <- list(); k <- 1L
  for (sex_value in c("Female", "Male")) {
    bundles <- readRDS(file.path(
      repo_root, "models", paste0("corrected_H_18_69_", sex_value, "_bundles.rds")
    ))
    for (imp in 1:5) {
      dat <- validation %>% filter(imputation == imp, as.character(sex) == sex_value)
      b <- bundles[[paste0("imp", imp, "_", sex_value)]]
      sc <- score_fit(b$model, b$family, b$train_data, dat)
      scored_h[[k]] <- bind_cols(dat, sc %>% rename(p_H = percentile, z_H = z))
      k <- k + 1L
    }
  }
  validation <- bind_rows(scored_h)
}
validation <- validation %>%
  mutate(low_p5_H = as.integer(p_H < .05), low_p10_H = as.integer(p_H < .10))

# -----------------------------------------------------------------------------
# Full-data age-only fits and BIC audit.
# -----------------------------------------------------------------------------
age_only_bundles <- list()
fit_audit <- list(); fa <- 1L
for (sex_value in c("Female", "Male")) {
  family_name <- if (sex_value == "Female") "BCCG" else "BCT"
  for (imp in 1:5) {
    dat <- development %>%
      filter(imputation == imp, sex == sex_value) %>%
      mutate(w_model = wtmec_correct / mean(wtmec_correct)) %>%
      droplevels()
    result <- fit_age_only_with_retries(dat, sex_value, family_name)
    fit <- result$fit
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      stop("Age-only full fit failed for ", sex_value, " imputation ", imp)
    }
    key <- paste0("imp", imp, "_", sex_value)
    age_only_bundles[[key]] <- list(
      imputation = imp, sex = sex_value, family = family_name,
      train_data = dat, model = fit
    )
    fit_audit[[fa]] <- tibble(
      sex = sex_value, imputation = imp, family = family_name, n = nrow(dat),
      converged = fit$converged, iterations = fit$iter, effective_df = fit$df.fit,
      AIC = AIC(fit), BIC = GAIC(fit, k = log(nrow(dat))), attempt = result$attempt
    )
    fa <- fa + 1L
  }
}
fit_audit <- bind_rows(fit_audit)
write_csv(fit_audit, file.path(out_dir, "age_only_full_fit_audit.csv"))

# -----------------------------------------------------------------------------
# PSU-grouped 5-fold OOF calibration. Fold allocation exactly mirrors the final
# H-model pipeline: set.seed(20260812 + 1000*imp + 17 for women, +0 for men).
# -----------------------------------------------------------------------------
oof_rows <- list(); oo <- 1L
for (sex_value in c("Female", "Male")) {
  family_name <- if (sex_value == "Female") "BCCG" else "BCT"
  for (imp in 1:5) {
    dat <- development %>% filter(imputation == imp, sex == sex_value) %>% droplevels()
    units <- unique(as.character(dat$psu_pool))
    set.seed(20260812 + 1000 * imp + ifelse(sex_value == "Female", 17, 0))
    fold_map <- tibble(unit = units, fold = sample(rep(1:5, length.out = length(units))))
    fold_assignment <- fold_map$fold[match(as.character(dat$psu_pool), fold_map$unit)]
    anchors <- make_age_only_anchors(dat, sex_value)
    for (fold in 1:5) {
      raw_train <- dat[fold_assignment != fold, ]
      test <- dat[fold_assignment == fold, ]
      train <- bind_rows(raw_train, anchors) %>%
        mutate(w_model = wtmec_correct / mean(wtmec_correct)) %>%
        droplevels()
      result <- fit_age_only_with_retries(train, sex_value, family_name)
      fit <- result$fit
      if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
        stop("Age-only OOF fit failed for ", sex_value, " imp ", imp, " fold ", fold)
      }
      sc <- score_fit(fit, family_name, train, test)
      oof_rows[[oo]] <- bind_cols(test, sc) %>% mutate(fold = fold)
      oo <- oo + 1L
    }
  }
}
oof <- bind_rows(oof_rows) %>%
  mutate(
    wtmec_pool = wtmec_correct,
    low_p5 = as.integer(percentile < .05), low_p10 = as.integer(percentile < .10),
    strata_pool = factor(strata_pool), psu_pool = factor(psu_pool)
  )

oof_single <- map_dfr(1:5, function(imp) {
  map_dfr(c("Female", "Male"), function(s) {
    dat <- oof %>% filter(imputation == imp, sex == s)
    des <- make_design(dat)
    mn <- svymean(~z + low_p5 + low_p10, des, na.rm = TRUE)
    tibble(
      imputation = imp, sex = s,
      mean_z = as.numeric(coef(mn)[1]), mean_z_se = as.numeric(SE(mn)[1]),
      p5 = as.numeric(coef(mn)[2]), p5_se = as.numeric(SE(mn)[2]),
      p10 = as.numeric(coef(mn)[3]), p10_se = as.numeric(SE(mn)[3]),
      z_sd = safe_sd(des, "z")
    )
  })
})

oof_cal <- bind_rows(
  oof_single %>% transmute(imputation, sex, metric = "mean_z", estimate = mean_z, standard_error = mean_z_se),
  oof_single %>% transmute(imputation, sex, metric = "P5", estimate = p5, standard_error = p5_se),
  oof_single %>% transmute(imputation, sex, metric = "P10", estimate = p10, standard_error = p10_se)
) %>%
  pool_table(c("sex", "metric")) %>%
  left_join(oof_single %>% group_by(sex) %>% summarise(z_sd = mean(z_sd), .groups = "drop"), by = "sex")
write_csv(oof_cal, file.path(out_dir, "age_only_oof_calibration.csv"))

# -----------------------------------------------------------------------------
# Frozen temporal application of the age-only models.
# -----------------------------------------------------------------------------
scored_age_only <- list(); so <- 1L
for (sex_value in c("Female", "Male")) {
  family_name <- if (sex_value == "Female") "BCCG" else "BCT"
  for (imp in 1:5) {
    dat <- validation %>% filter(imputation == imp, as.character(sex) == sex_value)
    b <- age_only_bundles[[paste0("imp", imp, "_", sex_value)]]
    sc <- score_fit(b$model, family_name, b$train_data, dat)
    scored_age_only[[so]] <- bind_cols(dat, sc %>% rename(p_A = percentile, z_A = z))
    so <- so + 1L
  }
}
validation_scored <- bind_rows(scored_age_only) %>%
  mutate(low_p5_A = as.integer(p_A < .05), low_p10_A = as.integer(p_A < .10))

# Overall frozen temporal calibration for both references.
temporal_single <- map_dfr(1:5, function(imp) {
  map_dfr(c("Female", "Male"), function(s) {
    dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s)
    des <- make_design(dat)
    map_dfr(c("Age-only", "Age+height"), function(ref) {
      z_var <- if (ref == "Age-only") "z_A" else "z_H"
      p5_var <- if (ref == "Age-only") "low_p5_A" else "low_p5_H"
      p10_var <- if (ref == "Age-only") "low_p10_A" else "low_p10_H"
      m <- svymean(as.formula(paste0("~", z_var, "+", p5_var, "+", p10_var)), des, na.rm = TRUE)
      tibble(
        imputation = imp, sex = s, reference = ref,
        mean_z = as.numeric(coef(m)[1]), mean_z_se = as.numeric(SE(m)[1]),
        p5 = as.numeric(coef(m)[2]), p5_se = as.numeric(SE(m)[2]),
        p10 = as.numeric(coef(m)[3]), p10_se = as.numeric(SE(m)[3]),
        z_sd = safe_sd(des, z_var)
      )
    })
  })
})

temporal_cal <- bind_rows(
  temporal_single %>% transmute(imputation, sex, reference, metric = "mean_z", estimate = mean_z, standard_error = mean_z_se),
  temporal_single %>% transmute(imputation, sex, reference, metric = "P5", estimate = p5, standard_error = p5_se),
  temporal_single %>% transmute(imputation, sex, reference, metric = "P10", estimate = p10, standard_error = p10_se)
) %>%
  pool_table(c("sex", "reference", "metric")) %>%
  left_join(temporal_single %>% group_by(sex, reference) %>% summarise(z_sd = mean(z_sd), .groups = "drop"),
            by = c("sex", "reference"))
write_csv(temporal_cal, file.path(out_dir, "age_only_vs_age_height_calibration.csv"))

# Height associations for continuous z and lower-tail status.
height_single <- list(); hs <- 1L
for (imp in 1:5) {
  for (s in c("Female", "Male")) {
    dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s) %>% droplevels()
    dat$height10 <- dat$height_cm / 10
    des <- make_design(dat)

    for (ref in c("Age-only", "Age+height")) {
      z_var <- if (ref == "Age-only") "z_A" else "z_H"
      p5_var <- if (ref == "Age-only") "low_p5_A" else "low_p5_H"
      p10_var <- if (ref == "Age-only") "low_p10_A" else "low_p10_H"

      fit_z <- svyglm(as.formula(paste0(z_var, " ~ height10 + ns(age,df=3) + race + cycle")), design = des)
      fit_p5 <- svyglm(as.formula(paste0(p5_var, " ~ height10 + ns(age,df=3) + race + cycle")),
                       design = des, family = quasibinomial())
      fit_p10 <- svyglm(as.formula(paste0(p10_var, " ~ height10 + ns(age,df=3) + race + cycle")),
                        design = des, family = quasibinomial())

      height_single[[hs]] <- tibble(
        imputation = imp, sex = s, reference = ref,
        z_beta = coef(fit_z)["height10"], z_se = SE(fit_z)["height10"],
        p5_log_or = coef(fit_p5)["height10"], p5_se = SE(fit_p5)["height10"],
        p10_log_or = coef(fit_p10)["height10"], p10_se = SE(fit_p10)["height10"]
      )
      hs <- hs + 1L
    }
  }
}
height_single <- bind_rows(height_single)

continuous <- height_single %>%
  transmute(imputation, sex, reference, metric = "continuous_z_beta_per_10cm",
            estimate = z_beta, standard_error = z_se) %>%
  pool_table(c("sex", "reference", "metric"))

p5 <- height_single %>%
  transmute(imputation, sex, reference, metric = "P5_OR_per_10cm",
            estimate = p5_log_or, standard_error = p5_se) %>%
  group_by(sex, reference, metric) %>%
  group_modify(~ {
    p <- rubin_pool(.x$estimate, .x$standard_error)
    p %>% mutate(
      estimate = exp(estimate),
      standard_error = NA_real_,
      conf_low = exp(conf_low), conf_high = exp(conf_high)
    )
  }) %>% ungroup()

p10 <- height_single %>%
  transmute(imputation, sex, reference, metric = "P10_OR_per_10cm",
            estimate = p10_log_or, standard_error = p10_se) %>%
  group_by(sex, reference, metric) %>%
  group_modify(~ {
    p <- rubin_pool(.x$estimate, .x$standard_error)
    p %>% mutate(
      estimate = exp(estimate),
      standard_error = NA_real_,
      conf_low = exp(conf_low), conf_high = exp(conf_high)
    )
  }) %>% ungroup()

height_results <- bind_rows(continuous, p5, p10)
write_csv(height_results, file.path(out_dir, "age_only_vs_age_height_temporal.csv"))

# Restricted-cubic-spline audit using the same four weighted height knots as the
# final validation pipeline (P5/P35/P65/P95). The audit is secondary to the
# linear per-10-cm estimates reported in the manuscript.
weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  vapply(p, function(q) x[which(cumsum(w) / sum(w) >= q)[1]], numeric(1))
}

rcs_basis <- function(height_cm, knots_cm) {
  x <- height_cm / 10; k <- knots_cm / 10; K <- length(k)
  tp <- function(v) pmax(v, 0)^3
  cols <- lapply(1:(K - 2), function(j) {
    (tp(x - k[j]) - tp(x - k[K - 1]) * (k[K] - k[j]) / (k[K] - k[K - 1]) +
       tp(x - k[K]) * (k[K - 1] - k[j]) / (k[K] - k[K - 1])) / (k[K] - k[1])^2
  })
  out <- data.frame(h_linear = x)
  for (j in seq_along(cols)) out[[paste0("h_nl", j)]] <- cols[[j]]
  out
}

pool_multi_test <- function(betas, covs, indices) {
  m <- length(betas)
  Q <- do.call(rbind, lapply(betas, function(x) x[indices]))
  qbar <- colMeans(Q)
  U <- Reduce("+", lapply(covs, function(x) x[indices, indices, drop = FALSE])) / m
  B <- if (m > 1) cov(Q) else matrix(0, length(indices), length(indices))
  Tmat <- U + (1 + 1 / m) * B
  inv <- try(solve(Tmat), silent = TRUE)
  if (inherits(inv, "try-error")) inv <- qr.solve(Tmat)
  stat <- as.numeric(t(qbar) %*% inv %*% qbar)
  c(statistic = stat, df = length(indices), p_value = pchisq(stat, length(indices), lower.tail = FALSE))
}

spline_rows <- list(); sr <- 1L
for (s in c("Female", "Male")) {
  ref_dat <- validation_scored %>% filter(imputation == 1, as.character(sex) == s)
  knots <- weighted_quantile(ref_dat$height_cm, ref_dat$wtmec_pool, c(.05, .35, .65, .95))
  for (ref in c("Age-only", "Age+height")) {
    for (metric in c("continuous_z", "P5", "P10")) {
      betas <- list(); covs <- list()
      for (imp in 1:5) {
        dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s) %>% droplevels()
        basis <- rcs_basis(dat$height_cm, knots)
        dat$h_linear <- basis$h_linear; dat$h_nl1 <- basis$h_nl1; dat$h_nl2 <- basis$h_nl2
        outcome <- if (ref == "Age-only") {
          switch(metric, continuous_z = "z_A", P5 = "low_p5_A", P10 = "low_p10_A")
        } else {
          switch(metric, continuous_z = "z_H", P5 = "low_p5_H", P10 = "low_p10_H")
        }
        form <- as.formula(paste0(outcome, " ~ h_linear + h_nl1 + h_nl2 + ns(age,df=3) + race + cycle"))
        fit <- if (metric == "continuous_z") {
          svyglm(form, design = make_design(dat))
        } else {
          svyglm(form, design = make_design(dat), family = quasibinomial())
        }
        idx <- match(c("h_linear", "h_nl1", "h_nl2"), names(coef(fit)))
        betas[[imp]] <- coef(fit)[idx]
        covs[[imp]] <- vcov(fit)[idx, idx, drop = FALSE]
      }
      global <- pool_multi_test(betas, covs, 1:3)
      nonlin <- pool_multi_test(betas, covs, 2:3)
      spline_rows[[sr]] <- tibble(
        sex = s, reference = ref, metric = metric,
        global_p = unname(global["p_value"]),
        nonlinearity_p = unname(nonlin["p_value"]),
        knots_cm = paste(round(knots, 1), collapse = ";")
      )
      sr <- sr + 1L
    }
  }
}
spline_audit <- bind_rows(spline_rows)
write_csv(spline_audit, file.path(out_dir, "age_only_vs_age_height_spline_audit.csv"))

# Reclassification between age-only and age+height lower-tail markers.
reclass_single <- map_dfr(1:5, function(imp) {
  map_dfr(c("Female", "Male"), function(s) {
    dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s)
    dat <- dat %>% mutate(
      p5_changed = as.integer(low_p5_A != low_p5_H),
      p10_changed = as.integer(low_p10_A != low_p10_H)
    )
    des <- make_design(dat)
    m <- svymean(~p5_changed + p10_changed, des, na.rm = TRUE)
    tibble(
      imputation = imp, sex = s,
      p5 = as.numeric(coef(m)[1]), p5_se = as.numeric(SE(m)[1]),
      p10 = as.numeric(coef(m)[2]), p10_se = as.numeric(SE(m)[2])
    )
  })
})
reclass <- bind_rows(
  reclass_single %>% transmute(imputation, sex, metric = "P5_reclassification", estimate = p5, standard_error = p5_se),
  reclass_single %>% transmute(imputation, sex, metric = "P10_reclassification", estimate = p10, standard_error = p10_se)
) %>%
  pool_table(c("sex", "metric")) %>%
  mutate(percent = 100 * estimate, percent_low = 100 * pmax(0, conf_low), percent_high = 100 * pmin(1, conf_high))
write_csv(reclass, file.path(out_dir, "age_only_vs_age_height_reclassification.csv"))

# Manuscript-level summary.
h_manifest <- read_csv(file.path(repo_root, "models", "model_manifest.csv"), show_col_types = FALSE)
bic_summary <- bind_rows(
  fit_audit %>% group_by(sex) %>% summarise(reference = "Age-only", summed_BIC = sum(BIC), .groups = "drop"),
  h_manifest %>% group_by(sex) %>% summarise(reference = "Age+height", summed_BIC = sum(BIC), .groups = "drop")
)

summary_out <- bind_rows(
  bic_summary %>% transmute(sex, reference, metric = "summed_BIC", estimate = summed_BIC, conf_low = NA_real_, conf_high = NA_real_),
  height_results %>% select(sex, reference, metric, estimate, conf_low, conf_high),
  reclass %>% transmute(sex, reference = "Age-only vs Age+height", metric, estimate = percent, conf_low = percent_low, conf_high = percent_high)
) %>% arrange(sex, metric, reference)
write_csv(summary_out, file.path(out_dir, "age_only_vs_age_height_summary.csv"))

message("Age-only versus age+height sensitivity analysis completed.")
message("Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
