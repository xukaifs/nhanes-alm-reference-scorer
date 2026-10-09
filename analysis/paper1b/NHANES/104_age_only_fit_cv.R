options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gamlss); library(gamlss.dist)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("Female", "Male")) {
  stop("Usage: Rscript 104_age_only_fit_cv.R <Female|Male>")
}
sex_value <- args[[1]]
family_name <- if (sex_value == "Female") "BCCG" else "BCT"
set.seed(20260820)

root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
frozen_dir <- file.path(root_dir, "First Paper", "Paper1_FINAL_FROZEN_20260814")
out_dir <- file.path(root_dir, "First Paper", "age_only_comparator_20260820")
cache_dir <- file.path(out_dir, "cache")
model_dir <- file.path(out_dir, "models")
log_dir <- file.path(out_dir, "logs")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

raw_path <- file.path(frozen_dir, "cache", "development_corrected_weight_narrow.csv.gz")
stopifnot(file.exists(raw_path))
log_path <- file.path(log_dir, paste0("age_only_", sex_value, ".log"))
log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(line, "\n"); cat(line, "\n", file = log_path, append = TRUE)
}

development <- read_csv(raw_path, show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    SEQN = as.numeric(SEQN), imputation = as.integer(imputation), cycle = as.character(cycle),
    examined = RIDSTATR == 2,
    sex = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ NA_character_),
    age = as.numeric(RIDAGEYR), pregnant = coalesce(RIAGENDR == 2 & RIDEXPRG == 1, FALSE),
    height_m = as.numeric(BMXHT) / 100, height_cm = as.numeric(BMXHT),
    bmi = as.numeric(BMXBMI), alm_kg = as.numeric(alm_g) / 1000,
    wtmec_correct = as.numeric(wtmec_correct),
    strata = as.character(SDMVSTRA), psu = as.character(SDMVPSU)
  ) %>%
  filter(examined, sex == sex_value, age >= 18, age <= 69, !pregnant,
         is.finite(height_m), height_m > 0, is.finite(alm_kg), alm_kg > 0,
         is.finite(wtmec_correct), wtmec_correct > 0, !is.na(strata), !is.na(psu)) %>%
  distinct(imputation, SEQN, .keep_all = TRUE) %>%
  mutate(
    bmi = if_else(is.finite(bmi), bmi, median(bmi[is.finite(bmi)], na.rm = TRUE)),
    age_band = cut(age, c(18,20,30,40,50,60,70), right = FALSE, include.lowest = TRUE,
                   labels = c("18-19","20-29","30-39","40-49","50-59","60-69")),
    strata_pool = interaction(cycle, strata, drop = TRUE),
    psu_pool = interaction(cycle, strata, psu, drop = TRUE)
  )

fit_selected <- function(train_data, mu_step = 1, algorithm = "RS", start_mode = "lm") {
  family_function <- get(family_name, envir = asNamespace("gamlss.dist"))
  formula_text <- if (sex_value == "Female") {
    "alm_kg ~ pb(age, df=3, inter=10)"
  } else {
    "alm_kg ~ pb(age, df=3)"
  }
  lm_start <- lm(alm_kg ~ age, data = train_data, weights = w_model)
  mu_start <- if (start_mode == "lm") pmax(as.numeric(predict(lm_start, train_data)), .5) else NULL
  method_call <- switch(algorithm, RS = quote(RS()), CG = quote(CG()), mixed = quote(mixed()))
  fit_call <- substitute(
    gamlss(as.formula(FORMULA_TEXT), sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
           family = family_function(), data = train_data, weights = w_model, mu.start = START,
           method = METHOD, control = gamlss.control(n.cyc = 500, trace = FALSE, mu.step = STEP)),
    list(FORMULA_TEXT = formula_text, START = mu_start, METHOD = method_call, STEP = mu_step)
  )
  fit <- eval(fit_call)
  fit$call$family <- as.call(list(as.name(family_name)))
  fit$call$data <- quote(train_data)
  fit
}

fit_with_retries <- function(train_data) {
  attempts <- tribble(
    ~mu_step, ~algorithm, ~start_mode,
    1.00, "RS", "lm", .50, "RS", "lm", .25, "RS", "lm", .10, "RS", "lm",
    .50, "mixed", "lm", .25, "mixed", "lm", .25, "CG", "lm", .50, "RS", "default"
  )
  last <- NULL
  for (i in seq_len(nrow(attempts))) {
    fit <- try(fit_selected(train_data, attempts$mu_step[i], attempts$algorithm[i], attempts$start_mode[i]),
               silent = TRUE)
    last <- fit
    if (!inherits(fit, "try-error") && isTRUE(fit$converged)) {
      return(list(fit = fit, attempt = i, mu_step = attempts$mu_step[i],
                  algorithm = attempts$algorithm[i], start_mode = attempts$start_mode[i]))
    }
  }
  list(fit = last, attempt = nrow(attempts), mu_step = NA_real_, algorithm = NA_character_,
       start_mode = NA_character_)
}

predict_parameters <- function(model, train_data, new_data) {
  new_data <- as.data.frame(new_data[, intersect(names(new_data), names(train_data)), drop = FALSE])
  result <- try(suppressWarnings(predictAll(model, newdata = new_data, data = as.data.frame(train_data),
                                             type = "response")), silent = TRUE)
  if (inherits(result, "try-error")) {
    result <- suppressWarnings(predictAll(model, newdata = new_data, type = "response"))
  }
  as.data.frame(result)
}

score_model <- function(model, train_data, new_data) {
  pars <- predict_parameters(model, train_data, new_data)
  pfun <- get(paste0("p", family_name), envir = asNamespace("gamlss.dist"))
  call_args <- list(q = new_data$alm_kg, mu = pars$mu, sigma = pars$sigma)
  if ("nu" %in% names(pars)) call_args$nu <- pars$nu
  if ("tau" %in% names(pars)) call_args$tau <- pars$tau
  prob <- as.numeric(do.call(pfun, call_args))
  prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
  tibble(percentile_age_only = prob, z_age_only = qnorm(prob))
}

make_anchors <- function(dat) {
  lm_fit <- lm(alm_kg ~ age, data = dat, weights = wtmec_correct)
  covs <- tibble(age = c(18, median(dat$age), 69))
  tibble(
    SEQN = -(1:nrow(covs)), imputation = unique(dat$imputation), cycle = "anchor", examined = TRUE,
    sex = sex_value, age = covs$age, pregnant = FALSE,
    height_m = median(dat$height_m), height_cm = 100 * median(dat$height_m),
    bmi = median(dat$bmi, na.rm = TRUE), alm_kg = pmax(as.numeric(predict(lm_fit, covs)), .5),
    wtmec_correct = mean(dat$wtmec_correct) * 1e-4, strata = "anchor",
    psu = paste0("anchor", 1:nrow(covs)),
    age_band = cut(covs$age, c(18,20,30,40,50,60,70), right = FALSE, include.lowest = TRUE,
                   labels = c("18-19","20-29","30-39","40-49","50-59","60-69")),
    strata_pool = factor("anchor"), psu_pool = factor(paste0("anchor", 1:nrow(covs)))
  )
}

# Full corrected-weight models. The family and age smooth are intentionally held
# identical to the final H model so this sensitivity isolates height conditioning.
bundle_path <- file.path(model_dir, paste0("age_only_18_69_", sex_value, "_bundles.rds"))
fit_audit_path <- file.path(cache_dir, paste0("age_only_full_fit_", sex_value, ".csv"))
bundles <- list(); full_rows <- list()
for (imp in 1:5) {
  dat <- development %>% filter(imputation == imp) %>%
    mutate(w_model = wtmec_correct / mean(wtmec_correct)) %>% droplevels()
  log_message("full imp", imp); start <- Sys.time(); result <- fit_with_retries(dat); fit <- result$fit
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
    stop("Full fit failed: ", sex_value, " imp", imp, "; ", substr(as.character(fit), 1, 500))
  }
  bundles[[paste0("imp", imp, "_", sex_value)]] <- list(
    imputation = imp, sex = sex_value, model_type = "age_only", family = family_name,
    age_df = 3, height_df = 0, sigma_spec = "constant", age_domain = c(18,69),
    train_data = dat, model = fit
  )
  full_rows[[imp]] <- tibble(
    sex = sex_value, imputation = imp, model = "Age-only", family = family_name,
    mu_formula = if_else(sex_value == "Female", "pb(age,df=3,inter=10)", "pb(age,df=3)"),
    n = nrow(dat), converged = fit$converged, iterations = fit$iter,
    effective_df = fit$df.fit, global_deviance = fit$G.deviance,
    AIC = AIC(fit), BIC = GAIC(fit, k = log(nrow(dat))), attempt = result$attempt,
    elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
    weight_version = "Corrected pooled development weights"
  )
}
saveRDS(bundles, bundle_path)
write_csv(bind_rows(full_rows), fit_audit_path)

# PSU-grouped OOF CV with the exact deterministic folds used by the final H model.
prediction_path <- file.path(cache_dir, paste0("age_only_OOF_", sex_value, ".csv.gz"))
audit_path <- file.path(cache_dir, paste0("age_only_CV_audit_", sex_value, ".csv"))
prediction_rows <- list(); audit_rows <- list(); pc <- ac <- 1L
for (imp in 1:5) {
  dat <- development %>% filter(imputation == imp) %>% droplevels()
  units <- unique(as.character(dat$psu_pool))
  set.seed(20260812 + 1000 * imp + ifelse(sex_value == "Female", 17, 0))
  fold_map <- tibble(unit = units, fold = sample(rep(1:5, length.out = length(units))))
  fold_assignment <- fold_map$fold[match(as.character(dat$psu_pool), fold_map$unit)]
  anchors <- make_anchors(dat)
  for (fold in 1:5) {
    raw_train <- dat[fold_assignment != fold, ]; test <- dat[fold_assignment == fold, ]
    train <- bind_rows(raw_train, anchors) %>%
      mutate(w_model = wtmec_correct / mean(wtmec_correct)) %>% droplevels()
    log_message("CV imp", imp, " fold", fold); start <- Sys.time()
    result <- fit_with_retries(train); fit <- result$fit
    if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
      stop("CV fit failed: ", sex_value, " imp", imp, " fold", fold)
    }
    scored <- score_model(fit, train, test)
    if (any(!is.finite(scored$z_age_only))) stop("Invalid CV prediction")
    prediction_rows[[pc]] <- bind_cols(test, scored) %>% mutate(fold = fold); pc <- pc + 1L
    audit_rows[[ac]] <- tibble(
      sex = sex_value, imputation = imp, fold = fold, train_n = nrow(raw_train), test_n = nrow(test),
      converged = fit$converged, iterations = fit$iter, attempt = result$attempt,
      elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
      valid_predictions = sum(is.finite(scored$z_age_only)), expected_predictions = nrow(test)
    ); ac <- ac + 1L
    write_csv(bind_rows(audit_rows), audit_path)
    rm(fit); gc(FALSE)
  }
}
predictions <- bind_rows(prediction_rows) %>%
  mutate(across(c(strata, psu, strata_pool, psu_pool), as.character),
         low_p5_age_only = as.numeric(percentile_age_only < .05),
         low_p10_age_only = as.numeric(percentile_age_only < .10))
write_csv(predictions, prediction_path)
write_csv(bind_rows(audit_rows), audit_path)
log_message("age-only full and CV complete")
