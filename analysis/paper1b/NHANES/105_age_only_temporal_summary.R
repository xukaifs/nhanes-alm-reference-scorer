options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gamlss); library(gamlss.dist); library(survey); library(splines)
})

set.seed(20260820)
root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
frozen_dir <- file.path(root_dir, "First Paper", "Paper1_FINAL_FROZEN_20260814")
out_dir <- file.path(root_dir, "First Paper", "age_only_comparator_20260820")
cache_dir <- file.path(out_dir, "cache")
model_dir <- file.path(out_dir, "models")
figure_dir <- file.path(out_dir, "figures")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  file.path(frozen_dir, "cache", "validation_corrected_H_18_59.rds"),
  file.path(model_dir, "age_only_18_69_Female_bundles.rds"),
  file.path(model_dir, "age_only_18_69_Male_bundles.rds"),
  file.path(cache_dir, "age_only_OOF_Female.csv.gz"),
  file.path(cache_dir, "age_only_OOF_Male.csv.gz")
)
stopifnot(all(file.exists(required)))

weighted_quantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]; o <- order(x); x <- x[o]; w <- w[o]
  vapply(p, function(q) x[which(cumsum(w) / sum(w) >= q)[1]], numeric(1))
}
make_design <- function(dat) {
  svydesign(ids = ~psu_pool, strata = ~strata_pool, weights = ~analysis_weight,
            nest = TRUE, data = dat)
}
rubin_pool <- function(q, se, null = 0) {
  ok <- is.finite(q) & is.finite(se) & se >= 0; q <- q[ok]; se <- se[ok]; m <- length(q)
  if (!m) return(tibble(m = 0, estimate = NA_real_, standard_error = NA_real_, df = NA_real_,
                        conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
  qbar <- mean(q); u <- mean(se^2); b <- if (m > 1) var(q) else 0
  tvar <- u + (1 + 1/m) * b; s <- sqrt(tvar)
  df <- if (m > 1 && b > 0) (m - 1) * (1 + u / ((1 + 1/m) * b))^2 else Inf
  crit <- qt(.975, df)
  tibble(m = m, estimate = qbar, standard_error = s, df = df,
         conf_low = qbar - crit * s, conf_high = qbar + crit * s,
         p_value = if_else(s > 0, 2 * pt(abs((qbar - null) / s), df = df, lower.tail = FALSE), NA_real_))
}
pool_table <- function(dat, groups) {
  dat %>% group_by(across(all_of(groups))) %>%
    group_modify(~rubin_pool(.x$estimate, .x$standard_error)) %>% ungroup()
}
safe_sd <- function(des, var) {
  x <- try(svyvar(as.formula(paste0("~", var)), des, na.rm = TRUE), silent = TRUE)
  if (inherits(x, "try-error")) return(NA_real_)
  sqrt(pmax(as.numeric(coef(x)[1]), 0))
}
pool_multi_test <- function(betas, covs, indices) {
  m <- length(betas); Q <- do.call(rbind, lapply(betas, function(x) x[indices])); qbar <- colMeans(Q)
  U <- Reduce("+", lapply(covs, function(x) x[indices, indices, drop = FALSE])) / m
  B <- if (m > 1) cov(Q) else matrix(0, length(indices), length(indices)); T <- U + (1 + 1/m) * B
  inv <- try(solve(T), silent = TRUE); if (inherits(inv, "try-error")) inv <- qr.solve(T)
  stat <- as.numeric(t(qbar) %*% inv %*% qbar)
  c(statistic = stat, df = length(indices), p_value = pchisq(stat, length(indices), lower.tail = FALSE))
}

predict_parameters <- function(bundle, new_data) {
  train <- as.data.frame(bundle$train_data)
  nd <- as.data.frame(new_data[, intersect(names(new_data), names(train)), drop = FALSE])
  ans <- try(suppressWarnings(predictAll(bundle$model, newdata = nd, data = train, type = "response")),
             silent = TRUE)
  if (inherits(ans, "try-error")) ans <- suppressWarnings(predictAll(bundle$model, newdata = nd,
                                                                       type = "response"))
  as.data.frame(ans)
}
score_bundle <- function(bundle, new_data) {
  pars <- predict_parameters(bundle, new_data)
  pfun <- get(paste0("p", bundle$family), asNamespace("gamlss.dist"))
  args <- list(q = new_data$alm_kg, mu = pars$mu, sigma = pars$sigma)
  if ("nu" %in% names(pars)) args$nu <- pars$nu
  if ("tau" %in% names(pars)) args$tau <- pars$tau
  p <- as.numeric(do.call(pfun, args)); p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  tibble(percentile_age_only = p, z_age_only = qnorm(p))
}

# Frozen temporal scoring without refitting or use of validation outcomes in model selection.
validation <- readRDS(file.path(frozen_dir, "cache", "validation_corrected_H_18_59.rds")) %>%
  mutate(sex = factor(as.character(sex), levels = c("Female", "Male")),
         race = factor(as.character(race)), cycle = factor(as.character(cycle)),
         strata_pool = factor(strata_pool), psu_pool = factor(psu_pool),
         analysis_weight = wtmec_pool)
score_rows <- list(); counter <- 1L
for (sex_value in c("Female", "Male")) {
  bundles <- readRDS(file.path(model_dir, paste0("age_only_18_69_", sex_value, "_bundles.rds")))
  for (imp in 1:5) {
    dat <- validation %>% filter(imputation == imp, as.character(sex) == sex_value)
    score_rows[[counter]] <- bind_cols(dat, score_bundle(bundles[[paste0("imp", imp, "_", sex_value)]], dat))
    counter <- counter + 1L
  }
}
validation_scored <- bind_rows(score_rows) %>% mutate(
  low_p5_age_only = as.numeric(percentile_age_only < .05),
  low_p10_age_only = as.numeric(percentile_age_only < .10)
)
saveRDS(validation_scored, file.path(cache_dir, "validation_age_only_vs_H.rds"))
write_csv(validation_scored %>% select(
  SEQN, imputation, cycle, sex, race, age, height_cm, height_m, bmi, alm_kg,
  analysis_weight, strata_pool, psu_pool,
  z_age_only, percentile_age_only, low_p5_age_only, low_p10_age_only,
  z_H = z_H_corrected, percentile_H = percentile_H_corrected,
  low_p5_H = low_p5_H_corrected, low_p10_H = low_p10_H_corrected
), file.path(cache_dir, "validation_age_only_vs_H.csv.gz"))

# Full-fit information criteria; no re-selection is performed.
age_fit <- map_dfr(c("Female", "Male"), ~read_csv(
  file.path(cache_dir, paste0("age_only_full_fit_", .x, ".csv")), show_col_types = FALSE))
h_fit <- map_dfr(c("Female", "Male"), ~read_csv(
  file.path(frozen_dir, "cache", paste0("corrected_H_full_fit_", .x, ".csv")), show_col_types = FALSE)) %>%
  mutate(model = "Age+height", mu_formula = if_else(sex == "Female",
    "pb(age,df=3,inter=10)+pb(height,df=3,inter=10)", "pb(age,df=3)+pb(height,df=3)"),
    weight_version = "Corrected pooled development weights")
fit_detail <- bind_rows(age_fit, h_fit %>% select(any_of(names(age_fit))))
fit_summary <- fit_detail %>% group_by(sex, model, family, mu_formula, weight_version) %>%
  summarise(imputations = n(), all_converged = all(converged), total_BIC = sum(BIC),
            mean_BIC = mean(BIC), mean_effective_df = mean(effective_df), .groups = "drop") %>%
  group_by(sex) %>% mutate(delta_total_BIC_vs_age_height = total_BIC - total_BIC[model == "Age+height"]) %>%
  ungroup()
write_csv(fit_detail, file.path(out_dir, "01_model_fit_detail.csv"))
write_csv(fit_summary, file.path(out_dir, "01b_model_fit_summary.csv"))

# Internal OOF calibration for the comparator and existing final H model.
age_oof <- map_dfr(c("Female", "Male"), ~read_csv(
  file.path(cache_dir, paste0("age_only_OOF_", .x, ".csv.gz")), show_col_types = FALSE, progress = FALSE)) %>%
  transmute(SEQN, imputation, sex, strata_pool, psu_pool, analysis_weight = wtmec_correct,
            model = "Age-only", z_score = z_age_only, low_p5 = low_p5_age_only,
            low_p10 = low_p10_age_only)
h_oof <- map_dfr(c("Female", "Male"), ~read_csv(
  file.path(frozen_dir, "cache", paste0("corrected_H_OOF_", .x, ".csv.gz")),
  show_col_types = FALSE, progress = FALSE)) %>%
  transmute(SEQN, imputation, sex, strata_pool, psu_pool, analysis_weight = wtmec_correct,
            model = "Age+height", z_score = z, low_p5, low_p10)
oof <- bind_rows(age_oof, h_oof) %>%
  mutate(strata_pool = factor(strata_pool), psu_pool = factor(psu_pool))
oof_single <- map_dfr(1:5, function(imp) map_dfr(c("Female", "Male"), function(s) map_dfr(
  c("Age-only", "Age+height"), function(mod) {
    dat <- oof %>% filter(imputation == imp, sex == s, model == mod) %>% droplevels()
    des <- make_design(dat); mn <- svymean(~z_score + low_p5 + low_p10, des, na.rm = TRUE)
    tibble(imputation = imp, sex = s, model = mod, metric = c("mean_z", "P5", "P10"),
           estimate = as.numeric(coef(mn)), standard_error = as.numeric(SE(mn)),
           z_sd = safe_sd(des, "z_score"), unweighted_n = nrow(dat))
  }
)))
oof_cal <- oof_single %>% pool_table(c("sex", "model", "metric")) %>%
  left_join(oof_single %>% group_by(sex, model, metric) %>% summarise(
    z_sd = mean(z_sd), unweighted_n = mean(unweighted_n), .groups = "drop"),
    by = c("sex", "model", "metric")) %>% mutate(
      percent = if_else(metric %in% c("P5", "P10"), 100 * estimate, NA_real_),
      percent_low = if_else(metric %in% c("P5", "P10"), 100 * pmax(conf_low, 0), NA_real_),
      percent_high = if_else(metric %in% c("P5", "P10"), 100 * pmin(conf_high, 1), NA_real_))
write_csv(oof_cal, file.path(out_dir, "02_internal_OOF_calibration.csv"))

# Temporal overall calibration for both frozen references.
temporal_long <- bind_rows(
  validation_scored %>% transmute(across(everything()), model = "Age-only",
                                  z_score = z_age_only, low_p5 = low_p5_age_only,
                                  low_p10 = low_p10_age_only),
  validation_scored %>% transmute(across(everything()), model = "Age+height",
                                  z_score = z_H_corrected, low_p5 = low_p5_H_corrected,
                                  low_p10 = low_p10_H_corrected)
)
temporal_single <- map_dfr(1:5, function(imp) map_dfr(c("Female", "Male"), function(s) map_dfr(
  c("Age-only", "Age+height"), function(mod) {
    dat <- temporal_long %>% filter(imputation == imp, as.character(sex) == s, model == mod) %>% droplevels()
    des <- make_design(dat); mn <- svymean(~z_score + low_p5 + low_p10, des, na.rm = TRUE)
    tibble(imputation = imp, sex = s, model = mod, metric = c("mean_z", "P5", "P10"),
           estimate = as.numeric(coef(mn)), standard_error = as.numeric(SE(mn)),
           z_sd = safe_sd(des, "z_score"), unweighted_n = nrow(dat))
  }
)))
temporal_cal <- temporal_single %>% pool_table(c("sex", "model", "metric")) %>%
  left_join(temporal_single %>% group_by(sex, model, metric) %>% summarise(
    z_sd = mean(z_sd), unweighted_n = mean(unweighted_n), .groups = "drop"),
    by = c("sex", "model", "metric")) %>% mutate(
      percent = if_else(metric %in% c("P5", "P10"), 100 * estimate, NA_real_),
      percent_low = if_else(metric %in% c("P5", "P10"), 100 * pmax(conf_low, 0), NA_real_),
      percent_high = if_else(metric %in% c("P5", "P10"), 100 * pmin(conf_high, 1), NA_real_))
write_csv(temporal_cal, file.path(out_dir, "03_temporal_calibration.csv"))

# Same adjusted linear and restricted-cubic-spline stature diagnostics for both models.
rcs_basis <- function(height_cm, knots_cm) {
  x <- height_cm / 10; k <- knots_cm / 10; K <- length(k); tp <- function(v) pmax(v, 0)^3
  cols <- lapply(1:(K - 2), function(j) {
    (tp(x - k[j]) - tp(x - k[K-1]) * (k[K] - k[j]) / (k[K] - k[K-1]) +
       tp(x - k[K]) * (k[K-1] - k[j]) / (k[K] - k[K-1])) / (k[K] - k[1])^2
  })
  out <- data.frame(h_linear = x)
  for (j in seq_along(cols)) out[[paste0("h_nl", j)]] <- cols[[j]]
  out
}
marginal_prediction <- function(fit, dat, height_value, knots, binary) {
  nd <- dat; nd$height_cm <- height_value; nd$height_m <- height_value / 100
  basis <- rcs_basis(nd$height_cm, knots)
  nd$h_linear <- basis$h_linear; nd$h_nl1 <- basis$h_nl1; nd$h_nl2 <- basis$h_nl2
  X <- model.matrix(delete.response(terms(fit)), nd); bn <- names(coef(fit)); miss <- setdiff(bn, colnames(X))
  if (length(miss)) X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
  X <- X[, bn, drop = FALSE]; beta <- coef(fit); V <- vcov(fit); w <- dat$analysis_weight
  eta <- as.numeric(X %*% beta)
  if (binary) {
    p <- plogis(eta); est <- weighted.mean(p, w); grad <- colSums(X * (p * (1-p) * w)) / sum(w)
  } else {
    xbar <- colSums(X * w) / sum(w); est <- sum(xbar * beta); grad <- xbar
  }
  se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
  c(estimate = est, standard_error = se)
}

height_tests <- list(); height_preds <- list(); ht <- hp <- 1L
for (s in c("Female", "Male")) {
  ref <- validation_scored %>% filter(imputation == 1, as.character(sex) == s)
  knots <- weighted_quantile(ref$height_cm, ref$analysis_weight, c(.05, .35, .65, .95))
  points <- weighted_quantile(ref$height_cm, ref$analysis_weight, c(.10, .25, .50, .75, .90))
  names(points) <- c("P10", "P25", "P50", "P75", "P90")
  for (mod in c("Age-only", "Age+height")) for (outcome in c("z", "P5", "P10")) {
    linear_q <- linear_se <- numeric(); spline_b <- list(); spline_v <- list(); pred_single <- list(); ps <- 1L
    for (imp in 1:5) {
      dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s) %>% droplevels()
      dat$outcome_value <- if (mod == "Age-only") {
        if (outcome == "z") dat$z_age_only else if (outcome == "P5") dat$low_p5_age_only else dat$low_p10_age_only
      } else {
        if (outcome == "z") dat$z_H_corrected else if (outcome == "P5") dat$low_p5_H_corrected else dat$low_p10_H_corrected
      }
      basis <- rcs_basis(dat$height_cm, knots)
      dat$h_linear <- basis$h_linear; dat$h_nl1 <- basis$h_nl1; dat$h_nl2 <- basis$h_nl2
      f_lin <- outcome_value ~ h_linear + ns(age, df = 3) + race + cycle
      f_spl <- outcome_value ~ h_linear + h_nl1 + h_nl2 + ns(age, df = 3) + race + cycle
      if (outcome == "z") {
        fit_lin <- svyglm(f_lin, design = make_design(dat)); fit_spl <- svyglm(f_spl, design = make_design(dat))
      } else {
        fit_lin <- svyglm(f_lin, design = make_design(dat), family = quasibinomial())
        fit_spl <- svyglm(f_spl, design = make_design(dat), family = quasibinomial())
      }
      linear_q[imp] <- coef(fit_lin)["h_linear"]; linear_se[imp] <- SE(fit_lin)["h_linear"]
      idx <- match(c("h_linear", "h_nl1", "h_nl2"), names(coef(fit_spl)))
      spline_b[[imp]] <- coef(fit_spl)[idx]; spline_v[[imp]] <- vcov(fit_spl)[idx, idx, drop = FALSE]
      for (pt in names(points)) {
        pr <- marginal_prediction(fit_spl, dat, points[[pt]], knots, binary = outcome != "z")
        pred_single[[ps]] <- tibble(imputation = imp, height_point = pt, height_cm = points[[pt]],
                                    estimate = pr["estimate"], standard_error = pr["standard_error"])
        ps <- ps + 1L
      }
    }
    lin <- rubin_pool(linear_q, linear_se); global <- pool_multi_test(spline_b, spline_v, 1:3)
    nonlin <- pool_multi_test(spline_b, spline_v, 2:3)
    height_tests[[ht]] <- lin %>% transmute(
      sex = s, model = mod, outcome = outcome, effect_scale = if_else(outcome == "z", "beta_per_10cm", "OR_per_10cm"),
      estimate = if_else(outcome == "z", estimate, exp(estimate)),
      conf_low = if_else(outcome == "z", conf_low, exp(conf_low)),
      conf_high = if_else(outcome == "z", conf_high, exp(conf_high)), p_value,
      spline_global_p = global["p_value"], spline_nonlinearity_p = nonlin["p_value"],
      knot_spec = paste(round(knots, 1), collapse = ";"))
    ht <- ht + 1L
    pp <- bind_rows(pred_single) %>% pool_table(c("height_point", "height_cm"))
    height_preds[[hp]] <- pp %>% transmute(
      sex = s, model = mod, outcome = outcome, height_point, height_cm,
      estimate, conf_low, conf_high, spline_global_p = global["p_value"],
      spline_nonlinearity_p = nonlin["p_value"])
    hp <- hp + 1L
  }
}
height_effects <- bind_rows(height_tests)
height_predictions <- bind_rows(height_preds)
write_csv(height_effects, file.path(out_dir, "04_temporal_height_effects.csv"))
write_csv(height_predictions, file.path(out_dir, "05_temporal_adjusted_height_predictions.csv"))

# Survey-weighted reclassification between age-only and age+height references.
reclass_single <- list(); rc <- 1L
for (imp in 1:5) for (s in c("Female", "Male")) for (tail in c("P5", "P10")) {
  dat <- validation_scored %>% filter(imputation == imp, as.character(sex) == s) %>% droplevels()
  age_low <- if (tail == "P5") dat$low_p5_age_only else dat$low_p10_age_only
  h_low <- if (tail == "P5") dat$low_p5_H_corrected else dat$low_p10_H_corrected
  dat <- dat %>% mutate(
    neither = as.numeric(age_low == 0 & h_low == 0), age_only_only = as.numeric(age_low == 1 & h_low == 0),
    age_height_only = as.numeric(age_low == 0 & h_low == 1), both = as.numeric(age_low == 1 & h_low == 1),
    discordant = as.numeric(age_low != h_low), abs_delta_z = abs(z_age_only - z_H_corrected)
  )
  des <- make_design(dat); m <- svymean(~neither + age_only_only + age_height_only + both + discordant + abs_delta_z,
                                        des, na.rm = TRUE)
  reclass_single[[rc]] <- tibble(
    imputation = imp, sex = s, tail = tail,
    metric = c("Neither low", "Age-only only low", "Age+height only low", "Both low",
               "Total reclassified", "Mean absolute delta z"),
    estimate = as.numeric(coef(m)), standard_error = as.numeric(SE(m)),
    unweighted_n = nrow(dat), z_correlation = cor(dat$z_age_only, dat$z_H_corrected, use = "complete.obs")
  ); rc <- rc + 1L
}
reclassification <- bind_rows(reclass_single) %>% pool_table(c("sex", "tail", "metric")) %>%
  left_join(bind_rows(reclass_single) %>% group_by(sex, tail, metric) %>% summarise(
    unweighted_n = mean(unweighted_n), z_correlation = mean(z_correlation), .groups = "drop"),
    by = c("sex", "tail", "metric")) %>% mutate(
      percent = if_else(metric != "Mean absolute delta z", 100 * estimate, NA_real_),
      percent_low = if_else(metric != "Mean absolute delta z", 100 * pmax(conf_low, 0), NA_real_),
      percent_high = if_else(metric != "Mean absolute delta z", 100 * pmin(conf_high, 1), NA_real_))
write_csv(reclassification, file.path(out_dir, "06_temporal_reclassification.csv"))

# Manuscript-ready compact table.
supp_height <- height_effects %>% transmute(
  sex, comparison = "Height association", metric = case_when(
    outcome == "z" ~ "Continuous z: beta per 10 cm",
    outcome == "P5" ~ "P5: odds ratio per 10 cm",
    TRUE ~ "P10: odds ratio per 10 cm"),
  model, estimate, conf_low, conf_high, p_value, spline_nonlinearity_p)
supp_cal <- temporal_cal %>% filter(metric %in% c("P5", "P10")) %>% transmute(
  sex, comparison = "Temporal prevalence", metric = paste0(metric, " weighted prevalence (%)"),
  model, estimate = percent, conf_low = percent_low, conf_high = percent_high,
  p_value, spline_nonlinearity_p = NA_real_)
supp_reclass <- reclassification %>% filter(metric == "Total reclassified") %>% transmute(
  sex, comparison = "Reclassification", metric = paste0(tail, " total reclassified (%)"),
  model = "Age-only vs age+height", estimate = percent, conf_low = percent_low,
  conf_high = percent_high, p_value, spline_nonlinearity_p = NA_real_)
supp_table <- bind_rows(supp_height, supp_cal, supp_reclass)
write_csv(supp_table, file.path(out_dir, "07_manuscript_Supplementary_Table_S15.csv"))

# Adjusted height curves for Supplementary Figure S3 (base R for dependency stability).
plot_data <- height_predictions %>% mutate(
  estimate_plot = if_else(outcome == "z", estimate, 100 * estimate),
  low_plot = if_else(outcome == "z", conf_low, 100 * pmax(conf_low, 0)),
  high_plot = if_else(outcome == "z", conf_high, 100 * pmin(conf_high, 1)))
plot_path <- file.path(figure_dir, "Supplementary_Figure_S3_age_only_vs_age_height.png")
png(plot_path, width = 2800, height = 2600, res = 320, bg = "white")
par(mfrow = c(3, 2), mar = c(4.2, 4.4, 3.0, 1.1), oma = c(1.0, 0.5, 3.0, 0.5),
    las = 1, family = "sans")
cols <- c("Age-only" = "#C75B39", "Age+height" = "#176B87")
fills <- c("Age-only" = adjustcolor("#C75B39", alpha.f = .14),
           "Age+height" = adjustcolor("#176B87", alpha.f = .14))
for (outcome in c("z", "P5", "P10")) for (s in c("Female", "Male")) {
  panel <- plot_data %>% filter(.data$outcome == .env$outcome, as.character(.data$sex) == .env$s)
  yr <- range(c(panel$low_plot, panel$high_plot), finite = TRUE)
  pad <- max(diff(yr) * .08, ifelse(outcome == "z", .05, .3)); yr <- yr + c(-pad, pad)
  plot(range(panel$height_cm), yr, type = "n", xlab = "Height (cm)",
       ylab = if (outcome == "z") "Adjusted mean z" else paste0("Adjusted ", outcome, " prevalence (%)"),
       main = s, bty = "l")
  abline(h = ifelse(outcome == "z", 0, ifelse(outcome == "P5", 5, 10)),
         col = "grey60", lty = 2, lwd = 1)
  for (mod in c("Age-only", "Age+height")) {
    z <- panel %>% filter(.data$model == .env$mod) %>% arrange(height_cm)
    polygon(c(z$height_cm, rev(z$height_cm)), c(z$low_plot, rev(z$high_plot)),
            col = fills[[mod]], border = NA)
    lines(z$height_cm, z$estimate_plot, col = cols[[mod]], lwd = 2.2)
    points(z$height_cm, z$estimate_plot, col = cols[[mod]], pch = 16, cex = .75)
  }
  if (outcome == "z" && s == "Female") {
    legend("topleft", legend = names(cols), col = cols, lwd = 2.2, pch = 16,
           bty = "n", cex = .85)
  }
}
mtext("Frozen age-only versus age-and-height ALM references", outer = TRUE,
      side = 3, line = 1.0, font = 2, cex = 1.15)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "99_session_info.txt"))
cat("Age-only comparator summary complete at", out_dir, "\n")
