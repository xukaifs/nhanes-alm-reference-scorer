# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({
  library(survey)
  library(gamlss)
  library(gamlss.dist)
})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
paper_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/First Paper/Final Paper/\u65b0\u5efa\u6587\u4ef6\u5939"
scorer_root <- file.path(paper_root, "nhanes-alm-reference-scorer-v2.1.0")

dat <- readRDS(file.path(root, "data", "knhanes_2008_2011_analysis_19_69.rds"))
dat$row_key <- seq_len(nrow(dat))
source(file.path(scorer_root, "R", "score_conditional_alm.R"), local = .GlobalEnv)
model_files <- file.path(scorer_root, "models", c(
  "corrected_H_18_69_Female_bundles.rds", "corrected_H_18_69_Male_bundles.rds"
))
bundles <- load_alm_reference(model_files)
score_input <- data.frame(
  row_key = dat$row_key, sex = dat$sex_label, age = dat$age,
  height_m = dat$height_m, alm_kg = dat$alm_raw_kg
)
scored <- score_conditional_alm(score_input, bundles = bundles, warn = FALSE)$pooled
dat$z_uncorrected <- scored$model_averaged_z[match(dat$row_key, scored$row_key)]
dat$height10 <- dat$height_cm / 10
dat$almi_uncorrected <- dat$alm_raw_kg / dat$height_m^2
dat$alm_bmi_uncorrected <- dat$alm_raw_kg / dat$bmi
dat$low_p5_uncorrected <- dat$z_uncorrected < qnorm(.05)
dat$low_p10_uncorrected <- dat$z_uncorrected < qnorm(.10)
dat$low_ewgsop2_uncorrected <- ifelse(dat$sex_label == "Male", dat$almi_uncorrected < 7,
                                      dat$almi_uncorrected < 5.5)
dat$low_fnih_uncorrected <- ifelse(dat$sex_label == "Male", dat$alm_bmi_uncorrected < .789,
                                   dat$alm_bmi_uncorrected < .512)

make_design <- function(x) svydesign(
  ids = ~psu_pool, strata = ~strata_pool, weights = ~wt_pool,
  nest = TRUE, data = x
)
design <- make_design(dat)

safe_prop <- function(d, variable) {
  ans <- svyciprop(as.formula(paste0("~", variable)), d, method = "logit", na.rm = TRUE)
  c(est = coef(ans)[1], lcl = confint(ans)[1], ucl = confint(ans)[2])
}

cal_rows <- list()
slope_rows <- list()
or_rows <- list()
outcomes <- c(low_p5_uncorrected = "Conditional P5", low_p10_uncorrected = "Conditional P10",
              low_ewgsop2_uncorrected = "EWGSOP2", low_fnih_uncorrected = "FNIH")
for (sex_value in c("Female", "Male")) {
  d <- subset(design, sex_label == sex_value)
  mz <- svymean(~z_uncorrected, d)
  p5 <- safe_prop(d, "low_p5_uncorrected")
  p10 <- safe_prop(d, "low_p10_uncorrected")
  cal_rows[[sex_value]] <- data.frame(
    sex = sex_value, n = sum(dat$sex_label == sex_value),
    mean_z = coef(mz)[1], mean_z_lcl = confint(mz)[1], mean_z_ucl = confint(mz)[2],
    sd_z = sqrt(as.numeric(svyvar(~z_uncorrected, d))[1]),
    p5_percent = 100 * p5[1], p5_lcl = 100 * p5[2], p5_ucl = 100 * p5[3],
    p10_percent = 100 * p10[1], p10_lcl = 100 * p10[2], p10_ucl = 100 * p10[3]
  )
  fit_slope <- svyglm(z_uncorrected ~ height10 + age + factor(year), design = d)
  b <- coef(fit_slope)["height10"]
  se <- sqrt(vcov(fit_slope)["height10", "height10"])
  slope_rows[[sex_value]] <- data.frame(
    sex = sex_value, slope_per_10cm = b, lcl = b - qnorm(.975) * se,
    ucl = b + qnorm(.975) * se, p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE)
  )
  for (outcome in names(outcomes)) {
    fit <- svyglm(as.formula(paste0(outcome, " ~ height10 + age + factor(year)")),
                  design = d, family = quasibinomial())
    b <- coef(fit)["height10"]
    se <- sqrt(vcov(fit)["height10", "height10"])
    or_rows[[length(or_rows) + 1L]] <- data.frame(
      sex = sex_value, definition = outcomes[[outcome]], or_per_10cm = exp(b),
      lcl = exp(b - qnorm(.975) * se), ucl = exp(b + qnorm(.975) * se),
      p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE)
    )
  }
}

write.csv(do.call(rbind, cal_rows), file.path(root, "tables", "08a_uncorrected_alm_calibration.csv"), row.names = FALSE)
write.csv(do.call(rbind, slope_rows), file.path(root, "tables", "08b_uncorrected_alm_z_height_slope.csv"), row.names = FALSE)
write.csv(do.call(rbind, or_rows), file.path(root, "tables", "08c_uncorrected_alm_height_or.csv"), row.names = FALSE)
print(do.call(rbind, cal_rows), row.names = FALSE)
print(do.call(rbind, slope_rows), row.names = FALSE)
print(do.call(rbind, or_rows), row.names = FALSE)
