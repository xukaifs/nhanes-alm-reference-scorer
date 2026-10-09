# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
tab <- file.path(root, "tables")

required <- c(
  "00_limb_composition_identity_audit.csv", "01_sample_flow.csv", "01b_overall_flow.csv",
  "02_qc_weighted_distributions.csv", "03_frozen_calibration.csv", "04_z_height_slopes.csv",
  "05_low_alm_height_or.csv", "06a_height_quintile_cutpoints.csv",
  "06b_height_quintile_prevalence.csv", "07_classification_discordance.csv",
  "08a_uncorrected_alm_calibration.csv", "08b_uncorrected_alm_z_height_slope.csv",
  "08c_uncorrected_alm_height_or.csv", "99_software_versions.csv"
)
stopifnot(all(file.exists(file.path(tab, required))))

flow <- read.csv(file.path(tab, "01_sample_flow.csv"))
stopifnot(
  identical(flow$year, 2008:2011),
  all(flow$dxa_n == flow$matched_all_n),
  sum(flow$final_common_analytic_19_69_n) == 16140,
  sum(flow$scorer_domain_valid_n) == 16146
)

identity <- read.csv(file.path(tab, "00_limb_composition_identity_audit.csv"))
stopifnot(
  max(identity$mean_abs_mass_minus_ln_ft, na.rm = TRUE) < 1e-4,
  min(identity$mean_abs_mass_minus_ln_ft_bmc, na.rm = TRUE) > 100
)

cal <- read.csv(file.path(tab, "03_frozen_calibration.csv"))
primary_cal <- subset(cal, scope == "Primary 19-69 pooled")
stopifnot(nrow(primary_cal) == 2, all(is.finite(primary_cal$weighted_mean_z)), all(is.finite(primary_cal$p5_percent)))

slopes <- read.csv(file.path(tab, "04_z_height_slopes.csv"))
primary_slopes <- subset(slopes, scope == "Primary 19-69 pooled")
stopifnot(
  nrow(primary_slopes) == 4,
  all(abs(subset(primary_slopes, score == "z_age_only")$estimate) > 0.5),
  all(abs(subset(primary_slopes, score == "z_conditional")$estimate) < 0.2)
)

ors <- read.csv(file.path(tab, "05_low_alm_height_or.csv"))
primary_ors <- subset(ors, scope == "Primary 19-69 pooled")
stopifnot(nrow(primary_ors) == 8, all(is.finite(primary_ors$or_per_10cm)))

prev <- read.csv(file.path(tab, "06b_height_quintile_prevalence.csv"))
stopifnot(nrow(prev) == 40, all(prev$prevalence_percent >= 0 & prev$prevalence_percent <= 100))

discordance <- read.csv(file.path(tab, "07_classification_discordance.csv"))
stopifnot(nrow(discordance) == 16, all(discordance$weighted_percent >= 0 & discordance$weighted_percent <= 100))

cat("All KNHANES output verification checks passed.\n")
