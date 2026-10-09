# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({
  library(haven)
  library(survey)
})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
raw_all <- file.path(root, "raw_all_sav")
raw_dxa <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/Knhanes\u6570\u636e\u5e93"
out <- file.path(root, "tables_official_weight")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

years <- 2008:2011
all_files <- setNames(file.path(raw_all, c("HN08_ALL.sav", "HN09_all.sav", "HN10_all.sav", "HN11_all.sav")), years)
dxa_files <- setNames(file.path(raw_dxa, sprintf("hn%02d_dxa.sav", 8:11)), years)
weight_map <- c(`2008` = "wt_ex1", `2009` = "wt_dw", `2010` = "wt_itvex", `2011` = "wt_ex1")
pool_numerator <- c(`2008` = 108, `2009` = 199, `2010` = 192, `2011` = 80)

read_selected <- function(path, wanted_lower) {
  header <- read_sav(path, n_max = 1)
  selected <- names(header)[tolower(names(header)) %in% wanted_lower]
  required_missing <- setdiff(wanted_lower, tolower(selected))
  required_missing <- setdiff(required_missing, c("wt_ex1", "wt_dw", "wt_itvex"))
  if (length(required_missing)) stop("Missing variables in ", basename(path), ": ", paste(required_missing, collapse = ", "))
  x <- read_sav(path, col_select = all_of(selected))
  names(x) <- tolower(names(x))
  as.data.frame(x)
}

clean_num <- function(x) suppressWarnings(as.numeric(zap_labels(zap_missing(x))))
not_special <- function(x) is.finite(x) & abs(x - 9999) > .02 & abs(x - 9999.9) > .02 &
  abs(x - 9999.99) > .02 & abs(x - 99999) > .02

all_wanted <- c("id", "psu", "sex", "age", "kstrata", "he_ht", "he_wt", "wt_ex1", "wt_dw", "wt_itvex")
dxa_wanted <- c("id", "dw_lrm_ln", "dw_rrm_ln", "dw_llg_ln", "dw_rlg_ln",
                "dw_lrm_bmc", "dw_rrm_bmc", "dw_llg_bmc", "dw_rlg_bmc")

parts <- list()
for (year in years) {
  a <- read_selected(all_files[[as.character(year)]], all_wanted)
  d <- read_selected(dxa_files[[as.character(year)]], dxa_wanted)
  a$id <- trimws(as.character(a$id)); d$id <- trimws(as.character(d$id))
  x <- merge(d, a, by = "id", all.x = TRUE, sort = FALSE)
  x$year <- year
  for (v in setdiff(names(x), c("id", "psu"))) x[[v]] <- clean_num(x[[v]])

  ln <- c("dw_lrm_ln", "dw_rrm_ln", "dw_llg_ln", "dw_rlg_ln")
  bmc <- c("dw_lrm_bmc", "dw_rrm_bmc", "dw_llg_bmc", "dw_rlg_bmc")
  ln_ok <- vapply(x[ln], function(z) not_special(z) & z > 0 & z < 20000, logical(nrow(x)))
  bmc_ok <- vapply(x[bmc], function(z) not_special(z) & z > 0 & z < 900, logical(nrow(x)))
  component_ok <- rowSums(ln_ok) == 4L & rowSums(bmc_ok) == 4L
  x$alm_soft_kg <- NA_real_
  x$alm_soft_kg[component_ok] <- (rowSums(x[component_ok, ln, drop = FALSE]) -
                                   rowSums(x[component_ok, bmc, drop = FALSE])) / 1000
  x$alm_corrected_kg <- x$alm_soft_kg * .946
  x$alm_ln_sum_corrected_kg <- NA_real_
  x$alm_ln_sum_corrected_kg[rowSums(ln_ok) == 4L] <- rowSums(x[rowSums(ln_ok) == 4L, ln, drop = FALSE]) / 1000 * .946
  x$height_cm <- ifelse(is.finite(x$he_ht) & x$he_ht >= 120 & x$he_ht <= 220, x$he_ht, NA_real_)
  x$weight_kg <- ifelse(is.finite(x$he_wt) & x$he_wt >= 20 & x$he_wt <= 300, x$he_wt, NA_real_)
  x$female <- as.numeric(x$sex == 2)
  x$wt_year <- x[[weight_map[[as.character(year)]]]]
  x$wt_pool <- x$wt_year * pool_numerator[[as.character(year)]] / 579
  x$strata_pool <- interaction(x$year, x$kstrata, drop = TRUE)
  x$psu_pool <- interaction(x$year, x$psu, drop = TRUE)
  parts[[as.character(year)]] <- x[, c("id", "year", "age", "female", "height_cm", "weight_kg",
    "alm_soft_kg", "alm_corrected_kg", "alm_ln_sum_corrected_kg", "wt_year", "wt_pool",
    "kstrata", "psu", "strata_pool", "psu_pool")]
  rm(a, d, x); gc()
}
dat <- do.call(rbind, parts)

criteria <- list(
  `Legacy subtract-BMC ALM valid` = with(dat, is.finite(age) & age >= 50 & is.finite(alm_corrected_kg)),
  `Task Force LN*0.946 ALM valid` = with(dat, is.finite(age) & age >= 50 & is.finite(alm_ln_sum_corrected_kg)),
  `Task Force sample + survey design` = with(dat, is.finite(age) & age >= 50 & is.finite(alm_ln_sum_corrected_kg) &
                                                 is.finite(wt_pool) & wt_pool > 0 & is.finite(kstrata) &
                                                 !is.na(psu) & nzchar(psu)),
  `Plus height and weight` = with(dat, is.finite(age) & age >= 50 & is.finite(alm_ln_sum_corrected_kg) &
                                      is.finite(height_cm) & is.finite(weight_kg)),
  `Complete-case survey design` = with(dat, is.finite(age) & age >= 50 & is.finite(alm_ln_sum_corrected_kg) &
                                  is.finite(height_cm) & is.finite(weight_kg) & is.finite(wt_pool) & wt_pool > 0 &
                                  is.finite(kstrata) & !is.na(psu) & nzchar(psu))
)
flow <- do.call(rbind, lapply(names(criteria), function(nm) {
  keep <- criteria[[nm]]
  data.frame(criterion = nm, n = sum(keep, na.rm = TRUE), female_n = sum(keep & dat$female, na.rm = TRUE),
             male_n = sum(keep & !dat$female, na.rm = TRUE))
}))
write.csv(flow, file.path(out, "08_taskforce_checksum_flow.csv"), row.names = FALSE)

keep <- criteria[["Task Force sample + survey design"]]
chk <- dat[keep, ]
design <- svydesign(ids = ~psu_pool, strata = ~strata_pool, weights = ~wt_pool, nest = TRUE, data = chk)
means <- svymean(~age + female + height_cm + weight_kg + alm_ln_sum_corrected_kg + alm_corrected_kg + alm_soft_kg,
                 design, na.rm = TRUE)
ci <- confint(means)
published <- c(age = 62.1, female = .539, height_cm = 159.6, weight_kg = 61.4,
               alm_ln_sum_corrected_kg = 17.2, alm_corrected_kg = NA_real_, alm_soft_kg = NA_real_)
checksum <- data.frame(
  metric = names(coef(means)), estimate = as.numeric(coef(means)),
  lcl = ci[, 1], ucl = ci[, 2], published_target = unname(published[names(coef(means))])
)
checksum$difference <- checksum$estimate - checksum$published_target
checksum$estimate[checksum$metric == "female"] <- 100 * checksum$estimate[checksum$metric == "female"]
checksum$lcl[checksum$metric == "female"] <- 100 * checksum$lcl[checksum$metric == "female"]
checksum$ucl[checksum$metric == "female"] <- 100 * checksum$ucl[checksum$metric == "female"]
checksum$published_target[checksum$metric == "female"] <- 100 * checksum$published_target[checksum$metric == "female"]
checksum$difference[checksum$metric == "female"] <- 100 * checksum$difference[checksum$metric == "female"]
write.csv(checksum, file.path(out, "09_taskforce_checksum_estimates.csv"), row.names = FALSE)

year_counts <- aggregate(id ~ year, data = chk, FUN = length)
names(year_counts)[2] <- "n"
write.csv(year_counts, file.path(out, "10_taskforce_checksum_year_counts.csv"), row.names = FALSE)
saveRDS(chk, file.path(root, "data", "knhanes_2008_2011_taskforce_checksum_age50plus.rds"), compress = "xz")

cat("Task Force checksum complete.\n")
print(flow, row.names = FALSE)
print(checksum, row.names = FALSE)
print(year_counts, row.names = FALSE)
