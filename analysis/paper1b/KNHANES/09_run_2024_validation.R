# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(gamlss)
  library(gamlss.dist)
  library(splines)
})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
path <- file.path(root, "raw_all_sav", "HN24_ALL.sav")
out <- file.path(root, "tables_2024")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
paper_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/First Paper/Final Paper/\u65b0\u5efa\u6587\u4ef6\u5939"
scorer_root <- file.path(paper_root, "nhanes-alm-reference-scorer-v2.1.0")
age_only_root <- file.path(paper_root, "age_only_comparator_20260820", "models")

wanted <- c("ID", "psu", "sex", "age", "kstrata", "wt_itvex", "HE_ht", "HE_wt", "HE_BMI",
  "DW_Lrm_LN", "DW_Rrm_LN", "DW_Llg_LN", "DW_Rlg_LN",
  "DW_Lrm_BMC", "DW_Rrm_BMC", "DW_Llg_BMC", "DW_Rlg_BMC", "DW_WBT_pFT", "DW_SBT_pFT",
  "GS_mea_r_1", "GS_mea_r_2", "GS_mea_l_1", "GS_mea_l_2")
x <- read_sav(path, col_select = all_of(wanted))
names(x) <- tolower(names(x))
x <- as.data.frame(x)
x$id <- trimws(as.character(x$id)); x$psu <- trimws(as.character(x$psu))
for (v in setdiff(names(x), c("id", "psu"))) {
  x[[v]] <- suppressWarnings(as.numeric(zap_labels(zap_missing(x[[v]]))))
}
not_special <- function(z) is.finite(z) & abs(z - 9999) > .02 & abs(z - 9999.9) > .02 &
  abs(z - 9999.99) > .02 & abs(z - 99999) > .02

ln <- c("dw_lrm_ln", "dw_rrm_ln", "dw_llg_ln", "dw_rlg_ln")
bmc <- c("dw_lrm_bmc", "dw_rrm_bmc", "dw_llg_bmc", "dw_rlg_bmc")
ln_ok <- vapply(x[ln], function(z) not_special(z) & z > 0 & z < 20000, logical(nrow(x)))
bmc_ok <- vapply(x[bmc], function(z) not_special(z) & z > 0 & z < 900, logical(nrow(x)))
x$valid_four_ln <- rowSums(ln_ok) == 4L
x$valid_four_bmc <- rowSums(bmc_ok) == 4L
x$alm_kg <- NA_real_
component_ok <- x$valid_four_ln
# The published 2024 Task Force construction is the sum of the four reported
# limb LN values, with no 0.946 correction. This reproduces its ALM mean exactly.
x$alm_kg[component_ok] <- rowSums(x[component_ok, ln, drop = FALSE]) / 1000
x$alm_kg[!is.finite(x$alm_kg) | x$alm_kg <= 5 | x$alm_kg >= 70] <- NA_real_
x$height_cm <- ifelse(is.finite(x$he_ht) & x$he_ht >= 120 & x$he_ht <= 220, x$he_ht, NA_real_)
x$height_m <- x$height_cm / 100
x$weight_kg <- ifelse(is.finite(x$he_wt) & x$he_wt >= 20 & x$he_wt <= 300, x$he_wt, NA_real_)
x$bmi <- ifelse(is.finite(x$he_bmi) & x$he_bmi >= 10 & x$he_bmi <= 80,
                x$he_bmi, x$weight_kg / x$height_m^2)
x$body_fat_pct <- ifelse(is.finite(x$dw_sbt_pft) & x$dw_sbt_pft >= 1 & x$dw_sbt_pft <= 80,
                         x$dw_sbt_pft, NA_real_)
x$sex_label <- ifelse(x$sex == 1, "Male", ifelse(x$sex == 2, "Female", NA_character_))
x$wt <- ifelse(is.finite(x$wt_itvex) & x$wt_itvex > 0, x$wt_itvex, NA_real_)

grip_vars <- c("gs_mea_r_1", "gs_mea_r_2", "gs_mea_l_1", "gs_mea_l_2")
for (v in grip_vars) x[[v]] <- ifelse(not_special(x[[v]]) & x[[v]] > 0 & x[[v]] < 100, x[[v]], NA_real_)
row_max <- function(z) apply(z, 1, function(v) if (all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE))
x$grip_right_best <- row_max(x[c("gs_mea_r_1", "gs_mea_r_2")])
x$grip_left_best <- row_max(x[c("gs_mea_l_1", "gs_mea_l_2")])
x$grip_combined <- ifelse(is.finite(x$grip_right_best) & is.finite(x$grip_left_best),
                          x$grip_right_best + x$grip_left_best, NA_real_)
x$grip_single_max <- row_max(x[grip_vars])

domain <- is.finite(x$age) & x$age >= 40 & x$age <= 69 & !is.na(x$sex_label) &
  is.finite(x$alm_kg) & is.finite(x$height_m)
x$row_key <- seq_len(nrow(x))
source(file.path(scorer_root, "R", "score_conditional_alm.R"), local = .GlobalEnv)
h_bundles <- load_alm_reference(file.path(scorer_root, "models",
  c("corrected_H_18_69_Female_bundles.rds", "corrected_H_18_69_Male_bundles.rds")))
age_bundles <- c(readRDS(file.path(age_only_root, "age_only_18_69_Female_bundles.rds")),
                 readRDS(file.path(age_only_root, "age_only_18_69_Male_bundles.rds")))
score_input <- data.frame(row_key = x$row_key[domain], sex = x$sex_label[domain], age = x$age[domain],
                          height_m = x$height_m[domain], alm_kg = x$alm_kg[domain])
h <- score_conditional_alm(score_input, bundles = h_bundles, warn = TRUE)$pooled
a <- score_conditional_alm(score_input, bundles = age_bundles, warn = FALSE)$pooled
x$z_conditional <- x$z_age_only <- NA_real_
x$scorer_caution <- FALSE
x$z_conditional[match(h$row_key, x$row_key)] <- h$model_averaged_z
x$z_age_only[match(a$row_key, x$row_key)] <- a$model_averaged_z
x$scorer_caution[match(h$row_key, x$row_key)] <- h$use_caution
x$almi <- x$alm_kg / x$height_m^2
x$alm_bmi <- x$alm_kg / x$bmi
x$low_p5 <- is.finite(x$z_conditional) & x$z_conditional < qnorm(.05)
x$low_p10 <- is.finite(x$z_conditional) & x$z_conditional < qnorm(.10)
x$low_ewgsop2 <- ifelse(x$sex_label == "Male", x$almi < 7,
                        ifelse(x$sex_label == "Female", x$almi < 5.5, NA))
x$low_fnih <- ifelse(x$sex_label == "Male", x$alm_bmi < .789,
                     ifelse(x$sex_label == "Female", x$alm_bmi < .512, NA))
x$height10 <- x$height_cm / 10
analytic <- domain & is.finite(x$z_conditional) & is.finite(x$z_age_only) & is.finite(x$bmi) &
  is.finite(x$wt) & x$wt > 0 & is.finite(x$kstrata) & !is.na(x$psu) & nzchar(x$psu)
dat <- x[analytic, ]
saveRDS(dat, file.path(root, "data", "knhanes_2024_analysis_40_69.rds"), compress = "xz")

flow <- data.frame(
  stage = c("2024 ALL records", "Age 40-69", "Four valid limb LN", "Four valid limb BMC",
            "ALM + height in scorer domain", "Final common analytic", "Combined grip measured", "Max single-hand measured",
            "Scorer recommended-range caution"),
  n = c(nrow(x), sum(x$age >= 40 & x$age <= 69, na.rm = TRUE), sum(x$valid_four_ln), sum(x$valid_four_bmc),
        sum(domain), nrow(dat), sum(is.finite(dat$grip_combined)), sum(is.finite(dat$grip_single_max)),
        sum(dat$scorer_caution, na.rm = TRUE)))
write.csv(flow, file.path(out, "01_sample_flow.csv"), row.names = FALSE)

make_design <- function(z) svydesign(ids = ~psu, strata = ~kstrata, weights = ~wt, nest = TRUE, data = z)
safe_prop <- function(design, variable) {
  values <- design$variables[[variable]]; values <- values[!is.na(values)]
  if (!length(values)) return(c(estimate = NA, lcl = NA, ucl = NA))
  if (all(!values)) return(c(estimate = 0, lcl = 0, ucl = 0))
  if (all(values)) return(c(estimate = 1, lcl = 1, ucl = 1))
  ans <- svyciprop(as.formula(paste0("~", variable)), design, method = "logit", na.rm = TRUE)
  c(estimate = coef(ans)[1], lcl = confint(ans)[1], ucl = confint(ans)[2])
}
extract_linear <- function(fit, term = "height10") {
  b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  c(estimate = b, se = se, lcl = b - qnorm(.975) * se, ucl = b + qnorm(.975) * se,
    p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}

scope_data <- list(`Primary 40-69` = dat, `Sensitivity 40-59` = dat[dat$age <= 59, ])
cal_rows <- slope_rows <- or_rows <- list()
outcomes <- c(low_p5 = "Conditional P5", low_p10 = "Conditional P10", low_ewgsop2 = "EWGSOP2", low_fnih = "FNIH")
for (scope in names(scope_data)) for (sex_value in c("Female", "Male")) {
  z <- scope_data[[scope]]; des <- subset(make_design(z), sex_label == sex_value)
  mz <- svymean(~z_conditional, des, na.rm = TRUE); p5 <- safe_prop(des, "low_p5"); p10 <- safe_prop(des, "low_p10")
  cal_rows[[length(cal_rows) + 1L]] <- data.frame(scope = scope, sex = sex_value,
    n = sum(z$sex_label == sex_value), weighted_mean_z = coef(mz)[1], mean_z_lcl = confint(mz)[1], mean_z_ucl = confint(mz)[2],
    weighted_sd_z = sqrt(as.numeric(svyvar(~z_conditional, des, na.rm = TRUE))[1]),
    p5_percent = 100*p5[1], p5_lcl = 100*p5[2], p5_ucl = 100*p5[3],
    p10_percent = 100*p10[1], p10_lcl = 100*p10[2], p10_ucl = 100*p10[3])
  for (score in c("z_age_only", "z_conditional")) {
    e <- extract_linear(svyglm(as.formula(paste0(score, " ~ height10 + age")), design = des))
    slope_rows[[length(slope_rows) + 1L]] <- data.frame(scope = scope, sex = sex_value, score = score,
      estimate = e[1], se = e[2], lcl = e[3], ucl = e[4], p_value = e[5])
  }
  for (outcome in names(outcomes)) {
    e <- extract_linear(svyglm(as.formula(paste0(outcome, " ~ height10 + age")), design = des, family = quasibinomial()))
    or_rows[[length(or_rows) + 1L]] <- data.frame(scope = scope, sex = sex_value, definition = outcomes[[outcome]],
      log_or = e[1], se = e[2], or_per_10cm = exp(e[1]), or_lcl = exp(e[3]), or_ucl = exp(e[4]), p_value = e[5])
  }
}
calibration <- do.call(rbind, cal_rows); slopes <- do.call(rbind, slope_rows); ors <- do.call(rbind, or_rows)
write.csv(calibration, file.path(out, "02_calibration.csv"), row.names = FALSE)
write.csv(slopes, file.path(out, "03_z_height_slopes.csv"), row.names = FALSE)
write.csv(ors, file.path(out, "04_low_alm_height_or.csv"), row.names = FALSE)

weighted_quantile <- function(v, w, probs) {
  ok <- is.finite(v) & is.finite(w) & w > 0; v <- v[ok]; w <- w[ok]
  ord <- order(v); v <- v[ord]; w <- w[ord]
  v[pmin(length(v), findInterval(probs, cumsum(w)/sum(w)) + 1L)]
}
dat$height_quintile <- NA_character_; cut_rows <- list()
for (sex_value in c("Female", "Male")) {
  idx <- dat$sex_label == sex_value; cuts <- weighted_quantile(dat$height_cm[idx], dat$wt[idx], c(.2,.4,.6,.8))
  dat$height_quintile[idx] <- as.character(cut(dat$height_cm[idx], c(-Inf,cuts,Inf), labels=paste0("Q",1:5), include.lowest=TRUE))
  cut_rows[[sex_value]] <- data.frame(sex=sex_value,boundary=c("Q20","Q40","Q60","Q80"),height_cm=cuts)
}
dat$height_quintile <- factor(dat$height_quintile, levels=paste0("Q",1:5)); desq <- make_design(dat)
prev_rows <- list()
for (sex_value in c("Female","Male")) for (q in paste0("Q",1:5)) for (outcome in names(outcomes)) {
  dz <- subset(desq, sex_label == sex_value & height_quintile == q); p <- safe_prop(dz,outcome)
  prev_rows[[length(prev_rows)+1L]] <- data.frame(sex=sex_value,height_quintile=q,definition=outcomes[[outcome]],
    n=sum(dat$sex_label==sex_value & dat$height_quintile==q), prevalence_percent=100*p[1],lcl=100*p[2],ucl=100*p[3])
}
write.csv(do.call(rbind,cut_rows),file.path(out,"05_height_quintile_cutpoints.csv"),row.names=FALSE)
write.csv(do.call(rbind,prev_rows),file.path(out,"06_height_quintile_prevalence.csv"),row.names=FALSE)

discordance_levels <- c("Neither low", "Conventional-only low", "Conditional-only low", "Both low")
add_discordance <- function(z, conventional) {
  z$discordance <- factor(ifelse(!z$low_p5 & !z[[conventional]], "Neither low",
    ifelse(!z$low_p5 & z[[conventional]], "Conventional-only low",
      ifelse(z$low_p5 & !z[[conventional]], "Conditional-only low", "Both low"))), levels=discordance_levels)
  z
}
disc_rows <- list()
for (sex_value in c("Female","Male")) for (comparison in c("FNIH","EWGSOP2")) {
  conventional <- if(comparison=="FNIH") "low_fnih" else "low_ewgsop2"
  z <- add_discordance(dat[dat$sex_label==sex_value,],conventional); dz <- make_design(z)
  for (g in discordance_levels) {
    dz$variables$group_indicator <- dz$variables$discordance == g; p <- safe_prop(dz,"group_indicator")
    disc_rows[[length(disc_rows)+1L]] <- data.frame(sex=sex_value,comparison=comparison,group=g,
      n=sum(z$discordance==g),weighted_percent=100*p[1],lcl=100*p[2],ucl=100*p[3])
  }
}
write.csv(do.call(rbind,disc_rows),file.path(out,"07_classification_discordance.csv"),row.names=FALSE)

# Grip-strength characterization: primary bilateral sum and maximum single-hand sensitivity.
marginal_vector <- function(fit, z, group_value) {
  nd <- z; nd$discordance <- factor(group_value, levels=levels(z$discordance))
  mm <- model.matrix(delete.response(terms(fit)), nd); bn <- names(coef(fit)); miss <- setdiff(bn,colnames(mm))
  if(length(miss)) mm <- cbind(mm,matrix(0,nrow(mm),length(miss),dimnames=list(NULL,miss)))
  mm <- mm[,bn,drop=FALSE]; as.numeric(colSums(mm*z$wt)/sum(z$wt))
}
model_specs <- list(Raw="1", M1="ns(age,df=3)+sex_label+height_m",
                    M2="ns(age,df=3)+sex_label+height_m+bmi",
                    M3="ns(age,df=3)+sex_label+height_m+bmi+body_fat_pct")
grip_rows <- flow_rows <- list()
for (grip in c("grip_combined","grip_single_max")) for (domain_name in c("Full sample","Normal BMI")) {
  for (comparison in c("FNIH","EWGSOP2")) for (model_name in names(model_specs)) {
    conventional <- if(comparison=="FNIH") "low_fnih" else "low_ewgsop2"
    z <- dat
    if(domain_name=="Normal BMI") z <- z[is.finite(z$bmi) & z$bmi>=18.5 & z$bmi<25,]
    z <- add_discordance(z,conventional)
    needed <- is.finite(z[[grip]]) & is.finite(z$age) & is.finite(z$height_m) & !is.na(z$sex_label)
    if(model_name %in% c("M2","M3")) needed <- needed & is.finite(z$bmi)
    if(model_name=="M3") needed <- needed & is.finite(z$body_fat_pct)
    z <- droplevels(z[needed,]); present <- discordance_levels[discordance_levels %in% as.character(unique(z$discordance))]
    z$discordance <- factor(as.character(z$discordance),levels=present)
    formula <- as.formula(paste0(grip," ~ discordance + ",model_specs[[model_name]]))
    fit <- svyglm(formula,design=make_design(z)); beta <- coef(fit); V <- vcov(fit)
    vectors <- setNames(lapply(present,function(g)marginal_vector(fit,z,g)),present); ref <- vectors[["Neither low"]]
    flow_rows[[length(flow_rows)+1L]] <- data.frame(grip=grip,domain=domain_name,comparison=comparison,model=model_name,
      analysis_n=nrow(z),neither_n=sum(z$discordance=="Neither low"),conventional_only_n=sum(z$discordance=="Conventional-only low"),
      conditional_only_n=sum(z$discordance=="Conditional-only low"),both_n=sum(z$discordance=="Both low"))
    for(g in present) {
      v <- vectors[[g]]; cvec <- v-ref; se_mean <- sqrt(as.numeric(t(v)%*%V%*%v)); se_diff <- sqrt(as.numeric(t(cvec)%*%V%*%cvec))
      zd <- z; zd$group_indicator <- zd$discordance==g; p <- safe_prop(make_design(zd),"group_indicator")
      grip_rows[[length(grip_rows)+1L]] <- data.frame(grip=grip,domain=domain_name,comparison=comparison,model=model_name,group=g,
        analysis_n=nrow(z),group_n=sum(z$discordance==g),weighted_percent=100*p[1],
        adjusted_mean=sum(v*beta),mean_lcl=sum(v*beta)-qnorm(.975)*se_mean,mean_ucl=sum(v*beta)+qnorm(.975)*se_mean,
        difference_vs_neither=sum(cvec*beta),diff_lcl=sum(cvec*beta)-qnorm(.975)*se_diff,diff_ucl=sum(cvec*beta)+qnorm(.975)*se_diff,
        p_value=ifelse(g=="Neither low",NA,2*pnorm(abs(sum(cvec*beta)/se_diff),lower.tail=FALSE)))
    }
  }
}
write.csv(do.call(rbind,grip_rows),file.path(out,"08_grip_discordance_models.csv"),row.names=FALSE)
write.csv(do.call(rbind,flow_rows),file.path(out,"09_grip_model_flow.csv"),row.names=FALSE)

# Direct era comparison restricted to the common 40-69 age band.
old <- readRDS(file.path(root,"data","knhanes_2008_2011_analysis_19_69_official_weight.rds"))
old <- old[old$age>=40,]
era_rows <- list()
for (era in c("2008-2011","2024")) for (sex_value in c("Female","Male")) {
  z <- if(era=="2008-2011") old else dat; dz <- subset(if(era=="2008-2011") {
    svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wt_pool,nest=TRUE,data=z)
  } else make_design(z), sex_label==sex_value)
  mz <- svymean(~z_conditional,dz,na.rm=TRUE); p5 <- safe_prop(dz,"low_p5")
  fit <- svyglm(z_conditional~I(height_cm/10)+age,design=dz); term <- "I(height_cm/10)"; e <- extract_linear(fit,term)
  era_rows[[length(era_rows)+1L]] <- data.frame(era=era,sex=sex_value,n=sum(z$sex_label==sex_value),
    mean_z=coef(mz)[1],sd_z=sqrt(as.numeric(svyvar(~z_conditional,dz,na.rm=TRUE))[1]),p5_percent=100*p5[1],
    conditional_z_slope=e[1],slope_lcl=e[3],slope_ucl=e[4])
}
write.csv(do.call(rbind,era_rows),file.path(out,"10_era_comparison_40_69.csv"),row.names=FALSE)

cat("KNHANES 2024 validation complete.\n")
print(flow,row.names=FALSE)
print(calibration,row.names=FALSE)
print(subset(slopes,scope=="Primary 40-69"),row.names=FALSE)
print(subset(ors,scope=="Primary 40-69"),row.names=FALSE)
print(do.call(rbind,era_rows),row.names=FALSE)
