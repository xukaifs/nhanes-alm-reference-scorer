# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({library(survey); library(splines)})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
out <- file.path(root, "tables_2024")
dat <- readRDS(file.path(root, "data", "knhanes_2024_analysis_40_69.rds"))
dat$sex_label <- factor(dat$sex_label, levels = c("Female", "Male"))
age_basis <- ns(dat$age, df=3)
dat$age_ns1 <- age_basis[,1]; dat$age_ns2 <- age_basis[,2]; dat$age_ns3 <- age_basis[,3]
make_design <- function(z) svydesign(ids=~psu, strata=~kstrata, weights=~wt, nest=TRUE, data=z)

model_specs <- list(
  M0 = "z_conditional",
  M1 = "z_conditional + age_ns1 + age_ns2 + age_ns3 + sex_label + height_m",
  M2 = "z_conditional + age_ns1 + age_ns2 + age_ns3 + sex_label + height_m + bmi",
  M3 = "z_conditional + age_ns1 + age_ns2 + age_ns3 + sex_label + height_m + bmi + body_fat_pct"
)
outcomes <- c(grip_combined="Bilateral combined grip", grip_single_max="Max single-hand grip")

coef_row <- function(fit, term, outcome, model, exposure="conditional z") {
  b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term,term])
  data.frame(outcome=outcomes[[outcome]], model=model, exposure=exposure,
    n=nrow(fit$survey.design$variables), beta_per_1z=b, se=se,
    lcl=b-qnorm(.975)*se, ucl=b+qnorm(.975)*se,
    p_value=2*pnorm(abs(b/se),lower.tail=FALSE), stringsAsFactors=FALSE)
}

rows <- list()
for(outcome in names(outcomes)) {
  models_to_run <- if(outcome=="grip_combined") names(model_specs) else c("M1","M2","M3")
  for(model in models_to_run) {
    required <- c(outcome,"z_conditional","age_ns1","age_ns2","age_ns3","sex_label","height_m","wt","psu","kstrata")
    if(model %in% c("M2","M3")) required <- c(required,"bmi")
    if(model=="M3") required <- c(required,"body_fat_pct")
    z <- dat[complete.cases(dat[required]),]
    cat("Fitting",outcome,model,"N=",nrow(z),"sex=",paste(names(table(z$sex_label)),table(z$sex_label),collapse=","),"\n")
    fit <- svyglm(as.formula(paste(outcome,"~",model_specs[[model]])),design=make_design(z))
    row <- coef_row(fit,"z_conditional",outcome,model)
    mg <- svymean(as.formula(paste0("~",outcome)),make_design(z),na.rm=TRUE)
    row$weighted_mean_grip <- coef(mg)[1]
    row$weighted_mean_grip_lcl <- confint(mg)[1]
    row$weighted_mean_grip_ucl <- confint(mg)[2]
    rows[[length(rows)+1L]] <- row
  }
}
main <- do.call(rbind,rows)
write.csv(main,file.path(out,"14_continuous_z_grip_models.csv"),row.names=FALSE)

# M3 sex interaction for both grip definitions. Sex-specific slopes are emitted only when P<0.05.
interaction_rows <- sex_rows <- list()
for(outcome in names(outcomes)) {
  required <- c(outcome,"z_conditional","age_ns1","age_ns2","age_ns3","sex_label","height_m","bmi","body_fat_pct","wt","psu","kstrata")
  z <- dat[complete.cases(dat[required]),]
  cat("Fitting interaction",outcome,"N=",nrow(z),"sex=",paste(names(table(z$sex_label)),table(z$sex_label),collapse=","),"\n")
  fit <- svyglm(as.formula(paste(outcome,
    "~ z_conditional*sex_label + age_ns1 + age_ns2 + age_ns3 + height_m + bmi + body_fat_pct")),design=make_design(z))
  term <- "z_conditional:sex_labelMale"; b <- coef(fit)[term]; se <- sqrt(vcov(fit)[term,term])
  p <- 2*pnorm(abs(b/se),lower.tail=FALSE)
  interaction_rows[[length(interaction_rows)+1L]] <- data.frame(outcome=outcomes[[outcome]],model="M3",
    n=nrow(z),interaction_beta=b,se=se,lcl=b-qnorm(.975)*se,ucl=b+qnorm(.975)*se,p_interaction=p)
  if(is.finite(p) && p<.05) {
    V <- vcov(fit); bw <- coef(fit)["z_conditional"]; sew <- sqrt(V["z_conditional","z_conditional"])
    bm <- bw+b; sem <- sqrt(V["z_conditional","z_conditional"]+V[term,term]+2*V["z_conditional",term])
    sex_rows[[length(sex_rows)+1L]] <- data.frame(outcome=outcomes[[outcome]],model="M3 interaction",sex=c("Female","Male"),n=c(sum(z$sex_label=="Female"),sum(z$sex_label=="Male")),
      beta_per_1z=c(bw,bm),se=c(sew,sem),lcl=c(bw-qnorm(.975)*sew,bm-qnorm(.975)*sem),
      ucl=c(bw+qnorm(.975)*sew,bm+qnorm(.975)*sem),
      p_value=c(2*pnorm(abs(bw/sew),lower.tail=FALSE),2*pnorm(abs(bm/sem),lower.tail=FALSE)))
  }
}
interaction_table <- do.call(rbind,interaction_rows)
write.csv(interaction_table,file.path(out,"15_z_sex_interaction.csv"),row.names=FALSE)
if(length(sex_rows)) write.csv(do.call(rbind,sex_rows),file.path(out,"15b_sex_specific_slopes_if_interaction.csv"),row.names=FALSE)

# Sex-specific survey-weighted centering; a pure location shift must not alter the M3 slope.
des_all <- make_design(dat)
means <- sapply(levels(dat$sex_label),function(s)unname(coef(svymean(~z_conditional,subset(des_all,sex_label==s),na.rm=TRUE))[1]))
dat$z_korea_centered <- dat$z_conditional-unname(means[as.character(dat$sex_label)])
z <- dat[complete.cases(dat[c("grip_combined","z_conditional","z_korea_centered","age_ns1","age_ns2","age_ns3","sex_label","height_m","bmi","body_fat_pct","wt","psu","kstrata")]),]
fit_raw <- svyglm(grip_combined~z_conditional+age_ns1+age_ns2+age_ns3+sex_label+height_m+bmi+body_fat_pct,design=make_design(z))
fit_ctr <- svyglm(grip_combined~z_korea_centered+age_ns1+age_ns2+age_ns3+sex_label+height_m+bmi+body_fat_pct,design=make_design(z))
raw <- coef_row(fit_raw,"z_conditional","grip_combined","M3","conditional z")
ctr <- coef_row(fit_ctr,"z_korea_centered","grip_combined","M3","sex-centered conditional z")
centered <- rbind(raw,ctr); centered$sex_weighted_mean_female <- means["Female"]
centered$sex_weighted_mean_male <- means["Male"]
centered$slope_difference_from_raw <- centered$beta_per_1z-raw$beta_per_1z
write.csv(centered,file.path(out,"16_korea_centered_z_sensitivity.csv"),row.names=FALSE)

# Primary nonlinear test: combined grip, M3 covariates. A natural/restricted cubic
# spline with 3 total df is decomposed into its linear term plus 2 nonlinear df.
z <- dat[complete.cases(dat[c("grip_combined","z_conditional","age_ns1","age_ns2","age_ns3","sex_label","height_m","bmi","body_fat_pct","wt","psu","kstrata")]),]
wquant <- function(v,w,p){o<-order(v);v<-v[o];w<-w[o];v[pmin(length(v),findInterval(p,cumsum(w)/sum(w))+1L)]}
internal <- wquant(z$z_conditional,z$wt,c(.35,.65)); boundary <- wquant(z$z_conditional,z$wt,c(.05,.95))
B <- ns(z$z_conditional,knots=internal,Boundary.knots=boundary,intercept=FALSE)
# Remove the intercept/linear subspace, leaving two independent nonlinear directions.
R <- qr.resid(qr(cbind(1,z$z_conditional)),B); qR <- qr(R,tol=1e-9)
if(qR$rank!=2L) stop("Unexpected nonlinear spline rank: ",qR$rank)
keep_cols <- qR$pivot[seq_len(qR$rank)]; z$z_nl1 <- R[,keep_cols[1]]; z$z_nl2 <- R[,keep_cols[2]]
fit_spline <- svyglm(grip_combined~z_conditional+z_nl1+z_nl2+age_ns1+age_ns2+age_ns3+sex_label+height_m+bmi+body_fat_pct,design=make_design(z))
overall <- regTermTest(fit_spline,~z_conditional+z_nl1+z_nl2,method="Wald")
nonlinear <- regTermTest(fit_spline,~z_nl1+z_nl2,method="Wald")
nonlinear_table <- data.frame(outcome="Bilateral combined grip",model="M3 restricted cubic spline",
  n=nrow(z),internal_knot_1=internal[1],internal_knot_2=internal[2],boundary_knot_1=boundary[1],boundary_knot_2=boundary[2],
  overall_wald=as.numeric(overall$Ftest),overall_df_num=as.numeric(overall$df),overall_df_den=as.numeric(overall$ddf),overall_p=as.numeric(overall$p),
  nonlinearity_wald=as.numeric(nonlinear$Ftest),nonlinearity_df_num=as.numeric(nonlinear$df),nonlinearity_df_den=as.numeric(nonlinear$ddf),nonlinearity_p=as.numeric(nonlinear$p))
write.csv(nonlinear_table,file.path(out,"17_continuous_z_nonlinearity_test.csv"),row.names=FALSE)

cat("Continuous conditional-z grip analyses complete.\n")
print(main,row.names=FALSE)
print(interaction_table,row.names=FALSE)
if(length(sex_rows)) print(do.call(rbind,sex_rows),row.names=FALSE)
print(centered,row.names=FALSE)
print(nonlinear_table,row.names=FALSE)
