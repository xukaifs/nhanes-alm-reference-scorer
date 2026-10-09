
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gamlss); library(gamlss.dist); library(survey)
})

set.seed(20260897)

root_dir <- normalizePath(".", winslash="/", mustWork=TRUE)
paper_dir <- file.path(root_dir, "First Paper")
out_dir <- file.path(paper_dir, "age_range_18_84_corrected_20260815")
table_dir <- file.path(out_dir, "tables")
model_dir <- file.path(out_dir, "models")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(model_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

# Prefer the final corrected-weight cache created during the frozen rerun.
raw_path <- file.path(paper_dir, "Paper1_FINAL_FROZEN_20260814", "cache",
                      "development_corrected_weight_narrow.csv.gz")
if (!file.exists(raw_path)) {
  stop("Missing corrected cache: ", raw_path,
       "\nRun 91a_make_corrected_weight_cache.py first.", call.=FALSE)
}

log_path <- file.path(log_dir, "corrected_18_84_model_fit.log")
log_message <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                 paste0(..., collapse=""))
  cat(line, "\n"); cat(line, "\n", file=log_path, append=TRUE)
}

development <- read_csv(raw_path, show_col_types=FALSE, progress=FALSE) %>%
  transmute(
    SEQN=as.numeric(SEQN), imputation=as.integer(imputation),
    cycle=as.character(cycle),
    sex=case_when(RIAGENDR==1~"Male", RIAGENDR==2~"Female", TRUE~NA_character_),
    age=as.numeric(RIDAGEYR),
    pregnant=coalesce(RIAGENDR==2 & RIDEXPRG==1, FALSE),
    height_m=as.numeric(BMXHT)/100,
    alm_kg=as.numeric(alm_g)/1000,
    wtmec_correct=as.numeric(wtmec_correct),
    strata=as.character(SDMVSTRA), psu=as.character(SDMVPSU)
  ) %>%
  filter(
    age>=18, age<=84, sex %in% c("Female","Male"), !pregnant,
    is.finite(height_m), height_m>0, is.finite(alm_kg), alm_kg>0,
    is.finite(wtmec_correct), wtmec_correct>0, !is.na(strata), !is.na(psu)
  ) %>%
  distinct(imputation, SEQN, .keep_all=TRUE) %>%
  mutate(
    age_band=cut(age, c(18,20,30,40,50,60,70,80,85), right=FALSE,
                 include.lowest=TRUE,
                 labels=c("18-19","20-29","30-39","40-49",
                          "50-59","60-69","70-79","80-84")),
    strata_pool=interaction(cycle,strata,drop=TRUE),
    psu_pool=interaction(cycle,strata,psu,drop=TRUE)
  )

get_training_set <- function(imp, sex_value) {
  development %>% filter(imputation==imp, sex==sex_value) %>%
    mutate(w_model=wtmec_correct/mean(wtmec_correct)) %>% droplevels()
}

sample_summary <- development %>% filter(imputation==1) %>%
  group_by(sex) %>% summarise(
    age_range="18-84", unweighted_n=n(),
    age_min=min(age), age_max=max(age), unique_cycles=n_distinct(cycle),
    unique_strata=n_distinct(strata_pool), unique_psu=n_distinct(psu_pool),
    .groups="drop")
write_csv(sample_summary, file.path(table_dir,"01_sample_summary.csv"))

sample_by_band <- development %>% filter(imputation==1) %>%
  count(sex,age_band,cycle,name="unweighted_n") %>%
  complete(sex,age_band,cycle,fill=list(unweighted_n=0L))
write_csv(sample_by_band,file.path(table_dir,"02_n_by_age_sex_cycle.csv"))

sigma_formula_from_spec <- function(spec) {
  switch(spec,
    constant=~1,
    age_df3=~pb(age,df=3),
    height_df3=~pb(height_m,df=3),
    age_height_df3=~pb(age,df=3)+pb(height_m,df=3),
    stop("Unknown sigma specification: ",spec))
}

fit_candidate <- function(train_data,family_name,age_df,sigma_spec="constant",
                          n_cycles=300, mu_step=1, algorithm="RS",
                          start_mode="default") {
  fam <- get(family_name,envir=asNamespace("gamlss.dist"))
  method_call <- switch(algorithm,RS=quote(RS()),CG=quote(CG()),mixed=quote(mixed()))
  mu_start <- NULL
  if (start_mode=="lm") {
    lm_fit <- lm(alm_kg~age+height_m,data=train_data,weights=w_model)
    mu_start <- pmax(as.numeric(predict(lm_fit,train_data)),0.5)
  }
  fit_call <- substitute(
    gamlss(as.formula(FORMULA),
      sigma.formula=sigma_formula_from_spec(SIGMA),
      nu.formula=~1,tau.formula=~1,
      family=fam(),data=train_data,weights=w_model,
      mu.start=START,method=METHOD,
      control=gamlss.control(n.cyc=NCYC,trace=FALSE,mu.step=STEP)),
    list(FORMULA=paste0("alm_kg ~ pb(age, df = ",age_df,
                        ") + pb(height_m, df = 3)"),
         SIGMA=sigma_spec,START=mu_start,METHOD=method_call,
         NCYC=n_cycles,STEP=mu_step))
  fit <- eval(fit_call)
  fit$call$family <- as.call(list(as.name(family_name)))
  fit$call$data <- quote(train_data)
  fit
}

fit_with_retries <- function(train_data,family_name,age_df,sigma_spec,n_cycles=300) {
  attempts <- tribble(
    ~mu_step,~algorithm,~start_mode,
    1.00,"RS","default",
    1.00,"RS","lm",
    0.50,"RS","lm",
    0.25,"RS","lm",
    0.10,"RS","lm",
    0.50,"mixed","lm",
    0.25,"mixed","lm",
    0.25,"CG","lm"
  )
  last <- NULL
  for(i in seq_len(nrow(attempts))) {
    fit <- try(fit_candidate(train_data,family_name,age_df,sigma_spec,n_cycles,
                             attempts$mu_step[i],attempts$algorithm[i],
                             attempts$start_mode[i]),silent=TRUE)
    last <- fit
    if(!inherits(fit,"try-error") && isTRUE(fit$converged))
      return(list(fit=fit,attempt=i))
  }
  list(fit=last,attempt=nrow(attempts))
}

diag_row <- function(fit,imp,sex_value,family_name,age_df,sigma_spec,n,attempt=NA_integer_) {
  if(inherits(fit,"try-error") || !isTRUE(fit$converged)) {
    return(tibble(imputation=imp,sex=sex_value,family=family_name,
      age_df=age_df,height_df=3,sigma_spec=sigma_spec,analysis_n=n,
      converged=FALSE,iterations=NA_real_,effective_df=NA_real_,
      global_deviance=NA_real_,AIC=NA_real_,BIC=Inf,attempt=attempt,
      error=substr(as.character(fit),1,500)))
  }
  tibble(imputation=imp,sex=sex_value,family=family_name,
    age_df=age_df,height_df=3,sigma_spec=sigma_spec,analysis_n=n,
    converged=TRUE,iterations=fit$iter,effective_df=fit$df.fit,
    global_deviance=fit$G.deviance,AIC=AIC(fit),
    BIC=GAIC(fit,k=log(n)),attempt=attempt,error=NA_character_)
}

# Stage A: EXACT original candidate grid.
rows <- list(); k <- 1L
for(sex_value in c("Female","Male")) for(imp in 1:5) {
  train <- get_training_set(imp,sex_value)
  for(family_name in c("NO","BCCG","BCPE","BCT")) for(age_df in c(3,4,5)) {
    log_message("Stage A ",sex_value," imp",imp," ",family_name," age_df=",age_df)
    rr <- fit_with_retries(train,family_name,age_df,"constant",300)
    rows[[k]] <- diag_row(rr$fit,imp,sex_value,family_name,age_df,
                          "constant",nrow(train),rr$attempt)
    k <- k+1L; rm(rr); gc(FALSE)
  }
}
stage_a <- bind_rows(rows)
write_csv(stage_a,file.path(table_dir,"04a_family_age_df_model_comparison.csv"))
stage_a_summary <- stage_a %>%
  group_by(sex,family,age_df,height_df,sigma_spec) %>%
  summarise(total_BIC=if(all(converged & is.finite(BIC)))sum(BIC) else Inf,
            mean_BIC=if(all(converged & is.finite(BIC)))mean(BIC) else Inf,
            converged_imputations=sum(converged & is.finite(BIC)),.groups="drop") %>%
  group_by(sex) %>% arrange(total_BIC,.by_group=TRUE) %>%
  mutate(BIC_rank=row_number(),delta_BIC=total_BIC-min(total_BIC)) %>% ungroup()
write_csv(stage_a_summary,file.path(table_dir,"04b_family_age_df_BIC_summary.csv"))
stage_a_selected <- stage_a_summary %>% filter(BIC_rank==1)
write_csv(stage_a_selected,file.path(table_dir,"04c_family_age_df_selected.csv"))
stopifnot(nrow(stage_a_selected)==2,all(is.finite(stage_a_selected$total_BIC)))

# Stage B: EXACT original sigma grid.
rows <- list(); k <- 1L
for(sex_value in c("Female","Male")) {
  sel <- stage_a_selected %>% filter(sex==sex_value)
  for(imp in 1:5) {
    train <- get_training_set(imp,sex_value)
    for(sigma_spec in c("constant","age_df3","height_df3","age_height_df3")) {
      log_message("Stage B ",sex_value," imp",imp," sigma=",sigma_spec)
      rr <- fit_with_retries(train,sel$family[[1]],sel$age_df[[1]],sigma_spec,350)
      rows[[k]] <- diag_row(rr$fit,imp,sex_value,sel$family[[1]],
                            sel$age_df[[1]],sigma_spec,nrow(train),rr$attempt)
      k <- k+1L; rm(rr); gc(FALSE)
    }
  }
}
stage_b <- bind_rows(rows)
write_csv(stage_b,file.path(table_dir,"05a_sigma_model_comparison.csv"))
stage_b_summary <- stage_b %>%
  group_by(sex,family,age_df,height_df,sigma_spec) %>%
  summarise(total_BIC=if(all(converged & is.finite(BIC)))sum(BIC) else Inf,
            mean_BIC=if(all(converged & is.finite(BIC)))mean(BIC) else Inf,
            converged_imputations=sum(converged & is.finite(BIC)),.groups="drop") %>%
  group_by(sex) %>% arrange(total_BIC,.by_group=TRUE) %>%
  mutate(BIC_rank=row_number(),delta_BIC=total_BIC-min(total_BIC)) %>% ungroup()
write_csv(stage_b_summary,file.path(table_dir,"05b_sigma_BIC_summary.csv"))
selected <- stage_b_summary %>% filter(BIC_rank==1)
write_csv(selected,file.path(table_dir,"05c_final_specification.csv"))
stopifnot(nrow(selected)==2,all(is.finite(selected$total_BIC)))

# Final full-data models
bundles <- list(); rows <- list(); k <- 1L
for(sex_value in c("Female","Male")) {
  sel <- selected %>% filter(sex==sex_value)
  for(imp in 1:5) {
    train <- get_training_set(imp,sex_value)
    log_message("Final fit ",sex_value," imp",imp)
    rr <- fit_with_retries(train,sel$family[[1]],sel$age_df[[1]],
                           sel$sigma_spec[[1]],500)
    fit <- rr$fit
    if(inherits(fit,"try-error") || !isTRUE(fit$converged))
      stop("Final fit failed: ",sex_value," imp",imp)
    bundles[[paste0("imp",imp,"_",sex_value)]] <- list(
      imputation=imp,sex=sex_value,family=sel$family[[1]],
      age_df=sel$age_df[[1]],height_df=3,sigma_spec=sel$sigma_spec[[1]],
      age_domain=c(18,84),train_data=train,model=fit)
    rows[[k]] <- diag_row(fit,imp,sex_value,sel$family[[1]],sel$age_df[[1]],
                          sel$sigma_spec[[1]],nrow(train),rr$attempt)
    k <- k+1L
  }
}
saveRDS(bundles,file.path(model_dir,"corrected_age18_84_model_bundles.rds"))
write_csv(bind_rows(rows),file.path(table_dir,"06_final_model_fit.csv"))
log_message("DONE corrected 18-84 model selection and final fitting")
print(selected)
