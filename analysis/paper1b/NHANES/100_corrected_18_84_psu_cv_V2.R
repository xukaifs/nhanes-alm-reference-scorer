options(stringsAsFactors=FALSE, survey.lonely.psu="adjust")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gamlss); library(gamlss.dist); library(survey)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1 || !args[1] %in% c("Female","Male"))
  stop("Usage: Rscript 100_corrected_18_84_psu_cv_V2.R <Female|Male>")
sex_value <- args[1]

find_root <- function(start) {
  d <- normalizePath(start,winslash="/",mustWork=TRUE)
  for(i in 0:10) {
    if(file.exists(file.path(d,"First Paper","raw_for_refit","merged_training_1999_2006.csv.gz")))
      return(d)
    p <- dirname(d); if(identical(p,d)) break; d <- p
  }
  stop("Could not find project root")
}
root_dir <- find_root(getwd())
paper_dir <- file.path(root_dir,"First Paper")
out_dir <- file.path(paper_dir,"age_range_18_84_corrected_20260815")
table_dir <- file.path(out_dir,"tables")
cache_dir <- file.path(out_dir,"cache")
log_dir <- file.path(out_dir,"logs")
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(log_dir,recursive=TRUE,showWarnings=FALSE)

log_path <- file.path(log_dir,paste0("corrected_18_84_CV_",sex_value,"_V2.log"))
logmsg <- function(...) {
  x <- paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S")," | ",paste0(...,collapse=""))
  cat(x,"\n"); cat(x,"\n",file=log_path,append=TRUE)
}

raw_path <- file.path(paper_dir,"Paper1_FINAL_FROZEN_20260814","cache",
                      "development_corrected_weight_narrow.csv.gz")
spec_path <- file.path(table_dir,"05c_final_specification.csv")
stopifnot(file.exists(raw_path),file.exists(spec_path))

selected <- read_csv(spec_path,show_col_types=FALSE) %>% filter(sex==sex_value)
stopifnot(nrow(selected)==1)
family_name <- selected$family[[1]]
age_df <- selected$age_df[[1]]
sigma_spec <- selected$sigma_spec[[1]]

logmsg("Selected spec: sex=",sex_value,", family=",family_name,
       ", age_df=",age_df,", height_df=3, sigma=",sigma_spec)

development <- read_csv(raw_path,show_col_types=FALSE,progress=FALSE) %>%
  transmute(
    SEQN=as.numeric(SEQN), imputation=as.integer(imputation), cycle=as.character(cycle),
    sex=case_when(RIAGENDR==1~"Male",RIAGENDR==2~"Female",TRUE~NA_character_),
    age=as.numeric(RIDAGEYR), pregnant=coalesce(RIAGENDR==2&RIDEXPRG==1,FALSE),
    height_m=as.numeric(BMXHT)/100, alm_kg=as.numeric(alm_g)/1000,
    wtmec_correct=as.numeric(wtmec_correct),
    strata=as.character(SDMVSTRA), psu=as.character(SDMVPSU)
  ) %>%
  filter(sex==sex_value,age>=18,age<=84,!pregnant,
         is.finite(height_m),height_m>0,is.finite(alm_kg),alm_kg>0,
         is.finite(wtmec_correct),wtmec_correct>0,!is.na(strata),!is.na(psu)) %>%
  distinct(imputation,SEQN,.keep_all=TRUE) %>%
  mutate(
    age_band=cut(age,c(18,20,30,40,50,60,70,80,85),right=FALSE,include.lowest=TRUE,
      labels=c("18-19","20-29","30-39","40-49","50-59","60-69","70-79","80-84")),
    strata_pool=interaction(cycle,strata,drop=TRUE),
    psu_pool=interaction(cycle,strata,psu,drop=TRUE)
  )

sigma_formula_from_spec <- function(spec) {
  switch(spec,
    constant=~1,
    age_df3=~pb(age,df=3),
    height_df3=~pb(height_m,df=3),
    age_height_df3=~pb(age,df=3)+pb(height_m,df=3),
    stop("Unknown sigma spec: ",spec))
}

fit_one <- function(train,mu_step,algorithm,start_mode,ncyc) {
  fam <- get(family_name,envir=asNamespace("gamlss.dist"))
  mu_start <- NULL
  if(start_mode=="lm") {
    lmfit <- lm(alm_kg~age+height_m,data=train,weights=w_model)
    mu_start <- pmax(as.numeric(predict(lmfit,train)),0.5)
  }
  method_call <- switch(algorithm,RS=quote(RS()),CG=quote(CG()),mixed=quote(mixed()))
  fcall <- substitute(
    gamlss(as.formula(FORM),
      sigma.formula=sigma_formula_from_spec(SIG),
      nu.formula=~1,tau.formula=~1,
      family=fam(),data=train,weights=w_model,mu.start=START,
      method=METHOD,
      control=gamlss.control(n.cyc=NCYC,trace=FALSE,mu.step=STEP)),
    list(FORM=paste0("alm_kg ~ pb(age, df=",age_df,") + pb(height_m, df=3)"),
         SIG=sigma_spec,START=mu_start,METHOD=method_call,NCYC=ncyc,STEP=mu_step))
  fit <- eval(fcall)
  fit$call$family <- as.call(list(as.name(family_name)))
  fit$call$data <- quote(train)
  fit
}

attempt_grid <- tribble(
  ~mu_step,~algorithm,~start_mode,~ncyc,
  1.00,"RS","lm",500,
  0.50,"RS","lm",500,
  0.25,"RS","lm",500,
  0.10,"RS","lm",500,
  0.50,"mixed","lm",500,
  0.25,"mixed","lm",500,
  0.25,"CG","lm",500,
  0.50,"RS","default",500,
  0.25,"RS","lm",1000,
  0.10,"RS","lm",1000,
  0.05,"RS","lm",1000,
  0.25,"mixed","lm",1000,
  0.10,"mixed","lm",1000,
  0.10,"CG","lm",1000
)

fit_retry <- function(train,imp,fold) {
  attempts_out <- list()
  last <- NULL
  for(i in seq_len(nrow(attempt_grid))) {
    a <- attempt_grid[i,]
    start <- Sys.time()
    fit <- try(fit_one(train,a$mu_step,a$algorithm,a$start_mode,a$ncyc),silent=TRUE)
    elapsed <- as.numeric(difftime(Sys.time(),start,units="secs"))
    last <- fit
    if(inherits(fit,"try-error")) {
      err <- substr(as.character(fit),1,600)
      conv <- FALSE; iter <- NA_real_; dev <- NA_real_
    } else {
      err <- NA_character_; conv <- isTRUE(fit$converged)
      iter <- fit$iter; dev <- fit$G.deviance
    }
    logmsg("imp",imp," fold",fold," attempt",i,
           " converged=",conv," algorithm=",a$algorithm,
           " step=",a$mu_step," start=",a$start_mode," ncyc=",a$ncyc,
           if(!is.na(err)) paste0(" error=",err) else "")
    attempts_out[[i]] <- tibble(
      sex=sex_value,imputation=imp,fold=fold,attempt=i,
      mu_step=a$mu_step,algorithm=a$algorithm,start_mode=a$start_mode,ncyc=a$ncyc,
      converged=conv,iterations=iter,global_deviance=dev,
      elapsed_seconds=elapsed,error=err
    )
    if(!inherits(fit,"try-error") && isTRUE(fit$converged))
      return(list(fit=fit,audit=bind_rows(attempts_out),attempt=i))
  }
  list(fit=last,audit=bind_rows(attempts_out),attempt=NA_integer_)
}

predict_parameters <- function(model,train,newdata) {
  x <- as.data.frame(newdata[,intersect(names(newdata),names(train)),drop=FALSE])
  p <- try(suppressWarnings(predictAll(model,newdata=x,data=as.data.frame(train),type="response")),silent=TRUE)
  if(inherits(p,"try-error"))
    p <- suppressWarnings(predictAll(model,newdata=x,type="response"))
  as.data.frame(p)
}
score_model <- function(model,train,newdata) {
  p <- predict_parameters(model,train,newdata)
  pf <- get(paste0("p",family_name),envir=asNamespace("gamlss.dist"))
  aa <- list(q=newdata$alm_kg,mu=p$mu,sigma=p$sigma)
  if("nu"%in%names(p))aa$nu<-p$nu
  if("tau"%in%names(p))aa$tau<-p$tau
  pr <- as.numeric(do.call(pf,aa))
  pr <- pmin(pmax(pr,1e-10),1-1e-10)
  tibble(percentile=pr,z=qnorm(pr))
}

make_anchors <- function(dat) {
  lmfit <- lm(alm_kg~age+height_m,data=dat,weights=wtmec_correct)
  covs <- tibble(
    age=c(18,84,median(dat$age),median(dat$age)),
    height_m=c(median(dat$height_m),median(dat$height_m),
               min(dat$height_m),max(dat$height_m))
  )
  tibble(
    SEQN=-(1:4),imputation=unique(dat$imputation),cycle="anchor",sex=sex_value,
    age=covs$age,pregnant=FALSE,height_m=covs$height_m,
    alm_kg=pmax(as.numeric(predict(lmfit,covs)),0.5),
    wtmec_correct=mean(dat$wtmec_correct)*1e-4,
    strata="anchor",psu=paste0("anchor",1:4),
    age_band=cut(covs$age,c(18,20,30,40,50,60,70,80,85),
                 right=FALSE,include.lowest=TRUE,
                 labels=c("18-19","20-29","30-39","40-49","50-59","60-69","70-79","80-84")),
    strata_pool=factor("anchor"),
    psu_pool=factor(paste0("anchor",1:4))
  )
}

fold_files <- c()
audit_all <- list()
for(imp in 1:5) {
  dat <- development %>% filter(imputation==imp) %>% droplevels()
  units <- unique(as.character(dat$psu_pool))
  set.seed(20260812+1000*imp+ifelse(sex_value=="Female",17,0))
  fold_map <- tibble(unit=units,fold=sample(rep(1:5,length.out=length(units))))
  fold_assignment <- fold_map$fold[match(as.character(dat$psu_pool),fold_map$unit)]
  anchors <- make_anchors(dat)

  for(fold in 1:5) {
    fold_file <- file.path(cache_dir,
      paste0("checkpoint_corrected_18_84_",sex_value,"_imp",imp,"_fold",fold,".csv.gz"))
    fold_audit <- file.path(cache_dir,
      paste0("checkpoint_corrected_18_84_",sex_value,"_imp",imp,"_fold",fold,"_attempts.csv"))

    if(file.exists(fold_file)) {
      logmsg("Using checkpoint imp",imp," fold",fold)
      fold_files <- c(fold_files,fold_file)
      if(file.exists(fold_audit)) audit_all[[length(audit_all)+1]] <- read_csv(fold_audit,show_col_types=FALSE)
      next
    }

    rawtrain <- dat[fold_assignment!=fold,]
    test <- dat[fold_assignment==fold,]
    train <- bind_rows(rawtrain,anchors) %>%
      mutate(w_model=wtmec_correct/mean(wtmec_correct)) %>% droplevels()

    logmsg("START imp",imp," fold",fold," train_n=",nrow(rawtrain)," test_n=",nrow(test))
    rr <- fit_retry(train,imp,fold)
    write_csv(rr$audit,fold_audit)
    audit_all[[length(audit_all)+1]] <- rr$audit

    fit <- rr$fit
    if(inherits(fit,"try-error") || !isTRUE(fit$converged)) {
      write_csv(bind_rows(audit_all),
        file.path(cache_dir,paste0("corrected_18_84_CV_attempts_",sex_value,".csv")))
      stop("All convergence attempts failed for ",sex_value," imp",imp," fold",fold,
           ". Inspect: ",fold_audit)
    }

    sc <- try(score_model(fit,train,test),silent=TRUE)
    if(inherits(sc,"try-error") || any(!is.finite(sc$z))) {
      stop("Prediction failed after converged fit: ",sex_value," imp",imp," fold",fold)
    }

    out <- bind_cols(test,sc) %>% mutate(fold=fold)
    write_csv(out,fold_file)
    fold_files <- c(fold_files,fold_file)
    rm(fit,rr,sc,out); gc(FALSE)
  }
}

pred <- map_dfr(fold_files,~read_csv(.x,show_col_types=FALSE,progress=FALSE)) %>%
  distinct(imputation,SEQN,.keep_all=TRUE) %>%
  mutate(across(c(strata,psu,strata_pool,psu_pool),as.character),
         low_p5=as.numeric(percentile<.05),low_p10=as.numeric(percentile<.10))

final_pred <- file.path(cache_dir,paste0("corrected_18_84_OOF_",sex_value,".csv.gz"))
final_audit <- file.path(cache_dir,paste0("corrected_18_84_CV_attempts_",sex_value,".csv"))
write_csv(pred,final_pred)
write_csv(bind_rows(audit_all),final_audit)

expected <- nrow(development)
if(nrow(pred)!=expected)
  stop("OOF row count mismatch: expected ",expected,", got ",nrow(pred))

logmsg("DONE. OOF rows=",nrow(pred))
cat("\nDONE ",sex_value,"\n",final_pred,"\n")
