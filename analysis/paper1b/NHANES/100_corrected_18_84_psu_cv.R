
options(stringsAsFactors=FALSE,survey.lonely.psu="adjust")
suppressPackageStartupMessages({
  library(readr);library(dplyr);library(tidyr);library(purrr)
  library(gamlss);library(gamlss.dist);library(survey)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1 || !args[1] %in% c("Female","Male"))
  stop("Usage: Rscript 100_corrected_18_84_psu_cv.R <Female|Male>")
sex_value <- args[1]

root_dir <- normalizePath(".",winslash="/",mustWork=TRUE)
paper_dir <- file.path(root_dir,"First Paper")
out_dir <- file.path(paper_dir,"age_range_18_84_corrected_20260815")
table_dir <- file.path(out_dir,"tables"); log_dir <- file.path(out_dir,"logs")
cache_dir <- file.path(out_dir,"cache")
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(log_dir,recursive=TRUE,showWarnings=FALSE)

raw_path <- file.path(paper_dir,"Paper1_FINAL_FROZEN_20260814","cache",
                      "development_corrected_weight_narrow.csv.gz")
selected <- read_csv(file.path(table_dir,"05c_final_specification.csv"),
                     show_col_types=FALSE) %>% filter(sex==sex_value)
stopifnot(nrow(selected)==1)
family_name <- selected$family[[1]]
age_df <- selected$age_df[[1]]
sigma_spec <- selected$sigma_spec[[1]]

development <- read_csv(raw_path,show_col_types=FALSE,progress=FALSE) %>%
  transmute(
    SEQN=as.numeric(SEQN),imputation=as.integer(imputation),cycle=as.character(cycle),
    sex=case_when(RIAGENDR==1~"Male",RIAGENDR==2~"Female",TRUE~NA_character_),
    age=as.numeric(RIDAGEYR),pregnant=coalesce(RIAGENDR==2&RIDEXPRG==1,FALSE),
    height_m=as.numeric(BMXHT)/100,alm_kg=as.numeric(alm_g)/1000,
    wtmec_correct=as.numeric(wtmec_correct),
    strata=as.character(SDMVSTRA),psu=as.character(SDMVPSU)
  ) %>%
  filter(sex==sex_value,age>=18,age<=84,!pregnant,
         is.finite(height_m),height_m>0,is.finite(alm_kg),alm_kg>0,
         is.finite(wtmec_correct),wtmec_correct>0,!is.na(strata),!is.na(psu)) %>%
  distinct(imputation,SEQN,.keep_all=TRUE) %>%
  mutate(
    age_band=cut(age,c(18,20,30,40,50,60,70,80,85),right=FALSE,
                 include.lowest=TRUE,
                 labels=c("18-19","20-29","30-39","40-49","50-59","60-69","70-79","80-84")),
    strata_pool=interaction(cycle,strata,drop=TRUE),
    psu_pool=interaction(cycle,strata,psu,drop=TRUE)
  )

sigma_formula <- function(spec)switch(spec,constant=~1,age_df3=~pb(age,df=3),
  height_df3=~pb(height_m,df=3),
  age_height_df3=~pb(age,df=3)+pb(height_m,df=3),stop("unknown sigma"))

fit_once <- function(train,mu_step=1,algorithm="RS") {
  fam <- get(family_name,envir=asNamespace("gamlss.dist"))
  lmfit <- lm(alm_kg~age+height_m,data=train,weights=w_model)
  start <- pmax(as.numeric(predict(lmfit,train)),0.5)
  method_call <- switch(algorithm,RS=quote(RS()),CG=quote(CG()),mixed=quote(mixed()))
  call <- substitute(
    gamlss(as.formula(FORMULA),sigma.formula=sigma_formula(SIG),
      nu.formula=~1,tau.formula=~1,family=fam(),data=train,weights=w_model,
      mu.start=START,method=METHOD,
      control=gamlss.control(n.cyc=500,trace=FALSE,mu.step=STEP)),
    list(FORMULA=paste0("alm_kg~pb(age,df=",age_df,")+pb(height_m,df=3)"),
         SIG=sigma_spec,START=start,METHOD=method_call,STEP=mu_step))
  fit <- eval(call);fit$call$family<-as.call(list(as.name(family_name)));fit$call$data<-quote(train);fit
}
fit_retry <- function(train) {
  specs <- tribble(~step,~alg,1,"RS",.5,"RS",.25,"RS",.1,"RS",.5,"mixed",.25,"mixed",.25,"CG")
  last<-NULL
  for(i in 1:nrow(specs)){
    f<-try(fit_once(train,specs$step[i],specs$alg[i]),silent=TRUE);last<-f
    if(!inherits(f,"try-error")&&isTRUE(f$converged))return(f)
  }; last
}
predpars <- function(fit,train,nd) {
  x<-as.data.frame(nd[,intersect(names(nd),names(train)),drop=FALSE])
  a<-try(suppressWarnings(predictAll(fit,newdata=x,data=as.data.frame(train),type="response")),silent=TRUE)
  if(inherits(a,"try-error"))a<-suppressWarnings(predictAll(fit,newdata=x,type="response"))
  as.data.frame(a)
}
score <- function(fit,train,nd) {
  p<-predpars(fit,train,nd);pf<-get(paste0("p",family_name),asNamespace("gamlss.dist"))
  a<-list(q=nd$alm_kg,mu=p$mu,sigma=p$sigma);if("nu"%in%names(p))a$nu<-p$nu;if("tau"%in%names(p))a$tau<-p$tau
  pr<-pmin(pmax(as.numeric(do.call(pf,a)),1e-10),1-1e-10)
  tibble(percentile=pr,z=qnorm(pr))
}
make_anchors <- function(dat){
  lmfit<-lm(alm_kg~age+height_m,data=dat,weights=wtmec_correct)
  covs<-tibble(age=c(18,84,median(dat$age),median(dat$age)),
               height_m=c(median(dat$height_m),median(dat$height_m),
                          min(dat$height_m),max(dat$height_m)))
  tibble(SEQN=-(1:4),imputation=unique(dat$imputation),cycle="anchor",sex=sex_value,
    age=covs$age,pregnant=FALSE,height_m=covs$height_m,
    alm_kg=pmax(as.numeric(predict(lmfit,covs)),.5),
    wtmec_correct=mean(dat$wtmec_correct)*1e-4,strata="anchor",
    psu=paste0("anchor",1:4),age_band=NA,
    strata_pool=factor("anchor"),psu_pool=factor(paste0("anchor",1:4)))
}

predrows<-list();aud<-list();pc<-ac<-1L
for(imp in 1:5){
  dat<-development%>%filter(imputation==imp)%>%droplevels()
  units<-unique(as.character(dat$psu_pool))
  # Exact original age-range CV seed convention.
  set.seed(20260812+1000*imp+ifelse(sex_value=="Female",17,0))
  fmap<-tibble(unit=units,fold=sample(rep(1:5,length.out=length(units))))
  foldid<-fmap$fold[match(as.character(dat$psu_pool),fmap$unit)]
  anchors<-make_anchors(dat)
  for(fold in 1:5){
    rawtrain<-dat[foldid!=fold,];test<-dat[foldid==fold,]
    train<-bind_rows(rawtrain,anchors)%>%
      mutate(w_model=wtmec_correct/mean(wtmec_correct))%>%droplevels()
    fit<-fit_retry(train)
    if(inherits(fit,"try-error")||!isTRUE(fit$converged))
      stop("CV failed ",sex_value," imp",imp," fold",fold)
    sc<-score(fit,train,test)
    if(any(!is.finite(sc$z)))stop("Invalid CV predictions")
    predrows[[pc]]<-bind_cols(test,sc)%>%mutate(fold=fold);pc<-pc+1L
    aud[[ac]]<-tibble(sex=sex_value,imputation=imp,fold=fold,
      train_n=nrow(rawtrain),test_n=nrow(test),converged=fit$converged,
      iterations=fit$iter);ac<-ac+1L
  }
}
pred<-bind_rows(predrows)%>%
  mutate(across(c(strata,psu,strata_pool,psu_pool),as.character),
         low_p5=as.numeric(percentile<.05),low_p10=as.numeric(percentile<.10))
write_csv(pred,file.path(cache_dir,paste0("corrected_18_84_OOF_",sex_value,".csv.gz")))
write_csv(bind_rows(aud),file.path(cache_dir,paste0("corrected_18_84_CV_audit_",sex_value,".csv")))
cat("DONE ",sex_value,"\n")
