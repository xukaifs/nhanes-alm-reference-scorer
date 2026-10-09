options(stringsAsFactors=FALSE,survey.lonely.psu="adjust")
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(purrr);library(gamlss);library(gamlss.dist)})

args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2||!args[1]%in%c("Female","Male")||!args[2]%in%c("Hcommon","HB"))
  stop("Usage: Rscript 94_corrected_HB_fit_cv.R <Female|Male> <Hcommon|HB>")
sex_value<-args[1];model_name<-args[2];family_name<-if(sex_value=="Female")"BCCG"else"BCT"
set.seed(20260814)
root_dir<-normalizePath(".",winslash="/",mustWork=TRUE)
out_dir<-file.path(root_dir,"First Paper","Paper1_FINAL_FROZEN_20260814")
cache_dir<-file.path(out_dir,"cache");model_dir<-file.path(out_dir,"models");log_dir<-file.path(out_dir,"logs")
walk(c(cache_dir,model_dir,log_dir),dir.create,recursive=TRUE,showWarnings=FALSE)
raw_path<-file.path(cache_dir,"development_corrected_weight_narrow.csv.gz")
tag<-paste(sex_value,model_name,sep="_");log_path<-file.path(log_dir,paste0("corrected_",tag,".log"))
log_message<-function(...){line<-paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S")," | ",paste0(...,collapse=""));cat(line,"\n");cat(line,"\n",file=log_path,append=TRUE)}

development<-read_csv(raw_path,show_col_types=FALSE,progress=FALSE)%>%transmute(
  SEQN=as.numeric(SEQN),imputation=as.integer(imputation),cycle=as.character(cycle),examined=RIDSTATR==2,
  sex=case_when(RIAGENDR==1~"Male",RIAGENDR==2~"Female",TRUE~NA_character_),
  age=as.numeric(RIDAGEYR),pregnant=coalesce(RIAGENDR==2&RIDEXPRG==1,FALSE),
  height_cm=as.numeric(BMXHT),height_m=as.numeric(BMXHT)/100,bmi=as.numeric(BMXBMI),bmi10=as.numeric(BMXBMI)/10,
  alm_kg=as.numeric(alm_g)/1000,wtmec_legacy=as.numeric(wtmec_pool),wtmec_correct=as.numeric(wtmec_correct),
  strata=as.character(SDMVSTRA),psu=as.character(SDMVPSU))%>%
  filter(examined,sex==sex_value,age>=18,age<=69,!pregnant,is.finite(height_m),height_m>0,
         is.finite(bmi),bmi>0,is.finite(alm_kg),alm_kg>0,is.finite(wtmec_correct),wtmec_correct>0,
         !is.na(strata),!is.na(psu))%>%distinct(imputation,SEQN,.keep_all=TRUE)%>%
  mutate(strata_pool=interaction(cycle,strata,drop=TRUE),psu_pool=interaction(cycle,strata,psu,drop=TRUE))

height_variable<-if(sex_value=="Female")"height_cm"else"height_m"
age_term<-if(sex_value=="Female")"pb(age,df=3,inter=10)"else"pb(age,df=3)"
height_term<-if(sex_value=="Female")paste0("pb(",height_variable,",df=3,inter=10)")else paste0("pb(",height_variable,",df=3)")
mu_formula<-as.formula(paste0("alm_kg~",age_term,"+",height_term,if(model_name=="HB")"+pb(bmi10,df=3,inter=10)"else""))
linear_formula<-as.formula(paste0("alm_kg~age+",height_variable,if(model_name=="HB")"+bmi10"else""))

fit_model<-function(train_data,mu_step=1,algorithm="RS",start_mode="lm"){
  family_function<-get(family_name,envir=asNamespace("gamlss.dist"));lm_start<-lm(linear_formula,data=train_data,weights=w_model)
  selected_start<-if(start_mode=="lm")pmax(as.numeric(predict(lm_start,train_data)),.5)else NULL
  method_call<-switch(algorithm,RS=quote(RS()),CG=quote(CG()),mixed=quote(mixed()))
  call<-substitute(gamlss(MU,sigma.formula=~1,nu.formula=~1,tau.formula=~1,family=family_function(),data=train_data,
    weights=w_model,mu.start=START,method=METHOD,control=gamlss.control(n.cyc=500,trace=FALSE,mu.step=STEP)),
    list(MU=mu_formula,START=selected_start,METHOD=method_call,STEP=mu_step))
  fit<-eval(call);fit$call$family<-as.call(list(as.name(family_name)));fit$call$data<-quote(train_data);fit
}
fit_with_retries<-function(train_data){
  attempts<-tribble(~mu_step,~algorithm,~start_mode,1,"RS","lm",.5,"RS","lm",.25,"RS","lm",.1,"RS","lm",
    .5,"mixed","lm",.25,"mixed","lm",.25,"CG","lm",.5,"RS","default")
  last<-NULL
  for(i in 1:nrow(attempts)){fit<-try(fit_model(train_data,attempts$mu_step[i],attempts$algorithm[i],attempts$start_mode[i]),silent=TRUE);last<-fit
    if(!inherits(fit,"try-error")&&isTRUE(fit$converged))return(list(fit=fit,attempt=i,algorithm=attempts$algorithm[i],mu_step=attempts$mu_step[i]))}
  list(fit=last,attempt=nrow(attempts),algorithm=NA_character_,mu_step=NA_real_)
}
predict_parameters<-function(model,train_data,new_data){nd<-as.data.frame(new_data[,intersect(names(new_data),names(train_data)),drop=FALSE]);
  ans<-try(suppressWarnings(predictAll(model,newdata=nd,data=as.data.frame(train_data),type="response")),silent=TRUE)
  if(inherits(ans,"try-error"))ans<-suppressWarnings(predictAll(model,newdata=nd,type="response"));as.data.frame(ans)}
score_model<-function(model,train_data,new_data){pars<-predict_parameters(model,train_data,new_data);pfun<-get(paste0("p",family_name),asNamespace("gamlss.dist"));
  a<-list(q=new_data$alm_kg,mu=pars$mu,sigma=pars$sigma);if("nu"%in%names(pars))a$nu<-pars$nu;if("tau"%in%names(pars))a$tau<-pars$tau
  p<-pmin(pmax(as.numeric(do.call(pfun,a)),1e-10),1-1e-10);tibble(percentile=p,z=qnorm(p))}

bundle_path<-file.path(model_dir,paste0("corrected_",tag,"_bundles.rds"));fit_path<-file.path(cache_dir,paste0("corrected_",tag,"_full_fit.csv"))
bundles<-list();full_rows<-list()
for(imp in 1:5){dat<-development%>%filter(imputation==imp)%>%mutate(w_model=wtmec_correct/mean(wtmec_correct))%>%droplevels();
  log_message("full imp",imp);start<-Sys.time();res<-fit_with_retries(dat);fit<-res$fit
  if(inherits(fit,"try-error")||!isTRUE(fit$converged))stop(tag," full fit failed imp",imp," ",substr(as.character(fit),1,400))
  bundles[[paste0("imp",imp)]]<-list(sex=sex_value,model=model_name,family=family_name,train_data=dat,fit=fit)
  full_rows[[imp]]<-tibble(sex=sex_value,model=model_name,family=family_name,imputation=imp,n=nrow(dat),converged=fit$converged,
    iterations=fit$iter,effective_df=fit$df.fit,global_deviance=fit$G.deviance,AIC=AIC(fit),BIC=GAIC(fit,k=log(nrow(dat))),
    attempt=res$attempt,elapsed_seconds=as.numeric(difftime(Sys.time(),start,units="secs")))}
saveRDS(bundles,bundle_path);write_csv(bind_rows(full_rows),fit_path)

make_anchors<-function(full_dat,raw_train){covs<-tibble(
  age=c(18,69,median(full_dat$age),median(full_dat$age),median(full_dat$age),median(full_dat$age)),
  height_cm=c(median(full_dat$height_cm),median(full_dat$height_cm),min(full_dat$height_cm),max(full_dat$height_cm),median(full_dat$height_cm),median(full_dat$height_cm)),
  bmi=c(median(full_dat$bmi),median(full_dat$bmi),median(full_dat$bmi),median(full_dat$bmi),min(full_dat$bmi),max(full_dat$bmi)))%>%
  mutate(height_m=height_cm/100,bmi10=bmi/10);lm_fit<-lm(linear_formula,data=raw_train,weights=wtmec_correct)
  tibble(SEQN=-(1:nrow(covs)),imputation=unique(full_dat$imputation),cycle="anchor",examined=TRUE,sex=sex_value,age=covs$age,pregnant=FALSE,
    height_cm=covs$height_cm,height_m=covs$height_m,bmi=covs$bmi,bmi10=covs$bmi10,alm_kg=pmax(as.numeric(predict(lm_fit,covs)),.5),
    wtmec_legacy=0,wtmec_correct=mean(raw_train$wtmec_correct)*1e-4,strata="anchor",psu=paste0("anchor",1:nrow(covs)),
    strata_pool=factor("anchor"),psu_pool=factor(paste0("anchor",1:nrow(covs))))}

pred_path<-file.path(cache_dir,paste0("corrected_",tag,"_OOF.csv.gz"));audit_path<-file.path(cache_dir,paste0("corrected_",tag,"_CV_audit.csv"))
pred_rows<-list();audit_rows<-list();pc<-ac<-1L
for(imp in 1:5){dat<-development%>%filter(imputation==imp)%>%droplevels();units<-unique(as.character(dat$psu_pool));
  set.seed(20260812+1000*imp+ifelse(sex_value=="Female",17,0));fold_map<-tibble(unit=units,fold=sample(rep(1:5,length.out=length(units))));
  fold_assignment<-fold_map$fold[match(as.character(dat$psu_pool),fold_map$unit)]
  for(fold in 1:5){raw_train<-dat[fold_assignment!=fold,];test<-dat[fold_assignment==fold,];anchors<-make_anchors(dat,raw_train)
    train<-bind_rows(raw_train,anchors)%>%mutate(w_model=wtmec_correct/mean(wtmec_correct))%>%droplevels();log_message("CV imp",imp," fold",fold)
    start<-Sys.time();res<-fit_with_retries(train);fit<-res$fit;if(inherits(fit,"try-error")||!isTRUE(fit$converged))stop(tag," CV failed")
    sc<-score_model(fit,train,test);if(any(!is.finite(sc$z)))stop(tag," invalid CV scores")
    pred_rows[[pc]]<-bind_cols(test,sc)%>%mutate(fold=fold);pc<-pc+1L
    audit_rows[[ac]]<-tibble(sex=sex_value,model=model_name,imputation=imp,fold=fold,train_n=nrow(raw_train),test_n=nrow(test),
      converged=fit$converged,iterations=fit$iter,attempt=res$attempt,elapsed_seconds=as.numeric(difftime(Sys.time(),start,units="secs")),
      valid_predictions=sum(is.finite(sc$z)),expected_predictions=nrow(test));ac<-ac+1L;write_csv(bind_rows(audit_rows),audit_path);rm(fit);gc(FALSE)}}
pred<-bind_rows(pred_rows)%>%mutate(across(c(strata,psu,strata_pool,psu_pool),as.character),low_p5=as.numeric(percentile<.05),low_p10=as.numeric(percentile<.10))
write_csv(pred,pred_path);write_csv(bind_rows(audit_rows),audit_path);log_message("complete")
