
options(stringsAsFactors=FALSE,survey.lonely.psu="adjust")
suppressPackageStartupMessages({
  library(readr);library(dplyr);library(tidyr);library(purrr)
  library(gamlss);library(gamlss.dist);library(survey);library(splines)
})

root_dir<-normalizePath(".",winslash="/",mustWork=TRUE)
paper_dir<-file.path(root_dir,"First Paper")
out_dir<-file.path(paper_dir,"age_range_18_84_corrected_20260815")
table_dir<-file.path(out_dir,"tables");cache_dir<-file.path(out_dir,"cache")
model_dir<-file.path(out_dir,"models")
dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)

make_design<-function(d)svydesign(ids=~psu_pool,strata=~strata_pool,
  weights=~wtmec_pool,nest=TRUE,data=d)

rubin<-function(q,se,null=0){
  ok<-is.finite(q)&is.finite(se)&se>=0;q<-q[ok];se<-se[ok];m<-length(q)
  if(!m)return(tibble(m=0,estimate=NA_real_,standard_error=NA_real_,df=NA_real_,
    conf_low=NA_real_,conf_high=NA_real_,p_value=NA_real_))
  qb<-mean(q);u<-mean(se^2);b<-if(m>1)var(q)else 0;t<-u+(1+1/m)*b;s<-sqrt(t)
  df<-if(m>1&&b>0)(m-1)*(1+u/((1+1/m)*b))^2 else Inf;cr<-qt(.975,df)
  tibble(m=m,estimate=qb,standard_error=s,df=df,
         conf_low=qb-cr*s,conf_high=qb+cr*s,
         p_value=2*pt(abs((qb-null)/s),df=df,lower.tail=FALSE))
}
pool_table<-function(d,g)d%>%group_by(across(all_of(g)))%>%
  group_modify(~rubin(.x$estimate,.x$standard_error))%>%ungroup()
safe_sd<-function(des,var){x<-svyvar(as.formula(paste0("~",var)),des,na.rm=TRUE);
  sqrt(pmax(as.numeric(coef(x)[1]),0))}

# ---------------- OOF calibration ----------------
oof<-map_dfr(c("Female","Male"),function(s)
  read_csv(file.path(cache_dir,paste0("corrected_18_84_OOF_",s,".csv.gz")),
           show_col_types=FALSE,progress=FALSE)%>%mutate(sex=s))%>%
  mutate(wtmec_pool=as.numeric(wtmec_correct),
         strata_pool=factor(strata_pool),psu_pool=factor(psu_pool),
         low_p5=as.numeric(percentile<.05),low_p10=as.numeric(percentile<.10),
         audit_band=cut(age,c(18,20,40,60,70,80,85),right=FALSE,include.lowest=TRUE,
           labels=c("18-19","20-39","40-59","60-69","70-79","80-84")))

single<-list();k<-1L
for(imp in 1:5)for(s in c("Female","Male")){
  dd<-oof%>%filter(imputation==imp,sex==s)
  for(band in c("Overall","18-19","20-39","40-59","60-69","70-79","80-84")){
    x<-if(band=="Overall")dd else dd%>%filter(as.character(audit_band)==band)
    if(!nrow(x))next
    des<-make_design(x)
    mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE)
    for(metric in c("mean_z","below_p5","below_p10")){
      term<-switch(metric,mean_z="z",below_p5="low_p5",below_p10="low_p10")
      single[[k]]<-tibble(imputation=imp,sex=s,age_band=band,metric=metric,
        estimate=as.numeric(coef(mn)[term]),standard_error=as.numeric(SE(mn)[term]),
        z_sd=safe_sd(des,"z"),unweighted_n=nrow(x));k<-k+1L
    }
  }
}
single<-bind_rows(single)
oofsum<-pool_table(single,c("sex","age_band","metric"))%>%
  left_join(single%>%group_by(sex,age_band,metric)%>%
    summarise(z_sd=mean(z_sd,na.rm=TRUE),unweighted_n=mean(unweighted_n),.groups="drop"),
    by=c("sex","age_band","metric"))%>%
  mutate(percent=if_else(metric%in%c("below_p5","below_p10"),100*estimate,NA_real_))
write_csv(oofsum,file.path(table_dir,"10_corrected_18_84_psu_cv_summary.csv"))

# ---------------- Frozen temporal validation ----------------
validation_paths<-file.path(root_dir,"新建文件夹",
  paste0("07_validation_reference_imputation",1:5,".csv.gz"))
stopifnot(all(file.exists(validation_paths)))
bundles84<-readRDS(file.path(model_dir,"corrected_age18_84_model_bundles.rds"))

predpars<-function(bundle,nd){
  tr<-as.data.frame(bundle$train_data)
  x<-as.data.frame(nd[,intersect(names(nd),names(tr)),drop=FALSE])
  p<-try(suppressWarnings(predictAll(bundle$model,newdata=x,data=tr,type="response")),silent=TRUE)
  if(inherits(p,"try-error"))p<-suppressWarnings(predictAll(bundle$model,newdata=x,type="response"))
  as.data.frame(p)
}
score_bundle<-function(bundle,nd){
  p<-predpars(bundle,nd);pf<-get(paste0("p",bundle$family),asNamespace("gamlss.dist"))
  a<-list(q=nd$alm_kg,mu=p$mu,sigma=p$sigma);if("nu"%in%names(p))a$nu<-p$nu;if("tau"%in%names(p))a$tau<-p$tau
  pr<-pmin(pmax(as.numeric(do.call(pf,a)),1e-10),1-1e-10)
  tibble(percentile_18_84=pr,z_18_84=qnorm(pr))
}

prepare_val<-function(path,imp){
  read_csv(path,show_col_types=FALSE,progress=FALSE)%>%
    transmute(SEQN=as.numeric(SEQN),imputation=imp,cycle=as.character(cycle),
      sex=as.character(sex),race=as.character(race),age=as.numeric(age),
      pregnant=as.logical(pregnant),height_cm=as.numeric(height_cm),
      height_m=as.numeric(height_m),bmi=as.numeric(bmi),alm_kg=as.numeric(alm_kg),
      wtmec_pool=as.numeric(wtmec_pool),strata=as.character(strata),psu=as.character(psu))%>%
    distinct(SEQN,.keep_all=TRUE)%>%
    filter(age>=18,age<=59,sex%in%c("Female","Male"),!coalesce(pregnant,FALSE),
      is.finite(height_m),height_m>0,is.finite(alm_kg),alm_kg>0,
      is.finite(wtmec_pool),wtmec_pool>0,!is.na(strata),!is.na(psu))%>%
    mutate(race=factor(race),cycle=factor(cycle),
      sex=factor(sex,levels=c("Female","Male")),
      strata_pool=factor(interaction(cycle,strata,drop=TRUE)),
      psu_pool=factor(interaction(cycle,strata,psu,drop=TRUE)))
}
vals<-map2(validation_paths,1:5,prepare_val)
vr<-list();k<-1L
for(imp in 1:5)for(s in c("Female","Male")){
  dd<-vals[[imp]]%>%filter(as.character(sex)==s)
  sc<-score_bundle(bundles84[[paste0("imp",imp,"_",s)]],dd)
  vr[[k]]<-bind_cols(dd,sc)%>%
    mutate(low_p5_18_84=as.numeric(percentile_18_84<.05),
           low_p10_18_84=as.numeric(percentile_18_84<.10));k<-k+1L
}
val84<-bind_rows(vr)
saveRDS(val84,file.path(cache_dir,"validation_corrected_18_84_18_59.rds"))

# Merge corrected 18-69 final validation scores.
val69path<-file.path(paper_dir,"Paper1_FINAL_FROZEN_20260814","cache",
                     "validation_corrected_H_18_59.rds")
stopifnot(file.exists(val69path))
val69<-readRDS(val69path)
# Detect corrected H score names robustly.
z69name<-intersect(c("z","z_H","z_18_69","z_corrected_H"),names(val69))[1]
p69name<-intersect(c("percentile","percentile_H","percentile_18_69","percentile_corrected_H"),names(val69))[1]
if(is.na(z69name)||is.na(p69name))stop("Could not identify corrected 18-69 score columns")
val69small<-val69%>%select(SEQN,imputation,all_of(z69name),all_of(p69name))%>%
  rename(z_18_69=all_of(z69name),percentile_18_69=all_of(p69name))
val<-val84%>%left_join(val69small,by=c("SEQN","imputation"))%>%
  mutate(low_p5_18_69=as.numeric(percentile_18_69<.05),
         low_p10_18_69=as.numeric(percentile_18_69<.10),
         abs_delta_z=abs(z_18_84-z_18_69),
         changed_p5=as.numeric(low_p5_18_84!=low_p5_18_69),
         changed_p10=as.numeric(low_p10_18_84!=low_p10_18_69))

# Temporal calibration.
cal_single<-list();k<-1L
for(imp in 1:5)for(s in c("Female","Male")){
  dd<-val%>%filter(imputation==imp,as.character(sex)==s);des<-make_design(dd)
  for(model in c("18-69","18-84")){
    zvar<-paste0("z_",gsub("-","_",model))
    p5<-paste0("low_p5_",gsub("-","_",model));p10<-paste0("low_p10_",gsub("-","_",model))
    mn<-svymean(as.formula(paste0("~",zvar,"+",p5,"+",p10)),des,na.rm=TRUE)
    for(metric in c("mean_z","below_p5","below_p10")){
      term<-switch(metric,mean_z=zvar,below_p5=p5,below_p10=p10)
      cal_single[[k]]<-tibble(imputation=imp,sex=s,model=model,metric=metric,
        estimate=as.numeric(coef(mn)[term]),standard_error=as.numeric(SE(mn)[term]),
        z_sd=safe_sd(des,zvar));k<-k+1L
    }
  }
}
cs<-bind_rows(cal_single)
cal<-pool_table(cs,c("sex","model","metric"))%>%
  left_join(cs%>%group_by(sex,model,metric)%>%summarise(z_sd=mean(z_sd),.groups="drop"),
            by=c("sex","model","metric"))%>%
  mutate(percent=if_else(metric%in%c("below_p5","below_p10"),100*estimate,NA_real_))
write_csv(cal,file.path(table_dir,"11_corrected_temporal_18_69_vs_18_84.csv"))

# Score perturbation/reclassification.
st_single<-list();k<-1L
for(imp in 1:5)for(s in c("Female","Male")){
  dd<-val%>%filter(imputation==imp,as.character(sex)==s);des<-make_design(dd)
  mm<-svymean(~abs_delta_z+changed_p5+changed_p10,des,na.rm=TRUE)
  for(v in c("abs_delta_z","changed_p5","changed_p10")){
    st_single[[k]]<-tibble(imputation=imp,sex=s,metric=v,
      estimate=as.numeric(coef(mm)[v]),standard_error=as.numeric(SE(mm)[v]));k<-k+1L
  }
}
stab<-pool_table(bind_rows(st_single),c("sex","metric"))%>%
  mutate(percent=if_else(grepl("changed",metric),100*estimate,NA_real_))
write_csv(stab,file.path(table_dir,"12_corrected_score_stability_18_69_vs_18_84.csv"))

# Height-related P5/P10 ORs for corrected 18-84 only.
or_single<-list();k<-1L
for(imp in 1:5)for(s in c("Female","Male")){
  dd<-val%>%filter(imputation==imp,as.character(sex)==s)
  for(outcome in c("low_p5_18_84","low_p10_18_84")){
    fit<-svyglm(as.formula(paste(outcome,
      "~ I(height_cm/10) + ns(age,df=3) + race + cycle")),
      design=make_design(dd),family=quasibinomial())
    cf<-coef(fit)["I(height_cm/10)"];se<-SE(fit)["I(height_cm/10)"]
    or_single[[k]]<-tibble(imputation=imp,sex=s,outcome=outcome,
      estimate=as.numeric(cf),standard_error=as.numeric(se));k<-k+1L
  }
}
orp<-pool_table(bind_rows(or_single),c("sex","outcome"))%>%
  mutate(OR=exp(estimate),OR_low=exp(conf_low),OR_high=exp(conf_high))
write_csv(orp,file.path(table_dir,"13_corrected_18_84_height_OR.csv"))

cat("\n=== Corrected 18-84 PSU-CV ===\n")
print(oofsum%>%filter(age_band%in%c("Overall","70-79","80-84")))
cat("\n=== Corrected temporal 18-69 vs 18-84 ===\n");print(cal)
cat("\n=== Corrected score stability ===\n");print(stab)
cat("\n=== Corrected 18-84 height OR ===\n");print(orp)
