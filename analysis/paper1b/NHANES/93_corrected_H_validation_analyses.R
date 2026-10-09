options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(purrr)
  library(gamlss); library(gamlss.dist); library(survey); library(splines)
})

set.seed(20260814)
root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
paper_dir <- file.path(root_dir, "First Paper")
out_dir <- file.path(paper_dir, "Paper1_FINAL_FROZEN_20260814")
cache_dir <- file.path(out_dir, "cache")
model_dir <- file.path(out_dir, "models")
old_validation_path <- file.path(paper_dir, "final_age_range_comparison_20260814", "cache",
                                 "validation_scores_three_scorers_18_59.rds")
stopifnot(file.exists(old_validation_path))

weighted_quantile <- function(x,w,p) {
  ok <- is.finite(x)&is.finite(w)&w>0; x<-x[ok];w<-w[ok];o<-order(x);x<-x[o];w<-w[o]
  vapply(p,function(q)x[which(cumsum(w)/sum(w)>=q)[1]],numeric(1))
}
make_design <- function(dat) svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wtmec_pool,
                                        nest=TRUE,data=dat)
rubin_pool <- function(q,se,null=0) {
  ok<-is.finite(q)&is.finite(se)&se>=0;q<-q[ok];se<-se[ok];m<-length(q)
  if(!m)return(tibble(m=0,estimate=NA_real_,standard_error=NA_real_,df=NA_real_,
                      conf_low=NA_real_,conf_high=NA_real_,p_value=NA_real_))
  qbar<-mean(q);u<-mean(se^2);b<-if(m>1)var(q)else 0;tvar<-u+(1+1/m)*b;s<-sqrt(tvar)
  df<-if(m>1&&b>0)(m-1)*(1+u/((1+1/m)*b))^2 else Inf;crit<-qt(.975,df)
  tibble(m=m,estimate=qbar,standard_error=s,df=df,conf_low=qbar-crit*s,conf_high=qbar+crit*s,
         p_value=2*pt(abs((qbar-null)/s),df=df,lower.tail=FALSE))
}
pool_table <- function(dat,groups) dat%>%group_by(across(all_of(groups)))%>%
  group_modify(~rubin_pool(.x$estimate,.x$standard_error))%>%ungroup()
safe_mean <- function(des,var) {
  x<-try(svymean(as.formula(paste0("~",var)),des,na.rm=TRUE),silent=TRUE)
  if(inherits(x,"try-error"))return(c(estimate=NA_real_,se=NA_real_))
  c(estimate=as.numeric(coef(x)[1]),se=as.numeric(SE(x)[1]))
}
safe_sd <- function(des,var) {
  x<-try(svyvar(as.formula(paste0("~",var)),des,na.rm=TRUE),silent=TRUE)
  if(inherits(x,"try-error"))return(NA_real_)
  sqrt(pmax(as.numeric(coef(x)[1]),0))
}

# Score the unchanged 18-69 H specification after correcting development weights.
validation <- readRDS(old_validation_path) %>% mutate(
  sex=factor(as.character(sex),levels=c("Female","Male")),race=factor(as.character(race)),
  cycle=factor(as.character(cycle)),strata_pool=factor(strata_pool),psu_pool=factor(psu_pool)
)

predict_parameters <- function(bundle,new_data) {
  train<-as.data.frame(bundle$train_data)
  nd<-as.data.frame(new_data[,intersect(names(new_data),names(train)),drop=FALSE])
  ans<-try(suppressWarnings(predictAll(bundle$model,newdata=nd,data=train,type="response")),silent=TRUE)
  if(inherits(ans,"try-error"))ans<-suppressWarnings(predictAll(bundle$model,newdata=nd,type="response"))
  as.data.frame(ans)
}
score_bundle <- function(bundle,new_data) {
  pars<-predict_parameters(bundle,new_data);pfun<-get(paste0("p",bundle$family),asNamespace("gamlss.dist"))
  args<-list(q=new_data$alm_kg,mu=pars$mu,sigma=pars$sigma)
  if("nu"%in%names(pars))args$nu<-pars$nu;if("tau"%in%names(pars))args$tau<-pars$tau
  p<-as.numeric(do.call(pfun,args));p<-pmin(pmax(p,1e-10),1-1e-10)
  tibble(percentile_H_corrected=p,z_H_corrected=qnorm(p))
}

score_rows<-list();counter<-1L
for(sex_value in c("Female","Male")){
  bundles<-readRDS(file.path(model_dir,paste0("corrected_H_18_69_",sex_value,"_bundles.rds")))
  for(imp in 1:5){
    dat<-validation%>%filter(imputation==imp,as.character(sex)==sex_value)
    sc<-score_bundle(bundles[[paste0("imp",imp,"_",sex_value)]],dat)
    score_rows[[counter]]<-bind_cols(dat,sc);counter<-counter+1L
  }
}
validation_new<-bind_rows(score_rows)%>%mutate(
  low_p5_H_corrected=as.numeric(percentile_H_corrected<.05),
  low_p10_H_corrected=as.numeric(percentile_H_corrected<.10),
  low_p5_H_legacy=low_p5_18_69,low_p10_H_legacy=low_p10_18_69
)
saveRDS(validation_new,file.path(cache_dir,"validation_corrected_H_18_59.rds"))
write_csv(validation_new%>%select(SEQN,imputation,cycle,sex,race,age,height_cm,height_m,bmi,alm_kg,almi,
  alm_bmi_ratio,body_fat_pct,grip_strength_kg,wtmec_pool,strata,psu,strata_pool,psu_pool,
  z_H_corrected,percentile_H_corrected,low_p5_H_corrected,low_p10_H_corrected,
  z_H_legacy=z_18_69,percentile_H_legacy=percentile_18_69,low_p5_H_legacy,low_p10_H_legacy,
  low_ewgsop2,low_fnih),file.path(cache_dir,"validation_corrected_H_18_59.csv.gz"))

# Old-vs-corrected impact audit in the same validation participants.
impact_single<-map_dfr(1:5,function(imp)map_dfr(c("Female","Male"),function(s){
  dat<-validation_new%>%filter(imputation==imp,as.character(sex)==s)%>%mutate(
    abs_delta_z=abs(z_H_corrected-z_18_69),p5_changed=as.numeric(low_p5_H_corrected!=low_p5_H_legacy),
    p10_changed=as.numeric(low_p10_H_corrected!=low_p10_H_legacy))
  des<-make_design(dat);m<-svymean(~abs_delta_z+p5_changed+p10_changed,des,na.rm=TRUE)
  tibble(imputation=imp,sex=s,mean_abs_delta_z=coef(m)[1],mean_abs_delta_z_se=SE(m)[1],
    p5_reclassification=coef(m)[2],p5_reclassification_se=SE(m)[2],
    p10_reclassification=coef(m)[3],p10_reclassification_se=SE(m)[3],
    weighted_correlation=cor(dat$z_H_corrected,dat$z_18_69,use="complete.obs"),n=nrow(dat))
}))
impact<-bind_rows(
  impact_single%>%transmute(imputation,sex,metric="mean_abs_delta_z",estimate=mean_abs_delta_z,standard_error=mean_abs_delta_z_se),
  impact_single%>%transmute(imputation,sex,metric="p5_reclassification",estimate=p5_reclassification,standard_error=p5_reclassification_se),
  impact_single%>%transmute(imputation,sex,metric="p10_reclassification",estimate=p10_reclassification,standard_error=p10_reclassification_se)
)%>%pool_table(c("sex","metric"))%>%left_join(impact_single%>%group_by(sex)%>%summarise(
  weighted_correlation=mean(weighted_correlation),unweighted_n=mean(n),.groups="drop"),by="sex")
write_csv(impact,file.path(cache_dir,"corrected_weight_impact_audit.csv"))

# Corrected-weight internal PSU-CV summary.
oof<-map_dfr(c("Female","Male"),~read_csv(file.path(cache_dir,paste0("corrected_H_OOF_",.x,".csv.gz")),
                                           show_col_types=FALSE,progress=FALSE))%>%
  mutate(strata_pool=factor(strata_pool),psu_pool=factor(psu_pool),wtmec_pool=wtmec_correct)
cv_single<-map_dfr(1:5,function(imp)map_dfr(c("Female","Male"),function(s){
  dat<-oof%>%filter(imputation==imp,sex==s);des<-make_design(dat)
  mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE)
  tibble(imputation=imp,sex=s,metric=c("mean_z","P5","P10"),estimate=as.numeric(coef(mn)),
         standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z"),unweighted_n=nrow(dat))
}))
cv_pooled<-cv_single%>%pool_table(c("sex","metric"))%>%left_join(cv_single%>%group_by(sex,metric)%>%
  summarise(z_sd=mean(z_sd),unweighted_n=mean(unweighted_n),.groups="drop"),by=c("sex","metric"))%>%
  mutate(percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),
         percent_low=if_else(metric%in%c("P5","P10"),100*pmax(conf_low,0),NA_real_),
         percent_high=if_else(metric%in%c("P5","P10"),100*pmin(conf_high,1),NA_real_),
         weight_version="Corrected pooled development weights")
write_csv(cv_pooled,file.path(out_dir,"04_internal_PSU_CV.csv"))

# Reusable H temporal calibration by arbitrary group.
calibration_by<-function(data,group_var,group_label){
  single<-map_dfr(1:5,function(imp){
    dat_imp<-data%>%filter(imputation==imp)
    groups<-if(group_var=="overall")"Overall" else unique(as.character(dat_imp[[group_var]]))
    map_dfr(c("Female","Male"),function(s)map_dfr(groups,function(g){
      dat<-dat_imp%>%filter(as.character(sex)==s)
      if(group_var!="overall")dat<-dat%>%filter(as.character(.data[[group_var]])==g)
      des<-make_design(dat);mn<-svymean(~z_H_corrected+low_p5_H_corrected+low_p10_H_corrected,des,na.rm=TRUE)
      tibble(imputation=imp,sex=s,group_type=group_label,group=g,metric=c("mean_z","P5","P10"),
             estimate=as.numeric(coef(mn)),standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z_H_corrected"),
             unweighted_n=nrow(dat))
    }))
  })
  single%>%pool_table(c("sex","group_type","group","metric"))%>%left_join(single%>%
    group_by(sex,group_type,group,metric)%>%summarise(z_sd=mean(z_sd),unweighted_n=mean(unweighted_n),.groups="drop"),
    by=c("sex","group_type","group","metric"))%>%mutate(
      percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),
      percent_low=if_else(metric%in%c("P5","P10"),100*pmax(0,conf_low),NA_real_),
      percent_high=if_else(metric%in%c("P5","P10"),100*pmin(1,conf_high),NA_real_))
}
temporal_overall<-calibration_by(validation_new,"overall","overall")
cycle_cal<-calibration_by(validation_new,"cycle","cycle")
race_cal<-calibration_by(validation_new,"race","race_ethnicity")
write_csv(temporal_overall,file.path(out_dir,"05_temporal_calibration.csv"))
write_csv(cycle_cal,file.path(out_dir,"06_cycle_calibration.csv"))
write_csv(race_cal,file.path(out_dir,"07_race_calibration.csv"))

# Restricted cubic spline basis with four weighted knots (P5/P35/P65/P95).
rcs_basis<-function(height_cm,knots_cm){
  x<-height_cm/10;k<-knots_cm/10;K<-length(k);tp<-function(v)pmax(v,0)^3
  cols<-lapply(1:(K-2),function(j){
    (tp(x-k[j])-tp(x-k[K-1])*(k[K]-k[j])/(k[K]-k[K-1])+
       tp(x-k[K])*(k[K-1]-k[j])/(k[K]-k[K-1]))/(k[K]-k[1])^2
  })
  out<-data.frame(h_linear=x);for(j in seq_along(cols))out[[paste0("h_nl",j)]]<-cols[[j]];out
}
pool_multi_test<-function(betas,covs,indices){
  m<-length(betas);Q<-do.call(rbind,lapply(betas,function(x)x[indices]));qbar<-colMeans(Q)
  U<-Reduce("+",lapply(covs,function(x)x[indices,indices,drop=FALSE]))/m
  B<-if(m>1)cov(Q)else matrix(0,length(indices),length(indices));T<-U+(1+1/m)*B
  inv<-try(solve(T),silent=TRUE);if(inherits(inv,"try-error"))inv<-qr.solve(T)
  stat<-as.numeric(t(qbar)%*%inv%*%qbar);c(statistic=stat,df=length(indices),p_value=pchisq(stat,length(indices),lower.tail=FALSE))
}
marginal_prediction<-function(fit,dat,height_value,knots){
  nd<-dat;nd$height_cm<-height_value;nd$height_m<-height_value/100
  basis<-rcs_basis(nd$height_cm,knots);nd$h_linear<-basis$h_linear;nd$h_nl1<-basis$h_nl1;nd$h_nl2<-basis$h_nl2
  X<-model.matrix(delete.response(terms(fit)),nd);bn<-names(coef(fit));miss<-setdiff(bn,colnames(X))
  if(length(miss))X<-cbind(X,matrix(0,nrow(X),length(miss),dimnames=list(NULL,miss)));X<-X[,bn,drop=FALSE]
  beta<-coef(fit);V<-vcov(fit);eta<-as.numeric(X%*%beta);p<-plogis(eta);w<-dat$wtmec_pool
  est<-weighted.mean(p,w);grad<-colSums(X*(p*(1-p)*w))/sum(w);se<-sqrt(as.numeric(t(grad)%*%V%*%grad))
  c(estimate=est,standard_error=se)
}

outcome_map<-c(conditional_P5="low_p5_H_corrected",conditional_P10="low_p10_H_corrected",
               EWGSOP2_low="low_ewgsop2",FNIH_low="low_fnih")
height_tests<-list();height_preds<-list();ht<-hp<-1L
for(s in c("Female","Male")){
  ref<-validation_new%>%filter(imputation==1,as.character(sex)==s)
  knots<-weighted_quantile(ref$height_cm,ref$wtmec_pool,c(.05,.35,.65,.95))
  points<-weighted_quantile(ref$height_cm,ref$wtmec_pool,c(.10,.25,.50,.75,.90));names(points)<-c("P10","P25","P50","P75","P90")
  for(label in names(outcome_map)){
    linear_q<-linear_se<-numeric();spline_b<-list();spline_v<-list();pred_single<-list();ps<-1L
    for(imp in 1:5){
      dat<-validation_new%>%filter(imputation==imp,as.character(sex)==s)%>%droplevels()
      basis<-rcs_basis(dat$height_cm,knots);dat$h_linear<-basis$h_linear;dat$h_nl1<-basis$h_nl1;dat$h_nl2<-basis$h_nl2
      outcome<-outcome_map[[label]]
      f_lin<-as.formula(paste0(outcome,"~h_linear+ns(age,df=3)+race+cycle"))
      f_spl<-as.formula(paste0(outcome,"~h_linear+h_nl1+h_nl2+ns(age,df=3)+race+cycle"))
      fit_lin<-svyglm(f_lin,design=make_design(dat),family=quasibinomial())
      fit_spl<-svyglm(f_spl,design=make_design(dat),family=quasibinomial())
      linear_q[imp]<-coef(fit_lin)["h_linear"];linear_se[imp]<-SE(fit_lin)["h_linear"]
      idx<-match(c("h_linear","h_nl1","h_nl2"),names(coef(fit_spl)))
      spline_b[[imp]]<-coef(fit_spl)[idx];spline_v[[imp]]<-vcov(fit_spl)[idx,idx,drop=FALSE]
      for(pt in names(points)){
        pr<-marginal_prediction(fit_spl,dat,points[[pt]],knots)
        pred_single[[ps]]<-tibble(imputation=imp,height_point=pt,height_cm=points[[pt]],
          estimate=pr["estimate"],standard_error=pr["standard_error"]);ps<-ps+1L
      }
    }
    lin<-rubin_pool(linear_q,linear_se);global<-pool_multi_test(spline_b,spline_v,1:3);nonlin<-pool_multi_test(spline_b,spline_v,2:3)
    height_tests[[ht]]<-lin%>%transmute(sex=s,outcome=label,record_type="model_test",height_point=NA_character_,
      height_cm=NA_real_,linear_log_OR_per10cm=estimate,linear_OR_per10cm=exp(estimate),
      linear_OR_low=exp(conf_low),linear_OR_high=exp(conf_high),linear_p=p_value,
      spline_global_p=global["p_value"],spline_nonlinearity_p=nonlin["p_value"],
      predicted_prevalence=NA_real_,predicted_low=NA_real_,predicted_high=NA_real_,knot_spec=paste(round(knots,1),collapse=";"));ht<-ht+1L
    pp<-bind_rows(pred_single)%>%pool_table(c("height_point","height_cm"))
    height_preds[[hp]]<-pp%>%transmute(sex=s,outcome=label,record_type="predicted_prevalence",height_point,height_cm,
      linear_log_OR_per10cm=NA_real_,linear_OR_per10cm=NA_real_,linear_OR_low=NA_real_,linear_OR_high=NA_real_,linear_p=NA_real_,
      spline_global_p=global["p_value"],spline_nonlinearity_p=nonlin["p_value"],
      predicted_prevalence=estimate,predicted_low=pmax(0,conf_low),predicted_high=pmin(1,conf_high),
      knot_spec=paste(round(knots,1),collapse=";"));hp<-hp+1L
  }
}
height_linearity<-bind_rows(bind_rows(height_tests),bind_rows(height_preds))
write_csv(height_linearity,file.path(out_dir,"08_height_linearity.csv"))

# Height-quintile classification drift for all four definitions.
cut_rows<-list();prev_single<-list();cc<-pc<-1L
for(s in c("Female","Male")){
  ref<-validation_new%>%filter(imputation==1,as.character(sex)==s);cuts<-weighted_quantile(ref$height_cm,ref$wtmec_pool,seq(0,1,.2))
  cuts[1]<--Inf;cuts[6]<-Inf;cut_rows[[cc]]<-tibble(sex=s,boundary=0:5,height_cm=cuts);cc<-cc+1L
  for(imp in 1:5){
    dat<-validation_new%>%filter(imputation==imp,as.character(sex)==s)%>%mutate(height_quintile=cut(height_cm,cuts,labels=paste0("Q",1:5),include.lowest=TRUE))
    for(q in paste0("Q",1:5))for(label in names(outcome_map)){
      sub<-dat%>%filter(as.character(height_quintile)==q);val<-safe_mean(make_design(sub),outcome_map[[label]])
      prev_single[[pc]]<-tibble(imputation=imp,sex=s,outcome=label,height_quintile=q,n=nrow(sub),estimate=val["estimate"],standard_error=val["se"]);pc<-pc+1L
    }
  }
}
height_prev<-bind_rows(prev_single)%>%pool_table(c("sex","outcome","height_quintile"))%>%
  left_join(bind_rows(prev_single)%>%group_by(sex,outcome,height_quintile)%>%summarise(unweighted_n=mean(n),.groups="drop"),
            by=c("sex","outcome","height_quintile"))%>%mutate(percent=100*estimate,percent_low=100*pmax(0,conf_low),percent_high=100*pmin(1,conf_high))
height_drift<-bind_rows(
  height_linearity%>%filter(record_type=="model_test")%>%mutate(section="continuous_height"),
  height_prev%>%mutate(section="height_quintile")
)
write_csv(height_drift,file.path(out_dir,"09_height_classification_drift.csv"))

# Full-sample and normal-BMI grip discordance analyses, 2011-2014, age 18-59.
discordance_levels<-c("Neither low","Conventional-only low","Conditional-only low","Both low")
add_discordance<-function(dat,conventional){
  cnd<-dat$low_p5_H_corrected;conv<-dat[[conventional]]
  dat$discordance<-factor(case_when(cnd==0&conv==0~"Neither low",cnd==0&conv==1~"Conventional-only low",
    cnd==1&conv==0~"Conditional-only low",cnd==1&conv==1~"Both low",TRUE~NA_character_),levels=discordance_levels);dat
}
marginal_group_vector<-function(fit,dat,group_value){
  nd<-dat;nd$discordance<-factor(group_value,levels=levels(dat$discordance));X<-model.matrix(delete.response(terms(fit)),nd)
  bn<-names(coef(fit));miss<-setdiff(bn,colnames(X));if(length(miss))X<-cbind(X,matrix(0,nrow(X),length(miss),dimnames=list(NULL,miss)))
  X<-X[,bn,drop=FALSE];as.numeric(colSums(X*dat$wtmec_pool)/sum(dat$wtmec_pool))
}
grip_models<-tribble(~model,~required,~adjustment,
  "Raw","base","1","M1","base","ns(age,df=3)+sex+height_m+race+cycle",
  "M2","bmi","ns(age,df=3)+sex+height_m+race+cycle+bmi",
  "M3","fat","ns(age,df=3)+sex+height_m+race+cycle+bmi+body_fat_pct")
grip_run<-function(domain){
  single<-list();flow<-list();counter<-fc<-1L
  for(imp in 1:5)for(comparison in c("FNIH","EWGSOP2")){
    conventional<-if(comparison=="FNIH")"low_fnih" else "low_ewgsop2"
    base<-validation_new%>%filter(imputation==imp,as.character(cycle)%in%c("2011-2012","2013-2014"))
    if(domain=="Normal BMI")base<-base%>%filter(is.finite(bmi),bmi>=18.5,bmi<25)
    base<-add_discordance(base,conventional)%>%filter(!is.na(discordance))
    for(i in seq_len(nrow(grip_models))){
      spec<-grip_models[i,];dat<-base%>%filter(is.finite(grip_strength_kg))
      if(spec$required%in%c("bmi","fat"))dat<-dat%>%filter(is.finite(bmi))
      if(spec$required=="fat")dat<-dat%>%filter(is.finite(body_fat_pct))
      dat<-dat%>%filter(is.finite(age),is.finite(height_m),!is.na(sex),!is.na(race),!is.na(cycle))
      present<-discordance_levels[discordance_levels%in%as.character(unique(dat$discordance))]
      if(!"Neither low"%in%present||length(present)<2)next
      dat<-dat%>%mutate(discordance=factor(as.character(discordance),levels=present))%>%droplevels()
      formula<-if(spec$model=="Raw")grip_strength_kg~discordance else as.formula(paste("grip_strength_kg~discordance+",spec$adjustment))
      fit<-svyglm(formula,design=make_design(dat));beta<-coef(fit);V<-vcov(fit)
      vectors<-setNames(map(present,~marginal_group_vector(fit,dat,.x)),present);ref<-vectors[["Neither low"]]
      for(g in present){
        v<-vectors[[g]];contrast<-v-ref;indicator<-as.numeric(as.character(dat$discordance)==g)
        des<-make_design(dat%>%mutate(group_indicator=indicator));pv<-safe_mean(des,"group_indicator")
        single[[counter]]<-tibble(imputation=imp,domain=domain,comparison=comparison,model=spec$model,group=g,
          analysis_n=nrow(dat),unweighted_group_n=sum(indicator),prevalence=pv["estimate"],prevalence_se=pv["se"],
          grip_mean=sum(v*beta),grip_mean_se=sqrt(as.numeric(t(v)%*%V%*%v)),
          difference=sum(contrast*beta),difference_se=sqrt(as.numeric(t(contrast)%*%V%*%contrast)));counter<-counter+1L
      }
      flow[[fc]]<-tibble(imputation=imp,domain=domain,comparison=comparison,model=spec$model,analysis_n=nrow(dat),
        neither_n=sum(as.character(dat$discordance)=="Neither low"),
        conventional_only_n=sum(as.character(dat$discordance)=="Conventional-only low"),
        conditional_only_n=sum(as.character(dat$discordance)=="Conditional-only low"),
        both_n=sum(as.character(dat$discordance)=="Both low"));fc<-fc+1L
    }
  }
  s<-bind_rows(single);sizes<-s%>%group_by(domain,comparison,model,group)%>%summarise(
    analysis_n=mean(analysis_n),unweighted_group_n=mean(unweighted_group_n),.groups="drop")
  pooled<-bind_rows(
    s%>%transmute(imputation,domain,comparison,model,group,measure="weighted_prevalence",estimate=prevalence,standard_error=prevalence_se),
    s%>%transmute(imputation,domain,comparison,model,group,measure="grip_mean_kg",estimate=grip_mean,standard_error=grip_mean_se),
    s%>%transmute(imputation,domain,comparison,model,group,measure="difference_vs_neither_kg",estimate=difference,standard_error=difference_se)
  )%>%pool_table(c("domain","comparison","model","group","measure"))%>%left_join(sizes,by=c("domain","comparison","model","group"))%>%
    mutate(percent=if_else(measure=="weighted_prevalence",100*estimate,NA_real_),
           percent_low=if_else(measure=="weighted_prevalence",100*pmax(0,conf_low),NA_real_),
           percent_high=if_else(measure=="weighted_prevalence",100*pmin(1,conf_high),NA_real_))
  list(results=pooled,flow=bind_rows(flow)%>%group_by(domain,comparison,model)%>%summarise(across(c(analysis_n,neither_n,
    conventional_only_n,conditional_only_n,both_n),mean),min_analysis_n=min(analysis_n),max_analysis_n=max(analysis_n),.groups="drop"))
}
grip_full<-grip_run("Full sample");grip_normal<-grip_run("Normal BMI")
write_csv(grip_full$results,file.path(out_dir,"10_fullsample_grip_discordance.csv"))
write_csv(grip_normal$results,file.path(out_dir,"11_normalBMI_grip_sensitivity.csv"))
write_csv(bind_rows(grip_full$flow,grip_normal$flow),file.path(cache_dir,"grip_sample_flow.csv"))

cat("Corrected H validation, height, calibration, and grip analyses completed at",out_dir,"\n")
