options(stringsAsFactors=FALSE,survey.lonely.psu="adjust")
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(purrr);library(gamlss);library(gamlss.dist);library(survey);library(splines)})
set.seed(20260814)
root_dir<-normalizePath(".",winslash="/",mustWork=TRUE);out_dir<-file.path(root_dir,"First Paper","Paper1_FINAL_FROZEN_20260814")
cache_dir<-file.path(out_dir,"cache");model_dir<-file.path(out_dir,"models");fig_dir<-file.path(out_dir,"figures");dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)

make_design<-function(dat)svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wtmec_pool,nest=TRUE,data=dat)
rubin<-function(q,se,null=0){ok<-is.finite(q)&is.finite(se)&se>=0;q<-q[ok];se<-se[ok];m<-length(q);if(!m)return(tibble(m=0,estimate=NA_real_,standard_error=NA_real_,df=NA_real_,conf_low=NA_real_,conf_high=NA_real_,p_value=NA_real_));
  qb<-mean(q);u<-mean(se^2);b<-if(m>1)var(q)else 0;t<-u+(1+1/m)*b;s<-sqrt(t);df<-if(m>1&&b>0)(m-1)*(1+u/((1+1/m)*b))^2 else Inf;cr<-qt(.975,df)
  tibble(m=m,estimate=qb,standard_error=s,df=df,conf_low=qb-cr*s,conf_high=qb+cr*s,p_value=2*pt(abs((qb-null)/s),df=df,lower.tail=FALSE))}
pool_table<-function(d,g)d%>%group_by(across(all_of(g)))%>%group_modify(~rubin(.x$estimate,.x$standard_error))%>%ungroup()
safe_sd<-function(des,var){x<-svyvar(as.formula(paste0("~",var)),des,na.rm=TRUE);sqrt(pmax(as.numeric(coef(x)[1]),0))}
wq<-function(x,w,p){ok<-is.finite(x)&is.finite(w)&w>0;x<-x[ok];w<-w[ok];o<-order(x);x<-x[o];w<-w[o];vapply(p,function(q)x[which(cumsum(w)/sum(w)>=q)[1]],numeric(1))}

# ---------------- Corrected development OOF ----------------
specs<-crossing(sex=c("Female","Male"),model=c("Hcommon","HB"))%>%mutate(path=file.path(cache_dir,paste0("corrected_",sex,"_",model,"_OOF.csv.gz")))
oof<-map_dfr(seq_len(nrow(specs)),function(i)read_csv(specs$path[i],show_col_types=FALSE,progress=FALSE)%>%mutate(sex=specs$sex[i],model=specs$model[i]))%>%
  mutate(wtmec_pool=wtmec_correct,strata_pool=factor(strata_pool),psu_pool=factor(psu_pool),cycle=factor(cycle),bmi_group=cut(bmi,c(-Inf,18.5,25,30,35,Inf),right=FALSE,
    labels=c("<18.5","18.5-24.9","25.0-29.9","30.0-34.9",">=35")))

overall_single<-map_dfr(1:5,function(imp)map_dfr(c("Female","Male"),function(s)map_dfr(c("Hcommon","HB"),function(mo){d<-oof%>%filter(imputation==imp,sex==s,model==mo);des<-make_design(d);mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE);tibble(imputation=imp,sex=s,model=mo,metric=c("mean_z","P5","P10"),estimate=as.numeric(coef(mn)),standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z"),n=nrow(d))})))
overall<-overall_single%>%pool_table(c("sex","model","metric"))%>%left_join(overall_single%>%group_by(sex,model,metric)%>%summarise(z_sd=mean(z_sd),unweighted_n=mean(n),.groups="drop"),by=c("sex","model","metric"))%>%
  mutate(percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),record_type="OOF_overall")

effects_single<-list();ec<-1L
for(imp in 1:5)for(s in c("Female","Male"))for(mo in c("Hcommon","HB")){
  d<-oof%>%filter(imputation==imp,sex==s,model==mo)%>%droplevels();des<-make_design(d)
  for(outcome in c("z","low_p5","low_p10"))for(exposure in c("height","bmi")){
    x<-if(exposure=="height")"I(height_cm/10)"else"I(bmi/5)";fam<-if(outcome=="z")gaussian()else quasibinomial()
    fit<-svyglm(as.formula(paste(outcome,"~",x,"+ns(age,df=3)+cycle")),design=des,family=fam);idx<-2
    effects_single[[ec]]<-tibble(imputation=imp,sex=s,model=mo,outcome=outcome,exposure=exposure,unit=ifelse(exposure=="height","10 cm","5 kg/m2"),
      estimate=coef(fit)[idx],standard_error=SE(fit)[idx]);ec<-ec+1L
  }
}
effects<-bind_rows(effects_single)%>%pool_table(c("sex","model","outcome","exposure","unit"))%>%mutate(
  effect_type=if_else(outcome=="z","beta","odds_ratio"),odds_ratio=if_else(outcome=="z",NA_real_,exp(estimate)),
  or_low=if_else(outcome=="z",NA_real_,exp(conf_low)),or_high=if_else(outcome=="z",NA_real_,exp(conf_high)),record_type="OOF_effect")

cut_map<-list();for(s in c("Female","Male")){ref<-oof%>%filter(imputation==1,sex==s,model=="Hcommon");cuts<-wq(ref$height_cm,ref$wtmec_pool,seq(0,1,.2));cuts[1]<--Inf;cuts[6]<-Inf;cut_map[[s]]<-cuts}
oof<-oof%>%rowwise()%>%ungroup();oof$height_quintile<-NA_character_
for(s in c("Female","Male"))oof$height_quintile[oof$sex==s]<-as.character(cut(oof$height_cm[oof$sex==s],cut_map[[s]],labels=paste0("Q",1:5),include.lowest=TRUE))
strat_single<-list();sc<-1L
for(imp in 1:5)for(s in c("Female","Male"))for(mo in c("Hcommon","HB"))for(strat in c("height_quintile","bmi_group")){
  groups<-if(strat=="height_quintile")paste0("Q",1:5)else levels(oof$bmi_group)
  for(g in groups){d<-oof%>%filter(imputation==imp,sex==s,model==mo,as.character(.data[[strat]])==g);des<-make_design(d);mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE)
    strat_single[[sc]]<-tibble(imputation=imp,sex=s,model=mo,stratifier=strat,group=g,metric=c("mean_z","P5","P10"),estimate=as.numeric(coef(mn)),standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z"),n=nrow(d));sc<-sc+1L}}
stratified<-bind_rows(strat_single)%>%pool_table(c("sex","model","stratifier","group","metric"))%>%left_join(bind_rows(strat_single)%>%group_by(sex,model,stratifier,group,metric)%>%summarise(z_sd=mean(z_sd),unweighted_n=mean(n),.groups="drop"),by=c("sex","model","stratifier","group","metric"))%>%
  mutate(percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),record_type="OOF_stratified")

fit_summary<-map_dfr(seq_len(nrow(specs)),function(i)read_csv(file.path(cache_dir,paste0("corrected_",specs$sex[i],"_",specs$model[i],"_full_fit.csv")),show_col_types=FALSE))%>%
  group_by(sex,model)%>%summarise(mean_BIC=mean(BIC),sum_BIC=sum(BIC),mean_AIC=mean(AIC),mean_edf=mean(effective_df),all_converged=all(converged),.groups="drop")%>%
  group_by(sex)%>%mutate(delta_mean_BIC_vs_H=mean_BIC-mean_BIC[model=="Hcommon"])%>%ungroup()%>%mutate(record_type="full_fit_BIC")
stability<-map_dfr(seq_len(nrow(specs)),function(i)read_csv(file.path(cache_dir,paste0("corrected_",specs$sex[i],"_",specs$model[i],"_CV_audit.csv")),show_col_types=FALSE))%>%
  group_by(sex,model)%>%summarise(folds=n(),all_converged=all(converged),valid_predictions=sum(valid_predictions),expected_predictions=sum(expected_predictions),
    fallback_folds=sum(attempt>1),max_seconds=max(elapsed_seconds),.groups="drop")%>%mutate(record_type="CV_stability")
development<-bind_rows(overall,effects,stratified,fit_summary,stability)%>%mutate(weight_version="Corrected pooled development weights")
write_csv(development,file.path(out_dir,"14_H_vs_HB_development.csv"))

# ---------------- Frozen temporal validation ----------------
base<-readRDS(file.path(cache_dir,"validation_corrected_H_18_59.rds"))%>%filter(is.finite(bmi),bmi>0)%>%mutate(
  sex=as.character(sex),race=factor(as.character(race)),cycle=factor(as.character(cycle)),strata_pool=factor(strata_pool),psu_pool=factor(psu_pool),
  bmi10=bmi/10,bmi_group=cut(bmi,c(-Inf,18.5,25,30,35,Inf),right=FALSE,labels=c("<18.5","18.5-24.9","25.0-29.9","30.0-34.9",">=35")))
predict_pars<-function(bundle,nd){
  tr<-as.data.frame(bundle$train_data);x<-as.data.frame(nd[,intersect(names(nd),names(tr)),drop=FALSE])
  # Frozen safe prediction domain: clamp predictors only to the observed
  # development support.  ALM (the response being scored) is never clamped.
  for(v in intersect(c("age","height_cm","height_m","bmi","bmi10"),names(x))){
    lo<-min(tr[[v]],na.rm=TRUE);hi<-max(tr[[v]],na.rm=TRUE);x[[v]]<-pmin(pmax(x[[v]],lo),hi)
  }
  a<-try(suppressWarnings(predictAll(bundle$fit,newdata=x,data=tr,type="response")),silent=TRUE)
  if(inherits(a,"try-error"))a<-suppressWarnings(predictAll(bundle$fit,newdata=x,type="response"));as.data.frame(a)
}
score<-function(bundle,nd){
  p<-predict_pars(bundle,nd);pf<-get(paste0("p",bundle$family),asNamespace("gamlss.dist"));valid<-is.finite(p$mu)&p$mu>0&is.finite(p$sigma)&p$sigma>0
  pr<-rep(NA_real_,nrow(nd));a<-list(q=nd$alm_kg[valid],mu=p$mu[valid],sigma=p$sigma[valid]);if("nu"%in%names(p))a$nu<-p$nu[valid];if("tau"%in%names(p))a$tau<-p$tau[valid]
  pr[valid]<-pmin(pmax(as.numeric(do.call(pf,a)),1e-10),1-1e-10)
  tibble(percentile=pr,z=qnorm(pr),low_p5=if_else(is.finite(pr),as.numeric(pr<.05),NA_real_),
         low_p10=if_else(is.finite(pr),as.numeric(pr<.10),NA_real_),prediction_valid=is.finite(pr))
}
temporal_rows<-list();support_rows<-list();tc<-sr<-1L
for(s in c("Female","Male"))for(mo in c("Hcommon","HB")){
  bundles<-readRDS(file.path(model_dir,paste0("corrected_",s,"_",mo,"_bundles.rds")))
  for(imp in 1:5){
    d<-base%>%filter(imputation==imp,sex==s);bundle<-bundles[[paste0("imp",imp)]];tr<-bundle$train_data
    vars<-intersect(c("age","height_cm","height_m","bmi","bmi10"),names(tr));outside<-rep(FALSE,nrow(d))
    for(v in vars)outside<-outside|d[[v]]<min(tr[[v]],na.rm=TRUE)|d[[v]]>max(tr[[v]],na.rm=TRUE)
    support_rows[[sr]]<-tibble(sex=s,model=mo,imputation=imp,n=nrow(d),outside_support_n=sum(outside),outside_support_pct=100*mean(outside))
    scored<-try(score(bundle,d),silent=TRUE)
    if(inherits(scored,"try-error"))stop("Temporal scoring failed for ",s," ",mo," imp",imp,": ",as.character(scored))
    support_rows[[sr]]$invalid_prediction_n<-sum(!scored$prediction_valid)
    support_rows[[sr]]$invalid_prediction_pct<-100*mean(!scored$prediction_valid)
    temporal_rows[[tc]]<-bind_cols(d,scored)%>%mutate(model=mo);tc<-tc+1L;sr<-sr+1L
  }}
write_csv(bind_rows(support_rows),file.path(cache_dir,"HB_temporal_safe_domain_audit.csv"))
tv_all<-bind_rows(temporal_rows)
valid_keys<-tv_all%>%group_by(imputation,sex,SEQN)%>%summarise(models=n_distinct(model),all_models_valid=all(prediction_valid),.groups="drop")%>%
  filter(models==2,all_models_valid)
tv<-tv_all%>%semi_join(valid_keys,by=c("imputation","sex","SEQN"))
tv$height_quintile<-NA_character_
for(s in c("Female","Male")){ref<-tv%>%filter(imputation==1,sex==s,model=="Hcommon");cuts<-wq(ref$height_cm,ref$wtmec_pool,seq(0,1,.2));cuts[1]<--Inf;cuts[6]<-Inf;tv$height_quintile[tv$sex==s]<-as.character(cut(tv$height_cm[tv$sex==s],cuts,labels=paste0("Q",1:5),include.lowest=TRUE))}

tv_over_single<-map_dfr(1:5,function(imp)map_dfr(c("Female","Male"),function(s)map_dfr(c("Hcommon","HB"),function(mo){d<-tv%>%filter(imputation==imp,sex==s,model==mo);des<-make_design(d);mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE);tibble(imputation=imp,sex=s,model=mo,metric=c("mean_z","P5","P10"),estimate=as.numeric(coef(mn)),standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z"),n=nrow(d))})))
tv_over<-tv_over_single%>%pool_table(c("sex","model","metric"))%>%left_join(tv_over_single%>%group_by(sex,model,metric)%>%summarise(z_sd=mean(z_sd),unweighted_n=mean(n),.groups="drop"),by=c("sex","model","metric"))%>%mutate(percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),record_type="temporal_overall")

tv_eff_single<-list();ec<-1L
for(imp in 1:5)for(s in c("Female","Male"))for(mo in c("Hcommon","HB")){
  d<-tv%>%filter(imputation==imp,sex==s,model==mo)%>%droplevels();des<-make_design(d)
  for(outcome in c("z","low_p5","low_p10"))for(exposure in c("height","bmi")){
    x<-if(exposure=="height")"I(height_cm/10)"else"I(bmi/5)";fam<-if(outcome=="z")gaussian()else quasibinomial();fit<-svyglm(as.formula(paste(outcome,"~",x,"+ns(age,df=3)+race+cycle")),design=des,family=fam)
    tv_eff_single[[ec]]<-tibble(imputation=imp,sex=s,model=mo,outcome=outcome,exposure=exposure,unit=ifelse(exposure=="height","10 cm","5 kg/m2"),estimate=coef(fit)[2],standard_error=SE(fit)[2]);ec<-ec+1L}}
tv_eff<-bind_rows(tv_eff_single)%>%pool_table(c("sex","model","outcome","exposure","unit"))%>%mutate(effect_type=if_else(outcome=="z","beta","odds_ratio"),odds_ratio=if_else(outcome=="z",NA_real_,exp(estimate)),or_low=if_else(outcome=="z",NA_real_,exp(conf_low)),or_high=if_else(outcome=="z",NA_real_,exp(conf_high)),record_type="temporal_effect")

tv_strat_single<-list();sc<-1L
for(imp in 1:5)for(s in c("Female","Male"))for(mo in c("Hcommon","HB"))for(strat in c("height_quintile","bmi_group")){
  groups<-if(strat=="height_quintile")paste0("Q",1:5)else levels(tv$bmi_group)
  for(g in groups){d<-tv%>%filter(imputation==imp,sex==s,model==mo,as.character(.data[[strat]])==g);des<-make_design(d);mn<-svymean(~z+low_p5+low_p10,des,na.rm=TRUE)
    tv_strat_single[[sc]]<-tibble(imputation=imp,sex=s,model=mo,stratifier=strat,group=g,metric=c("mean_z","P5","P10"),estimate=as.numeric(coef(mn)),standard_error=as.numeric(SE(mn)),z_sd=safe_sd(des,"z"),n=nrow(d));sc<-sc+1L}}
tv_strat<-bind_rows(tv_strat_single)%>%pool_table(c("sex","model","stratifier","group","metric"))%>%left_join(bind_rows(tv_strat_single)%>%group_by(sex,model,stratifier,group,metric)%>%summarise(z_sd=mean(z_sd),unweighted_n=mean(n),.groups="drop"),by=c("sex","model","stratifier","group","metric"))%>%mutate(percent=if_else(metric%in%c("P5","P10"),100*estimate,NA_real_),record_type="temporal_stratified")

# Adjusted temporal BMI spline curves at age 45, median height, first cycle.
spline_single<-list();sp<-1L
for(imp in 1:5)for(s in c("Female","Male"))for(mo in c("Hcommon","HB")){
  d<-tv%>%filter(imputation==imp,sex==s,model==mo)%>%droplevels()
  bb<-ns(d$bmi,df=3);ab<-ns(d$age,df=3);hb<-ns(d$height_m,df=3)
  d$bmi_s1<-bb[,1];d$bmi_s2<-bb[,2];d$bmi_s3<-bb[,3];d$age_s1<-ab[,1];d$age_s2<-ab[,2];d$age_s3<-ab[,3]
  d$height_s1<-hb[,1];d$height_s2<-hb[,2];d$height_s3<-hb[,3]
  fit<-svyglm(z~bmi_s1+bmi_s2+bmi_s3+age_s1+age_s2+age_s3+height_s1+height_s2+height_s3+cycle,design=make_design(d))
  bmi_grid<-seq(15,50,1);age_grid<-rep(45,length(bmi_grid));height_grid<-rep(weighted.mean(d$height_m,d$wtmec_pool),length(bmi_grid))
  bp<-predict(bb,bmi_grid);ap<-predict(ab,age_grid);hp<-predict(hb,height_grid)
  grid<-data.frame(bmi_s1=bp[,1],bmi_s2=bp[,2],bmi_s3=bp[,3],age_s1=ap[,1],age_s2=ap[,2],age_s3=ap[,3],
    height_s1=hp[,1],height_s2=hp[,2],height_s3=hp[,3],cycle=factor(levels(d$cycle)[1],levels=levels(d$cycle)))
  pr<-predict(fit,newdata=grid,se.fit=TRUE);spline_single[[sp]]<-tibble(imputation=imp,sex=s,model=mo,bmi=bmi_grid,estimate=as.numeric(pr),standard_error=as.numeric(SE(pr)));sp<-sp+1L}
spline_curve<-bind_rows(spline_single)%>%pool_table(c("sex","model","bmi"));write_csv(spline_curve,file.path(cache_dir,"HB_temporal_BMI_spline.csv"))

temporal_out<-bind_rows(tv_over,tv_eff,tv_strat)%>%mutate(weight_version="Corrected development weights; frozen application")
write_csv(temporal_out,file.path(out_dir,"15_H_vs_HB_temporal.csv"))
saveRDS(tv,file.path(cache_dir,"H_HB_temporal_scores.rds"))

png(file.path(fig_dir,"HB_temporal_adjusted_BMI_spline.png"),width=3000,height=1500,res=250)
par(mfrow=c(1,2),mar=c(5,5,3,1));for(s in c("Female","Male")){plot(NA,xlim=c(15,50),ylim=range(spline_curve$conf_low,spline_curve$conf_high,na.rm=TRUE),xlab="BMI (kg/m2)",ylab="Adjusted temporal z",main=s);abline(h=0,lty=2,col="grey50");for(mo in c("Hcommon","HB")){d<-spline_curve%>%filter(sex==s,model==mo);col<-if(mo=="Hcommon")"#C94C4C"else"#0072B2";lines(d$bmi,d$estimate,col=col,lwd=3)};legend("topleft",c("H-common","HB"),col=c("#C94C4C","#0072B2"),lwd=3,bty="n")};dev.off()
cat("Corrected H/HB development and frozen temporal validation completed\n")
