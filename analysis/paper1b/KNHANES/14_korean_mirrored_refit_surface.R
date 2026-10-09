# -*- coding: UTF-8 -*-
options(stringsAsFactors=FALSE, survey.lonely.psu="adjust")
invisible(Sys.setlocale("LC_ALL","English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy","C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({library(gamlss);library(gamlss.dist);library(survey);library(ggplot2)})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
out <- file.path(root,"us_korea_surface"); fig_dir <- file.path(out,"figures"); model_dir <- file.path(out,"models")
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE);dir.create(model_dir,recursive=TRUE,showWarnings=FALSE)
paper_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/First Paper/Final Paper/\u65b0\u5efa\u6587\u4ef6\u5939"
scorer_root <- file.path(paper_root,"nhanes-alm-reference-scorer-v2.1.0")
source(file.path(scorer_root,"R","score_conditional_alm.R"),local=.GlobalEnv)
check_alm_reference_environment(strict=TRUE)
us_bundles <- load_alm_reference(file.path(scorer_root,"models",c(
  "corrected_H_18_69_Female_bundles.rds","corrected_H_18_69_Male_bundles.rds")))

dat <- readRDS(file.path(root,"data","knhanes_2008_2011_analysis_19_69_official_weight.rds"))
dat <- dat[is.finite(dat$age)&dat$age>=19&dat$age<=69&is.finite(dat$height_m)&is.finite(dat$alm_kg)&
  is.finite(dat$wt_pool)&dat$wt_pool>0&!is.na(dat$psu_pool)&!is.na(dat$strata_pool),]
dat$height10 <- dat$height_cm/10

weighted_quantile <- function(x,w,p){ok<-is.finite(x)&is.finite(w)&w>0;x<-x[ok];w<-w[ok];o<-order(x);x<-x[o];w<-w[o];x[pmin(length(x),findInterval(p,cumsum(w)/sum(w))+1L)]}
safe_prop <- function(des,var){v<-des$variables[[var]];if(all(!v))return(c(est=0,lcl=0,ucl=0));if(all(v))return(c(est=1,lcl=1,ucl=1));
  a<-svyciprop(as.formula(paste0("~",var)),des,method="logit",na.rm=TRUE);c(est=coef(a)[1],lcl=confint(a)[1],ucl=confint(a)[2])}

fit_one <- function(train,sex_value){
  family_name <- if(sex_value=="Female")"BCCG" else "BCT"
  family_function <- get(family_name,envir=asNamespace("gamlss.dist"))
  formula_text <- if(sex_value=="Female")
    "alm_kg ~ pb(age, df=3, inter=10) + pb(height_m, df=3, inter=10)" else
    "alm_kg ~ pb(age, df=3) + pb(height_m, df=3)"
  train$w_model <- train$wt_pool/mean(train$wt_pool)
  lm_start <- lm(alm_kg~age+height_m,data=train,weights=w_model)
  start <- pmax(as.numeric(predict(lm_start,train)),.5)
  attempts <- list(
    list(step=1,start=start),list(step=.5,start=start),list(step=.25,start=start),
    list(step=.1,start=start),list(step=.5,start=NULL),list(step=.25,start=NULL))
  last <- NULL
  for(i in seq_along(attempts)){
    a<-attempts[[i]]
    fit<-try(gamlss(as.formula(formula_text),sigma.formula=~1,nu.formula=~1,tau.formula=~1,
      family=family_function(),data=train,weights=w_model,mu.start=a$start,
      control=gamlss.control(n.cyc=500,trace=FALSE,mu.step=a$step)),silent=TRUE)
    last<-fit
    if(!inherits(fit,"try-error")&&isTRUE(fit$converged)){
      fit$call$formula <- as.formula(formula_text)
      fit$call$family <- as.call(list(as.name(family_name)))
      fit$call$data <- quote(train_data)
      return(list(fit=fit,train=train,attempt=i,family=family_name,formula=formula_text))
    }
  }
  stop("Korean mirrored fit failed for ",sex_value,": ",substr(as.character(last),1,500))
}

local <- list(); fit_audit <- list(); calibration <- list(); slope_rows <- list()
for(sex_value in c("Female","Male")){
  train <- dat[dat$sex_label==sex_value,]
  ans <- fit_one(train,sex_value);fit<-ans$fit;train<-ans$train
  bundle <- list(imputation=1L,sex=sex_value,family=ans$family,age_df=3,height_df=3,
    sigma_spec="constant",age_domain=c(19,69),train_data=train,model=fit)
  local[[sex_value]] <- bundle
  saveRDS(bundle,file.path(model_dir,paste0("Korean_mirrored_H_19_69_",sex_value,".rds")),compress="xz")
  pars <- predict_reference_parameters(bundle,train[,c("age","height_m")])
  p <- pmin(pmax(reference_cdf(bundle,train$alm_kg,pars),1e-10),1-1e-10);train$z_local<-qnorm(p)
  train$low_local_p5<-train$z_local<qnorm(.05);train$low_local_p10<-train$z_local<qnorm(.10)
  des<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wt_pool,nest=TRUE,data=train)
  mz<-svymean(~z_local,des);vz<-svyvar(~z_local,des);p5<-safe_prop(des,"low_local_p5");p10<-safe_prop(des,"low_local_p10")
  zc<-train$z_local;ww<-train$wt_pool;zm<-weighted.mean(zc,ww);zsd<-sqrt(weighted.mean((zc-zm)^2,ww))
  skew<-weighted.mean(((zc-zm)/zsd)^3,ww);kurt<-weighted.mean(((zc-zm)/zsd)^4,ww)-3
  calibration[[sex_value]]<-data.frame(sex=sex_value,n=nrow(train),weighted_mean_z=coef(mz)[1],mean_lcl=confint(mz)[1],mean_ucl=confint(mz)[2],
    weighted_sd_z=sqrt(as.numeric(vz)[1]),p5_percent=100*p5[1],p5_lcl=100*p5[2],p5_ucl=100*p5[3],
    p10_percent=100*p10[1],p10_lcl=100*p10[2],p10_ucl=100*p10[3],weighted_skewness=skew,weighted_excess_kurtosis=kurt)
  sf<-svyglm(z_local~height10+age,design=des);cf<-coef(sf);V<-vcov(sf)
  for(term in c("height10","age")){b<-cf[term];se<-sqrt(V[term,term]);slope_rows[[length(slope_rows)+1L]]<-data.frame(sex=sex_value,term=term,
    estimate=b,se=se,lcl=b-qnorm(.975)*se,ucl=b+qnorm(.975)*se,p_value=2*pnorm(abs(b/se),lower.tail=FALSE))}
  fit_audit[[sex_value]]<-data.frame(sex=sex_value,family=ans$family,mu_formula=ans$formula,sigma_formula="~1",nu_formula="~1",tau_formula="~1",
    n=nrow(train),weight_mean=mean(train$w_model),converged=fit$converged,iterations=fit$iter,effective_df=fit$df.fit,
    global_deviance=fit$G.deviance,AIC=AIC(fit),BIC=GAIC(fit,k=log(nrow(train))),fit_attempt=ans$attempt)
}
write.csv(do.call(rbind,fit_audit),file.path(out,"01_local_model_fit_audit.csv"),row.names=FALSE)
write.csv(do.call(rbind,calibration),file.path(out,"02_local_z_calibration.csv"),row.names=FALSE)
write.csv(do.call(rbind,slope_rows),file.path(out,"03_local_z_residual_slopes.csv"),row.names=FALSE)

# Primary domain: overlap of weighted 1st-99th percentile height ranges.
# Sensitivity: overlap of weighted 5th-95th ranges and frozen recommended-use limits.
range_rows<-list();domains<-list();rep_heights<-list()
for(sex_value in c("Female","Male")){
  ub<-us_bundles[[paste0("imp1_",sex_value)]];u<-ub$train_data;u<-u[u$age>=19&u$age<=69,]
  k<-local[[sex_value]]$train_data
  probs<-c(.01,.05,.25,.50,.75,.95,.99)
  uq<-weighted_quantile(u$height_m,u$w_model,probs);kq<-weighted_quantile(k$height_m,k$wt_pool,probs)
  rec<-alm_reference_domain[alm_reference_domain$sex==sex_value,]
  primary_lo<-ceiling(max(uq[1],kq[1])*100)/100;primary_hi<-floor(min(uq[7],kq[7])*100)/100
  strict_lo<-ceiling(max(uq[2],kq[2],rec$height_recommended_min_m)*100)/100
  strict_hi<-floor(min(uq[6],kq[6],rec$height_recommended_max_m)*100)/100
  domains[[sex_value]]<-list(primary=c(primary_lo,primary_hi),strict=c(strict_lo,strict_hi))
  # Equal population weight for representative shared-domain height percentiles.
  uv<-u$height_m[u$height_m>=primary_lo&u$height_m<=primary_hi];uw<-u$w_model[u$height_m>=primary_lo&u$height_m<=primary_hi];uw<-uw/sum(uw)*.5
  kv<-k$height_m[k$height_m>=primary_lo&k$height_m<=primary_hi];kw<-k$wt_pool[k$height_m>=primary_lo&k$height_m<=primary_hi];kw<-kw/sum(kw)*.5
  rh<-weighted_quantile(c(uv,kv),c(uw,kw),c(.25,.50,.75));rep_heights[[sex_value]]<-rh
  range_rows[[sex_value]]<-data.frame(sex=sex_value,probability=probs,us_height_m=uq,korea_height_m=kq,
    primary_common_min_m=primary_lo,primary_common_max_m=primary_hi,strict_common_min_m=strict_lo,strict_common_max_m=strict_hi,
    representative_p25_m=rh[1],representative_p50_m=rh[2],representative_p75_m=rh[3])
}
write.csv(do.call(rbind,range_rows),file.path(out,"04_common_height_domains.csv"),row.names=FALSE)

predict_local_centiles <- function(bundle,newdata,probs){
  pars<-predict_reference_parameters(bundle,newdata[,c("age","height_m")]);out<-newdata
  for(p in probs){nm<-paste0("alm_p",round(100*p));out[[nm]]<-reference_quantile(bundle,p,pars)}
  out
}
probs<-c(.05,.10,.50,.90,.95);surface_rows<-list()
for(sex_value in c("Female","Male")){
  dm<-domains[[sex_value]]$primary
  g<-expand.grid(sex=sex_value,age=19:69,height_m=seq(dm[1],dm[2],by=.01),KEEP.OUT.ATTRS=FALSE)
  us<-reference_centiles(g,probabilities=probs,bundles=us_bundles,warn=FALSE)
  kr<-predict_local_centiles(local[[sex_value]],g,probs)
  outg<-g;outg$height_cm<-100*outg$height_m
  for(pp in c(5,10,50,90,95)){
    outg[[paste0("us_p",pp)]]<-us[[paste0("alm_p",pp)]];outg[[paste0("korea_p",pp)]]<-kr[[paste0("alm_p",pp)]]
    outg[[paste0("delta_p",pp,"_kg")]]<-outg[[paste0("korea_p",pp)]]-outg[[paste0("us_p",pp)]]
    outg[[paste0("delta_p",pp,"_pct")]]<-100*outg[[paste0("delta_p",pp,"_kg")]]/outg[[paste0("us_p",pp)]]
  }
  outg$us_width_p95_p5<-outg$us_p95-outg$us_p5;outg$korea_width_p95_p5<-outg$korea_p95-outg$korea_p5
  outg$width_p95_p5_ratio<-outg$korea_width_p95_p5/outg$us_width_p95_p5
  outg$us_width_p90_p10<-outg$us_p90-outg$us_p10;outg$korea_width_p90_p10<-outg$korea_p90-outg$korea_p10
  outg$width_p90_p10_ratio<-outg$korea_width_p90_p10/outg$us_width_p90_p10
  surface_rows[[sex_value]]<-outg
}
surface<-do.call(rbind,surface_rows);write.csv(surface,file.path(out,"05_common_grid_centiles_and_differences.csv"),row.names=FALSE)

summarize_surface <- function(z,sex_value,domain_name){
  d<-z$delta_p50_kg;dp<-z$delta_p50_pct;r90<-z$width_p95_p5_ratio;r80<-z$width_p90_p10_ratio
  lmfit<-lm(korea_p50~us_p50,data=z);resid_offset<-d-mean(d);trend<-lm(d~I((age-19)/10)+I((height_cm-mean(height_cm))/10),data=z)
  tc<-summary(trend)$coefficients
  data.frame(sex=sex_value,domain=domain_name,n_grid=nrow(z),delta_p50_mean_kg=mean(d),delta_p50_sd_kg=sd(d),
    delta_p50_q25_kg=quantile(d,.25),delta_p50_median_kg=median(d),delta_p50_q75_kg=quantile(d,.75),delta_p50_min_kg=min(d),delta_p50_max_kg=max(d),
    delta_p50_mean_pct=mean(dp),width_p95_p5_ratio_mean=mean(r90),width_p95_p5_ratio_median=median(r90),
    width_p90_p10_ratio_mean=mean(r80),width_p90_p10_ratio_median=median(r80),
    p50_correlation=cor(z$us_p50,z$korea_p50),p50_r_squared=cor(z$us_p50,z$korea_p50)^2,
    korea_on_us_intercept=coef(lmfit)[1],korea_on_us_slope=coef(lmfit)[2],korea_on_us_r_squared=summary(lmfit)$r.squared,
    offset_adjusted_rmse_kg=sqrt(mean(resid_offset^2)),offset_adjusted_max_abs_kg=max(abs(resid_offset)),
    delta_trend_per_10y=tc[2,1],delta_trend_10y_lcl=tc[2,1]-qnorm(.975)*tc[2,2],delta_trend_10y_ucl=tc[2,1]+qnorm(.975)*tc[2,2],
    delta_trend_per_10cm=tc[3,1],delta_trend_10cm_lcl=tc[3,1]-qnorm(.975)*tc[3,2],delta_trend_10cm_ucl=tc[3,1]+qnorm(.975)*tc[3,2])
}
summary_rows<-list()
for(sex_value in c("Female","Male")){
  z<-surface[surface$sex==sex_value,];summary_rows[[length(summary_rows)+1L]]<-summarize_surface(z,sex_value,"Primary common 1st-99th")
  ds<-domains[[sex_value]]$strict;zs<-z[z$height_m>=ds[1]&z$height_m<=ds[2],]
  summary_rows[[length(summary_rows)+1L]]<-summarize_surface(zs,sex_value,"Strict recommended/common 5th-95th")
}
surface_summary<-do.call(rbind,summary_rows);write.csv(surface_summary,file.path(out,"06_location_scale_geometry_summary.csv"),row.names=FALSE)

surface$age_band<-cut(surface$age,c(19,30,40,50,60,70),right=FALSE,labels=c("19-29","30-39","40-49","50-59","60-69"))
age_rows<-list()
for(sex_value in c("Female","Male"))for(band in levels(surface$age_band)){
  z<-surface[surface$sex==sex_value&surface$age_band==band,]
  age_rows[[length(age_rows)+1L]]<-data.frame(sex=sex_value,age_band=band,n_grid=nrow(z),
    delta_p50_mean_kg=mean(z$delta_p50_kg),delta_p50_mean_pct=mean(z$delta_p50_pct),
    delta_p5_mean_kg=mean(z$delta_p5_kg),delta_p5_mean_pct=mean(z$delta_p5_pct),
    delta_p10_mean_kg=mean(z$delta_p10_kg),delta_p10_mean_pct=mean(z$delta_p10_pct),
    width_p95_p5_ratio_mean=mean(z$width_p95_p5_ratio),width_p90_p10_ratio_mean=mean(z$width_p90_p10_ratio))
}
write.csv(do.call(rbind,age_rows),file.path(out,"07_age_pattern_summary.csv"),row.names=FALSE)

# Representative age curves.
curve_rows<-list()
for(sex_value in c("Female","Male")){
  h<-rep_heights[[sex_value]]
  g<-expand.grid(sex=sex_value,age=19:69,height_m=h,KEEP.OUT.ATTRS=FALSE)
  height_labels<-sprintf(c("P25: %.1f cm","P50: %.1f cm","P75: %.1f cm"),100*h)
  g$height_level=factor(rep(height_labels,each=51),levels=height_labels)
  us<-reference_centiles(g,probabilities=c(.05,.50),bundles=us_bundles,warn=FALSE);kr<-predict_local_centiles(local[[sex_value]],g,c(.05,.50))
  for(i in seq_len(nrow(g)))for(pp in c(5,50))for(ref in c("U.S.","Korea")){
    val<-if(ref=="U.S.")us[[paste0("alm_p",pp)]][i] else kr[[paste0("alm_p",pp)]][i]
    curve_rows[[length(curve_rows)+1L]]<-data.frame(sex=sex_value,age=g$age[i],height_m=g$height_m[i],height_level=g$height_level[i],
      centile=paste0("P",pp),reference=ref,alm_kg=val)
  }
}
curves<-do.call(rbind,curve_rows);write.csv(curves,file.path(out,"08_representative_age_curves.csv"),row.names=FALSE)
for(sex_value in c("Female","Male")){
  p<-ggplot(curves[curves$sex==sex_value,],aes(age,alm_kg,color=reference,linetype=height_level))+geom_line(linewidth=.9)+
    facet_wrap(~centile,ncol=1,scales="free_y")+scale_color_manual(values=c("U.S."="#2166AC","Korea"="#B2182B"))+
    labs(title=paste0(sex_value,": U.S. and Korean conditional ALM curves"),x="Age (years)",y="Conditional ALM (kg)",color="Reference",linetype="Shared height")+
    theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir,paste0("Figure_representative_curves_",tolower(sex_value),".png")),p,width=7.2,height=7,dpi=300)
  ggsave(file.path(fig_dir,paste0("Figure_representative_curves_",tolower(sex_value),".pdf")),p,width=7.2,height=7)
}

# P50 difference heatmaps (kg primary; percent supplementary).
for(sex_value in c("Female","Male")){
  z<-surface[surface$sex==sex_value,]
  for(metric in c("delta_p50_kg","delta_p50_pct")){
    unit<-if(metric=="delta_p50_kg")"kg" else "%"
    p<-ggplot(z,aes(height_cm,age,fill=.data[[metric]]))+geom_tile()+
      scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0)+
      labs(title=paste0(sex_value,": Korean minus U.S. conditional median"),subtitle=paste0("Difference in ",unit),
        x="Height (cm)",y="Age (years)",fill=paste0("Korea−U.S.\n(",unit,")"))+theme_bw(base_size=11)+theme(plot.title=element_text(face="bold"))
    stem<-paste0("Figure_delta_P50_heatmap_",tolower(sex_value),if(metric=="delta_p50_pct")"_percent" else "_kg")
    ggsave(file.path(fig_dir,paste0(stem,".png")),p,width=7.2,height=5.5,dpi=300)
    ggsave(file.path(fig_dir,paste0(stem,".pdf")),p,width=7.2,height=5.5)
  }
}

cat("Korean mirrored refit and U.S.-Korea surface comparison complete.\n")
print(do.call(rbind,fit_audit),row.names=FALSE);print(do.call(rbind,calibration),row.names=FALSE)
print(surface_summary,row.names=FALSE)
