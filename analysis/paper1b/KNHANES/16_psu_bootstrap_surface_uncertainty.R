# -*- coding: UTF-8 -*-
options(stringsAsFactors=FALSE, survey.lonely.psu="adjust")
invisible(Sys.setlocale("LC_ALL","English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy","C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({library(gamlss);library(gamlss.dist);library(survey);library(parallel);library(ggplot2)})

args <- commandArgs(trailingOnly=TRUE)
B <- if(length(args)>=1)as.integer(args[1]) else 500L
n_workers <- if(length(args)>=2)as.integer(args[2]) else min(4L,max(1L,detectCores()-1L))
seed <- if(length(args)>=3)as.integer(args[3]) else 260905L
stopifnot(is.finite(B),B>=20L,is.finite(n_workers),n_workers>=1L)

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
surface_dir <- file.path(root,"us_korea_surface")
out <- file.path(surface_dir,"bootstrap_uncertainty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
paper_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/First Paper/Final Paper/\u65b0\u5efa\u6587\u4ef6\u5939"
scorer_root <- file.path(paper_root,"nhanes-alm-reference-scorer-v2.1.0")
source(file.path(scorer_root,"R","score_conditional_alm.R"),local=.GlobalEnv)
check_alm_reference_environment(strict=TRUE)

dat <- readRDS(file.path(root,"data","knhanes_2008_2011_analysis_19_69_official_weight.rds"))
dat <- dat[is.finite(dat$age)&dat$age>=19&dat$age<=69&is.finite(dat$height_m)&is.finite(dat$alm_kg)&
  is.finite(dat$wt_pool)&dat$wt_pool>0&!is.na(dat$psu_pool)&!is.na(dat$strata_pool),]
dat <- droplevels(dat)

# Fixed U.S. predictions and common grid from the locked point-estimate analysis.
fixed <- read.csv(file.path(surface_dir,"05_common_grid_centiles_and_differences.csv"),check.names=FALSE)
domains <- unique(read.csv(file.path(surface_dir,"04_common_height_domains.csv"),check.names=FALSE)[,
  c("sex","primary_common_min_m","primary_common_max_m","strict_common_min_m","strict_common_max_m")])
point <- read.csv(file.path(surface_dir,"06_location_scale_geometry_summary.csv"),check.names=FALSE)
local <- list(
  Female=readRDS(file.path(surface_dir,"models","Korean_mirrored_H_19_69_Female.rds")),
  Male=readRDS(file.path(surface_dir,"models","Korean_mirrored_H_19_69_Male.rds")))
sex_index <- lapply(c("Female","Male"),function(s)which(dat$sex_label==s));names(sex_index)<-c("Female","Male")
fixed_grid <- lapply(c("Female","Male"),function(s)fixed[fixed$sex==s,c("sex","age","height_m","height_cm","us_p5","us_p10","us_p50","us_p90","us_p95")])
names(fixed_grid)<-c("Female","Male")
curve_source<-read.csv(file.path(surface_dir,"08_representative_age_curves.csv"),check.names=FALSE)
curve_us<-curve_source[curve_source$centile=="P50"&curve_source$reference=="U.S.",c("sex","age","height_m","height_level","alm_kg")]
names(curve_us)[5]<-"us_p50"
curve_kr<-curve_source[curve_source$centile=="P50"&curve_source$reference=="Korea",c("sex","age","height_m","height_level","alm_kg")]
names(curve_kr)[5]<-"korea_point_p50"
curve_source<-merge(curve_us,curve_kr,by=c("sex","age","height_m","height_level"),sort=FALSE)
curve_source$height_cm<-100*curve_source$height_m
curve_grid<-split(curve_source,curve_source$sex)
strict_ranges <- setNames(lapply(c("Female","Male"),function(s){z<-domains[domains$sex==s,];c(z$strict_common_min_m,z$strict_common_max_m)}),c("Female","Male"))
start_mu <- lapply(c("Female","Male"),function(s){
  idx<-sex_index[[s]];pars<-predict_reference_parameters(local[[s]],dat[idx,c("age","height_m")]);as.numeric(pars$mu)
});names(start_mu)<-c("Female","Male")
start_const <- lapply(c("Female","Male"),function(s){
  fit<-local[[s]]$model
  c(sigma=unname(fit$sigma.fv[1]),nu=unname(fit$nu.fv[1]),tau=if(s=="Male")unname(fit$tau.fv[1]) else NA_real_)
});names(start_const)<-c("Female","Male")

# Generate survey-package stratified PSU bootstrap replicate weights once. The same
# replicate is used for women and men, preserving the original pooled survey structure.
design <- svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wt_pool,nest=TRUE,data=dat)
set.seed(seed)
rep_design <- as.svrepdesign(design,type="bootstrap",replicates=B,mse=TRUE)
rep_weights <- weights(rep_design,type="analysis")
stopifnot(identical(dim(rep_weights),c(nrow(dat),B)))

fit_boot <- function(sex_value,w_all){
  idx<-sex_index[[sex_value]];w<-w_all[idx];keep<-is.finite(w)&w>0
  train_data<-dat[idx[keep],];train_data$w_model<-w[keep]/mean(w[keep])
  family_name<-if(sex_value=="Female")"BCCG" else "BCT"
  family_function<-get(family_name,envir=asNamespace("gamlss.dist"))
  formula_text<-if(sex_value=="Female")
    "alm_kg ~ pb(age, df=3, inter=10) + pb(height_m, df=3, inter=10)" else
    "alm_kg ~ pb(age, df=3) + pb(height_m, df=3)"
  package_fit<-function(fit,attempt,converged,algorithm){
    fit$call$formula<-as.formula(formula_text);fit$call$family<-as.call(list(as.name(family_name)));fit$call$data<-quote(train_data)
    bundle<-list(imputation=1L,sex=sex_value,family=family_name,age_df=3,height_df=3,sigma_spec="constant",
      age_domain=c(19,69),train_data=train_data,model=fit)
    list(bundle=bundle,attempt=attempt,n=nrow(train_data),iterations=fit$iter,converged=converged,algorithm=algorithm,
      deviance=fit$G.deviance,sigma=unname(fit$sigma.fv[1]),nu=unname(fit$nu.fv[1]),
      tau=if(sex_value=="Male")unname(fit$tau.fv[1]) else NA_real_)
  }
  accept_fit<-function(fit,attempt,algorithm){
    if(!inherits(fit,"try-error")&&isTRUE(fit$converged)){
      return(package_fit(fit,attempt,TRUE,algorithm))
    }
    NULL
  }

  # Primary optimizer: the original RS algorithm used for the locked model.
  first_args<-list(formula=as.formula(formula_text),sigma.formula=~1,nu.formula=~1,tau.formula=~1,
    family=family_function(),data=train_data,weights=train_data$w_model,
    mu.start=start_mu[[sex_value]][keep],sigma.start=rep(start_const[[sex_value]]["sigma"],nrow(train_data)),
    nu.start=rep(start_const[[sex_value]]["nu"],nrow(train_data)),
    control=gamlss.control(n.cyc=30L,trace=FALSE))
  if(sex_value=="Male")first_args$tau.start<-rep(start_const[[sex_value]]["tau"],nrow(train_data))
  fit<-try(do.call(gamlss,first_args),silent=TRUE)
  accepted<-accept_fit(fit,1L,"RS");if(!is.null(accepted))return(accepted)
  if(inherits(fit,"try-error"))stop(substr(as.character(fit),1,300))

  # Numerical rescue for the same likelihood/model: continue from the RS state,
  # then alternate one RS cycle with CG cycles while damping all parameter updates.
  fit2<-try(gamlss(as.formula(formula_text),sigma.formula=~1,nu.formula=~1,tau.formula=~1,
    family=family_function(),data=train_data,weights=w_model,start.from=fit,method=mixed(1,40),
    control=gamlss.control(trace=FALSE,c.crit=.01,mu.step=.5,sigma.step=.5,nu.step=.5,tau.step=.5)),silent=TRUE)
  accepted<-accept_fit(fit2,2L,"mixed_RS_CG_1_40");if(!is.null(accepted))return(accepted)
  last<-if(inherits(fit2,"try-error"))fit else fit2

  fit3<-try(gamlss(as.formula(formula_text),sigma.formula=~1,nu.formula=~1,tau.formula=~1,
    family=family_function(),data=train_data,weights=w_model,start.from=last,method=mixed(2,80),
    control=gamlss.control(trace=FALSE,c.crit=.01,mu.step=.25,sigma.step=.25,nu.step=.25,tau.step=.25)),silent=TRUE)
  accepted<-accept_fit(fit3,3L,"mixed_RS_CG_2_80");if(!is.null(accepted))return(accepted)
  last<-if(inherits(fit3,"try-error"))last else fit3
  package_fit(last,3L,FALSE,"mixed_RS_CG_2_80_nonconverged")
}

predict_direct <- function(bundle,g,probabilities){
  # Exact direct evaluation of the two fitted pb() functions. This is algebraically
  # identical to predict.gamlss(newdata=...), but avoids predict.gamlss's expensive
  # augmented-data refit for every bootstrap model and every distribution parameter.
  fit<-bundle$model
  mu_eta<-fit$mu.coefficients[1]+fit$mu.coefficients[2]*g$age+fit$mu.coefficients[3]*g$height_m+
    getSmo(fit,parameter="mu",which=1)$fun(g$age)+getSmo(fit,parameter="mu",which=2)$fun(g$height_m)
  pars<-list(mu=as.numeric(mu_eta),sigma=rep(fit$sigma.fv[1],nrow(g)),nu=rep(fit$nu.fv[1],nrow(g)))
  if(bundle$sex=="Male")pars$tau<-rep(fit$tau.fv[1],nrow(g))
  out<-lapply(probabilities,function(p)reference_quantile(bundle,p,pars));names(out)<-paste0("p",round(100*probabilities));out
}

summarize_boot <- function(bundle,sex_value,domain_name){
  g<-fixed_grid[[sex_value]]
  if(domain_name=="Strict recommended/common 5th-95th"){
    rr<-strict_ranges[[sex_value]];g<-g[g$height_m>=rr[1]&g$height_m<=rr[2],]
  }
  pred<-predict_direct(bundle,g,c(.05,.10,.50,.90,.95));kp5<-pred$p5;kp10<-pred$p10;kp50<-pred$p50;kp90<-pred$p90;kp95<-pred$p95
  delta<-kp50-g$us_p50;delta_pct<-100*delta/g$us_p50
  width_ratio<-(kp95-kp5)/(g$us_p95-g$us_p5)
  width_ratio_80<-(kp90-kp10)/(g$us_p90-g$us_p10)
  fit_geometry<-lm(kp50~g$us_p50)
  fit_trend<-lm(delta~I((age-19)/10)+I((height_cm-mean(height_cm))/10),data=g)
  data.frame(sex=sex_value,domain=domain_name,
    delta_p50_mean_kg=mean(delta),delta_p50_mean_pct=mean(delta_pct),
    width_p95_p5_ratio_mean=mean(width_ratio),width_p90_p10_ratio_mean=mean(width_ratio_80),p50_r_squared=cor(kp50,g$us_p50)^2,
    korea_on_us_slope=unname(coef(fit_geometry)[2]),offset_adjusted_rmse_kg=sqrt(mean((delta-mean(delta))^2)),
    delta_trend_per_10y=unname(coef(fit_trend)[2]),delta_trend_per_10cm=unname(coef(fit_trend)[3]))
}

curve_boot <- function(bundle,sex_value){
  g<-curve_grid[[sex_value]]
  p50<-predict_direct(bundle,g,.50)$p50
  data.frame(sex=sex_value,age=g$age,height_m=g$height_m,height_cm=g$height_cm,height_level=g$height_level,
    us_p50=g$us_p50,korea_point_p50=g$korea_point_p50,korea_boot_p50=p50)
}

run_task <- function(task){
  i<-task$replicate;sex_value<-task$sex;metrics<-list();curve<-NULL;started<-proc.time()[3]
  ans<-try(fit_boot(sex_value,rep_weights[,i]),silent=TRUE)
  if(inherits(ans,"try-error")){
    diag<-data.frame(replicate=i,sex=sex_value,status="FAIL",attempt=NA_integer_,algorithm="ERROR",n_positive=NA_integer_,iterations=NA_integer_,
      deviance=NA_real_,sigma=NA_real_,nu=NA_real_,tau=NA_real_,elapsed_sec=proc.time()[3]-started,error=substr(as.character(ans),1,300))
  }else if(!isTRUE(ans$converged)){
    diag<-data.frame(replicate=i,sex=sex_value,status="FAIL",attempt=ans$attempt,algorithm=ans$algorithm,n_positive=ans$n,iterations=ans$iterations,
      deviance=ans$deviance,sigma=ans$sigma,nu=ans$nu,tau=ans$tau,elapsed_sec=proc.time()[3]-started,error="Did not satisfy GAMLSS convergence criterion")
  }else{
    for(dm in c("Primary common 1st-99th","Strict recommended/common 5th-95th")){
      sm<-try(summarize_boot(ans$bundle,sex_value,dm),silent=TRUE)
      if(!inherits(sm,"try-error")){sm$replicate<-i;metrics[[length(metrics)+1L]]<-sm}
    }
    curve<-try(curve_boot(ans$bundle,sex_value),silent=TRUE)
    if(!inherits(curve,"try-error"))curve$replicate<-i else curve<-NULL
    diag<-data.frame(replicate=i,sex=sex_value,status="SUCCESS",attempt=ans$attempt,algorithm=ans$algorithm,n_positive=ans$n,iterations=ans$iterations,
      deviance=ans$deviance,sigma=ans$sigma,nu=ans$nu,tau=ans$tau,elapsed_sec=proc.time()[3]-started,error="")
  }
  list(metrics=if(length(metrics))do.call(rbind,metrics) else NULL,curve=curve,diagnostics=diag)
}

cat(sprintf("Starting %d stratified PSU bootstrap replicates with %d workers.\n",B,n_workers))
cl<-makeCluster(n_workers,outfile="")
on.exit(try(stopCluster(cl),silent=TRUE),add=TRUE)
clusterEvalQ(cl,{
  .libPaths(c("C:/Users/Public/CodexRLib42Copy","C:/Program Files/R/R-4.2.1/library"))
  suppressPackageStartupMessages({library(gamlss);library(gamlss.dist)})
  NULL
})
clusterExport(cl,c("scorer_root"),envir=environment())
clusterEvalQ(cl,{invisible(source(file.path(scorer_root,"R","score_conditional_alm.R"),local=.GlobalEnv));NULL})
clusterExport(cl,c("dat","sex_index","fixed_grid","curve_grid","strict_ranges","start_mu","rep_weights",
  "start_const","fit_boot","predict_direct","summarize_boot","curve_boot","run_task"),envir=environment())

tasks<-lapply(seq_len(B),function(i)list(list(replicate=i,sex="Female"),list(replicate=i,sex="Male")))
tasks<-unlist(tasks,recursive=FALSE)
script_version<-"2026-09-07-v2";checkpoint_path<-file.path(out,"bootstrap_checkpoint.rds")
all_results<-vector("list",length(tasks));batch_size<-n_workers;start_task<-1L
if(file.exists(checkpoint_path)){
  cp<-try(readRDS(checkpoint_path),silent=TRUE)
  if(!inherits(cp,"try-error")&&is.list(cp)&&is.list(cp$meta)&&identical(cp$meta$B,B)&&identical(cp$meta$seed,seed)&&
      identical(cp$meta$script_version,script_version)&&length(cp$results)<=length(tasks)){
    all_results[seq_along(cp$results)]<-cp$results;start_task<-length(cp$results)+1L
    cat(sprintf("Resuming verified checkpoint at sex-specific fit %d/%d.\n",start_task,length(tasks)))
  }
}
if(start_task<=length(tasks))for(lo in seq.int(start_task,length(tasks),by=batch_size)){
  hi<-min(length(tasks),lo+batch_size-1L);ids<-lo:hi
  batch<-parLapplyLB(cl,tasks[ids],run_task)
  all_results[ids]<-batch
  diag_so_far<-do.call(rbind,lapply(all_results[seq_len(hi)],function(x)if(is.null(x))NULL else x$diagnostics))
  saveRDS(list(meta=list(B=B,seed=seed,script_version=script_version),results=all_results[seq_len(hi)]),checkpoint_path,compress=FALSE)
  cat(sprintf("Completed %d/%d sex-specific fits (%d/%d replicates); successful fits %d/%d.\n",hi,2L*B,floor(hi/2),B,sum(diag_so_far$status=="SUCCESS"),nrow(diag_so_far)))
  flush.console()
}
stopCluster(cl)

metrics<-do.call(rbind,lapply(all_results,function(x)x$metrics))
curve_reps<-do.call(rbind,lapply(all_results,function(x)x$curve))
diagnostics<-do.call(rbind,lapply(all_results,function(x)x$diagnostics))
metrics<-metrics[,c("replicate","sex","domain","delta_p50_mean_kg","delta_p50_mean_pct","width_p95_p5_ratio_mean","width_p90_p10_ratio_mean",
  "p50_r_squared","korea_on_us_slope","offset_adjusted_rmse_kg","delta_trend_per_10y","delta_trend_per_10cm")]
write.csv(metrics,file.path(out,"10_bootstrap_surface_metrics_replicates.csv"),row.names=FALSE)
write.csv(diagnostics,file.path(out,"12_bootstrap_fit_diagnostics.csv"),row.names=FALSE)
write.csv(curve_reps,file.path(out,"15_bootstrap_representative_p50_replicates.csv"),row.names=FALSE)

metric_names<-c("delta_p50_mean_kg","delta_p50_mean_pct","width_p95_p5_ratio_mean","width_p90_p10_ratio_mean","p50_r_squared",
  "korea_on_us_slope","offset_adjusted_rmse_kg","delta_trend_per_10y","delta_trend_per_10cm")
ci_rows<-list()
for(sex_value in c("Female","Male"))for(dm in c("Primary common 1st-99th","Strict recommended/common 5th-95th"))for(m in metric_names){
  x<-metrics[metrics$sex==sex_value&metrics$domain==dm,m];x<-x[is.finite(x)]
  pe<-point[point$sex==sex_value&point$domain==dm,m]
  qq<-quantile(x,c(.025,.975),names=FALSE,type=6)
  ci_rows[[length(ci_rows)+1L]]<-data.frame(sex=sex_value,domain=dm,metric=m,point_estimate=pe,
    bootstrap_replicates=length(x),bootstrap_mean=mean(x),bootstrap_se=sd(x),bootstrap_bias=mean(x)-pe,
    percentile_95_lcl=qq[1],percentile_95_ucl=qq[2])
}
ci<-do.call(rbind,ci_rows)
write.csv(ci,file.path(out,"11_bootstrap_surface_metrics_ci.csv"),row.names=FALSE)

band_rows<-lapply(split(curve_reps,interaction(curve_reps$sex,curve_reps$height_level,curve_reps$age,drop=TRUE)),function(z){
  qq<-quantile(z$korea_boot_p50,c(.025,.50,.975),names=FALSE,type=6)
  data.frame(sex=z$sex[1],age=z$age[1],height_m=z$height_m[1],height_cm=z$height_cm[1],height_level=z$height_level[1],
    us_p50=z$us_p50[1],korea_point_p50=z$korea_point_p50[1],bootstrap_replicates=nrow(z),
    korea_boot_median=qq[2],korea_p50_lcl=qq[1],korea_p50_ucl=qq[3])
})
bands<-do.call(rbind,band_rows);bands<-bands[order(bands$sex,bands$height_cm,bands$age),]
write.csv(bands,file.path(out,"16_bootstrap_representative_p50_bands.csv"),row.names=FALSE)
fig_dir<-file.path(out,"figures");dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
for(sex_value in c("Female","Male")){
  z<-bands[bands$sex==sex_value,]
  p<-ggplot(z,aes(age,korea_point_p50))+geom_ribbon(aes(ymin=korea_p50_lcl,ymax=korea_p50_ucl),fill="#D6604D",alpha=.24)+
    geom_line(aes(y=korea_point_p50,color="Korean local"),linewidth=.9)+geom_line(aes(y=us_p50,color="Frozen U.S."),linewidth=.9)+
    facet_wrap(~height_level,nrow=1,scales="free_y")+scale_color_manual(values=c("Korean local"="#B2182B","Frozen U.S."="#2166AC"))+
    labs(title=paste0(sex_value,": conditional median ALM with PSU-bootstrap band"),subtitle="Shading: Korean local pointwise 95% percentile bootstrap CI",
      x="Age (years)",y="Conditional P50 ALM (kg)",color=NULL)+theme_bw(base_size=11)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
  ggsave(file.path(fig_dir,paste0("Figure_bootstrap_P50_band_",tolower(sex_value),".png")),p,width=10,height=4.8,dpi=300)
  ggsave(file.path(fig_dir,paste0("Figure_bootstrap_P50_band_",tolower(sex_value),".pdf")),p,width=10,height=4.8)
}

# Confirm that the fast direct smoother evaluation reproduces every locked point
# estimate before accepting bootstrap intervals based on it.
direct_point<-do.call(rbind,lapply(c("Female","Male"),function(s)do.call(rbind,lapply(
  c("Primary common 1st-99th","Strict recommended/common 5th-95th"),function(dm)summarize_boot(local[[s]],s,dm)))))
direct_compare<-merge(direct_point,point[,c("sex","domain",metric_names)],by=c("sex","domain"),suffixes=c("_direct","_locked"))
direct_diffs<-unlist(lapply(metric_names,function(m)abs(direct_compare[[paste0(m,"_direct")]]-direct_compare[[paste0(m,"_locked")]])))
max_direct_difference<-max(direct_diffs)

psu_by_stratum<-aggregate(psu_pool~strata_pool,dat,function(x)length(unique(x)))
singleton<-psu_by_stratum$strata_pool[psu_by_stratum$psu_pool==1]
audit<-data.frame(
  item=c("script version","R version","survey version","gamlss version","gamlss.dist version","bootstrap type","optimizer","seed","requested replicates","workers",
    "analytic N","strata","PSUs","singleton strata","singleton rows","singleton weighted percent","successful female fits","successful male fits"),
  value=c(script_version,as.character(getRversion()),as.character(packageVersion("survey")),as.character(packageVersion("gamlss")),as.character(packageVersion("gamlss.dist")),
    "survey::as.svrepdesign(type='bootstrap'); PSU resampling within strata","RS c.crit=.001; nonconverged fits rescued by mixed RS/CG c.crit=.01",
    seed,B,n_workers,nrow(dat),nrow(psu_by_stratum),length(unique(dat$psu_pool)),
    length(singleton),sum(dat$strata_pool%in%singleton),100*sum(dat$wt_pool[dat$strata_pool%in%singleton])/sum(dat$wt_pool),
    sum(diagnostics$sex=="Female"&diagnostics$status=="SUCCESS"),sum(diagnostics$sex=="Male"&diagnostics$status=="SUCCESS")))
write.csv(audit,file.path(out,"13_bootstrap_design_and_run_audit.csv"),row.names=FALSE)

success_rate<-with(diagnostics,tapply(status=="SUCCESS",sex,mean))
curve_ok<-nrow(bands)==306L&&all(is.finite(bands$korea_p50_lcl))&&all(bands$korea_p50_lcl<=bands$korea_boot_median&bands$korea_boot_median<=bands$korea_p50_ucl)
verification<-data.frame(check=c("Requested replicate count","Fit success rate >=95%","All nine metrics finite","Every CI ordered and contains finite endpoints",
  "Point estimates match locked surface summary","Representative P50 bands complete"),
  status=c(if(nrow(diagnostics)==2*B)"PASS" else "FAIL",if(all(success_rate>=.95))"PASS" else "FAIL",
    if(all(vapply(metrics[,metric_names,drop=FALSE],function(x)all(is.finite(x)),logical(1))))"PASS" else "FAIL",if(all(is.finite(ci$percentile_95_lcl)&is.finite(ci$percentile_95_ucl)&ci$percentile_95_lcl<=ci$percentile_95_ucl))"PASS" else "FAIL",
    if(all(is.finite(ci$point_estimate))&&max_direct_difference<1e-8)"PASS" else "FAIL",if(curve_ok)"PASS" else "FAIL"),
  detail=c(paste0(nrow(diagnostics)/2,"/",B),paste(names(success_rate),sprintf("%.1f%%",100*success_rate),collapse="; "),
    paste0(nrow(metrics)," metric rows"),paste0(nrow(ci)," intervals"),paste0("maximum direct-vs-locked difference=",format(max_direct_difference,scientific=TRUE)),
    paste0(nrow(bands)," age-height band rows")))
write.csv(verification,file.path(out,"14_bootstrap_verification.csv"),row.names=FALSE)
print(ci,row.names=FALSE);print(verification,row.names=FALSE)
if(any(verification$status!="PASS"))stop("Bootstrap verification failed")
cat("PSU bootstrap surface uncertainty analysis complete.\n")
