options(stringsAsFactors=FALSE,survey.lonely.psu="adjust")
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(purrr);library(gamlss);library(gamlss.dist);library(survey)})
root_dir<-normalizePath(".",winslash="/",mustWork=TRUE);out_dir<-file.path(root_dir,"First Paper","Paper1_FINAL_FROZEN_20260814")
cache_dir<-file.path(out_dir,"cache");model_dir<-file.path(out_dir,"models")
validation<-readRDS(file.path(cache_dir,"validation_corrected_H_18_59.rds"))%>%mutate(strata_pool=factor(strata_pool),psu_pool=factor(psu_pool))
wq<-function(x,w,p){o<-order(x);x<-x[o];w<-w[o];vapply(p,function(q)x[which(cumsum(w)/sum(w)>=q)[1]],numeric(1))}
make_design<-function(d)svydesign(ids=~psu_pool,strata=~strata_pool,weights=~wtmec_pool,nest=TRUE,data=d)
rubin<-function(q,se){m<-length(q);u<-mean(se^2);b<-if(m>1)var(q)else 0;t<-u+(1+1/m)*b;s<-sqrt(t);df<-if(b>0)(m-1)*(1+u/((1+1/m)*b))^2 else Inf;cr<-qt(.975,df);c(estimate=mean(q),se=s,low=mean(q)-cr*s,high=mean(q)+cr*s)}

fit_anchor<-function(train,s,family_name){
  formula_text<-if(s=="Female")"alm_kg~pb(age,df=3,inter=10)+pb(height_m,df=3,inter=10)"else"alm_kg~pb(age,df=3)+pb(height_m,df=3)"
  fam<-get(family_name,asNamespace("gamlss.dist"));lmfit<-lm(alm_kg~age+height_m,data=train,weights=w_model);start<-pmax(predict(lmfit,train),.5)
  fit<-gamlss(as.formula(formula_text),sigma.formula=~1,nu.formula=~1,tau.formula=~1,family=fam(),data=train,weights=w_model,
    mu.start=start,control=gamlss.control(n.cyc=500,trace=FALSE));fit$call$family<-as.call(list(as.name(family_name)));fit$call$data<-quote(train);fit
}
predict_pars<-function(model,train,nd){x<-as.data.frame(nd[,intersect(names(nd),names(train)),drop=FALSE]);a<-try(suppressWarnings(predictAll(model,newdata=x,data=as.data.frame(train),type="response")),silent=TRUE);if(inherits(a,"try-error"))a<-suppressWarnings(predictAll(model,newdata=x,type="response"));as.data.frame(a)}
score<-function(model,train,family_name,nd){p<-predict_pars(model,train,nd);pf<-get(paste0("p",family_name),asNamespace("gamlss.dist"));a<-list(q=nd$alm_kg,mu=p$mu,sigma=p$sigma);if("nu"%in%names(p))a$nu<-p$nu;if("tau"%in%names(p))a$tau<-p$tau;pr<-pmin(pmax(as.numeric(do.call(pf,a)),1e-10),1-1e-10);tibble(z_anchor=qnorm(pr),p_anchor=pr)}
quant<-function(model,train,family_name,nd,prob){p<-predict_pars(model,train,nd);qf<-get(paste0("q",family_name),asNamespace("gamlss.dist"));a<-list(p=prob,mu=p$mu,sigma=p$sigma);if("nu"%in%names(p))a$nu<-p$nu;if("tau"%in%names(p))a$tau<-p$tau;as.numeric(do.call(qf,a))}

person_single<-list();curve_rows<-list();ps<-cr<-1L
for(s in c("Female","Male")){
  family_name<-if(s=="Female")"BCCG"else"BCT";bundles<-readRDS(file.path(model_dir,paste0("corrected_H_18_69_",s,"_bundles.rds")))
  for(imp in 1:5){b<-bundles[[paste0("imp",imp,"_",s)]];raw<-b$train_data
    lmfit<-lm(alm_kg~age+height_m,data=raw,weights=wtmec_correct);covs<-tibble(age=c(18,69,median(raw$age),median(raw$age)),height_m=c(median(raw$height_m),median(raw$height_m),min(raw$height_m),max(raw$height_m)))
    anchors<-raw[rep(1,4),];anchors$SEQN<- -(1:4);anchors$age<-covs$age;anchors$height_m<-covs$height_m;anchors$height_cm<-100*covs$height_m
    anchors$alm_kg<-pmax(as.numeric(predict(lmfit,covs)),.5);anchors$wtmec_correct<-mean(raw$wtmec_correct)*1e-4;anchors$wtmec_legacy<-0
    anchors$cycle<-"anchor";anchors$strata<-"anchor";anchors$psu<-paste0("anchor",1:4);anchors$strata_pool<-factor("anchor");anchors$psu_pool<-factor(paste0("anchor",1:4))
    train<-bind_rows(raw,anchors)%>%mutate(w_model=wtmec_correct/mean(wtmec_correct))%>%droplevels();fit<-fit_anchor(train,s,family_name)
    d<-validation%>%filter(imputation==imp,as.character(sex)==s);sc<-score(fit,train,family_name,d);tmp<-bind_cols(d,sc)%>%mutate(abs_delta=abs(z_anchor-z_H_corrected),p5_changed=as.numeric((p_anchor<.05)!=low_p5_H_corrected),p10_changed=as.numeric((p_anchor<.10)!=low_p10_H_corrected))
    des<-make_design(tmp);mn<-svymean(~abs_delta+p5_changed+p10_changed,des,na.rm=TRUE)
    person_single[[ps]]<-tibble(imputation=imp,sex=s,mean_abs_delta_z=coef(mn)[1],mean_abs_delta_z_se=SE(mn)[1],p5_reclass=coef(mn)[2],p5_reclass_se=SE(mn)[2],p10_reclass=coef(mn)[3],p10_reclass_se=SE(mn)[3],max_abs_delta_z=max(tmp$abs_delta));ps<-ps+1L
    hs<-wq(raw$height_m,raw$wtmec_correct,c(.05,.5,.95));grid<-expand_grid(age=18:69,height_m=hs)%>%mutate(height_cm=100*height_m,alm_kg=1)
    for(prob in c(.05,.10)){q0<-quant(b$model,raw,family_name,grid,prob);q1<-quant(fit,train,family_name,grid,prob);curve_rows[[cr]]<-tibble(imputation=imp,sex=s,centile=prob,mean_abs_curve_delta_kg=mean(abs(q1-q0)),max_abs_curve_delta_kg=max(abs(q1-q0)));cr<-cr+1L}
  }
}
p<-bind_rows(person_single);person<-bind_rows(
  p%>%transmute(imputation,sex,metric="mean_abs_delta_z",estimate=mean_abs_delta_z,se=mean_abs_delta_z_se),
  p%>%transmute(imputation,sex,metric="P5_reclassification",estimate=p5_reclass,se=p5_reclass_se),
  p%>%transmute(imputation,sex,metric="P10_reclassification",estimate=p10_reclass,se=p10_reclass_se))%>%group_by(sex,metric)%>%group_modify(~{x<-rubin(.x$estimate,.x$se);tibble(estimate=x[1],standard_error=x[2],conf_low=x[3],conf_high=x[4])})%>%ungroup()%>%left_join(p%>%group_by(sex)%>%summarise(max_abs_delta_z=max(max_abs_delta_z),.groups="drop"),by="sex")%>%mutate(record_type="participant")
curves<-bind_rows(curve_rows)%>%group_by(sex,centile)%>%summarise(mean_abs_curve_delta_kg=mean(mean_abs_curve_delta_kg),max_abs_curve_delta_kg=max(max_abs_curve_delta_kg),.groups="drop")%>%mutate(record_type="curve")
write_csv(bind_rows(person,curves),file.path(cache_dir,"synthetic_anchor_audit.csv"))
cat("Synthetic anchor audit completed\n")
