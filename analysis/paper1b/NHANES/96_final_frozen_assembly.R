options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr);library(purrr)})
root_dir<-normalizePath(".",winslash="/",mustWork=TRUE);paper_dir<-file.path(root_dir,"First Paper")
out_dir<-file.path(paper_dir,"Paper1_FINAL_FROZEN_20260814");cache_dir<-file.path(out_dir,"cache")

# 01. Unified, machine-readable sample flow.
legacy_flow<-read_csv(file.path(paper_dir,"final_age_range_comparison_20260814","01_sample_flow_age_range.csv"),show_col_types=FALSE)%>%
  filter(age_range=="18-69",sex%in%c("Female","Male"))
dev_flow<-bind_rows(
  legacy_flow%>%transmute(analysis_component="Primary H development",population="NHANES 1999-2006; age 18-69",sex,
    comparison=NA_character_,model="Corrected-weight H",stage="Reference eligible",unweighted_n=reference_eligible_n,notes="Eligibility before valid ALM"),
  legacy_flow%>%transmute(analysis_component="Primary H development",population="NHANES 1999-2006; age 18-69",sex,
    comparison=NA_character_,model="Corrected-weight H",stage="Valid ALM/reference model",unweighted_n=valid_alm_n,notes="Five completed data sets; n shown per completed set")
)
val<-readRDS(file.path(cache_dir,"validation_corrected_H_18_59.rds"))
validation_flow<-val%>%filter(imputation==1)%>%count(sex,name="unweighted_n")%>%transmute(
  analysis_component="Primary H temporal validation",population="NHANES 2011-2018; age 18-59",sex=as.character(sex),
  comparison=NA_character_,model="Corrected-weight H",stage="Valid DXA and H score",unweighted_n,
  notes="Validation weights were already WTMEC2YR/4")
grip_flow<-read_csv(file.path(cache_dir,"grip_sample_flow.csv"),show_col_types=FALSE)%>%
  pivot_longer(c(neither_n,conventional_only_n,conditional_only_n,both_n),names_to="stage",values_to="unweighted_n")%>%
  transmute(analysis_component=if_else(domain=="Full sample","Full-sample grip","Normal-BMI grip sensitivity"),
    population="NHANES 2011-2014; age 18-59",sex="Both",comparison,model,stage,unweighted_n,
    notes=paste0("Model analysis n=",round(analysis_n)))
support<-read_csv(file.path(cache_dir,"HB_temporal_safe_domain_audit.csv"),show_col_types=FALSE)%>%
  group_by(sex,model)%>%summarise(n=mean(n),outside=mean(outside_support_n),invalid=mean(invalid_prediction_n),.groups="drop")
hb_flow<-support%>%transmute(analysis_component="Frozen H/HB temporal comparison",population="NHANES 2011-2018; age 18-59; BMI-complete",
  sex,comparison="H-common vs HB",model,stage="Scored before common-valid restriction",unweighted_n=n,
  notes=paste0("Outside observed predictor support: ",round(outside,1),"; invalid predictions: ",round(invalid,1)))
sample_flow<-bind_rows(dev_flow,validation_flow,grip_flow,hb_flow)%>%mutate(weight_version="Corrected/final unless explicitly stated")
write_csv(sample_flow,file.path(out_dir,"01_sample_flow.csv"))

# 02. Final weight audit, preserving both the discovered legacy field and the corrected pooled weight.
raw<-read_csv(file.path(cache_dir,"development_corrected_weight_narrow.csv.gz"),show_col_types=FALSE,progress=FALSE)%>%
  transmute(SEQN=as.numeric(SEQN),imputation=as.integer(imputation),cycle=as.character(cycle),examined=RIDSTATR==2,
    sex=case_when(RIAGENDR==1~"Male",RIAGENDR==2~"Female",TRUE~NA_character_),age=as.numeric(RIDAGEYR),
    pregnant=coalesce(RIAGENDR==2&RIDEXPRG==1,FALSE),height_cm=as.numeric(BMXHT),alm_kg=as.numeric(alm_g)/1000,
    source=if_else(cycle%in%c("1999-2000","2001-2002"),"WTMEC4YR","WTMEC2YR"),
    source_weight=if_else(source=="WTMEC4YR",as.numeric(WTMEC4YR),as.numeric(WTMEC2YR)),
    multiplier=if_else(source=="WTMEC4YR",.5,.25),legacy_weight=as.numeric(wtmec_pool),corrected_weight=as.numeric(wtmec_correct),
    strata=as.character(SDMVSTRA),psu=as.character(SDMVPSU))%>%
  filter(imputation==1,examined,sex%in%c("Female","Male"),age>=18,age<=69,!pregnant,is.finite(height_cm),height_cm>0,
    is.finite(alm_kg),alm_kg>0,is.finite(corrected_weight),corrected_weight>0,!is.na(strata),!is.na(psu))%>%
  group_by(sex)%>%mutate(centered_weight=corrected_weight/mean(corrected_weight))%>%ungroup()
dev_weight<-raw%>%group_by(sex,cycle,source,multiplier)%>%summarise(n=n(),source_weight_sum=sum(source_weight),legacy_sum_weight=sum(legacy_weight),
  sum_weight=sum(corrected_weight),centered_mean=mean(centered_weight),centered_min=min(centered_weight),centered_max=max(centered_weight),
  corrected_max_abs_error=max(abs(corrected_weight-source_weight*multiplier)),legacy_max_abs_error=max(abs(legacy_weight-source_weight*multiplier)),.groups="drop")%>%
  mutate(period="Development",weight_used_final="corrected_weight",gamlss_centering="corrected_weight / within-sex mean",
    survey_weight="corrected_weight",survey_strata="interaction(cycle, SDMVSTRA)",survey_psu="interaction(cycle, SDMVSTRA, SDMVPSU)",
    corrected_rule_verified=corrected_max_abs_error<1e-8,legacy_rule_verified=legacy_max_abs_error<1e-8)
overall<-raw%>%group_by(sex)%>%summarise(cycle="Overall 1999-2006",source="cycle-specific",multiplier=NA_real_,n=n(),source_weight_sum=sum(source_weight),
  legacy_sum_weight=sum(legacy_weight),sum_weight=sum(corrected_weight),centered_mean=mean(centered_weight),centered_min=min(centered_weight),centered_max=max(centered_weight),
  corrected_max_abs_error=max(abs(corrected_weight-source_weight*multiplier)),legacy_max_abs_error=max(abs(legacy_weight-source_weight*multiplier)),.groups="drop")%>%
  mutate(period="Development",weight_used_final="corrected_weight",gamlss_centering="corrected_weight / within-sex mean",survey_weight="corrected_weight",
    survey_strata="interaction(cycle, SDMVSTRA)",survey_psu="interaction(cycle, SDMVSTRA, SDMVPSU)",corrected_rule_verified=TRUE,legacy_rule_verified=FALSE)
weight_audit<-bind_rows(dev_weight,overall)%>%select(period,sex,cycle,source,multiplier,n,source_weight_sum,legacy_sum_weight,sum_weight,
  centered_mean,centered_min,centered_max,corrected_max_abs_error,legacy_max_abs_error,corrected_rule_verified,legacy_rule_verified,
  weight_used_final,gamlss_centering,survey_weight,survey_strata,survey_psu)%>%arrange(sex,cycle)
write_csv(weight_audit,file.path(out_dir,"02_weight_audit.csv"))

# 03. Ensure the frozen specification records the numerical implementation and correction.
spec<-read_csv(file.path(out_dir,"03_primary_H_model_specification.csv"),show_col_types=FALSE)%>%mutate(
  development_weight="1999-2002 WTMEC4YR x 0.5; 2003-2006 WTMEC2YR x 0.25",
  internal_basis_intervals=if_else(sex=="Female","10 (numerical stabilization; df remains 3)","GAMLSS default; df=3"),
  weight_audit_status="PASS after corrected-weight refit",frozen_parameter_date="2026-08-14")
write_csv(spec,file.path(out_dir,"03_primary_H_model_specification.csv"))

# 12-13 were closed before the weight error was discovered.  Keep them only as
# transparent qualitative audit trails; they are not sources for final primary numbers.
age_sens<-read_csv(file.path(paper_dir,"final_age_range_comparison_20260814","24_final_age_range_decision_table.csv"),show_col_types=FALSE)%>%
  mutate(weight_version="Legacy unscaled development field",final_use="Qualitative age-range decision only; excluded from corrected primary numeric claims")
write_csv(age_sens,file.path(out_dir,"12_age_range_18_84_sensitivity.csv"))
refine<-read_csv(file.path(out_dir,"13_female_P5_refinement_summary.csv"),show_col_types=FALSE)%>%mutate(
  weight_version="Legacy unscaled development field",final_use="Qualitative model-refinement audit only; no candidate replaced corrected primary H")
write_csv(refine,file.path(out_dir,"13_female_P5_refinement_summary.csv"))

# File-level freeze manifest.
files<-list.files(out_dir,pattern="^[0-9]{2}_.*\\.csv$",full.names=TRUE)
manifest<-map_dfr(files,function(p)tibble(file=basename(p),bytes=file.info(p)$size,rows=nrow(read_csv(p,show_col_types=FALSE,progress=FALSE)),
  modified=as.character(file.info(p)$mtime)))
write_csv(manifest,file.path(cache_dir,"final_csv_manifest.csv"))
cat("Final assembly and corrected weight audit completed\n")
