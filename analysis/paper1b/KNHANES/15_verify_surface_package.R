# -*- coding: UTF-8 -*-
options(stringsAsFactors=FALSE)
root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
out <- file.path(root,"us_korea_surface")

audit <- read.csv(file.path(out,"01_local_model_fit_audit.csv"),check.names=FALSE)
cal <- read.csv(file.path(out,"02_local_z_calibration.csv"),check.names=FALSE)
slopes <- read.csv(file.path(out,"03_local_z_residual_slopes.csv"),check.names=FALSE)
domains <- read.csv(file.path(out,"04_common_height_domains.csv"),check.names=FALSE)
grid <- read.csv(file.path(out,"05_common_grid_centiles_and_differences.csv"),check.names=FALSE)
geom <- read.csv(file.path(out,"06_location_scale_geometry_summary.csv"),check.names=FALSE)
ages <- read.csv(file.path(out,"07_age_pattern_summary.csv"),check.names=FALSE)
curves <- read.csv(file.path(out,"08_representative_age_curves.csv"),check.names=FALSE)

checks <- list()
add_check <- function(name,ok,detail){
  checks[[length(checks)+1L]] <<- data.frame(check=name,status=if(isTRUE(ok))"PASS" else "FAIL",detail=as.character(detail))
}

female_formula <- "alm_kg ~ pb(age, df=3, inter=10) + pb(height_m, df=3, inter=10)"
male_formula <- "alm_kg ~ pb(age, df=3) + pb(height_m, df=3)"
add_check("Exact mirrored model specifications",
  identical(audit$mu_formula[audit$sex=="Female"],female_formula) &&
    identical(audit$mu_formula[audit$sex=="Male"],male_formula) &&
    all(audit$sigma_formula=="~1") && all(audit$nu_formula=="~1") && all(audit$tau_formula=="~1") &&
    identical(audit$family, c("BCCG","BCT")),
  paste(audit$sex,audit$family,audit$mu_formula,collapse="; "))
add_check("Model convergence and sample size",all(audit$converged) && identical(audit$n,c(9240L,6986L)),
  paste0("converged=",paste(audit$converged,collapse=","),"; N=",paste(audit$n,collapse=",")))
add_check("Normalized fitting weights",all(abs(audit$weight_mean-1)<1e-10),paste(audit$weight_mean,collapse=","))
add_check("Local z mean and SD calibration",all(abs(cal$weighted_mean_z)<.01) && all(cal$weighted_sd_z>.98 & cal$weighted_sd_z<1.02),
  paste(cal$sex,sprintf("mean=%.5f SD=%.5f",cal$weighted_mean_z,cal$weighted_sd_z),collapse="; "))
add_check("Local z tail calibration",all(cal$p5_percent>3 & cal$p5_percent<7) && all(cal$p10_percent>8 & cal$p10_percent<12),
  paste(cal$sex,sprintf("P5=%.2f%% P10=%.2f%%",cal$p5_percent,cal$p10_percent),collapse="; "))
add_check("Residual age-height independence",all(slopes$p_value>.05),
  paste(slopes$sex,slopes$term,sprintf("b=%.5f p=%.3f",slopes$estimate,slopes$p_value),collapse="; "))
add_check("Common domains valid",all(domains$primary_common_min_m<domains$primary_common_max_m) &&
    all(domains$strict_common_min_m<domains$strict_common_max_m) &&
    all(domains$strict_common_min_m>=domains$primary_common_min_m) &&
    all(domains$strict_common_max_m<=domains$primary_common_max_m),
  paste(unique(domains$sex),collapse=","))
add_check("Common grid complete",nrow(grid)==2703 && !anyNA(grid),paste0("rows=",nrow(grid),"; missing=",sum(is.na(grid))))
for(ref in c("us","korea")){
  cc <- grid[,paste0(ref,"_p",c(5,10,50,90,95))]
  add_check(paste0(if(ref=="us")"U.S." else "Korean"," centiles ordered"),
    all(apply(cc,1,function(x)all(diff(x)>0))),paste0("rows checked=",nrow(cc)))
}
add_check("Surface summaries complete",nrow(geom)==4 && all(table(geom$sex)==2) && all(geom$p50_r_squared>=0 & geom$p50_r_squared<=1),
  paste(geom$sex,geom$domain,collapse="; "))
add_check("Scale ratios valid",all(geom$width_p95_p5_ratio_mean>0 & geom$width_p95_p5_ratio_mean<1) &&
    all(geom$width_p90_p10_ratio_mean>0 & geom$width_p90_p10_ratio_mean<1),
  paste(geom$sex,sprintf("P95-P5 ratio=%.3f",geom$width_p95_p5_ratio_mean),collapse="; "))
add_check("Age bands complete",nrow(ages)==10 && all(table(ages$sex)==5),paste0("rows=",nrow(ages)))
add_check("Representative curves complete",nrow(curves)==1224 && all(table(curves$sex)==612) &&
    length(unique(curves$height_level[curves$sex=="Female"]))==3 && length(unique(curves$height_level[curves$sex=="Male"]))==3,
  paste0("rows=",nrow(curves),"; labels=",paste(unique(curves$height_level),collapse=", ")))

expected_figures <- c(
  paste0("Figure_representative_curves_",rep(c("female","male"),each=2),rep(c(".png",".pdf"),2)),
  unlist(lapply(c("female","male"),function(s)unlist(lapply(c("_kg","_percent"),function(m)
    paste0("Figure_delta_P50_heatmap_",s,m,c(".png",".pdf")))))))
figure_paths <- file.path(out,"figures",expected_figures)
add_check("All figures generated",all(file.exists(figure_paths)) && all(file.info(figure_paths)$size>1000),
  paste0(sum(file.exists(figure_paths)),"/",length(figure_paths)," files present"))

result <- do.call(rbind,checks)
write.csv(result,file.path(out,"09_surface_package_verification.csv"),row.names=FALSE)
print(result,row.names=FALSE)
if(any(result$status!="PASS"))stop("Surface package verification failed: ",paste(result$check[result$status!="PASS"],collapse=", "))
cat("All surface package checks passed.\n")
