# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
ow <- file.path(root,"tables_official_weight"); y24 <- file.path(root,"tables_2024")
checks <- list(); add <- function(name,ok,detail="") checks[[length(checks)+1L]] <<- data.frame(check=name,pass=isTRUE(ok),detail=detail)

coef <- read.csv(file.path(ow,"00_official_pooling_coefficients.csv"))
add("Official pooling coefficients sum to one",abs(sum(coef$coefficient)-1)<1e-12,sprintf("sum=%.12f",sum(coef$coefficient)))
add("Official pooling numerators are 108/199/192/80",identical(coef$survey_areas,c(108L,199L,192L,80L)),paste(coef$survey_areas,collapse=","))

old <- readRDS(file.path(root,"data","knhanes_2008_2011_analysis_19_69_official_weight.rds"))
new <- readRDS(file.path(root,"data","knhanes_2024_analysis_40_69.rds"))
add("2008-2011 final N",nrow(old)==16226,paste("N=",nrow(old)))
add("2024 final N",nrow(new)==2513,paste("N=",nrow(new)))
add("2008-2011 age domain",all(old$age>=19&old$age<=69),sprintf("range=%d-%d",min(old$age),max(old$age)))
add("2024 age domain",all(new$age>=40&new$age<=69),sprintf("range=%d-%d",min(new$age),max(new$age)))
add("2024 ALM equals sum of limb LN",max(abs(new$alm_kg-rowSums(new[c("dw_lrm_ln","dw_rrm_ln","dw_llg_ln","dw_rlg_ln")])/1000))<1e-10)

flow_old <- read.csv(file.path(ow,"08_taskforce_checksum_flow.csv"))
add("Published 2008-2011 checksum N reproduced",flow_old$n[flow_old$criterion=="Task Force LN*0.946 ALM valid"]==9219)
est_old <- read.csv(file.path(ow,"09_taskforce_checksum_estimates.csv"))
target_old <- subset(est_old,!is.na(published_target))
add("Published 2008-2011 checksum estimates reproduce within 0.15 reported units",
    all(abs(target_old$estimate-target_old$published_target)<=.15),
    paste(target_old$metric,round(target_old$estimate,2),collapse="; "))

est24 <- read.csv(file.path(y24,"12_taskforce_checksum_estimates.csv"))
target24 <- subset(est24,!is.na(published_target))
flow24 <- read.csv(file.path(y24,"11_taskforce_checksum_flow.csv"))
add("Published 2024 checksum N reproduced",tail(flow24$n,1)==2481,paste("N=",tail(flow24$n,1)))
add("Published 2024 checksum estimates reproduce within 0.15 reported units",
    all(abs(target24$estimate-target24$published_target)<=.15),
    paste(target24$metric,round(target24$estimate,2),collapse="; "))

for(file in c(file.path(ow,"01_calibration.csv"),file.path(ow,"02_z_height_slopes.csv"),file.path(ow,"03_low_alm_height_or.csv"),
              file.path(y24,"02_calibration.csv"),file.path(y24,"03_z_height_slopes.csv"),file.path(y24,"04_low_alm_height_or.csv"),
              file.path(y24,"08_grip_discordance_models.csv"))) {
  z <- read.csv(file); numeric <- vapply(z,is.numeric,logical(1));
  add(paste("No infinite numeric estimates:",basename(file)),!any(is.infinite(as.matrix(z[numeric]))),paste("rows=",nrow(z)))
}
cal <- read.csv(file.path(ow,"01_calibration.csv")); cal24 <- read.csv(file.path(y24,"02_calibration.csv"))
add("All prevalence estimates are within 0-100",all(cal$p5_percent>=0&cal$p5_percent<=100&cal$p10_percent>=0&cal$p10_percent<=100)&
      all(cal24$p5_percent>=0&cal24$p5_percent<=100&cal24$p10_percent>=0&cal24$p10_percent<=100))
or_cols <- c("or_per_10cm","or_lcl","or_ucl")
ors <- rbind(read.csv(file.path(ow,"03_low_alm_height_or.csv"))[or_cols],
             read.csv(file.path(y24,"04_low_alm_height_or.csv"))[or_cols])
add("All odds ratios and confidence limits are positive",all(ors$or_per_10cm>0&ors$or_lcl>0&ors$or_ucl>0))

result <- do.call(rbind,checks); write.csv(result,file.path(root,"FINAL_VERIFICATION.csv"),row.names=FALSE)
if(!all(result$pass)){print(result,row.names=FALSE);stop("Final verification failed")}
cat("ALL FINAL VERIFICATION CHECKS PASSED (",nrow(result)," checks).\n",sep="")
print(result,row.names=FALSE)
