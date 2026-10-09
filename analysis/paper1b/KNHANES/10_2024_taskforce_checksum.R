# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({library(haven); library(survey)})
root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
path <- file.path(root,"raw_all_sav","HN24_ALL.sav"); out <- file.path(root,"tables_2024")
wanted <- c("ID","psu","sex","age","kstrata","wt_itvex","HE_ht","HE_wt","DXW_ex","GS_SP_DXA",
  "DW_Lrm_LN","DW_Rrm_LN","DW_Llg_LN","DW_Rlg_LN","DW_Lrm_BMC","DW_Rrm_BMC","DW_Llg_BMC","DW_Rlg_BMC",
  "DW_WBT_pFT","DW_SBT_pFT","GS_mea_r_1","GS_mea_r_2","GS_mea_l_1","GS_mea_l_2")
x <- as.data.frame(read_sav(path,col_select=all_of(wanted))); names(x) <- tolower(names(x))
x$psu <- trimws(as.character(x$psu))
for(v in setdiff(names(x),c("id","psu"))) x[[v]] <- suppressWarnings(as.numeric(zap_labels(zap_missing(x[[v]]))))
valid <- function(z,lo,hi) is.finite(z)&z>lo&z<hi&abs(z-9999)>.02&abs(z-9999.9)>.02&abs(z-9999.99)>.02
ln <- c("dw_lrm_ln","dw_rrm_ln","dw_llg_ln","dw_rlg_ln"); bmc <- c("dw_lrm_bmc","dw_rrm_bmc","dw_llg_bmc","dw_rlg_bmc")
lnok <- vapply(x[ln],valid,logical(nrow(x)),lo=0,hi=20000); bmcok <- vapply(x[bmc],valid,logical(nrow(x)),lo=0,hi=900)
ok <- rowSums(lnok)==4L & rowSums(bmcok)==4L
x$alm_kg <- NA_real_; x$alm_kg[ok] <- (rowSums(x[ok,ln,drop=FALSE])-rowSums(x[ok,bmc,drop=FALSE]))/1000
x$alm_ln_kg <- NA_real_; x$alm_ln_kg[rowSums(lnok)==4L] <- rowSums(x[rowSums(lnok)==4L,ln,drop=FALSE])/1000
x$height_cm <- ifelse(valid(x$he_ht,120,220),x$he_ht,NA_real_); x$weight_kg <- ifelse(valid(x$he_wt,20,300),x$he_wt,NA_real_)
x$female <- as.numeric(x$sex==2); x$body_fat_pct <- ifelse(valid(x$dw_wbt_pft,1,80),x$dw_wbt_pft,NA_real_)
x$subtotal_fat_pct <- ifelse(valid(x$dw_sbt_pft,1,80),x$dw_sbt_pft,NA_real_)
grip <- c("gs_mea_r_1","gs_mea_r_2","gs_mea_l_1","gs_mea_l_2")
for(v in grip)x[[v]]<-ifelse(valid(x[[v]],0,100),x[[v]],NA_real_)
x$grip_single_max <- apply(x[grip],1,function(z)if(all(!is.finite(z)))NA_real_ else max(z,na.rm=TRUE))
x$wt <- ifelse(is.finite(x$wt_itvex)&x$wt_itvex>0,x$wt_itvex,NA_real_)
base_keep <- is.finite(x$age)&x$age>=50&is.finite(x$alm_ln_kg)&is.finite(x$wt)&x$wt>0&is.finite(x$kstrata)&!is.na(x$psu)&nzchar(x$psu)
keep <- base_keep & is.finite(x$grip_single_max)
d <- x[keep,]; des <- svydesign(ids=~psu,strata=~kstrata,weights=~wt,nest=TRUE,data=d)
m <- svymean(~age+female+height_cm+weight_kg+alm_ln_kg+alm_kg+body_fat_pct+subtotal_fat_pct+grip_single_max,des,na.rm=TRUE); ci <- confint(m)
targets <- c(age=62.9,female=.519,height_cm=162.7,weight_kg=64.3,alm_ln_kg=17.0,alm_kg=NA,body_fat_pct=NA,subtotal_fat_pct=36.3,grip_single_max=31.1)
res <- data.frame(metric=names(coef(m)),estimate=as.numeric(coef(m)),lcl=ci[,1],ucl=ci[,2],published_target=targets[names(coef(m))])
for(col in c("estimate","lcl","ucl","published_target"))res[res$metric=="female",col]<-100*res[res$metric=="female",col]
res$difference <- res$estimate-res$published_target
write.csv(data.frame(stage=c("Age >=50 with valid 2024 ALM and survey design","Plus height and weight","Max single-hand grip measured (Task Force sample)"),
  n=c(sum(base_keep),sum(base_keep&is.finite(x$height_cm)&is.finite(x$weight_kg)),sum(keep))),file.path(out,"11_taskforce_checksum_flow.csv"),row.names=FALSE)
write.csv(res,file.path(out,"12_taskforce_checksum_estimates.csv"),row.names=FALSE)
selection <- rbind(
  data.frame(variable="DXW_ex",value=names(table(d$dxw_ex,useNA="ifany")),n=as.integer(table(d$dxw_ex,useNA="ifany"))),
  data.frame(variable="GS_SP_DXA",value=names(table(d$gs_sp_dxa,useNA="ifany")),n=as.integer(table(d$gs_sp_dxa,useNA="ifany")))
)
write.csv(selection,file.path(out,"13_taskforce_selection_flags.csv"),row.names=FALSE)
cat("2024 Task Force checksum complete; N =",nrow(d),"\n"); print(res,row.names=FALSE); print(selection,row.names=FALSE)
