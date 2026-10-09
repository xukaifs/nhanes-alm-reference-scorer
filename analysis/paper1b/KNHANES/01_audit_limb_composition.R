# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
dxa_root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/Knhanes\u6570\u636e\u5e93"
dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)

limbs <- c("Lrm", "Rrm", "Llg", "Rlg")
needed <- c("ID", unlist(lapply(limbs, function(x) paste0("DW_", x, c("_LN", "_BMC", "_FT", "_MS")))))

clean_numeric <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x) | x >= 9999] <- NA_real_
  x
}

rows <- list()
for (year in 2008:2011) {
  path <- file.path(dxa_root, sprintf("hn%02d_dxa.sav", year - 2000))
  x <- read_sav(path, col_select = any_of(needed))
  for (v in setdiff(names(x), "ID")) x[[v]] <- clean_numeric(x[[v]])
  for (limb in limbs) {
    ln <- x[[paste0("DW_", limb, "_LN")]]
    bmc <- x[[paste0("DW_", limb, "_BMC")]]
    ft <- x[[paste0("DW_", limb, "_FT")]]
    ms <- x[[paste0("DW_", limb, "_MS")]]
    complete <- is.finite(ln) & is.finite(bmc) & is.finite(ft) & is.finite(ms)
    rows[[length(rows) + 1L]] <- data.frame(
      year = year,
      limb = limb,
      n_complete = sum(complete),
      mean_ln_g = mean(ln[complete]),
      mean_bmc_g = mean(bmc[complete]),
      mean_abs_mass_minus_ln_ft_bmc = mean(abs(ms[complete] - ln[complete] - ft[complete] - bmc[complete])),
      mean_abs_mass_minus_ln_ft = mean(abs(ms[complete] - ln[complete] - ft[complete])),
      cor_mass_ln_ft_bmc = cor(ms[complete], ln[complete] + ft[complete] + bmc[complete]),
      cor_mass_ln_ft = cor(ms[complete], ln[complete] + ft[complete]),
      stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(root, "tables", "00_limb_composition_identity_audit.csv"), row.names = FALSE)
print(out, row.names = FALSE)
