# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE)
invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages(library(haven))

root <- "D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904"
path <- file.path(root, "raw_all_sav", "HN24_ALL.sav")
out <- file.path(root, "metadata")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
x <- read_sav(path, n_max = 1)
meta <- data.frame(variable = names(x), label = vapply(x, function(z) {
  lab <- attr(z, "label"); if (is.null(lab)) "" else as.character(lab)
}, character(1)), class = vapply(x, function(z) class(z)[1], character(1)))
write.csv(meta, file.path(out, "2024_all_variables.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pat <- paste(c("^id$", "^sex$", "^age$", "psu", "strata", "wt_", "he_ht", "he_wt", "he_bmi",
               "dxa", "^dw_", "lean", "bmc", "\uadfc\uc721", "\uc81c\uc9c0\ubc29", "\uc545\ub825", "grip", "^gs_"), collapse = "|")
candidates <- meta[grepl(pat, meta$variable, ignore.case = TRUE) |
                     grepl(pat, meta$label, ignore.case = TRUE), ]
write.csv(candidates, file.path(out, "2024_candidate_variables.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(candidates, row.names = FALSE)
