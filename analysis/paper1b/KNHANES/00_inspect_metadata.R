# -*- coding: UTF-8 -*-
options(stringsAsFactors = FALSE)

.libPaths(c("C:/Users/Public/CodexRLib42Copy", "C:/Program Files/R/R-4.2.1/library"))
suppressPackageStartupMessages(library(haven))

invisible(Sys.setlocale("LC_ALL", "English_United States.utf8"))
native_path <- identity

root <- native_path("D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/R\u5206\u6790\u7ed3\u679c/KNHANES_external_validation_20260904")
raw_all <- file.path(root, "raw_all_sav")
raw_dxa <- native_path("D:/DXA\u8eab\u4f53\u6210\u5206\u5206\u6790/Knhanes\u6570\u636e\u5e93")
dir.create(file.path(root, "metadata"), showWarnings = FALSE, recursive = TRUE)

years <- 2008:2011
all_files <- setNames(
  file.path(raw_all, c("HN08_ALL.sav", "HN09_all.sav", "HN10_all.sav", "HN11_all.sav")),
  years
)
dxa_files <- setNames(file.path(raw_dxa, sprintf("hn%02d_dxa.sav", 8:11)), years)

patterns <- paste(
  c(
    "age", "sex", "gender", "height", "weight", "body mass", "bmi",
    "psu", "strat", "weight", "wt_", "ID", "identifier",
    "arm", "leg", "lean", "fat-free", "bone mineral", "BMC",
    "\uc5f0\ub839", "\uc131\ubcc4", "\uc2e0\uc7a5", "\uccb4\uc911", "\uccb4\uc9c8\ub7c9", "\uac00\uc911\uce58",
    "\uce35\ud654", "\uc9c0\ubc29", "\uc81c\uc9c0\ubc29", "\uace8\uc5fc", "\ud314", "\ub2e4\ub9ac"
  ),
  collapse = "|"
)

inspect_one <- function(path, year, kind) {
  x <- read_sav(path, n_max = 1)
  labels <- vapply(x, function(z) {
    value <- attr(z, "label", exact = TRUE)
    if (is.null(value)) "" else as.character(value)
  }, character(1))
  meta <- data.frame(
    year = year,
    kind = kind,
    variable = names(x),
    label = labels,
    class = vapply(x, function(z) paste(class(z), collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(
    meta,
    file.path(root, "metadata", sprintf("%s_%s_variables.csv", year, kind)),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  meta[grepl(patterns, paste(meta$variable, meta$label), ignore.case = TRUE), ]
}

hits <- list()
for (year in years) {
  hits[[paste0(year, "_all")]] <- inspect_one(all_files[[as.character(year)]], year, "all")
  hits[[paste0(year, "_dxa")]] <- inspect_one(dxa_files[[as.character(year)]], year, "dxa")
}
hits <- do.call(rbind, hits)
write.csv(
  hits,
  file.path(root, "metadata", "candidate_variables.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

print(hits, row.names = FALSE)
