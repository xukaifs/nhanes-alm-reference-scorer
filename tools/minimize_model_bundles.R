options(stringsAsFactors = FALSE)

source_files <- file.path("models", c(
  "corrected_H_18_69_Female_bundles.rds",
  "corrected_H_18_69_Male_bundles.rds"
))
stopifnot(all(file.exists(source_files)))
out_dir <- file.path("models", "public")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
keep <- c("alm_kg", "age", "height_m", "w_model")

for (path in source_files) {
  bundles <- readRDS(path)
  for (i in seq_along(bundles)) {
    td <- as.data.frame(bundles[[i]]$train_data)
    stopifnot(all(keep %in% names(td)))
    bundles[[i]]$train_data <- td[, keep, drop = FALSE]
  }
  target <- file.path(out_dir, sub("_bundles\\.rds$", "_public_bundles.rds", basename(path)))
  saveRDS(bundles, target, version = 3)
  message("Wrote: ", target)
}

message("After minimization: point load_alm_reference(model_path=...) to the two public files, rerun tests, and regenerate SHA-256 hashes before release.")
