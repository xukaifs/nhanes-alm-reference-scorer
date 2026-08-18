options(stringsAsFactors = FALSE)
source("R/score_conditional_alm.R", encoding = "UTF-8")
models <- load_alm_reference()

# Conservative release grid using the recommended stature guardrails.
# Heights are 1-cm increments; ages span the final 18-69 domain.
female <- expand.grid(
  sex = "Female", age = 18:69, height_m = seq(1.46, 1.78, by = 0.01),
  stringsAsFactors = FALSE
)
male <- expand.grid(
  sex = "Male", age = 18:69, height_m = seq(1.57, 1.94, by = 0.01),
  stringsAsFactors = FALSE
)
grid <- rbind(female, male)
result <- reference_centiles(
  grid,
  probabilities = c(0.025,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.975),
  bundles = models,
  warn = FALSE
)
result$height_cm <- 100 * result$height_m
result <- result[, c("sex","age","height_cm", grep("^alm_p", names(result), value=TRUE))]
utils::write.csv(result, gzfile("data/age_height_reference_grid_v2.csv.gz"), row.names=FALSE)
message("Wrote data/age_height_reference_grid_v2.csv.gz with ", nrow(result), " rows.")
