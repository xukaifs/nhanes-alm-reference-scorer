source("R/score_conditional_alm.R", encoding = "UTF-8")

check_alm_reference_environment(strict = FALSE)
reference_models <- load_alm_reference()

example_patients <- data.frame(
  case = c("Female example", "Male example"),
  sex = c("Female", "Male"),
  age = c(45, 60),
  height_m = c(1.60, 1.75),
  alm_kg = c(15.0, 22.0)
)

example_scores <- score_conditional_alm(
  example_patients,
  bundles = reference_models
)
print(example_scores$pooled)

example_centiles <- reference_centiles(
  data.frame(
    sex = c("Female", "Male"),
    age = c(45, 60),
    height_m = c(1.60, 1.75)
  ),
  bundles = reference_models
)
print(example_centiles)
