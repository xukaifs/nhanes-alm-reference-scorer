source("R/score_conditional_alm.R", encoding = "UTF-8")

# This reports a warning, rather than stopping, if the installed package
# versions differ from the versions used for release validation. For an exact
# reproducibility check, run environment/verify_environment.R separately.
environment_status <- check_alm_reference_environment(strict = FALSE)

reference_models <- load_alm_reference(
  "models/revised_primary_model_bundles.rds"
)

example_patients <- data.frame(
  case = c(
    "Identical ALM: shorter stature",
    "Identical ALM: taller stature",
    "Near-identical ALMI: younger/taller",
    "Near-identical ALMI: older/shorter"
  ),
  sex = c("Female", "Female", "Male", "Male"),
  age = c(45, 45, 38, 58),
  height_m = c(1.55, 1.75, 1.790, 1.596),
  alm_kg = c(15.0, 15.0, 22.40, 17.82)
)

example_scores <- score_conditional_alm(
  example_patients,
  bundles = reference_models
)

print(example_scores$pooled)

# Rows 1-2 hold age and ALM constant while height changes.
# Rows 3-4 have nearly identical ALMI but different age and height.

# Create selected centiles for one age-height combination.
reference_example <- reference_centiles(
  data.frame(sex = "Female", age = 45, height_m = 1.60),
  bundles = reference_models
)
print(reference_example)
