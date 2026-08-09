options(stringsAsFactors = FALSE)

if (!file.exists(file.path("R", "score_conditional_alm.R"))) {
  stop("Run this test from the repository root.", call. = FALSE)
}

source(file.path("environment", "verify_environment.R"), encoding = "UTF-8")
source(file.path("R", "score_conditional_alm.R"), encoding = "UTF-8")

bundles <- load_alm_reference(
  file.path("models", "revised_primary_model_bundles.rds")
)
validate_reference_bundles(bundles)

minimal_training_columns <- c(
  "alm_kg", "age", "height_m", "w_model"
)
stopifnot(all(vapply(
  bundles,
  function(bundle) {
    identical(names(bundle$train_data), minimal_training_columns)
  },
  logical(1)
)))

saved_examples <- read.csv(
  file.path("data", "clinical_calculation_examples.csv"),
  check.names = FALSE
)
expected_case_labels <- c(
  "Identical ALM: shorter stature",
  "Identical ALM: taller stature",
  "Near-identical ALMI: younger/taller",
  "Near-identical ALMI: older/shorter"
)
expected_example_input <- data.frame(
  sex = c("Female", "Female", "Male", "Male"),
  age = c(45, 45, 38, 58),
  height_m = c(1.55, 1.75, 1.790, 1.596),
  alm_kg = c(15.0, 15.0, 22.40, 17.82)
)
stopifnot(
  nrow(saved_examples) == 4,
  identical(as.character(saved_examples$case), expected_case_labels),
  identical(as.character(saved_examples$sex), expected_example_input$sex),
  max(abs(saved_examples$age - expected_example_input$age)) < 1e-12,
  max(abs(
    saved_examples$height_cm / 100 -
      expected_example_input$height_m
  )) < 1e-12,
  max(abs(
    saved_examples$alm_kg -
      expected_example_input$alm_kg
  )) < 1e-12
)
example_input <- data.frame(
  sex = saved_examples$sex,
  age = saved_examples$age,
  height_m = saved_examples$height_cm / 100,
  alm_kg = saved_examples$alm_kg
)

scoring_warnings <- character()
example_scores <- withCallingHandlers(
  score_conditional_alm(
    example_input,
    bundles = bundles,
    warn = FALSE
  )$pooled,
  warning = function(warning_condition) {
    scoring_warnings <<- c(
      scoring_warnings,
      conditionMessage(warning_condition)
    )
    invokeRestart("muffleWarning")
  }
)
stopifnot(length(scoring_warnings) == 0)
stopifnot(all(
  abs(
    example_scores$model_averaged_z -
      saved_examples$mean_z
  ) < 1e-12
))
stopifnot(all(
  abs(
    example_scores$model_averaged_percentile_pct -
      saved_examples$mean_percentile_pct
  ) < 1e-12
))
stopifnot(all(
  abs(
    example_scores$almi_kg_m2 -
      saved_examples$almi_kg_m2
  ) < 1e-12
))
stopifnot(identical(
  example_scores$use_caution,
  saved_examples$use_caution
))
stopifnot(
  example_scores$model_averaged_z[3] < -1.645,
  example_scores$model_averaged_z[4] > -1.645
)

# Regression tests for one-person and single-sex scoring. Earlier release
# candidates only exercised mixed-sex batches, which could hide empty-branch
# length errors in the other sex's five model bundles.
single_female <- score_conditional_alm(
  example_input[1, , drop = FALSE],
  bundles = bundles,
  warn = FALSE
)
single_male <- score_conditional_alm(
  example_input[3, , drop = FALSE],
  bundles = bundles,
  warn = FALSE
)
stopifnot(
  nrow(single_female$pooled) == 1,
  nrow(single_female$by_imputation) == 5,
  is.finite(single_female$pooled$model_averaged_z),
  nrow(single_male$pooled) == 1,
  nrow(single_male$by_imputation) == 5,
  is.finite(single_male$pooled$model_averaged_z),
  abs(
    single_female$pooled$model_averaged_z -
      example_scores$model_averaged_z[1]
  ) < 1e-12,
  abs(
    single_male$pooled$model_averaged_z -
      example_scores$model_averaged_z[3]
  ) < 1e-12
)

invalid_score_input <- data.frame(
  sex = c("Other", "Female", "Female", "Female", NA),
  age = c(45, 60, 45, 45, 45),
  height_m = c(1.70, 1.60, 2.20, 1.60, 1.70),
  alm_kg = c(20, 15, 15, 0, 20)
)
invalid_scores <- score_conditional_alm(
  invalid_score_input,
  bundles = bundles,
  warn = FALSE
)
stopifnot(
  nrow(invalid_scores$pooled) == nrow(invalid_score_input),
  all(!invalid_scores$pooled$score_valid),
  all(is.na(invalid_scores$pooled$model_averaged_z)),
  all(is.na(invalid_scores$pooled$model_averaged_percentile))
)

# Numeric factors must be interpreted by their labels rather than their
# internal factor codes.
factor_input <- data.frame(
  sex = factor("Female"),
  age = factor("45"),
  height_m = factor("1.55"),
  alm_kg = factor("15")
)
factor_score <- score_conditional_alm(
  factor_input,
  bundles = bundles,
  warn = FALSE
)$pooled
stopifnot(
  factor_score$age == 45,
  factor_score$height_m == 1.55,
  factor_score$alm_kg == 15,
  abs(
    factor_score$model_averaged_z -
      example_scores$model_averaged_z[1]
  ) < 1e-12
)

empty_score <- score_conditional_alm(
  example_input[0, ],
  bundles = bundles,
  warn = FALSE
)
stopifnot(
  nrow(empty_score$pooled) == 0,
  nrow(empty_score$by_imputation) == 0
)

boundary_input <- data.frame(
  sex = c("Female", "Female", "Other", "Female"),
  age = c(45, 60, 45, 45),
  height_m = c(1.60, 1.60, 1.70, 1.40)
)
boundary_result <- suppressWarnings(reference_centiles(
  boundary_input,
  bundles = bundles
))
stopifnot(
  nrow(boundary_result) == nrow(boundary_input),
  identical(boundary_result$sex, boundary_input$sex),
  identical(
    boundary_result$centile_valid,
    c(TRUE, FALSE, FALSE, TRUE)
  ),
  identical(
    boundary_result$use_caution,
    c(FALSE, FALSE, FALSE, TRUE)
  ),
  all(is.na(boundary_result$alm_p50[c(2, 3)])),
  is.finite(boundary_result$alm_p50[1]),
  is.finite(boundary_result$alm_p50[4])
)

outside_result <- suppressWarnings(reference_centiles(
  data.frame(sex = "Female", age = 45, height_m = 2.20),
  bundles = bundles
))
stopifnot(
  !outside_result$centile_valid,
  is.na(outside_result$alm_p50)
)

bad_probability <- try(
  reference_centiles(
    data.frame(sex = "Female", age = 45, height_m = 1.60),
    probabilities = c(0.10, 1),
    bundles = bundles
  ),
  silent = TRUE
)
stopifnot(inherits(bad_probability, "try-error"))

bad_bundles <- try(
  score_conditional_alm(
    example_input,
    bundles = bundles[-1],
    warn = FALSE
  ),
  silent = TRUE
)
stopifnot(inherits(bad_bundles, "try-error"))

saved_grid <- read.csv(
  gzfile(file.path("data", "age_height_reference_grid.csv.gz")),
  check.names = FALSE
)
grid_input <- data.frame(
  sex = saved_grid$sex,
  age = saved_grid$age,
  height_m = saved_grid$height_cm / 100
)
recalculated_grid <- reference_centiles(
  grid_input,
  bundles = bundles,
  warn = FALSE
)
centile_columns <- grep("^alm_p", names(saved_grid), value = TRUE)
stopifnot(
  nrow(saved_grid) == 2840,
  nrow(recalculated_grid) == nrow(saved_grid),
  max(abs(
    as.matrix(recalculated_grid[centile_columns]) -
      as.matrix(saved_grid[centile_columns])
  )) < 1e-12
)

coefficient_components <- c(
  "mu.coefficients", "sigma.coefficients",
  "nu.coefficients", "tau.coefficients"
)
recalculated_coefficients <- do.call(
  rbind,
  lapply(seq_along(bundles), function(i) {
    bundle <- bundles[[i]]
    do.call(rbind, lapply(coefficient_components, function(component) {
      values <- bundle$model[[component]]
      if (is.null(values)) return(NULL)
      data.frame(
        model_key = names(bundles)[i],
        imputation = bundle$imputation,
        sex = bundle$sex,
        family = bundle$family,
        parameter_component = component,
        coefficient = names(values),
        value = as.numeric(values),
        stringsAsFactors = FALSE
      )
    }))
  })
)
saved_coefficients <- read.csv(
  file.path("models", "coefficient_snapshot.csv"),
  check.names = FALSE
)
key <- c(
  "model_key", "imputation", "sex", "family",
  "parameter_component", "coefficient"
)
coefficient_comparison <- merge(
  recalculated_coefficients,
  saved_coefficients,
  by = key,
  suffixes = c("_recalculated", "_saved")
)
stopifnot(
  nrow(recalculated_coefficients) == 55,
  nrow(saved_coefficients) == 55,
  nrow(coefficient_comparison) == 55,
  max(abs(
    coefficient_comparison$value_recalculated -
      coefficient_comparison$value_saved
  )) < 1e-12
)

message("All release tests passed.")
