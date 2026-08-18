options(stringsAsFactors = FALSE)

suppressWarnings(
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })
)

# Scoring domain for the final corrected-weight 18-69 models.
# Age bounds match the final development domain. Height and ALM guardrails are
# deliberately conservative: they retain the v1.0.0 20-59 observed and central
# reference ranges, which are guaranteed to lie within the expanded 18-69
# development sample. They prevent unsupported extrapolation until a fully
# regenerated 18-69 empirical boundary table is archived.
alm_reference_domain <- data.frame(
  sex = c("Female", "Male"),
  age_min = c(18, 18),
  age_max = c(69, 69),
  height_observed_min_m = c(1.316, 1.304),
  height_observed_max_m = c(1.868, 2.041),
  height_recommended_min_m = c(1.458, 1.570),
  height_recommended_max_m = c(1.781, 1.943),
  alm_observed_min_kg = c(7.8996, 11.2186),
  alm_observed_max_kg = c(45.1642, 64.6173),
  alm_recommended_min_kg = c(10.7093, 16.3025),
  alm_recommended_max_kg = c(33.3594, 46.2880),
  stringsAsFactors = FALSE
)

alm_reference_required_versions <- c(
  gamlss = "5.4-22",
  gamlss.dist = "6.1-1",
  gamlss.data = "6.0-7",
  nlme = "3.1-157"
)

check_alm_reference_environment <- function(strict = FALSE) {
  installed <- vapply(
    names(alm_reference_required_versions),
    function(package) {
      if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageDescription(package)$Version)
    },
    character(1)
  )
  mismatch <- is.na(installed) |
    installed != alm_reference_required_versions
  result <- data.frame(
    package = names(alm_reference_required_versions),
    required_version = unname(alm_reference_required_versions),
    installed_version = unname(installed),
    matches = !mismatch,
    stringsAsFactors = FALSE
  )
  if (any(mismatch)) {
    message <- paste0(
      "The frozen scorer was validated with: ",
      paste0(
        result$package,
        " ",
        result$required_version,
        collapse = ", "
      ),
      ". See environment/verify_environment.R."
    )
    if (strict) {
      stop(message, call. = FALSE)
    } else {
      warning(message, call. = FALSE)
    }
  }
  result
}

validate_reference_bundles <- function(bundles) {
  if (!is.list(bundles) || length(bundles) != 10) {
    stop(
      "bundles must contain exactly 10 fitted objects: five imputations ",
      "for each sex.",
      call. = FALSE
    )
  }
  required_training_columns <- c(
    "alm_kg", "age", "height_m", "w_model"
  )
  metadata <- lapply(seq_along(bundles), function(i) {
    bundle <- bundles[[i]]
    required_components <- c(
      "imputation", "sex", "family", "train_data", "model"
    )
    if (!is.list(bundle) ||
        !all(required_components %in% names(bundle))) {
      stop("Bundle ", i, " is incomplete.", call. = FALSE)
    }
    if (!all(required_training_columns %in% names(bundle$train_data))) {
      stop(
        "Bundle ", i, " lacks required training columns: ",
        paste(
          setdiff(
            required_training_columns,
            names(bundle$train_data)
          ),
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    if (!isTRUE(bundle$model$converged)) {
      stop("Bundle ", i, " is not marked as converged.", call. = FALSE)
    }
    data.frame(
      imputation = as.integer(bundle$imputation),
      sex = as.character(bundle$sex),
      family = as.character(bundle$family),
      stringsAsFactors = FALSE
    )
  })
  metadata <- do.call(rbind, metadata)
  expected <- expand.grid(
    imputation = 1:5,
    sex = c("Female", "Male"),
    stringsAsFactors = FALSE
  )
  observed_keys <- paste(metadata$imputation, metadata$sex)
  expected_keys <- paste(expected$imputation, expected$sex)
  if (!setequal(observed_keys, expected_keys) ||
      anyDuplicated(observed_keys)) {
    stop(
      "bundles must contain one model for every imputation (1-5) and sex.",
      call. = FALSE
    )
  }
  expected_family <- ifelse(
    metadata$sex == "Male", "BCT", "BCCG"
  )
  if (any(metadata$family != expected_family)) {
    stop(
      "Unexpected distribution family in bundles; expected male BCT and ",
      "female BCCG.",
      call. = FALSE
    )
  }
  invisible(metadata)
}

find_reference_model_files <- function(model_path = NULL) {
  if (!is.null(model_path)) {
    model_path <- as.character(model_path)
    if (length(model_path) == 1L && dir.exists(model_path)) {
      model_path <- file.path(
        model_path,
        c(
          "corrected_H_18_69_Female_bundles.rds",
          "corrected_H_18_69_Male_bundles.rds"
        )
      )
    }
    if (!all(file.exists(model_path))) {
      stop("One or more supplied reference-model files do not exist.", call. = FALSE)
    }
    return(normalizePath(model_path, winslash = "/", mustWork = TRUE))
  }

  defaults <- file.path(
    "models",
    c(
      "corrected_H_18_69_Female_bundles.rds",
      "corrected_H_18_69_Male_bundles.rds"
    )
  )
  if (!all(file.exists(defaults))) {
    stop(
      "Final 18-69 reference-model files were not found. Run from the ",
      "repository root or supply model_path explicitly.",
      call. = FALSE
    )
  }
  normalizePath(defaults, winslash = "/", mustWork = TRUE)
}

load_alm_reference <- function(model_path = NULL) {
  files <- find_reference_model_files(model_path)
  objects <- lapply(files, readRDS)
  bundles <- do.call(c, objects)
  names(bundles) <- vapply(
    bundles,
    function(bundle) paste0("imp", bundle$imputation, "_", bundle$sex),
    character(1)
  )
  validate_reference_bundles(bundles)
  bundles
}

run_safe_prediction <- function(expression) {
  withCallingHandlers(
    expression,
    warning = function(warning_condition) {
      normalized_message <- gsub(
        "\\s+",
        " ",
        trimws(conditionMessage(warning_condition))
      )
      known_warning <- startsWith(
        normalized_message,
        paste(
          "There is a discrepancy between the original and the re-fit",
          "used to achieve 'safe' predictions"
        )
      )
      if (known_warning) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

predict_reference_parameters <- function(bundle, new_data) {
  ans <- try(
    run_safe_prediction(
      predictAll(
        bundle$model,
        newdata = new_data,
        data = bundle$train_data,
        type = "response"
      )
    ),
    silent = TRUE
  )
  if (inherits(ans, "try-error")) {
    ans <- run_safe_prediction(
      predictAll(
        bundle$model,
        newdata = new_data,
        type = "response"
      )
    )
  }
  as.data.frame(ans)
}

reference_distribution_function <- function(prefix, bundle) {
  get(
    paste0(prefix, bundle$family),
    envir = asNamespace("gamlss.dist")
  )
}

distribution_arguments <- function(parameters) {
  output <- list(mu = parameters$mu, sigma = parameters$sigma)
  if ("nu" %in% names(parameters)) output$nu <- parameters$nu
  if ("tau" %in% names(parameters)) output$tau <- parameters$tau
  output
}

reference_cdf <- function(bundle, alm_kg, parameters) {
  arguments <- c(list(q = alm_kg), distribution_arguments(parameters))
  as.numeric(
    do.call(reference_distribution_function("p", bundle), arguments)
  )
}

reference_quantile <- function(bundle, probability, parameters) {
  arguments <- c(list(p = probability), distribution_arguments(parameters))
  as.numeric(
    do.call(reference_distribution_function("q", bundle), arguments)
  )
}

validate_alm_input <- function(input) {
  required <- c("sex", "age", "height_m", "alm_kg")
  missing_columns <- setdiff(required, names(input))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  input <- as.data.frame(input)
  input$.row_id <- seq_len(nrow(input))
  input$sex <- as.character(input$sex)
  # Character conversion first prevents numeric factors from being scored
  # by their internal level codes.
  input$age <- suppressWarnings(as.numeric(as.character(input$age)))
  input$height_m <- suppressWarnings(
    as.numeric(as.character(input$height_m))
  )
  input$alm_kg <- suppressWarnings(as.numeric(as.character(input$alm_kg)))
  input$almi_kg_m2 <- ifelse(
    is.finite(input$height_m) & input$height_m > 0,
    input$alm_kg / input$height_m^2,
    NA_real_
  )

  domain <- alm_reference_domain[
    match(input$sex, alm_reference_domain$sex),
    ,
    drop = FALSE
  ]
  input$valid_sex <- input$sex %in% alm_reference_domain$sex
  input$valid_age <- input$valid_sex &
    is.finite(input$age) &
    input$age >= domain$age_min &
    input$age <= domain$age_max
  input$height_within_observed_range <- input$valid_sex &
    is.finite(input$height_m) &
    input$height_m >= domain$height_observed_min_m &
    input$height_m <= domain$height_observed_max_m
  input$height_within_recommended_range <- input$valid_sex &
    is.finite(input$height_m) &
    input$height_m >= domain$height_recommended_min_m &
    input$height_m <= domain$height_recommended_max_m
  input$alm_positive <- is.finite(input$alm_kg) & input$alm_kg > 0
  input$alm_within_observed_range <- input$valid_sex &
    is.finite(input$alm_kg) &
    input$alm_kg >= domain$alm_observed_min_kg &
    input$alm_kg <= domain$alm_observed_max_kg
  input$alm_within_recommended_range <- input$valid_sex &
    is.finite(input$alm_kg) &
    input$alm_kg >= domain$alm_recommended_min_kg &
    input$alm_kg <= domain$alm_recommended_max_kg
  input$score_valid <- input$valid_sex & input$valid_age &
    input$height_within_observed_range & input$alm_positive
  input$use_caution <- input$score_valid &
    (!input$height_within_recommended_range |
       !input$alm_within_recommended_range)
  input
}

# Score age- and height-conditioned ALM.
#
# Required input columns:
#   sex: "Male" or "Female"
#   age: completed years (18-69)
#   height_m: standing height in metres
#   alm_kg: DXA appendicular lean mass in kilograms
#
# Model-averaged values are arithmetic means across the five
# imputation-specific reference models:
#   model_averaged_z = mean(z_m)
#   model_averaged_percentile = mean(percentile_m)
# These two quantities are averaged separately; neither is reconstructed by
# applying qnorm() or pnorm() to the other. They are not Rubin-pooled
# quantities. Percentiles are returned on both the 0-1 and 0-100 scales.
score_conditional_alm <- function(
    new_data,
    bundles = load_alm_reference(),
    warn = TRUE
) {
  validate_reference_bundles(bundles)
  input <- validate_alm_input(new_data)

  if (warn && any(!input$score_valid)) {
    warning(
      sum(!input$score_valid),
      " row(s) had invalid inputs or were outside the observed age/height ",
      "domain; scores were returned as NA.",
      call. = FALSE
    )
  }
  if (warn && any(input$use_caution, na.rm = TRUE)) {
    warning(
      sum(input$use_caution, na.rm = TRUE),
      " row(s) were scoreable but outside the recommended-use range; ",
      "interpret as extrapolation-prone.",
      call. = FALSE
    )
  }

  if (nrow(input) == 0) {
    input$model_averaged_percentile <- numeric()
    input$model_averaged_percentile_pct <- numeric()
    input$model_averaged_z <- numeric()
    input$mean_percentile <- numeric()
    input$mean_percentile_pct <- numeric()
    input$mean_z <- numeric()
    input$.row_id <- NULL
    return(list(
      pooled = input,
      by_imputation = data.frame(
        .row_id = integer(),
        imputation = integer(),
        sex = character(),
        percentile = numeric(),
        percentile_pct = numeric(),
        z = numeric(),
        model_key = character(),
        stringsAsFactors = FALSE
      ),
      domain = alm_reference_domain
    ))
  }

  scored <- vector("list", length(bundles))
  for (i in seq_along(bundles)) {
    bundle <- bundles[[i]]
    same_sex <- !is.na(input$sex) & input$sex == bundle$sex
    eligible <- input$score_valid & same_sex
    percentile <- rep(NA_real_, nrow(input))

    if (any(eligible)) {
      prediction_data <- input[
        eligible,
        c("age", "height_m", "alm_kg"),
        drop = FALSE
      ]
      parameters <- predict_reference_parameters(bundle, prediction_data)
      probability <- reference_cdf(
        bundle,
        prediction_data$alm_kg,
        parameters
      )
      percentile[eligible] <- pmin(pmax(probability, 1e-10), 1 - 1e-10)
    }

    model_key <- names(bundles)[i]
    if (is.null(model_key) || is.na(model_key) || !nzchar(model_key)) {
      model_key <- paste0(
        "imp", bundle$imputation, "_", bundle$sex
      )
    }
    n_same_sex <- sum(same_sex)
    scored[[i]] <- data.frame(
      .row_id = input$.row_id[same_sex],
      imputation = rep.int(as.integer(bundle$imputation), n_same_sex),
      sex = input$sex[same_sex],
      percentile = percentile[same_sex],
      percentile_pct = 100 * percentile[same_sex],
      z = qnorm(percentile[same_sex]),
      model_key = rep.int(model_key, n_same_sex),
      stringsAsFactors = FALSE
    )
  }

  by_imputation <- do.call(rbind, scored)
  if (nrow(by_imputation) == 0) {
    input$model_averaged_percentile <- NA_real_
    input$model_averaged_percentile_pct <- NA_real_
    input$model_averaged_z <- NA_real_
    input$mean_percentile <- input$model_averaged_percentile
    input$mean_percentile_pct <- input$model_averaged_percentile_pct
    input$mean_z <- input$model_averaged_z
    input$.row_id <- NULL
    return(list(
      pooled = input,
      by_imputation = by_imputation,
      domain = alm_reference_domain
    ))
  }
  pooled <- aggregate(
    cbind(percentile, percentile_pct, z) ~ .row_id + sex,
    data = by_imputation,
    FUN = function(x) {
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    },
    na.action = na.pass
  )
  names(pooled)[names(pooled) == "percentile"] <-
    "model_averaged_percentile"
  names(pooled)[names(pooled) == "percentile_pct"] <-
    "model_averaged_percentile_pct"
  names(pooled)[names(pooled) == "z"] <- "model_averaged_z"

  pooled <- merge(
    input,
    pooled,
    by = c(".row_id", "sex"),
    all.x = TRUE,
    sort = FALSE
  )
  pooled <- pooled[order(pooled$.row_id), ]
  pooled$.row_id <- NULL
  # Backward-compatible aliases retained for existing analysis scripts.
  pooled$mean_percentile <- pooled$model_averaged_percentile
  pooled$mean_percentile_pct <-
    pooled$model_averaged_percentile_pct
  pooled$mean_z <- pooled$model_averaged_z

  list(
    pooled = pooled,
    by_imputation = by_imputation,
    domain = alm_reference_domain
  )
}

# Return sex-specific conditional ALM centiles for an age-height grid.
reference_centiles <- function(
    new_data,
    probabilities = c(0.025, 0.05, 0.10, 0.25, 0.50,
                      0.75, 0.90, 0.95, 0.975),
    bundles = load_alm_reference(),
    warn = TRUE
) {
  validate_reference_bundles(bundles)
  required <- c("sex", "age", "height_m")
  missing_columns <- setdiff(required, names(new_data))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) == 0 ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1) ||
      anyDuplicated(probabilities)) {
    stop(
      "probabilities must be unique, finite values strictly between 0 and 1.",
      call. = FALSE
    )
  }

  grid <- as.data.frame(new_data)
  grid$.row_id <- seq_len(nrow(grid))
  grid$sex <- as.character(grid$sex)
  grid$age <- suppressWarnings(as.numeric(as.character(grid$age)))
  grid$height_m <- suppressWarnings(
    as.numeric(as.character(grid$height_m))
  )

  domain <- alm_reference_domain[
    match(grid$sex, alm_reference_domain$sex),
    ,
    drop = FALSE
  ]
  grid$valid_sex <- grid$sex %in% alm_reference_domain$sex
  grid$valid_age <- grid$valid_sex &
    is.finite(grid$age) &
    grid$age >= domain$age_min &
    grid$age <= domain$age_max
  grid$height_within_observed_range <- grid$valid_sex &
    is.finite(grid$height_m) &
    grid$height_m >= domain$height_observed_min_m &
    grid$height_m <= domain$height_observed_max_m
  grid$height_within_recommended_range <- grid$valid_sex &
    is.finite(grid$height_m) &
    grid$height_m >= domain$height_recommended_min_m &
    grid$height_m <= domain$height_recommended_max_m
  grid$centile_valid <- grid$valid_sex &
    grid$valid_age &
    grid$height_within_observed_range
  grid$use_caution <- grid$centile_valid &
    !grid$height_within_recommended_range

  if (warn && any(!grid$centile_valid)) {
    warning(
      sum(!grid$centile_valid),
      " row(s) had an invalid sex/age or were outside the observed height ",
      "domain; centiles were returned as NA.",
      call. = FALSE
    )
  }
  if (warn && any(grid$use_caution)) {
    warning(
      sum(grid$use_caution),
      " row(s) were within the observed height domain but outside the ",
      "recommended-use range; interpret as extrapolation-prone.",
      call. = FALSE
    )
  }

  probability_labels <- vapply(
    100 * probabilities,
    function(value) {
      if (abs(value - round(value)) < 1e-10) {
        as.character(round(value))
      } else {
        gsub(
          "\\.",
          "_",
          format(value, trim = TRUE, scientific = FALSE)
        )
      }
    },
    character(1)
  )
  probability_names <- paste0("alm_p", probability_labels)
  if (nrow(grid) == 0) {
    for (column in probability_names) grid[[column]] <- numeric()
    grid$.row_id <- NULL
    return(grid)
  }

  centile_rows <- list()
  counter <- 1L
  for (i in seq_along(bundles)) {
    bundle <- bundles[[i]]
    eligible <- grid$centile_valid &
      !is.na(grid$sex) &
      grid$sex == bundle$sex
    if (!any(eligible)) next
    prediction_data <- grid[eligible, c("age", "height_m"), drop = FALSE]
    parameters <- predict_reference_parameters(bundle, prediction_data)
    for (probability in probabilities) {
      centile_rows[[counter]] <- data.frame(
        .row_id = grid$.row_id[eligible],
        imputation = bundle$imputation,
        sex = bundle$sex,
        probability = probability,
        alm_kg = reference_quantile(bundle, probability, parameters),
        stringsAsFactors = FALSE
      )
      counter <- counter + 1L
    }
  }

  if (length(centile_rows) == 0) {
    for (column in probability_names) grid[[column]] <- NA_real_
    grid$.row_id <- NULL
    return(grid)
  }
  by_imputation <- do.call(rbind, centile_rows)
  # Quantile curves are averaged across the five model-specific quantiles at
  # each probability and age-height point. This is model averaging, not
  # Rubin pooling.
  pooled <- aggregate(
    alm_kg ~ .row_id + sex + probability,
    data = by_imputation,
    FUN = mean
  )
  wide <- reshape(
    pooled,
    idvar = c(".row_id", "sex"),
    timevar = "probability",
    direction = "wide"
  )
  centile_columns <- grep("^alm_kg\\.", names(wide))
  names(wide)[centile_columns] <- probability_names
  wide$sex <- NULL
  output <- merge(
    grid,
    wide,
    by = ".row_id",
    all.x = TRUE,
    sort = FALSE
  )
  for (column in probability_names) {
    if (!column %in% names(output)) output[[column]] <- NA_real_
  }
  output <- output[order(output$.row_id), ]
  output$.row_id <- NULL
  output
}
