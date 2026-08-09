options(stringsAsFactors = FALSE)

required_r <- "4.2.1"
required_packages <- c(
  gamlss = "5.4-22",
  gamlss.dist = "6.1-1",
  gamlss.data = "6.0-7",
  nlme = "3.1-157"
)

installed <- vapply(
  names(required_packages),
  function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageDescription(package)$Version)
  },
  character(1)
)

result <- data.frame(
  component = c("R", names(required_packages)),
  required_version = c(required_r, unname(required_packages)),
  installed_version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    unname(installed)
  ),
  stringsAsFactors = FALSE
)
result$matches <- result$required_version == result$installed_version
result$matches[is.na(result$matches)] <- FALSE

print(result, row.names = FALSE)
if (any(!result$matches)) {
  stop(
    "Environment mismatch. Reproduce the validated versions before using ",
    "the frozen fitted objects.",
    call. = FALSE
  )
}
message("Environment verification passed.")
