# Age- and height-conditioned ALM reference scorer

This repository contains the frozen research implementation of sex-specific,
age- and height-conditioned reference centiles for DXA-derived appendicular
lean mass (ALM) in U.S. adults aged 20–59 years.

## Frozen scoring definition

Five NHANES DXA completed-dataset models were fitted separately. For an
individual, the formal summaries are:

```text
model_averaged_z = mean(z_m)
model_averaged_percentile = mean(percentile_m)
```

These quantities are averaged separately on two nonlinear scales and should
be interpreted as distinct model-averaged summaries. Do not replace the
formal z score with `qnorm(model_averaged_percentile)` and do not replace the
formal percentile with `pnorm(model_averaged_z)`.

In the 10,830-participant temporal-validation cohort, the numerical audit
found no P5 classification discordance and one P10 discordance for a
participant at the threshold. Typical numerical differences between the two
representations were negligible.

## Files

- `R/score_conditional_alm.R`: scoring and centile functions.
- `models/revised_primary_model_bundles.rds`: frozen fitted model objects.
- `models/model_manifest.csv`: family, formula, domain, and model metadata.
- `models/coefficient_snapshot.csv`: coefficient snapshot for audit.
- `models/model_provenance.csv`: full-source and minimized-public hashes.
- `data/age_height_reference_grid.csv.gz`: age–height reference grid.
- `data/valid_input_ranges.csv`: observed and recommended-use boundaries.
- `data/clinical_calculation_examples.csv`: worked cases.
- `examples/example_usage.R`: minimal R example.
- `environment/session_info.txt`: software and package versions.
- `environment/verify_environment.R`: exact-version guard.
- `environment/package_versions.csv`: machine-readable version requirements.
- `renv.lock`: reproducible R package environment specification.
- `tests/run_release_tests.R`: scoring, grid, single-record, invalid-input,
  factor-input, boundary, and coefficient regression tests.
- `CITATION.cff`: machine-readable software citation metadata.
- `LICENSE` and `DATA_LICENSE`: code and data/model license terms.
- `RELEASE_FILE_SHA256.csv`: SHA-256 manifest for all other release files.

The public model bundles retain only `alm_kg`, `age`, `height_m`, and
`w_model` in their explicit training-data component. Direct NHANES
participant identifiers, race/ethnicity, and raw examination-weight columns
were removed from that component without changing any reference score or
centile.

## Validated software environment

The frozen scorer was tested with R 4.2.1, `gamlss` 5.4-22,
`gamlss.dist` 6.1-1, `gamlss.data` 6.0-7, and `nlme` 3.1-157. Run:

```r
source("environment/verify_environment.R")
```

for exact release validation. Routine users can run
`check_alm_reference_environment(strict = FALSE)` after sourcing the scorer;
version differences are then reported as a warning. The prediction helper
filters only the known GAMLSS safe-prediction discrepancy warning that was
quantified in the release audit; all other warnings remain visible.

Users of `renv` can restore the recorded package versions from the repository
root with `renv::restore()`. R 4.2.1 itself must be installed separately.

## Minimal use

To run the bundled worked example from the repository root, use:

```r
source("examples/example_usage.R")
```

The example script sources `R/score_conditional_alm.R` itself, so the two
scripts do not need to be sourced consecutively. For custom data, use:

```r
source("R/score_conditional_alm.R")
check_alm_reference_environment(strict = FALSE)
models <- readRDS("models/revised_primary_model_bundles.rds")

patients <- data.frame(
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

score_conditional_alm(patients, bundles = models)$pooled
```

The scorer accepts single records, single-sex batches, and mixed-sex
batches. Rows with invalid sex, age, height, or ALM are retained and return
`NA` scores together with the input-validity and range flags.

Rounded display values:

| Contrast | Sex | Age | Height, cm | ALM, kg | ALMI, kg/m² | Model-averaged z | Model-averaged percentile, % |
|---|---|---:|---:|---:|---:|---:|---:|
| Identical ALM: shorter stature | Female | 45 | 155.0 | 15.00 | 6.24 | -0.346 | 36.5 |
| Identical ALM: taller stature | Female | 45 | 175.0 | 15.00 | 4.90 | -1.915 | 2.8 |
| Near-identical ALMI: younger/taller | Male | 38 | 179.0 | 22.40 | 6.99 | -1.655 | 4.9 |
| Near-identical ALMI: older/shorter | Male | 58 | 159.6 | 17.82 | 7.00 | -0.880 | 18.9 |

The first pair fixes age and ALM while varying height. The second pair has
nearly identical ALMI but different age and height. They are hypothetical
algorithm examples, not identifiable participants or diagnostic
classifications. Full-precision outputs are stored in
`data/clinical_calculation_examples.csv`.

## Intended use and boundaries

- Research standardization only; this is not a stand-alone diagnostic tool.
- Valid ages: 20–59 years.
- Do not extrapolate beyond the sex-specific observed height ranges.
- Recommended height limits are the empirical (unweighted) 1st–99th
  percentiles of development-sample height.
- Recommended ALM limits are the empirical 0.5th–99.5th percentiles after
  stacking the five completed development datasets.
- Interpret in-range extreme stature and ALM with caution.
- The models were developed from NHANES Hologic DXA data.
- Recalibration should be assessed before use with other manufacturers,
  scanner models, software versions, countries, or clinical populations.
- Race/ethnicity subgroup offsets were substantial; a common score must not
  be described as equally calibrated in all subgroups.

## Reproducibility note

Displayed centile curves are arithmetic means of the five model-specific
quantiles at the same probability and age–height point. Survey estimates in
the manuscript use their stated complex-survey and multiple-imputation
procedures; the individual scoring averages above are not Rubin-pooled
quantities.

## Citation and archive

Use the metadata in `CITATION.cff` when citing this software. For an archived
release, use the version-specific DOI displayed on the corresponding Zenodo
record. The accompanying manuscript should also be cited when available.

- GitHub repository:
  <https://github.com/xukaifs/nhanes-alm-reference-scorer>
- Version: `1.0.0`
- Release date: 2026-08-06

## License

Code in `R/`, `examples/`, `tests/`, and the R scripts in `environment/` is
licensed under the MIT License; see `LICENSE`. The fitted model objects,
reference grids, worked-example outputs, documentation, and other data/model
artifacts are licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0); see `DATA_LICENSE`. NHANES source data are
not redistributed by this repository and remain subject to their original
public-use terms.
