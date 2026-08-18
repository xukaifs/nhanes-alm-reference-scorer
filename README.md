# Age- and height-conditioned ALM reference scorer

This repository contains the frozen **Version 2.0.0** research implementation
of sex-specific, age- and height-conditioned reference scoring for DXA-derived
appendicular lean mass (ALM) in U.S. adults aged **18–69 years**.

Version 2.0.0 corresponds to the final corrected-weight primary H model used in
the revised manuscript. It supersedes v1.0.0 for the current manuscript while
leaving the archived v1.0.0 release unchanged for reproducibility.

## Final model definition

- Development source: NHANES 1999–2006 completed DXA datasets.
- Development n: 16,018 (7,890 women; 8,128 men).
- Five completed DXA datasets are fitted separately for each sex.
- Women: BCCG; age df=3; height df=3; constant sigma.
- Men: BCT; age df=3; height df=3; constant sigma.
- Primary reference domain: ages 18–69 years.
- Fitting uses the corrected pooled MEC examination weights, centered within sex
  for numerical fitting.
- The public scoring target is absolute ALM conditional on sex, age, and stature;
  BMI is not a conditioning variable in the primary H model.

The supplied full-fit audit files show convergence for all five female BCCG
models and all five male BCT models.

## Frozen scoring definition

For each person, the scorer obtains a percentile and z score from each of the
five imputation-specific models of the same sex, then reports arithmetic model
averages separately on the two scales:

```text
model_averaged_z = mean(z_m)
model_averaged_percentile = mean(percentile_m)
```

These are distinct model-averaged summaries. Do not replace the formal z score
with `qnorm(model_averaged_percentile)` or the formal percentile with
`pnorm(model_averaged_z)`.

## Repository contents

- `R/score_conditional_alm.R` — scoring and centile functions.
- `models/corrected_H_18_69_Female_bundles.rds` — five final female models.
- `models/corrected_H_18_69_Male_bundles.rds` — five final male models.
- `models/model_manifest.csv` — final model specifications and fit diagnostics.
- `models/model_provenance.csv` — hashes of the frozen model artifacts.
- `models/audit/` — supplied full-fit audit tables.
- `logs/` — primary H and H/HB sensitivity fit/CV logs.
- `data/valid_input_ranges.csv` — conservative release guardrails.
- `examples/example_usage.R` — minimal scoring example.
- `tools/build_reference_grid.R` — regenerate an age-height centile grid directly
  from the v2.0.0 models.
- `tests/run_release_tests.R` — structural and live-scoring checks.
- `environment/` and `renv.lock` — recorded v1/final-project R environment.
- `CITATION.cff`, `LICENSE`, `DATA_LICENSE`, `PROVENANCE.md` — citation,
  licensing, and provenance metadata.
- `RELEASE_FILE_SHA256.csv` — SHA-256 manifest for release files.

### Why the v1.0.0 reference grid is not copied forward

The v1.0.0 precomputed grid and fixed worked-example outputs were based on the
older 20–59-year frozen models. They are intentionally **not** carried into
v2.0.0. A new grid should be generated from the final v2.0.0 models with
`tools/build_reference_grid.R`; this avoids mixing old derived numbers with new
model objects.

## Minimal use

Run from the repository root:

```r
source("R/score_conditional_alm.R")
check_alm_reference_environment(strict = FALSE)
models <- load_alm_reference()

patients <- data.frame(
  sex = c("Female", "Male"),
  age = c(45, 60),
  height_m = c(1.60, 1.75),
  alm_kg = c(15.0, 22.0)
)

score_conditional_alm(patients, bundles = models)$pooled
```

Create centiles at selected age-height points:

```r
reference_centiles(
  data.frame(
    sex = c("Female", "Male"),
    age = c(45, 60),
    height_m = c(1.60, 1.75)
  ),
  bundles = models
)
```

## Intended use and boundaries

- Research standardization only; this is not a stand-alone diagnostic tool.
- Final model age domain: 18–69 years.
- The release uses conservative height guardrails retained from v1.0.0; these
  ranges lie within the expanded development sample and reduce extrapolation.
- Values within the guard range but outside the recommended central range are
  scoreable but flagged for caution.
- ALM must be positive. The stored ALM caution range is conservative and is not
  a diagnostic boundary.
- The models were developed from U.S. NHANES Hologic DXA data.
- Calibration should be assessed before use with other scanner systems,
  software versions, countries, or clinical populations.
- Subgroup calibration offsets should not be interpreted as fixed biological
  norms.
- Conditional P5/P10 are population-reference positions, not outcome-derived
  sarcopenia decision limits.

## Reproducibility and public-use data note

The v2.0.0 RDS files in this assembly are the exact final corrected-weight fit
objects supplied from the analysis. Unlike the minimized v1.0.0 public bundle,
they were **not post-processed** in this environment. Their explicit training
frames retain public-use NHANES analysis fields used during fitting. See
`PROVENANCE.md` before public redistribution. The scoring results themselves do
not require users to access raw NHANES source files separately.

The recorded environment is R 4.2.1 with `gamlss` 5.4-22, `gamlss.dist` 6.1-1,
`gamlss.data` 6.0-7, and `nlme` 3.1-157. Run
`environment/verify_environment.R` for an exact version check.

## Citation and archive

Use `CITATION.cff` when citing the software. Create the GitHub tag/release
`v2.0.0`, allow Zenodo to archive that release, and then use the
version-specific Zenodo DOI for the manuscript's Code availability statement.

- Repository: https://github.com/xukaifs/nhanes-alm-reference-scorer
- Version: `2.0.0`
- Release date: 2026-08-18

## License

Code is released under the MIT License. Repository-authored documentation and
model metadata are covered by `DATA_LICENSE`. Original NHANES public-use terms
are not replaced by the repository license; see `PROVENANCE.md`.
