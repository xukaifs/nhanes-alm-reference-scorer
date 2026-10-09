# Age- and height-conditioned ALM reference scorer — v2.2.0

**Release date: 2026-10-09**

This repository accompanies Paper 1B on DXA-derived appendicular lean mass (ALM) normalization, frozen temporal validation, and cross-population replication in U.S. NHANES and Korean KNHANES.

## What v2.2.0 adds

Version 2.2.0 extends the existing public repository without changing the primary frozen U.S. age- and height-conditioned scorer introduced in v2.0.0. It also retains the age-only comparator reproducibility resources added in v2.1.0.

New Paper 1B reproducibility resources include:

- frozen temporal validation in NHANES 2011–2018;
- external replication in KNHANES 2008–2011 and KNHANES 2024;
- FNIH ALM/BMI and EWGSOP2 ALM/height² comparator analyses;
- cross-population calibration summaries and Korean mirrored-refit analyses;
- grip-strength contextual analyses;
- aggregate manuscript-verification outputs and audit files.

The central methodological distinction is between **structural transportability** (whether stature-related behavior of the frozen model replicates) and **percentile calibration/interchangeability** (whether U.S.-derived percentile positions can be transferred directly to Korea).

## Version lineage

- **v2.0.0** — final frozen 18–69-year age+height scorer and model objects.
- **v2.1.0** — age-only versus age+height comparator reproducibility resources.
- **v2.2.0** — Paper 1B temporal validation, Korean external replication, calibration, mirrored-refit, and grip-strength resources.

## Primary U.S. reference specification

- Development data: NHANES 1999–2006 completed DXA datasets.
- Development sample: 16,018 adults (7,890 women; 8,128 men), ages 18–69 years.
- Women: BCCG; age df=3; height df=3; constant sigma.
- Men: BCT; age df=3; height df=3; constant sigma.
- Five completed DXA datasets are fitted separately by sex.
- The scoring target is absolute ALM conditional on sex, age, and stature; BMI is not a conditioning variable in the primary model.
- Conditional P5/P10 are population-reference positions, not stand-alone sarcopenia diagnoses or outcome-derived decision limits.

## Repository layout

```text
R/                              scoring functions
models/                         frozen model files, metadata, provenance, and audit tables
analysis/paper1b/NHANES/        archived U.S. analysis scripts
analysis/paper1b/KNHANES/       archived Korean analysis scripts
analysis/reference/             controlled age-only comparator materials
results/NHANES/                 aggregate U.S. outputs
results/KNHANES/                aggregate Korean outputs
environment/                    recorded R environment information
examples/                       minimal scorer usage example
tests/                          scorer/release checks
tools/                          helper scripts
docs/                           data access and reproducibility notes
```

## Important: preserve the existing frozen model files

The public GitHub `main` branch already contains the two final 18–69-year frozen GAMLSS model bundles. They are unchanged in v2.2.0 and should **not** be deleted or replaced when merging this update package:

- `models/corrected_H_18_69_Female_bundles.rds`
- `models/corrected_H_18_69_Male_bundles.rds`

Expected SHA-256 values recorded in `models/model_provenance.csv`:

- Female: `B2FA1530C23DC16C2F7D48D4A56D533C3B6B49CDC76E6B45D74BC56008E095BD`
- Male: `A11A01B43D5C4C4922566990144612D814927D1BC6518ED4DCA7082F4FE6C8BA`

This handoff ZIP does not duplicate those large binary files. Merge the extracted contents into the existing repository and preserve the current `models/*.rds` files.

## Reproducibility

The analysis scripts are archived from the analysis workspace to preserve provenance. Several scripts contain historical local Windows paths and R-library paths. They are not silently rewritten here because doing so would create unvalidated code changes. Configure local paths before rerunning. See `docs/REPRODUCIBILITY.md`.

Participant-level NHANES and KNHANES data are not redistributed. Users should obtain the public-use source data from the official NCHS and KDCA portals under their respective terms of use. Only source code, model metadata, and aggregate/verification outputs are included here.

## Minimal scorer use

With the existing verified model bundles retained under `models/`:

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

## Intended use

This software is for research standardization and reproducibility. It is not a stand-alone diagnostic tool. Calibration should be assessed before use with other scanner systems, software versions, countries, or clinical populations.

## Citation and archive

Repository: https://github.com/xukaifs/nhanes-alm-reference-scorer

Version: `2.2.0`

Release date: **2026-10-09**

After the GitHub `v2.2.0` release is archived by Zenodo, use the **version-specific Zenodo DOI** in the Paper 1B manuscript. Do not insert a DOI until Zenodo has created the v2.2.0 record.

## License

Code is released under the MIT License. Repository-authored documentation and model metadata are covered by `DATA_LICENSE`. Original NHANES and KNHANES data remain subject to their source terms and are not redistributed here.
