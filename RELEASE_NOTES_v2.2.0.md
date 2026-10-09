# v2.2.0 — Paper 1B cross-population validation and calibration release

Release date: 2026-10-09

## Scope

Version 2.2.0 extends the frozen age- and height-conditioned ALM reference framework with the reproducibility materials used in Paper 1B. The primary U.S. scorer and the final 18–69-year frozen model objects are unchanged from v2.0.0; the age-only comparator resources introduced in v2.1.0 are retained.

## New in v2.2.0

- Frozen temporal validation in NHANES 2011–2018.
- External replication in KNHANES 2008–2011.
- Contemporary replication in KNHANES 2024.
- FNIH ALM/BMI and EWGSOP2 ALM/height² comparator analyses.
- Cross-population calibration analyses using the frozen U.S. reference.
- Locked-architecture Korean mirrored-refit summaries.
- Grip-strength contextual analyses.
- Aggregate result and verification files used to audit the Paper 1B manuscript.

## Interpretation boundary

This release distinguishes **structural transportability** from **percentile interchangeability**. Replication of reduced residual stature dependence does not imply that U.S.-derived percentile positions can be used directly as Korean population reference positions without local calibration.

## Data and frozen model objects

No participant-level NHANES or KNHANES data are redistributed. The existing GitHub repository already contains the verified final model bundles:

- `models/corrected_H_18_69_Female_bundles.rds`
- `models/corrected_H_18_69_Male_bundles.rds`

They are unchanged in v2.2.0 and should be preserved when merging this update package.
