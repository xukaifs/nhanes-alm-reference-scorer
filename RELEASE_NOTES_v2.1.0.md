# v2.1.0 — Age-only comparator reproducibility update

Release date: 2026-08-20

## Changes from v2.0.0

- Added `analysis/age_only_vs_age_height_comparator.R`.
- Added `analysis/README.md` describing the validation-data input contract.
- Added `analysis/expected/age_only_vs_age_height_expected.csv` containing the
  manuscript-level age-only versus age+height audit targets.
- Updated release metadata and integrity hashes.

## Unchanged from v2.0.0

- Primary age+height model specification.
- Female and male frozen RDS model objects.
- Public scoring functions and validity boundaries.
- Model manifest and primary fit audits.

This is therefore a reproducibility/documentation update for a new manuscript
sensitivity analysis, not a new primary scoring model. v2.0.0 remains archived
and should not be edited or deleted.
