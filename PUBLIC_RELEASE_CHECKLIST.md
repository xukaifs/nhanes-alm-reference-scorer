# Public release checklist for v2.1.0

Before creating the GitHub `v2.1.0` tag:

1. Extract this archive into the existing `nhanes-alm-reference-scorer` repository.
2. Confirm the primary v2.0.0 model RDS files and scorer were not altered.
3. In R 4.2.1 with the recorded packages restored, run `source("tests/run_release_tests.R")`.
4. For formal manuscript reproduction, run `analysis/age_only_vs_age_height_comparator.R`
   with the archived analysis-ready NHANES 2011-2018 validation file and compare
   the generated summary with `analysis/expected/age_only_vs_age_height_expected.csv`.
5. Confirm `CITATION.cff` says 2.1.0 and 2026-08-20.
6. Commit the repository state.
7. Create GitHub Release/tag `v2.1.0`. Do not edit or delete v1.0.0 or v2.0.0.
8. Let the existing Zenodo integration archive the new GitHub release.
9. Insert the new version-specific Zenodo DOI into the manuscript Code availability statement.
