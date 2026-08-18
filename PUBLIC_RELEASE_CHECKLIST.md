# Public release checklist for v2.0.0

Before creating the GitHub `v2.0.0` tag:

1. Extract this archive into the existing `nhanes-alm-reference-scorer` repository, replacing the current main-branch contents intentionally.
2. In R 4.2.1 with the recorded packages restored, run `source("tests/run_release_tests.R")`.
3. Optional but recommended for data minimization: run `source("tools/minimize_model_bundles.R")`, validate the resulting public RDS files, update the scorer paths if those files will replace the full final objects, rerun tests, and regenerate SHA-256 hashes.
4. If you want a precomputed v2 grid, run `source("tools/build_reference_grid.R")`; do not reuse the v1.0.0 grid.
5. Confirm `CITATION.cff` says 2.0.0 and 2026-08-18.
6. Commit the repository state.
7. Create GitHub Release/tag `v2.0.0`. Do not edit or delete v1.0.0.
8. Let the existing Zenodo integration archive the new GitHub release.
9. Insert the new version-specific Zenodo DOI into the manuscript Code availability statement.
