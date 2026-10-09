# Pre-release checklist for v2.2.0 (2026-10-09)

- [ ] Merge this package into the existing repository; do not replace the repository with an empty/new tree.
- [ ] Preserve `models/corrected_H_18_69_Female_bundles.rds` and `models/corrected_H_18_69_Male_bundles.rds` already present on `main`.
- [ ] Verify the two model hashes against `models/model_provenance.csv`.
- [ ] Confirm the previous public release is v2.1.0 and no remote tag/release already uses `v2.2.0`.
- [ ] Confirm `VERSION` = 2.2.0.
- [ ] Confirm `CITATION.cff` version = 2.2.0 and date-released = 2026-10-09.
- [ ] Review `README.md`, `RELEASE_NOTES_v2.2.0.md`, and `docs/REPRODUCIBILITY.md`.
- [ ] Confirm no participant-level NHANES or KNHANES data are included.
- [ ] Commit and push the merged repository state.
- [ ] Create GitHub tag/release `v2.2.0` titled `Age- and Height-Conditioned ALM Reference Scorer v2.2.0 — Cross-Population Validation and Calibration`.
- [ ] Allow Zenodo to archive the GitHub release as a new version.
- [ ] Copy the new version-specific Zenodo DOI into the final Paper 1B Data availability statement.
