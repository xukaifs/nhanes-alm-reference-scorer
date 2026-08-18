# Model provenance and release status

Version 2.0.0 uses the final corrected-weight age-height (H) models developed
from NHANES 1999–2006 completed DXA datasets for adults aged 18–69 years.

## Primary frozen objects

- `models/corrected_H_18_69_Female_bundles.rds`: five female BCCG models.
- `models/corrected_H_18_69_Male_bundles.rds`: five male BCT models.

The accompanying full-fit audit tables report n=7,890 per female imputation and
n=8,128 per male imputation, with all ten final fits converged. The primary
models use age df=3, height df=3, and constant sigma.

## Corrected development weighting

The final analysis used cycle-correct pooled MEC examination weights and
centered those weights within sex for GAMLSS fitting. This v2.0.0 release is
therefore not numerically interchangeable with the historical v1.0.0 20–59
model release.

## Exact-object preservation

The RDS files in this assembly are copied byte-for-byte from the final supplied
analysis objects. They were not re-fitted, coefficient-edited, or serialized
again during package assembly. Their SHA-256 hashes are recorded in
`models/model_provenance.csv` and `RELEASE_FILE_SHA256.csv`.

## Public-use training-frame note

Unlike the deliberately minimized v1.0.0 public bundle, these supplied final
RDS objects retain the analysis `train_data` frame used by GAMLSS safe
prediction. Inspection confirms that this frame contains public-use NHANES
analysis fields such as `SEQN`, `cycle`, `strata`, and `psu`, in addition to
ALM, age, height, and model weights. No claim is made that v2.0.0 has been
post-fit data-minimized.

If a stricter minimized public artifact is desired before tagging the GitHub
release, run `tools/minimize_model_bundles.R` in the recorded R environment and
then rerun `tests/run_release_tests.R` and regenerate the SHA-256 manifest.
Because package assembly here does not include an R runtime, the original final
RDS bytes are preserved rather than altered without validation.

## Derived files

The v1.0.0 20–59 precomputed reference grid and worked-example outputs are not
carried forward. `tools/build_reference_grid.R` regenerates a v2-compatible grid
directly from the frozen 18–69 models.
