# Model provenance and data minimization

The fitted model objects were derived from the public-use NHANES 1999-2006
DXA completed datasets. The statistical fits were not changed for this
release.

Before public packaging, the explicit `train_data` component in every bundle
was reduced to the four fields required by the GAMLSS safe-prediction path:

- `alm_kg`
- `age`
- `height_m`
- `w_model`

The release bundles do not retain `SEQN`, race/ethnicity, or the raw
examination-weight column in `train_data`. The fitted GAMLSS objects still
contain the response and normalized case weights required to reproduce safe
predictions; they do not contain direct personal identifiers.

Numerical regression tests require:

- identical four-case worked-example scores, including the same-ALM female
  pair and near-identical-ALMI male pair;
- identical 2,840-row age-height reference grid;
- identical 55-row coefficient snapshot;
- no known safe-prediction warning exposed to the user;
- invalid centile requests retained as rows with `NA` centiles.
