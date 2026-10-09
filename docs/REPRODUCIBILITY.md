# Reproducibility notes

## 1. Frozen U.S. scorer

Paper 1B uses the final corrected-weight 18–69-year age+height reference. The two verified binary model bundles are required for live scoring. Their expected filenames and hashes are listed in the top-level README and `models/model_provenance.csv`.

## 2. Archived analysis scripts

Scripts in `analysis/paper1b/` are preserved from the working analysis directories. Several contain historical absolute Windows paths and custom `.libPaths()` settings. Those paths must be replaced or parameterized locally before rerunning. They are preserved rather than silently edited so that the public archive matches the supplied analysis source.

## 3. Korean analyses

For KNHANES 2008–2011, the official pooled-weight outputs under `results/KNHANES/tables_official_weight/` are the primary historical Korean results. Equal-year-weight and uncorrected-ALM outputs are sensitivity/audit materials and should not replace the primary results. KNHANES 2024 outputs are under `results/KNHANES/tables_2024/`. Cross-population mirrored-refit/surface outputs are under `results/KNHANES/us_korea_surface/`.

## 4. U.S. analyses

Final frozen U.S. validation outputs are under `results/NHANES/Paper1_FINAL_FROZEN_20260814/`. The controlled age-only comparison is under `results/NHANES/age_only_comparator_20260820/`.

## 5. Privacy and redistribution

No raw participant-level files should be committed. Before future releases, inspect any RDS/RData/cache objects for embedded training or participant-level frames.

## 6. Environment

The archived scorer environment records R 4.2.1 with `gamlss` 5.4-22, `gamlss.dist` 6.1-1, `gamlss.data` 6.0-7, and `nlme` 3.1-157. See `environment/` and `renv.lock`.
