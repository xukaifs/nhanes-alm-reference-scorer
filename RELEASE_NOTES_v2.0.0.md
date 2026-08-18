# v2.0.0 — Final corrected-weight 18–69-year reference models

Release date: 2026-08-18

## Major changes from v1.0.0

- Primary age domain expanded from 20–59 to 18–69 years.
- Final development sample: 16,018 adults (7,890 women; 8,128 men).
- Replaced the historical model objects with the final corrected pooled-weight H fits.
- Women remain BCCG; men remain BCT; age and height df=3; sigma constant.
- Model files are split by sex, each containing five completed-dataset fits.
- Removed stale v1.0.0 precomputed centile grid and fixed worked-example outputs.
- Added a v2 reference-grid regeneration script.
- Added final fit-audit tables and primary/sensitivity fit logs.
- Updated scorer age validity to 18–69 years.
- Preserved conservative stature guardrails from the prior release to avoid unsupported extrapolation.

This release supersedes v1.0.0 for the current manuscript. The v1.0.0 tag and Zenodo record should remain unchanged.
