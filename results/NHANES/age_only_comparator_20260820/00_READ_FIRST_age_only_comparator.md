# Age-only versus age+height ALM reference sensitivity analysis

Date: 2026-08-20  
Status: completed, frozen, and integrated into the revised manuscript

## Question

Does an age-only ALM reference retain clinically important stature dependence that is reduced by explicitly conditioning on height?

This analysis was designed to isolate the height term. It did **not** reselect a family or tune a model against the 2011–2018 cohort.

## Fixed comparison

- Development: NHANES 1999–2006, ages 18–69.
- Temporal evaluation: NHANES 2011–2018, ages 18–59.
- Sex-specific family, age smooth, scale/shape structure, corrected pooled development weights, five completed DXA datasets, and deterministic PSU folds were identical between models.
- Age-only removed height from the final location predictor; age+height was the frozen primary H reference.
- All five full fits and all 25 PSU-CV fits per sex converged.

## Main findings

### 1. Overall calibration alone did not reveal the problem

Out-of-fold age-only calibration was superficially good:

| Sex | Mean z | SD(z) | P5 | P10 |
|---|---:|---:|---:|---:|
| Women | −0.0005 | 1.0015 | 4.63% | 9.47% |
| Men | −0.0001 | 1.0014 | 4.87% | 10.06% |

Therefore, near-nominal overall P5/P10 coverage does not establish stature stability.

### 2. Age-only retained pronounced temporal stature gradients

Adjusted height effects per 10-cm greater stature:

| Sex | Outcome | Age-only | Age+height |
|---|---|---:|---:|
| Women | Continuous z beta | +0.644 (0.602 to 0.687) | +0.018 (−0.028 to 0.063) |
| Women | P5 OR | 0.171 (0.137 to 0.214) | 0.619 (0.477 to 0.804) |
| Women | P10 OR | 0.166 (0.133 to 0.206) | 0.895 (0.746 to 1.073) |
| Men | Continuous z beta | +0.708 (0.666 to 0.750) | −0.052 (−0.103 to −0.002) |
| Men | P5 OR | 0.163 (0.123 to 0.215) | 1.014 (0.814 to 1.263) |
| Men | P10 OR | 0.190 (0.161 to 0.224) | 1.047 (0.875 to 1.253) |

Explicit height conditioning reduced the absolute continuous height gradient by 97.3% in women and 92.6% in men.

### 3. Fit favored the age+height reference

Summed BIC across the five completed datasets was higher for age-only by 7,912.8 points in women and 14,869.1 points in men.

### 4. The reference choice materially changed lower-tail assignment

Survey-weighted temporal reclassification:

| Sex | P5 | P10 |
|---|---:|---:|
| Women | 2.80% (2.38% to 3.22%) | 6.33% (5.54% to 7.11%) |
| Men | 4.42% (3.73% to 5.10%) | 7.73% (6.64% to 8.82%) |

The score correlations were 0.909 in women and 0.848 in men; mean absolute score differences were 0.350 and 0.446 z, respectively.

## Interpretation for the paper

This is a strong sensitivity analysis because it directly answers the contribution question raised by prior age-specific centile work:

- Age-only centiles describe life-course distributions.
- Age+height centiles address stature-conditioned reference position.
- An age-only model can look well calibrated in the aggregate while still assigning systematically lower reference positions to shorter adults.
- The result supports height conditioning without claiming perfect invariance: a residual female P5 gradient and a small male continuous-score gradient remain.

Recommended manuscript placement: one concise paragraph in Results and Discussion, with the full comparison in Supplementary Table S15 and Supplementary Figure S3.

## Files

- `01_model_fit_detail.csv`: all five full model fits per sex and model.
- `01b_model_fit_summary.csv`: convergence, BIC, and model specifications.
- `02_internal_OOF_calibration.csv`: pooled PSU out-of-fold calibration.
- `03_temporal_calibration.csv`: frozen temporal calibration.
- `04_temporal_height_effects.csv`: adjusted linear and spline height diagnostics.
- `05_temporal_adjusted_height_predictions.csv`: marginal predictions used in the figure.
- `06_temporal_reclassification.csv`: P5/P10 cross-classification and score differences.
- `07_manuscript_Supplementary_Table_S15.csv`: manuscript-ready compact table.
- `figures/Supplementary_Figure_S3_age_only_vs_age_height.png`: manuscript figure.
- `99_session_info.txt`: computational environment.

## Manuscript integration and QA

The revised manuscript contains the new comparator in the Abstract, Introduction, Methods, Results, Discussion, Conclusions, Supplementary Table S15, and Supplementary Figure S3. Wang et al. 2024 was inserted as reference 14 and subsequent references/citations were renumbered. Supplementary Table S6 was also updated to the final corrected-weight continuous-height estimates.

Package QA passed: 19 tables, 7 inline figures, 38 numbered references, successful reopen after save, and correct end-of-document order (`S15 caption → S15 table → note → S3 caption → S3 figure`). A 22-page A4 fallback render was visually inspected after the packaged LibreOffice renderer was unavailable on this host; no clipping, overlap, missing figures, or table overflow was observed.
