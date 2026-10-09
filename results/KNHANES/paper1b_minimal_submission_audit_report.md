# Paper 1B minimal submission audit

**Analysis status: PASS**

This audit used only the locked analysis-ready KNHANES objects, their existing survey-design variables, and locked result tables. It did not rescore ALM, re-estimate the U.S. reference, refit a Korean reference, or run any bootstrap analysis.

## 1. KNHANES 2024 P10 calibration

| Sex | Below U.S. P10, weighted % (95% CI) | Raw events / N | P5 consistency check |
|---|---:|---:|---:|
| Female | 39.02 (95% CI 35.80 to 42.33) | 577 / 1505 | 24.23% (locked 24.23%) |
| Male | 31.30 (95% CI 28.08 to 34.70) | 311 / 1008 | 18.07% (locked 18.07%) |

P10 used the prespecified threshold z < -1.281551566. P5 was recomputed only as a consistency check.

## 2. KNHANES 2024 age-only z slope

| Sex | Beta per 10-cm greater height | SE | 95% CI | P value |
|---|---:|---:|---:|---:|
| Female | 0.643 | 0.043 | 0.560 to 0.727 | <0.001 |
| Male | 0.726 | 0.043 | 0.642 to 0.809 | <0.001 |

Models used the locked 2024 analytic sample and survey design, with age adjustment and height expressed per 10 cm.

## 3. FNIH rare-event and near-separation audit

### Primary survey-logistic models

| Cohort | Sex | Log-odds coefficient | SE | OR per 10 cm (95% CI) | Converged | Finite coefficient/SE | Numerical warning | Zero-event quintile |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 2008-2011 (19-69) | Female | -3.252 | 0.146 | 0.039 (95% CI 0.029 to 0.052) | Yes | Yes | No | No |
| 2008-2011 (19-69) | Male | -3.040 | 0.144 | 0.048 (95% CI 0.036 to 0.063) | Yes | Yes | No | Yes |
| 2024 (40-69) | Female | -3.484 | 0.243 | 0.031 (95% CI 0.019 to 0.049) | Yes | Yes | No | No |
| 2024 (40-69) | Male | -2.886 | 0.285 | 0.056 (95% CI 0.032 to 0.098) | Yes | Yes | No | No |

### Sex-specific weighted height quintiles

| Cohort | Sex | Quintile | Raw FNIH events / N | Raw % | Survey-weighted % (95% CI) |
|---|---|---:|---:|---:|---:|
| 2008-2011 (19-69) | Female | Q1 | 518 / 2130 | 24.319 | 23.882 (95% CI 21.355 to 26.607) |
| 2008-2011 (19-69) | Female | Q2 | 70 / 1919 | 3.648 | 3.804 (95% CI 2.880 to 5.010) |
| 2008-2011 (19-69) | Female | Q3 | 25 / 1849 | 1.352 | 1.332 (95% CI 0.778 to 2.274) |
| 2008-2011 (19-69) | Female | Q4 | 4 / 1731 | 0.231 | 0.199 (95% CI 0.065 to 0.602) |
| 2008-2011 (19-69) | Female | Q5 | 1 / 1611 | 0.062 | 0.053 (95% CI 0.007 to 0.376) |
| 2008-2011 (19-69) | Male | Q1 | 453 / 1729 | 26.200 | 23.258 (95% CI 20.914 to 25.779) |
| 2008-2011 (19-69) | Male | Q2 | 80 / 1476 | 5.420 | 5.085 (95% CI 3.782 to 6.806) |
| 2008-2011 (19-69) | Male | Q3 | 28 / 1362 | 2.056 | 2.293 (95% CI 1.501 to 3.489) |
| 2008-2011 (19-69) | Male | Q4 | 0 / 1238 | 0.000 | 0.000 (95% CI 0.000 to 0.000) |
| 2008-2011 (19-69) | Male | Q5 | 3 / 1181 | 0.254 | 0.254 (95% CI 0.062 to 1.039) |
| 2024 (40-69) | Female | Q1 | 151 / 329 | 45.897 | 45.775 (95% CI 40.084 to 51.578) |
| 2024 (40-69) | Female | Q2 | 45 / 312 | 14.423 | 14.076 (95% CI 9.887 to 19.652) |
| 2024 (40-69) | Female | Q3 | 11 / 297 | 3.704 | 3.163 (95% CI 1.780 to 5.562) |
| 2024 (40-69) | Female | Q4 | 6 / 294 | 2.041 | 2.221 (95% CI 0.964 to 5.031) |
| 2024 (40-69) | Female | Q5 | 1 / 273 | 0.366 | 0.389 (95% CI 0.053 to 2.765) |
| 2024 (40-69) | Male | Q1 | 119 / 227 | 52.423 | 50.435 (95% CI 42.278 to 58.568) |
| 2024 (40-69) | Male | Q2 | 48 / 204 | 23.529 | 20.417 (95% CI 15.307 to 26.696) |
| 2024 (40-69) | Male | Q3 | 14 / 202 | 6.931 | 6.973 (95% CI 4.041 to 11.769) |
| 2024 (40-69) | Male | Q4 | 6 / 196 | 3.061 | 4.998 (95% CI 1.698 to 13.810) |
| 2024 (40-69) | Male | Q5 | 2 / 179 | 1.117 | 1.048 (95% CI 0.247 to 4.329) |

All eight locked sex-by-year FNIH ORs for 2008-2011 remained below 1: yes (female 0.027-0.051; male 0.031-0.066).
Sparse-event conclusion: one zero-event cell was present (2008-2011 men, Q4), while the adjacent highest quintile contained 3/1181 events. All four continuous-height survey-logistic models converged with finite coefficients and SEs and produced no numerical warnings. The isolated zero cell therefore documents sparsity but does not show complete separation of the prespecified continuous-height model.
Penalized/Firth logistic regression was not used to replace the prespecified complex-survey model. Any trigger and the reason for not forcing an incompatible sensitivity are recorded in the CSV audit table.

## 4. Conditional P5 x EWGSOP2 nesting audit

| Cohort | Sex | Classification cell | Raw N | Survey-weighted % (95% CI) |
|---|---|---|---:|---:|
| 2008-2011 (19-69) | Female | Neither low | 6504 | 68.837 (95% CI 67.326 to 70.309) |
| 2008-2011 (19-69) | Female | EWGSOP2 only | 1870 | 20.894 (95% CI 19.791 to 22.041) |
| 2008-2011 (19-69) | Female | Conditional P5 only | 0 | 0.000 (95% CI 0.000 to 0.000) |
| 2008-2011 (19-69) | Female | Both low | 866 | 10.269 (95% CI 9.347 to 11.271) |
| 2008-2011 (19-69) | Male | Neither low | 5721 | 82.825 (95% CI 81.610 to 83.976) |
| 2008-2011 (19-69) | Male | EWGSOP2 only | 579 | 6.968 (95% CI 6.341 to 7.652) |
| 2008-2011 (19-69) | Male | Conditional P5 only | 0 | 0.000 (95% CI 0.000 to 0.000) |
| 2008-2011 (19-69) | Male | Both low | 686 | 10.207 (95% CI 9.254 to 11.245) |
| 2024 (40-69) | Female | Neither low | 702 | 47.163 (95% CI 43.996 to 50.354) |
| 2024 (40-69) | Female | EWGSOP2 only | 440 | 28.610 (95% CI 26.300 to 31.036) |
| 2024 (40-69) | Female | Conditional P5 only | 0 | 0.000 (95% CI 0.000 to 0.000) |
| 2024 (40-69) | Female | Both low | 363 | 24.227 (95% CI 21.581 to 27.086) |
| 2024 (40-69) | Male | Neither low | 702 | 70.745 (95% CI 67.455 to 73.831) |
| 2024 (40-69) | Male | EWGSOP2 only | 129 | 11.187 (95% CI 9.248 to 13.473) |
| 2024 (40-69) | Male | Conditional P5 only | 0 | 0.000 (95% CI 0.000 to 0.000) |
| 2024 (40-69) | Male | Both low | 177 | 18.068 (95% CI 15.540 to 20.905) |

All participants classified below the frozen conditional P5 also met the EWGSOP2 low-ALMI criterion in this sample.

Nesting is a classification relationship in these samples and is not evidence of diagnostic superiority.

## Consistency checks

- Locked sample sizes reproduced: yes.
- Locked 2024 P5 fractions reproduced within 0.05 percentage points: yes.
- Locked 2024 conditional-z stature slopes reproduced from the existing result table: yes.
- Refitted FNIH ORs matched the locked main estimates: yes.
- Conditional-P5-only raw count was zero in every cohort-sex table: yes.

## Effect on Paper 1B conclusions

No audit result requires modification of the existing Paper 1B main conclusions.

**Main conclusions unchanged; manuscript can proceed to figure/manuscript finalization.**
