# KNHANES frozen external validation: corrected final results

Run date: 2026-09-05  
Primary external validation: KNHANES 2008–2011, age 19–69  
Direct-comparison sensitivity: KNHANES 2008–2011, age 19–59  
Contemporary replication: KNHANES 2024, age 40–69

## Bottom line

After correcting both the official pooled weights and the operational ALM definition, the principal finding remains strong but the calibration estimates change materially.

1. **The age+height structure transports across country and DXA era.** In 2008–2011, age-only z slopes of +0.739 SD/10 cm in women and +0.744 in men fell to +0.021 and −0.079 after height conditioning. In 2024, +0.643/+0.726 fell to −0.119/−0.069.
2. **Absolute U.S.-reference calibration does not transport.** In 2008–2011, mean z was −0.697 in women and −0.644 in men, with 10.27% and 10.21% below the U.S.-derived P5. In 2024, the shift was larger: mean z −1.059/−0.871 and P5 fractions 24.23%/18.07%.
3. **FNIH's inverse height gradient replicated in both Korean cohorts.** ORs per 10 cm were 0.039/0.048 in 2008–2011 and 0.031/0.056 in 2024 (women/men), whereas conditional-P5 ORs were much closer to 1.
4. **The functional signal also replicated in 2024.** In the fully adjusted model, the conditional-P5-only group had 3.18 kg lower bilateral combined grip than the neither-low group (95% CI −4.24 to −2.12); the corresponding max-single-hand difference was −1.63 kg (−2.20 to −1.06).

The defensible interpretation is therefore **successful structural external validation with incomplete distributional calibration**. Fractions below a U.S.-derived P5 must be described as reference-position fractions, not Korean disease prevalence.

## Technical corrections that supersede the first pass

The earlier `/4` weighting and `(LN−BMC)×0.946` analysis is superseded.

- Official 2008–2011 pooled coefficients are annual weight × `108/579`, `199/579`, `192/579`, and `80/579` for 2008–2011.
- For 2008–2011, operational ALM is `0.946 × sum(four limb LN)`. BMC is not subtracted.
- For 2024, operational ALM is `sum(four limb LN)`, without the 0.946 factor and without BMC subtraction.
- The 2024 body-fat covariate uses `DW_SBT_pFT`, and max-single-hand grip is the maximum of the four recorded trials.

These choices are empirically anchored by the published Task Force checksum. For 2008–2011 age ≥50, `LN×0.946` reproduces exactly N=9,219 and gives age 62.06, women 53.89%, height 159.66 cm, weight 61.38 kg, and ALM 17.23 kg, matching the published 62.1, 53.9%, 159.6, 61.4, and 17.2. The erroneous subtract-BMC construction gives ALM 16.29 kg. For 2024 age ≥50 with grip, the checksum also reproduces exactly N=2,481; `sum(LN)`, subtotal fat percentage, and max-single-hand grip reproduce the published means (17.02 kg, 36.24%, and 31.15 kg versus 17.0, 36.3%, and 31.1).

## KNHANES 2008–2011 sample and primary calibration

All 21,303 DXA IDs matched their annual ALL file. The final 19–69 analytic sample was 16,226 (9,240 women; 6,986 men): 3,041 in 2008, 6,102 in 2009, 5,103 in 2010, and 1,980 in 2011.

| Sex | N | Mean z (95% CI) | SD(z) | Below U.S. P5, % (95% CI) | Below U.S. P10, % |
|---|---:|---:|---:|---:|---:|
| Women | 9,240 | −0.697 (−0.724, −0.670) | 0.743 | 10.27 (9.35, 11.27) | 21.24 |
| Men | 6,986 | −0.644 (−0.672, −0.616) | 0.786 | 10.21 (9.25, 11.25) | 20.97 |

The official pooling rule changed results only slightly relative to equal-year pooling: mean z by −0.011 in each sex, P5 by +0.50 percentage points in women and +0.44 in men, and did not alter the substantive conclusions.

## Height dependence, 2008–2011 age 19–69

Survey-weighted models adjust for age and survey year.

| Sex | Score/definition | Effect per 10 cm (95% CI) |
|---|---|---:|
| Women | Age-only z slope | +0.739 (+0.703, +0.776) |
| Women | Age+height z slope | +0.021 (−0.018, +0.061) |
| Men | Age-only z slope | +0.744 (+0.713, +0.775) |
| Men | Age+height z slope | −0.079 (−0.120, −0.038) |
| Women | Conditional P5 OR | 0.848 (0.713, 1.008) |
| Men | Conditional P5 OR | 1.144 (0.970, 1.348) |
| Women | FNIH OR | **0.039 (0.029, 0.052)** |
| Men | FNIH OR | **0.048 (0.036, 0.063)** |

FNIH remained strongly inverse in every year (women 0.027–0.051; men 0.031–0.066). Across height quintiles, FNIH fell from 23.88% to 0.05% in women and 23.26% to 0.25% in men. Conditional P5 varied far less: 8.21% to 12.59% and 8.51% to 9.91%, respectively.

## Classification discordance, 2008–2011

For conditional P5 versus FNIH, weighted group proportions were:

| Sex | Neither | FNIH only | Conditional only | Both |
|---|---:|---:|---:|---:|
| Women | 84.67% | 5.06% | 9.42% | 0.85% |
| Men | 84.40% | 5.39% | 9.27% | 0.94% |

Conditional-only participants were taller and leaner by BMI than FNIH-only participants. Conditional P5 was completely nested within EWGSOP2 in both sexes; there were no conditional-only participants in that comparison.

## Prespecified sensitivities

The 19–59 direct-comparison analysis included 13,098 participants. Mean z was −0.739/−0.656 and the P5 fraction 11.25%/10.47% (women/men). Conditional z slopes were +0.024 and −0.065; conditional-P5 ORs 0.843 and 1.133; FNIH ORs 0.035 and 0.042. Conclusions were unchanged.

The 60–69 subset included 3,128 participants. Conditional slopes were −0.002 in women and −0.195 in men, versus age-only slopes +0.754/+0.753. This supports substantial attenuation but requires a caveat about residual continuous-score dependence in older Korean men.

Restricting to the scorer's recommended height range excluded 345 participants (214 women, 131 men). Mean z, SD, P5, conditional slope, and FNIH OR were nearly unchanged: conditional slopes +0.020/−0.052 and FNIH ORs 0.027/0.050.

Without the mandated 0.946 legacy correction, 2008–2011 mean z becomes −0.371/−0.249 and P5 fractions 3.98%/4.61%. Height conclusions remain essentially unchanged. This confirms that the correction affects location calibration but not the normalization comparison; the uncorrected version is sensitivity-only.

## KNHANES 2024 contemporary replication, age 40–69

The final common sample was 2,513 (1,505 women; 1,008 men). Combined bilateral grip was available for 2,317 and max-single-hand grip for 2,375.

| Sex | Mean z (95% CI) | SD(z) | Below U.S. P5, % (95% CI) | Age-only slope/10 cm | Conditional slope/10 cm |
|---|---:|---:|---:|---:|---:|
| Women | −1.059 (−1.112, −1.005) | 0.825 | 24.23 (21.58, 27.09) | +0.643 | −0.119 (−0.215, −0.024) |
| Men | −0.871 (−0.933, −0.809) | 0.826 | 18.07 (15.54, 20.91) | +0.726 | −0.069 (−0.173, +0.035) |

FNIH ORs per 10 cm were 0.031 (0.019, 0.049) in women and 0.056 (0.032, 0.098) in men. Conditional-P5 ORs were 1.245 (0.933, 1.661) and 1.277 (0.914, 1.785). Thus the conventional height bias replicated, while the conditional binary classification was statistically compatible with no height trend. The female continuous conditional score retained a modest inverse slope.

In the common 40–69 age band, the contemporary distribution was lower than in 2008–2011: female mean z −1.059 versus −0.550 and male −0.871 versus −0.628. This is consistent with the Task Force report's finding of lower relative muscle indices in 2024 and argues against treating the Korean location shift as only an old-scanner artifact.

## 2024 grip characterization

Fully adjusted M3 models include spline age, sex, height, BMI, and subtotal body-fat percentage. Differences are versus the neither-low group.

| Comparison/group | Combined grip difference, kg (95% CI) | Max-single-hand difference, kg (95% CI) |
|---|---:|---:|
| FNIH only | −2.82 (−4.32, −1.32) | −1.64 (−2.35, −0.93) |
| Conditional P5 only | **−3.18 (−4.24, −2.12)** | **−1.63 (−2.20, −1.06)** |
| Both low | −4.50 (−7.22, −1.78) | −2.29 (−3.58, −1.00) |
| EWGSOP2 only | −1.61 (−2.64, −0.58) | −0.95 (−1.43, −0.47) |
| Conditional P5 + EWGSOP2 both low | −4.16 (−5.37, −2.95) | −2.20 (−2.82, −1.58) |

There was no conditional-only group versus EWGSOP2 because conditional P5 was nested within EWGSOP2. Results are association/functional characterization, not prospective diagnostic validation.

## Verification and file status

All 20 automated final checks passed: pooling coefficients, sample sizes, age domains, ALM identity, both published checksum panels, finite estimates, valid prevalence ranges, and positive OR confidence limits. The frozen v2.1.0 scorer and its exact R package environment had already passed its release tests.

Use `tables_official_weight` for the corrected 2008–2011 results and `tables_2024` for the contemporary replication. The old first-pass numbers in prior chat/output history are superseded.
