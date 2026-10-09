# KNHANES mirrored-reference PSU bootstrap: final summary

## Analysis status

The prespecified design-consistent bootstrap is complete. All automated verification checks passed. The analysis used 500 stratified PSU bootstrap replicates, refit the locked Korean mirrored models without model reselection (women: BCCG; men: BCT), evaluated predictions on the fixed U.S.–Korea common age–height grid, and compared them with the fixed frozen U.S. reference.

- Women: 500/500 successful fits (100.0%).
- Men: 493/500 successful fits (98.6%). Seven male BCT replicates did not meet the convergence criterion after the prespecified RS/CG rescue sequence; they were retained in the diagnostics as failures and were not replaced.
- All nine requested metrics were finite among successful replicates.
- All 36 percentile confidence intervals were finite and correctly ordered.
- Directly recomputed point estimates matched the locked surface summary (maximum absolute difference: 3.87e-13).
- Representative-height P50 bands were complete (306 age–height rows), and both figures passed visual inspection.

## Primary common domain (1st–99th percentile overlap)

Values are point estimate (95% percentile bootstrap CI). Delta P50 is Korea minus the frozen U.S. reference.

| Metric | Women | Men |
|---|---:|---:|
| Mean delta P50, kg | -1.821 (-1.876, -1.759) | -2.022 (-2.115, -1.907) |
| Mean delta P50, % | -10.889 (-11.217, -10.528) | -8.009 (-8.372, -7.568) |
| (P95-P5) width ratio, Korea/U.S. | 0.548 (0.535, 0.561) | 0.679 (0.661, 0.693) |
| (P90-P10) width ratio, Korea/U.S. | 0.560 (0.547, 0.572) | 0.694 (0.676, 0.708) |
| P50 surface R-squared | 0.895 (0.878, 0.912) | 0.989 (0.972, 0.992) |
| Korea-on-U.S. P50 surface slope | 0.893 (0.852, 0.936) | 0.864 (0.830, 0.908) |
| Offset-adjusted RMSE, kg | 0.548 (0.504, 0.595) | 0.450 (0.406, 0.567) |
| Delta P50 trend per 10 years, kg | +0.352 (+0.311, +0.390) | +0.035 (-0.011, +0.091) |
| Delta P50 trend per 10 cm, kg | -0.182 (-0.276, -0.080) | -0.421 (-0.547, -0.250) |

## Strict common domain (5th–95th percentile overlap)

| Metric | Women | Men |
|---|---:|---:|
| Mean delta P50, kg | -1.851 (-1.905, -1.787) | -2.133 (-2.213, -2.041) |
| Mean delta P50, % | -10.986 (-11.308, -10.609) | -8.410 (-8.730, -8.051) |
| (P95-P5) width ratio, Korea/U.S. | 0.548 (0.535, 0.560) | 0.676 (0.659, 0.690) |
| (P90-P10) width ratio, Korea/U.S. | 0.559 (0.546, 0.571) | 0.691 (0.673, 0.705) |
| P50 surface R-squared | 0.780 (0.746, 0.811) | 0.984 (0.960, 0.988) |
| Korea-on-U.S. P50 surface slope | 0.851 (0.808, 0.900) | 0.841 (0.806, 0.878) |
| Offset-adjusted RMSE, kg | 0.540 (0.495, 0.591) | 0.367 (0.322, 0.460) |
| Delta P50 trend per 10 years, kg | +0.352 (+0.311, +0.390) | +0.035 (-0.011, +0.091) |
| Delta P50 trend per 10 cm, kg | -0.207 (-0.309, -0.089) | -0.471 (-0.614, -0.319) |

## Interpretation

The Korean conditional median ALM surface is clearly lower than the frozen U.S. reference in both sexes, by about 1.82 kg (10.89%) in women and 2.02 kg (8.01%) in men in the primary common domain. The Korean conditional distribution is also substantially narrower: the P95-P5 width is approximately 55% of the U.S. width in women and 68% in men.

The surface differences are not consistent with a simple constant offset. In women, the Korea-minus-U.S. P50 difference becomes less negative with age but more negative with increasing height; both trends have intervals excluding zero. In men, there is little evidence of a systematic age gradient because its interval crosses zero, whereas the increasingly negative difference with height is clear. Slopes below one and nonzero offset-adjusted RMSE further support sex-specific nonconstant geometry despite high overall surface correspondence.

Results in the strict 5th–95th percentile overlap domain are materially similar, supporting robustness to boundary-grid choices.

## Statistical scope and limitations

- The frozen U.S. reference was treated as fixed by design. Confidence intervals quantify Korean survey sampling and Korean mirrored-refit uncertainty, not uncertainty from re-estimating the U.S. reference.
- PSU resampling was performed within 97 survey strata. Six singleton strata (191 participants; 1.304% weighted) remained fixed under the survey bootstrap implementation.
- Confidence intervals are unadjusted 2.5th–97.5th percentile bootstrap intervals based on the successful sex-specific fits (women: 500; men: 493).
- The analysis does not treat grid cells as independent observations and does not add cell-level P values.
- No family search, model selection, local 2024 reference, 19–59 local curve, BMI-conditioned reference, centering sensitivity, or additional functional analysis was performed.

## Manuscript-ready core statements

In women, the Korean conditional P50 was on average 1.82 kg lower than the frozen U.S. reference (95% CI 1.76 to 1.88 kg lower), corresponding to a 10.89% deficit (95% CI 10.53% to 11.22%). In men, the corresponding difference was 2.02 kg lower (95% CI 1.91 to 2.12 kg lower), or 8.01% (95% CI 7.57% to 8.37%).

The Korea-minus-U.S. P50 difference varied with age and height in women (+0.352 kg per 10 years and -0.182 kg per 10 cm), while in men there was little evidence of an age gradient (+0.035 kg per 10 years; 95% CI -0.011 to +0.091) but a clear height gradient (-0.421 kg per 10 cm; 95% CI -0.547 to -0.250).

These results support cross-population transport of the broad age–height conditioning structure while demonstrating population-specific location, scale, and residual geometry differences that preclude interpreting the Korean surface as a simple constant shift of the U.S. reference.
