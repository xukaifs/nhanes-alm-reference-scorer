# KNHANES 2024 continuous conditional ALM z and grip strength

Population: KNHANES 2024, age 40–69 years.  
Primary outcome: bilateral combined grip (best right + best left).  
Sensitivity outcome: maximum single-hand grip across four trials.  
Survey design: KNHANES strata, PSU, and `wt_itvex`.

## Main models

| Outcome | Model | N | β per +1 conditional z, kg | 95% CI | P |
|---|---|---:|---:|---:|---:|
| Bilateral combined grip | M0: unadjusted | 2,317 | +5.89 | +4.95 to +6.82 | <0.001 |
| Bilateral combined grip | M1: age spline + sex + height | 2,317 | +4.33 | +3.78 to +4.88 | <0.001 |
| Bilateral combined grip | **M2: M1 + BMI (primary)** | 2,317 | **+5.59** | **+4.74 to +6.45** | **<0.001** |
| Bilateral combined grip | M3: M2 + subtotal body-fat % | 2,317 | +5.50 | +4.28 to +6.71 | <0.001 |
| Max single-hand grip | M1 | 2,375 | +2.34 | +2.08 to +2.60 | <0.001 |
| Max single-hand grip | M2 | 2,375 | +2.95 | +2.51 to +3.40 | <0.001 |
| Max single-hand grip | M3 | 2,375 | +2.97 | +2.32 to +3.62 | <0.001 |

Weighted mean grip was 64.34 kg for the bilateral sum and 33.49 kg for max-single-hand grip.

## Sex interaction

The prespecified M3 interaction was significant for both outcomes:

| Outcome | P interaction | Women β per +1 z (95% CI) | Men β per +1 z (95% CI) |
|---|---:|---:|---:|
| Bilateral combined grip | 0.0018 | +4.42 (+3.17 to +5.67) | +6.26 (+4.86 to +7.66) |
| Max single-hand grip | 0.0026 | +2.44 (+1.77 to +3.10) | +3.35 (+2.62 to +4.08) |

Therefore, the direction is highly consistent in both sexes, but the magnitude is larger in men. A pooled result can be shown for comparability, but it should be accompanied by the interaction and sex-specific estimates.

## Location-shift robustness

Sex-specific weighted mean conditional z values were −1.059 in women and −0.871 in men. Subtracting these Korean means and repeating the combined-grip M3 model changed the slope from 5.497828 to 5.497828 kg (difference <6×10⁻¹⁵). This confirms algebraically and empirically that the Korean location shift changes the score origin, not its grip-strength gradient.

## Nonlinearity

For combined grip under M3 adjustment, the restricted/natural cubic spline association was strongly significant overall (3-df Wald P=7.0×10⁻¹⁵), while the two-degree-of-freedom nonlinearity test was not significant (P=0.608). The linear per-1-z summary is therefore adequate; no spline figure is indicated.

## Suggested result sentence

In KNHANES 2024 adults aged 40–69 years, each 1-unit higher conditional ALM z-score was associated with 5.59 kg higher bilateral combined grip strength after adjustment for spline age, sex, height, and BMI (95% CI 4.74–6.45; P<0.001). The association remained after additional adjustment for body-fat percentage (5.50 kg; 95% CI 4.28–6.71; P<0.001), was linear (P for nonlinearity=0.608), and was stronger in men than women (P for interaction=0.0018).

This completes the prespecified KNHANES continuous-z functional replication. No additional KNHANES analyses are recommended at this stage.
