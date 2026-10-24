# Validation against MATISSE (MATLAB)

Same inputs on both sides: the 8 synthetic example datasets of matisseR (`inst/extdata`, 1 x 1; MATLAB reads
them with MATISSE's `csv2rates`), zero cells filled, New APC Model by WLS with overdispersion (MATISSE defaults),
and CAPRICORN on filled rates for three comparisons: populations B vs A and C vs A (1 x 1), and B, C, D vs A
(chunked to 5 x 5).

| Quantity | File |
|---|---|
| Estimable functions LAC, CAC, FTT, FCP, LD (log scale, as MATISSE `ef`) | `apc_ef.csv` |
| Key parameters: intercept, LAT, net drift, CAT, quadratic age/period/cohort | `apc_key.csv` |
| CAPRICORN global PH tests, homogeneity tests, composite tests | `cap_global.csv`, `cap_homogeneity.csv`, `cap_composite.csv` |
| CAPRICORN AICc table | `cap_aic.csv` |
| CAPRICORN contrasts of estimable functions vs the reference stratum | `cap_ef.csv` |

Steps:

1. R (package root): `Rscript validation/validate_r.R` writes `validation/r/` (included in the repository).
2. MATLAB, with MATISSE on the path, from this folder: `validate_matisse` writes `validation/matisse/`.
3. R (package root): `Rscript validation/compare.R` prints, for each quantity, the number of matched values,
   the largest absolute and relative differences, and `agree` when the relative difference is below 1e-6.

Status: R side run; the MATLAB side has not been run yet (`validate_matisse.m` is untested). Small
differences in the CAPRICORN Poisson fits can come from the iteration stopping rules (MATISSE: changes
below 1e-10; R `glm.fit`: relative deviance change below 1e-8).
