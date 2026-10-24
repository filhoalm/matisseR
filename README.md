# matisseR

An R port of methods in **MATISSE**, the MATLAB toolbox for age-period-cohort analysis of vital rates
by Philip S. Rosenberg (US National Cancer Institute), described in Miranda Filho & Rosenberg,
*Front Oncol* 2024 ([doi:10.3389/fonc.2023.1332429](https://doi.org/10.3389/fonc.2023.1332429)).

- **Lexis diagrams of rates**: rates objects, aggregation into larger cells, zero filling, interpolation
  of 5-year to single-year ages (`i5`), NCI APC Web Tool CSV input/output, Lexis plots.
- **The New Age-Period-Cohort Model** and its estimable functions: longitudinal and cross-sectional age
  curves, fitted temporal trends, fitted cohort pattern, local drifts and net drift.
- **CAPRICORN**, hypothesis-based comparative APC analysis: joint fits of K Lexis diagrams, proportional
  models PH-L, PH-T, PH-X, PH-A vs NPH, global, homogeneity and composite tests, AICc, estimable functions
  and rate ratios by stratum.
- **Social generations** (beyond MATISSE): fitted cohort pattern at any reference age, mean rates and rate
  ratios between generations and versus proxy parents, sites combined, early- vs late-onset local drifts,
  random-effects pooling of rate ratios across strata.

Eight **synthetic** example datasets are included (`matisse_example()`): single years of age 35-84 by calendar
years 1991-2019, with Poisson counts drawn from known age-period-cohort surfaces shaped like common cancer patterns
([data-raw/make_examples.R](data-raw/make_examples.R)). They are not registry data. For your own data, read an NCI
APC Web Tool CSV with `csv2rates()` or build a rates object from counts and person-years with `rates_make()`.

## Install

```r
# install.packages("remotes")
remotes::install_github("filhoalm/matisseR")
```

Dependencies: `data.table` (ggplot2, patchwork and scales for `lexis_plot`).

## Example

```r
library(matisseR)
matisse_example()                                  # index of the example datasets

R <- rates_fill(matisse_example(1))                # synthetic melanoma-like, female, population A, 1 x 1
M <- apc_fit(R)                                    # New APC Model, WLS (MATISSE default)
M$net_drift                                        # % per year
apc_ef(M, "ld")                                    # local drifts by age, with 95% CIs
apc_efs(M)                                         # LAC, CAC, FTT, FCP and LD

R5 <- rates_chunk(matisse_example(1), 5, 5)        # 5 x 5 cells, most recent periods kept
A5 <- rates_chunk(matisse_example(1), age_block = 5)    # 5-year age groups by single years
R1 <- rates_i5(A5, "cubic")                        # back to single-year ages (1 x 1 cells)

Rs <- lapply(matisse_example(c(2, 1)), rates_fill) # populations B and A (B differs in its period trend)
names(Rs) <- c("B", "A")
S <- capricorn(Rs)                                 # add overdispersion = TRUE for quasi-Poisson tests
S$tests$global; S$aic                              # PH-X fits best, as built in
cap_ef(S, "fcp")                                   # fitted cohort patterns and their ratio

M5 <- apc_fit(rates_fill(R5))                      # social generations, 5 x 5
f  <- apc_fcp(M5, age = 60)                        # fitted cohort pattern at age 60
fcp_generations(f)$table                           # mean rate by generation
fcp_contrast(f)                                    # Silent vs Greatest, Boomers vs Silent, ..., vs proxy parents
ld_onset(M5, cut = 50)                             # local drifts, early vs late onset
```

## Functions and their MATISSE sources

| R | MATISSE | |
|---|---|---|
| `csv2rates`, `rates2csv` | `utilities/csv2rates.m`, `@rates/private/RATES2CSV.m` | NCI APC Web Tool CSV |
| `matisse_example` | (`data/ratesdata.m`) | 8 synthetic example datasets (MATISSE ships SEER data instead) |
| `rates_make`, `rates_long`, `rates_uni` | `@rates`, `VEC.m`, `UNI.m` | rates objects |
| `rates_chunk` | `@rates/private/CHUNK.m` | larger cells, leftovers dropped from the start (LTF) |
| `rates_fill` | `@rates/private/FILL.m` | zero cells, adaptive fill value (`'adp'`) |
| `rates_i5`, `pchip` | `@rates/private/I5.m`, MATLAB `pchip` | 5-year to single-year ages |
| `lexis_plot` | `@rates/private/LEXIS.m` | Lexis diagram |
| `apc_design` | `@rates/private/GEODE.m` | APC design matrix |
| `apc_fit` | `@apc/private/APC.m` | WLS (default) or Poisson fit, overdispersion |
| `apc_ef`, `apc_efs`, `apc_contrast`, `apc_ref` | `utilities/macaroons.m` | estimable functions, reference cells |
| `capricorn`, `cap_map`, `cap_fit` | `@capricorn/private/CAPRICORN.m` | comparative APC, tests, AICc |
| `cap_ef` | `@capricorn/private/EF.m` | estimable functions by stratum and contrasts |

Beyond MATISSE (social generations, *JAMA Netw Open* 2024): `apc_fcp`, `apc_ef_vcov`, `social_generations`,
`generation_contrasts`, `fcp_generations`, `fcp_contrast`, `fcp_sum`, `ld_onset`, `pool_rr`.

Beyond MATISSE (fitted cohort share of a cause among all deaths, direct and ratio estimators; cause deaths per
all-cause death for `capricorn`): `fcp_share`, `rates_share`.

Beyond MATISSE (cohort-axis lag of a target population behind a reference: 2-strata comparative fits over a grid
of 1-year shifts of the reference, by default with the target cohort curvature = kappa x the reference's; best shift
and approximate profile-deviance interval): `cap_shift`.

Differences from MATISSE: covariances are computed by QR on unit-scaled columns (stable at 1 x 1
resolution; same values), CAPRICORN fits start from the NPH fit, and `capricorn(overdispersion = TRUE)`
adds quasi-Poisson tests and QAICc (MATISSE uses Poisson tests).

## Validation

1. **Internal checks** (`tests/testthat`, run with `devtools::test()`), on the synthetic examples at 1 x 1 and
   5 x 5: CSV reading and round trip; `chunk`, `fill`, `pchip` and `i5` invariants (totals preserved);
   Poisson fits recover the true net drift of each surface; AICc picks the PH model that holds by construction;
   the design spans the full APC model and WLS and Poisson fits equal age + period + cohort factor models;
   each estimable function equals the fitted rates along its reference line with the other deviations
   removed; local drifts equal the least-squares slope of fitted log rates over periods; CAPRICORN NPH
   equals separate APC fits, and each PH model equals the factor model with shared or stratum-specific
   age, period and cohort effects; stratum estimable functions and ratios match the single-model ones;
   the fitted cohort share of all causes among all causes is exactly 1 (also with zero cells), the direct
   share equals the FCP of `rates_share()`, the ratio share and its variance follow the delta method, the
   warning ignores cohorts projected to the reference age, and `capricorn()` on shares (NPH) equals separate
   Poisson fits; `cap_shift()` finds lag 0 for a target equal to the reference and +10 for a copy shifted by 10
   years (each model) and recovers lag 10 and kappa 0.6 for a target with 0.6 x the reference's cohort curvature,
   its profile deviances equal the factor models, and its interval follows the profile rule.
2. **Against MATISSE (MATLAB)** ([validation/](validation)): `validate_r.R` and `validate_matisse.m` write
   the same quantities for the 8 examples and three CAPRICORN comparisons; `compare.R` reports the
   differences. R results are in `validation/r/`; MATISSE results are pending (MATLAB needed).

## Citation

- Miranda Filho A, Rosenberg PS. Advances in statistical methods for cancer surveillance research: an
  age-period-cohort perspective. *Front Oncol* 2024;13:1332429. doi:10.3389/fonc.2023.1332429
- Comparative APC (CAPRICORN): Rosenberg PS et al., *BMC Med Res Methodol* 2023 (PMID 37853346).
- Social generations: Rosenberg PS, Miranda Filho A, *JAMA Netw Open* 2024
  ([doi:10.1001/jamanetworkopen.2024.15731](https://doi.org/10.1001/jamanetworkopen.2024.15731)).
- Rosenberg PS, Check DP, Anderson WF. A web tool for age-period-cohort analysis of cancer incidence and
  mortality rates. *Cancer Epidemiol Biomarkers Prev* 2014;23:2296-302.

## License

MIT, including the synthetic example datasets.
