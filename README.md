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

The 16 MATISSE example datasets (SEER incidence, single years of age 35-84 from 1992) are included.

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

R <- rates_fill(matisse_example(16))               # melanoma, female, non-Hispanic White, 1 x 1
M <- apc_fit(R)                                    # New APC Model, WLS (MATISSE default)
M$net_drift                                        # % per year
apc_ef(M, "ld")                                    # local drifts by age, with 95% CIs
apc_efs(M)                                         # LAC, CAC, FTT, FCP and LD

R5 <- rates_chunk(matisse_example(16), 5, 5)       # 5 x 5 cells, most recent periods kept
R1 <- rates_i5(R5 |> rates_chunk(1, 1), "cubic")   # single-year ages from 5-year groups

Rs <- lapply(matisse_example(9:10), rates_fill)    # thyroid, male: NHB vs NHW
names(Rs) <- c("NHB", "NHW")
S <- capricorn(Rs)                                 # add overdispersion = TRUE for quasi-Poisson tests
S$tests$global; S$aic
cap_ef(S, "fcp")                                   # fitted cohort patterns and their ratio
```

## Functions and their MATISSE sources

| R | MATISSE | |
|---|---|---|
| `csv2rates`, `rates2csv` | `utilities/csv2rates.m`, `@rates/private/RATES2CSV.m` | NCI APC Web Tool CSV |
| `matisse_example` | `data/ratesdata.m` | example datasets 1-16 |
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

Differences from MATISSE: covariances are computed by QR on unit-scaled columns (stable at 1 x 1
resolution; same values), CAPRICORN fits start from the NPH fit, and `capricorn(overdispersion = TRUE)`
adds quasi-Poisson tests and QAICc (MATISSE uses Poisson tests).

## Validation

1. **Internal checks** (`tests/testthat`, run with `devtools::test()`), on the MATISSE examples at 1 x 1 and
   5 x 5: CSV reading and round trip; `chunk`, `fill`, `pchip` and `i5` invariants (totals preserved);
   the design spans the full APC model and WLS and Poisson fits equal age + period + cohort factor models;
   each estimable function equals the fitted rates along its reference line with the other deviations
   removed; local drifts equal the least-squares slope of fitted log rates over periods; CAPRICORN NPH
   equals separate APC fits, and each PH model equals the factor model with shared or stratum-specific
   age, period and cohort effects; stratum estimable functions and ratios match the single-model ones.
2. **Against MATISSE (MATLAB)** ([validation/](validation)): `validate_r.R` and `validate_matisse.m` write
   the same quantities for the 16 examples and three CAPRICORN comparisons; `compare.R` reports the
   differences. R results are in `validation/r/`; MATISSE results are pending (MATLAB needed).

## Citation

- Miranda Filho A, Rosenberg PS. Advances in statistical methods for cancer surveillance research: an
  age-period-cohort perspective. *Front Oncol* 2024;13:1332429. doi:10.3389/fonc.2023.1332429
- Comparative APC (CAPRICORN): Rosenberg PS et al., *BMC Med Res Methodol* 2023 (PMID 37853346).
- Rosenberg PS, Check DP, Anderson WF. A web tool for age-period-cohort analysis of cancer incidence and
  mortality rates. *Cancer Epidemiol Biomarkers Prev* 2014;23:2296-302.

## License

MIT (code). The example datasets are redistributed from the MATISSE toolbox.
