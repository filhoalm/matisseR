#' matisseR: age-period-cohort analysis of vital rates, ported from MATISSE
#'
#' An R port of methods in the MATISSE MATLAB toolbox (Philip S. Rosenberg, US National Cancer
#' Institute): Lexis diagrams of rates (`rates_*`, [csv2rates()], [rates_i5()]), the New
#' Age-Period-Cohort Model and its estimable functions ([apc_fit()], [apc_ef()]) and hypothesis-based
#' comparative APC analysis ([capricorn()], [cap_ef()]). The MATISSE example datasets are included
#' ([matisse_example()]).
#'
#' @references Miranda Filho A, Rosenberg PS. Advances in statistical methods for cancer surveillance
#'   research: an age-period-cohort perspective. Front Oncol 2024;13:1332429. doi:10.3389/fonc.2023.1332429
#' @keywords internal
#' @import data.table
#' @importFrom stats approx glm.control glm.fit lm.wfit pchisq poisson qnorm setNames
#' @importFrom utils head tail
"_PACKAGE"

utils::globalVariables(c(".", "age", "year", "cases", "py", "i", "j", "age_lo", "age_mid", "per_lo", "per_mid",
                         "events", "offset", "coh_mid", "rate", "number", "AICc", "DEV", "delta", "model", "df"))
