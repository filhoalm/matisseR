## Fitted cohort share of a cause among all deaths: the New APC model fitted to cause-specific deaths with all-cause
## deaths as the offset (direct), or the ratio of the fitted cohort patterns of two rate models (ratio).

#' Cause-specific deaths per all-cause death as a rates object
#'
#' `Events` are the deaths from the cause and `Offset` the all-cause deaths in the same Lexis cell, both after
#' [rates_fill()] (zero or missing counts get MATISSE's adaptive fill value), so that [apc_fit()] models the share of
#' deaths due to the cause. [capricorn()] uses only `Events` and `Offset`, so K cause groups compared on this scale need
#' no wrapper (see [fcp_share()]).
#'
#' @param R_cause,R_all rates objects over the same Lexis diagram: deaths from the cause and from all causes, with
#'   person-years as offsets. These must be present ([rates_fill()] stops on missing offsets) but do not enter the share.
#' @return rates object (`Events` cause deaths, `Offset` all-cause deaths) with `n_filled` and `fill_value` of the cause.
#' @examples
#' # synthetic, female, population A: melanoma-like (1) and all cancers (7), 5 x 5
#' ex <- lapply(matisse_example(c(1, 7)), rates_chunk, age_block = 5, per_block = 5)
#' rates_uni(rates_share(ex[[1]], ex[[2]]))        # melanoma per cancer case
#' @export
rates_share <- function(R_cause, R_all) {
  stopifnot(identical(R_cause$ages, R_all$ages), identical(R_cause$periods, R_all$periods),
            all(R_cause$Events <= R_all$Events * (1 + 1e-8), na.rm = TRUE))
  R <- rates_fill(R_cause); R$Offset <- rates_fill(R_all)$Events
  R$event_label <- paste(R$event_label, "per all-cause death")
  R
}

#' Fitted cohort share of a cause among all deaths
#'
#' Share of the deaths at age `ref_age` due to a cause, by birth cohort, from the New APC model, with two estimators:
#' * `direct` (primary): [apc_fit()] on [rates_share()] (cause deaths with all-cause deaths as the offset; WLS by
#'   default) and the fitted cohort pattern at `ref_age` ([apc_fcp()] with `per = 1`), with the fit's own covariance.
#'   Given the all-cause deaths, cause deaths are binomial, and the Poisson (WLS) variance of the log share exceeds the
#'   binomial one by 1 / (1 - share). With a dispersion of 1 (`overdispersion = FALSE`, or deviance / df <= 1) the
#'   confidence intervals are therefore conservative (standard errors at most 1.4 times too wide below a share of 0.5).
#'   With overdispersion (the default; deviance / df is usually well above 1 for deaths by cause) the scale
#'   max(1, deviance / df) absorbs this factor and the intervals are roughly calibrated: in simulations on Danish
#'   female cancer deaths with extra-binomial variation they were, if anything, slightly too narrow. A log-linear model
#'   of a proportion is adequate only well below 1. A warning is given if a fitted share exceeds `warn_above` in the
#'   Lexis cells (`max_cell`) or in the cohort pattern of a cohort observed within 5 years of `ref_age`
#'   (`max_fcp_observed`); cohorts projected to `ref_age` (`max_fcp_all`) can exceed it, or 1, without a warning, and
#'   generation means built on them should be read with care. `frac_half` is the fraction of cells whose observed
#'   share exceeds 0.5.
#' * `ratio` (secondary check): FCP of the cause over FCP of all causes at `ref_age`, from two rate models (deaths per
#'   person-year), with delta-method variance treating the two fits as independent. Cause deaths are a subset of all
#'   deaths, so the two FCPs are positively correlated and this variance is conservative.
#'
#' Zero or missing death counts are filled as in [rates_fill()]; the direct estimator on all causes over all causes is
#' exactly 1 for every cohort. Ratios of shares between generations: `fcp_contrast(x$fcp$direct)`.
#'
#' CAPRICORN on the share scale: [capricorn()] takes rates objects whose offsets are all-cause deaths directly (it uses
#' only `Events` and `Offset`), e.g. for cancer and cardiovascular deaths `Rc`, `Rv` and all deaths `Ra`,
#' `capricorn(lapply(list(cancer = Rc, cvd = Rv), rates_share, R_all = Ra), overdispersion = TRUE)`; the stratum
#' estimable functions are then shares (`cap_ef(S, "fcp", per = 1)`), and in each cell a ratio of shares is a ratio of
#' cause-specific rates. The fitted ratios differ slightly from those of [capricorn()] on the rates (deaths per
#' person-year): the all-cause term, not exactly an APC surface, does not cancel between strata fitted with different
#' weights. [cap_ef()] reports at the central age group of the diagram (52.5 for ages 20-84), not at `ref_age`, unless
#' `ref` names another age midpoint.
#'
#' @param R_cause,R_all rates objects over the same Lexis diagram (as from [rates_make()] and [rates_chunk()]): deaths
#'   from the cause and from all causes, with person-years as offsets (used by the ratio estimator only).
#' @param ref_age reference age (see [apc_fcp()]); `NULL` for the central age group.
#' @param estimator `direct`, `ratio` or both.
#' @param gens generations ([social_generations()]).
#' @param alpha 1 - confidence level.
#' @param min_cohorts fewest cohorts for a generation mean ([fcp_generations()]).
#' @param warn_above warn if a fitted share in the Lexis cells or in the cohort pattern of an observed cohort exceeds
#'   this value.
#' @param ... passed to [apc_fit()] (`method`, `overdispersion`).
#' @return list with
#'   * `share`: data.table `estimator`, `ref_age`, `cohort` (midpoint), `share`, `lo`, `hi`, `se_log` (standard error
#'     of the log share);
#'   * `generations`: data.table `estimator`, `generation`, `n` (cohorts), `mean` share, `se`, `lo`, `hi` (CI on the
#'     log scale);
#'   * `check`: data.table `estimator`, `max_cell` (largest fitted share in the Lexis cells; direct estimator only),
#'     `max_fcp_observed` (largest share in the cohort pattern among the cohorts with an observed cell at ages within
#'     5 years of `ref_age`, e.g. 45-49 or 50-54 for 50), `max_fcp_all` (all cohorts, including those projected to
#'     `ref_age`), `frac_half`, `n_filled` (cells without cause deaths, filled);
#'   * `fcp`: FCP lists (`x`, `f`, `V`, `age`) by estimator, for [fcp_generations()] and [fcp_contrast()];
#'   * `models`: the fits (`direct`; `cause` and `all` for the ratio).
#' @examples
#' # synthetic, female, population A: melanoma-like (1) and all cancers (7), 5 x 5
#' ex <- lapply(matisse_example(c(1, 7)), rates_chunk, age_block = 5, per_block = 5)
#' Ra <- ex[[2]]
#' s <- fcp_share(ex[[1]], Ra, ref_age = 60)       # melanoma share of all cancers at 60
#' s$share[cohort >= 1930]
#' s$generations
#' fcp_contrast(s$fcp$direct)                      # ratios of shares between generations
#' fcp_share(Ra, Ra, 60, "direct")$share$share     # 1, with a warning (share above 0.95)
#' @export
fcp_share <- function(R_cause, R_all, ref_age = NULL, estimator = c("direct", "ratio"), gens = social_generations(),
                      alpha = 0.05, min_cohorts = 2, warn_above = 0.95, ...) {
  estimator <- match.arg(estimator, several.ok = TRUE)
  Rs <- rates_share(R_cause, R_all); q <- qnorm(1 - alpha / 2)
  M <- fcp <- list()
  if ("direct" %in% estimator) { M$direct <- apc_fit(Rs, ...); fcp$direct <- apc_fcp(M$direct, ref_age, per = 1) }
  if ("ratio" %in% estimator) {
    M$cause <- apc_fit(rates_fill(R_cause), ...); M$all <- apc_fit(rates_fill(R_all), ...)
    fc <- apc_fcp(M$cause, ref_age, per = 1); fa <- apc_fcp(M$all, ref_age, per = 1); f <- fc$f / fa$f
    fcp$ratio <- list(x = fc$x, f = f, V = (fc$V / tcrossprod(fc$f) + fa$V / tcrossprod(fa$f)) * tcrossprod(f), age = fc$age)
  }
  est <- names(fcp); G <- M[[1]]$G; a <- rep(G$age, G$P)
  near <- abs(a - fcp[[1]]$age) <= 5 & !is.na(c(R_all$Events))                 # observed cells within 5 years of ref_age
  obs <- round(fcp[[1]]$x, 6) %in% round((rep(G$per, each = G$A) - a)[near], 6)
  cell <- if ("direct" %in% est) max(exp(G$X %*% M$direct$B)) else NA_real_
  chk <- data.table(estimator = est, max_cell = unname(c(direct = cell, ratio = NA_real_)[est]),
                    max_fcp_observed = vapply(est, function(e) max(fcp[[e]]$f[obs]), 0, USE.NAMES = FALSE),
                    max_fcp_all = vapply(est, function(e) max(fcp[[e]]$f), 0, USE.NAMES = FALSE),
                    frac_half = mean(R_cause$Events / R_all$Events > 0.5, na.rm = TRUE), n_filled = Rs$n_filled)
  mx <- pmax(chk$max_cell, chk$max_fcp_observed, na.rm = TRUE)
  for (i in which(mx > warn_above)) warning(sprintf("fitted share above %g (%s estimator, max %.3f)", warn_above, est[i], mx[i]), call. = FALSE)
  share <- rbindlist(lapply(est, function(e) { f <- fcp[[e]]; s <- sqrt(diag(f$V)) / f$f
    data.table(estimator = e, ref_age = f$age, cohort = f$x, share = f$f, lo = f$f * exp(-q * s), hi = f$f * exp(q * s), se_log = s) }))
  gm <- rbindlist(lapply(est, function(e) { g <- fcp_generations(fcp[[e]], gens, min_cohorts)$table; s <- g$se / g$mean
    data.table(estimator = e, g, lo = g$mean * exp(-q * s), hi = g$mean * exp(q * s)) }))
  list(share = share, generations = gm, check = chk, fcp = fcp, models = M)
}
