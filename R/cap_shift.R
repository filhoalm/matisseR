## Cohort-axis lag of a target population relative to a reference (beyond MATISSE). The reference Lexis diagram is
## moved along the cohort axis by d years (years relabelled y -> y + d, age groups fixed) and fitted jointly with the
## target over one fixed window of target years by a 2-strata comparative model (Poisson fits as cap_fit); the lag is
## the shift that minimises the model's deviance, with an approximate profile-deviance interval.

#' Cohort-axis lag of a target population relative to a reference
#'
#' For each shift `d` the reference years are relabelled `y -> y + d` (age groups fixed, so its cohorts and periods
#' both move by `d` and reference cohort `c - d` lies beside target cohort `c`). Target and shifted reference are
#' restricted to the same window of target years, built identically as 5 x 5 Lexis diagrams ([rates_make()],
#' [rates_chunk()] with `per_block = 5`, [rates_fill()]) and fitted jointly, target first, by the comparative model and
#' by NPH (Poisson fits as [cap_fit()], the model started from the NPH fit as in [capricorn()]). A positive lag means
#' that the target follows the reference: target cohort `c` resembles reference cohort `c - lag`.
#'
#' **Window.** The same for every shift: target years complete for all `ages` whose shifted reference year is complete
#' for every shift in the grid (longest run of consecutive years, the most recent if tied, within `years` if given),
#' trimmed from the start to whole 5-year periods (as [rates_chunk()]); at least 3 periods are needed, so a wide grid
#' shortens the window.
#'
#' **Models** (blocks of [cap_blocks()]):
#' * `CohScale` (default): target cohort curvature (QuadCoh, HiOrdCoh) = `kappa` x reference cohort curvature;
#'   intercept, age effects, net drift and period curvature free. `kappa` is profiled at each shift ([optimize()] on
#'   the log scale over 0.1-10) and counts as one parameter. The lag aligns the shape of the two cohort patterns
#'   whatever their amplitude (e.g. the strength of a tobacco epidemic) and whatever the differences in age slope, net
#'   drift and period effects.
#' * `CohCurv`: `CohScale` with `kappa` = 1, i.e. cohort curvature shared (the factor model
#'   `s:(age + period) + cohort`). Unbiased only if the two cohort curvatures have the same amplitude: in simulations
#'   a target curvature 0.6 (1.5) x the reference's biased the lag by about +1.3 (-0.5) years.
#' * `PHL` ([cap_map()]): age (LAT, QuadAge, HiOrdAge) and period curvature (QuadPer, HiOrdPer) shared; intercept,
#'   net drift and cohort curvature free. A pure cohort lag is absorbed by the free cohort effects, so the PH-L
#'   profile is driven by the alignment of the period curvature (moved by `d` with the cohorts), not of the cohort
#'   patterns: with common period effects it is drawn towards 0, without period effects it is flat.
#' * `PHT` ([cap_map()]): net drift, period and cohort curvature shared; age effects free.
#'
#' **Profile.** `stat = "lrt"`: deviance of the model minus that of NPH at each shift (the CAPRICORN global statistic;
#' it removes the reference's own lack of fit, which changes with the reference years used); `"dev"`: deviance of the
#' model. `phi` = max(1, deviance / df) of the NPH fit at the grid minimum (1 if `overdispersion = FALSE`). `misfit` =
#' lrt_min / (df_nph - df), lrt_min the minimum of the lrt profile refined as the best shift below (at least 0): of the
#' order of `phi` when the model holds, much larger when the two cohort patterns differ in ways the model does not
#' absorb (shape; amplitude too for `CohCurv`).
#'
#' **Interval.** Best shift: grid minimum refined by the vertex of the parabola through it and its two neighbours (not
#' refined at the grid edge). Grid shifts with `z` = (D(d) - D at the grid minimum) / phi <= qchisq(level, 1) are
#' inside. If the neighbour of the grid minimum on one side is outside, the limit on that side is vertex -/+
#' sqrt(qchisq(level, 1) phi / a) on the parabola (curvature `a`), so the interval need not contain the grid minimum
#' (`min_out`); otherwise it is interpolated linearly between the outermost grid point inside and the next one. Flags:
#' `sub_grid` (no other grid point inside), `disjoint` (grid points between the limits outside), `edge_lo` / `edge_hi`
#' (the limit is the grid edge, the interval is open on that side).
#' The interval is approximate and covers sampling error only. It is too narrow under overdispersion (`phi` scales it
#' but coverage was about 0.8-0.9 in simulations with phi near 4), and it does not cover structural differences
#' between the two populations (e.g. cohort patterns of different shape, or of different amplitude for `CohCurv`),
#' which bias the lag without widening the interval: check `misfit` and compare models.
#'
#' @param target,ref annual data (data.frame or data.table), one row per 5-year age group and calendar year: `age`
#'   (group index as in [rates_make()]: 1 = 0-4, 2 = 5-9, ..., 18 = 85+), `year` (single calendar year, so that
#'   1-year shifts are possible), `cases` (events) and `py` (person-years), all numeric; rows with the same age and
#'   year are summed.
#' @param shifts candidate lags in whole years (at least 3), e.g. `-10:30`.
#' @param ages consecutive age-group indices to compare (at least 3), e.g. `7:17` for 30-84.
#' @param model `CohScale`, `CohCurv`, `PHL` or `PHT`.
#' @param stat profile statistic: `lrt` (model minus NPH deviance) or `dev` (model deviance).
#' @param years optional target years to which the window is restricted.
#' @param overdispersion scale the profile by the NPH dispersion `phi`.
#' @param level confidence level of the profile interval, in (0, 1).
#' @return list with
#'   * `summary`: one-row data.table `model`, `stat`, `shift` (refined best lag), `lo`, `hi`, `kappa` (`CohScale`, at
#'     the grid minimum; NA for the other models), `misfit`, `phi`, `edge_lo`, `edge_hi`, `sub_grid`, `disjoint`,
#'     `min_out` (grid minimum outside the interval), `grid_from`, `grid_to`, `year_from`, `year_to` (window),
#'     `periods`, `coh_from`, `coh_to` (cohort midpoints of the evaluation diagram), `n_failed` (shifts whose fits
#'     failed, left out);
#'   * `profile`: data.table by shift: `dev`, `df` (model), `kappa`, `dev_nph`, `df_nph`, `lrt`, `z` of the chosen
#'     statistic, `inside`;
#'   * `window` (first and last target year), `phi`, and `Rs`: the target and shifted reference rates objects at the
#'     grid minimum (for [capricorn()] or [cap_ef()]); the reference periods are labelled in target years (reference
#'     years + shift).
#' @examples
#' \donttest{
#' library(data.table)
#' set.seed(1)
#' ref <- CJ(age = 7:17, year = 1950:2019)
#' ref[, `:=`(c = year + 0.5 - (5 * age - 2.5), py = 1e6)]
#' ref[, cases := rpois(.N, py * exp(-12 + 0.09 * (5 * age - 2.5) + 0.6 * exp(-((c - 1925) / 15)^2) +
#'                                     0.15 * sin((year - 1950) / 7)))]
#' tgt <- copy(ref)[, year := year + 10]          # the target follows the reference by 10 years
#' s <- cap_shift(tgt, ref, shifts = 0:20, ages = 7:17)
#' s$summary
#' }
#' @importFrom stats qchisq optimize
#' @export
cap_shift <- function(target, ref, shifts, ages, model = c("CohScale", "CohCurv", "PHL", "PHT"),
                      stat = c("lrt", "dev"), years = NULL, overdispersion = TRUE, level = 0.95) {
  model <- match.arg(model); stat <- match.arg(stat)
  if (!is.numeric(ages) || length(ages) < 3 || anyNA(ages) || any(diff(ages) != 1) || any(ages < 1 | ages > 18))
    stop("ages must be at least 3 consecutive age-group indices in 1-18")
  if (!is.numeric(shifts) || anyNA(shifts) || any(shifts %% 1 != 0) || length(unique(shifts)) < 3)
    stop("shifts must be at least 3 distinct whole numbers of years")
  stopifnot("level must be a number in (0, 1)" = is.numeric(level) && length(level) == 1 && level > 0 && level < 1)
  shifts <- sort(unique(shifts))
  long <- function(x, nm) {
    if (!is.data.frame(x)) stop(nm, " must be a data.frame or data.table")
    miss <- setdiff(c("age", "year", "cases", "py"), names(x))
    if (length(miss)) stop(nm, " lacks column(s) ", paste(miss, collapse = ", "))
    x <- as.data.table(x)[, .(age, year, cases, py)]
    if (!all(vapply(x, is.numeric, TRUE))) stop(nm, ": age, year, cases and py must be numeric")
    if (any(x$year %% 1 != 0, na.rm = TRUE)) stop(nm, ": year must be single calendar years")
    miss <- setdiff(ages, x$age)
    if (length(miss)) stop(nm, " lacks age group(s) ", paste(miss, collapse = ", "), " (indices as in rates_make())")
    x[age %in% ages]
  }
  tg <- long(target, "target"); rf <- long(ref, "ref")
  full <- function(x) x[, .(ok = !anyNA(cases) && !anyNA(py) && all(cases >= 0) && all(py > 0)), by = .(age, year)
                        ][(ok), .N, by = year][N == length(ages), year]
  yr <- full(rf); ok <- sort(Reduce(intersect, lapply(shifts, function(d) yr + d), full(tg)))
  if (!is.null(years)) ok <- intersect(ok, years)
  run <- if (length(ok)) { r <- cumsum(c(1, diff(ok) != 1)); n <- tabulate(r); ok[r == max(which(n == max(n)))] } else ok
  np <- length(run) %/% 5
  if (np < 3) stop(sprintf(paste("cap_shift(): %d consecutive target year(s)%s have the shifted reference for every",
                                 "shift in %g to %g; at least 15 (3 five-year periods) are needed: narrow the shift grid",
                                 "or extend the series"), length(run),
                           if (length(run)) sprintf(" (%g-%g)", min(run), max(run)) else "", min(shifts), max(shifts)))
  win <- tail(run, 5 * np)
  mk <- function(x) rates_fill(rates_chunk(rates_make(x[year %in% win], ages), per_block = 5))
  mv <- function(d) rf[, .(age, year = year + d, cases, py)]
  Rt <- mk(tg); G <- apc_design(Rt); X <- G$X; N <- length(Rt$Events)
  I <- diag(ncol(X)); Pt <- cap_blocks(G)
  m0 <- cap_map(G, "NPH")
  m1 <- if (model %in% c("CohScale", "CohCurv")) list(E = I[, unlist(Pt[c(6, 9)]), drop = FALSE],
                                                      F = I[, unlist(Pt[c(1:5, 7, 8)]), drop = FALSE]) else cap_map(G, model)
  fit <- function(Y, O, m, k = 1, eta0 = NULL) {                  # as cap_fit, shared columns x k in the target rows
    Z <- cbind(kronecker(c(k, 1), X %*% m$E), kronecker(diag(2), X %*% m$F)); Z <- sweep(Z, 2, sqrt(colSums(Z^2)), "/")
    g <- function(e) tryCatch(suppressWarnings(glm.fit(Z, Y, etastart = e, offset = log(O), family = poisson(),
                                                       control = glm.control(maxit = 100))), error = function(err) NULL)
    f <- g(eta0); if (!isTRUE(f$converged) && !is.null(eta0)) f <- g(NULL)
    if (isTRUE(f$converged) && f$rank == ncol(Z)) f
  }
  dv <- function(f) if (is.null(f)) NA_real_ else f$deviance
  P <- rbindlist(lapply(shifts, function(d) {
    Rr <- mk(mv(d)); Y <- c(Rt$Events, Rr$Events); O <- c(Rt$Offset, Rr$Offset)
    f0 <- fit(Y, O, m0); e0 <- f0$linear.predictors
    k <- if (model == "CohScale" && !is.null(f0))
      exp(optimize(function(l) dv(fit(Y, O, m1, exp(l), e0)), log(c(0.1, 10)))$minimum) else 1
    data.table(shift = d, dev = if (is.null(f0)) NA_real_ else dv(fit(Y, O, m1, k, e0)),
               df = ncol(m1$E) + 2 * ncol(m1$F) + (model == "CohScale"), kappa = if (model == "CohScale") k else NA_real_,
               dev_nph = dv(f0), df_nph = 2 * ncol(X))
  }))
  P[, lrt := dev - dev_nph]
  nf <- sum(is.na(P$dev))
  if (nf) warning("cap_shift(): fits failed at shift(s) ", paste(P[is.na(dev), shift], collapse = ", "),
                  "; left out of the profile", call. = FALSE)
  Q <- P[!is.na(dev)]
  if (nrow(Q) < 3) stop("cap_shift(): fewer than 3 shifts could be fitted")
  im <- which.min(Q[[stat]])
  phi <- if (overdispersion) max(1, Q$dev_nph[im] / (2 * N - Q$df_nph[im])) else 1
  ci <- shift_ci(Q$shift, Q[[stat]], phi, level)
  P[!is.na(dev), `:=`(z = ci$z, inside = ci$inside)]
  sm <- data.table(model, stat, shift = ci$best, lo = ci$lo, hi = ci$hi, kappa = Q$kappa[im],
                   misfit = max(0, shift_ci(Q$shift, Q$lrt, 1)$vmin) / (Q$df_nph[im] - Q$df[im]), phi, edge_lo = ci$edge_lo, edge_hi = ci$edge_hi,
                   sub_grid = ci$sub_grid, disjoint = ci$disjoint, min_out = ci$min_out, grid_from = min(shifts),
                   grid_to = max(shifts), year_from = min(win), year_to = max(win), periods = np, coh_from = min(G$coh),
                   coh_to = max(G$coh), n_failed = nf)
  list(summary = sm, profile = P[], window = range(win), phi = phi,
       Rs = list(target = Rt, ref = mk(mv(Q$shift[im]))))
}

## Refined minimum and profile interval of a deviance profile v over shifts s (sorted), scaled by phi. Grid points are
## judged against the observed grid minimum; a limit next to it is taken on the parabola around the vertex (a chord
## from the minimum would cross too early), further out linearly between grid points.
shift_ci <- function(s, v, phi, level = 0.95) {
  n <- length(s); i <- which.min(v); best <- s[i]; vmin <- v[i]; a <- 0
  if (i > 1 && i < n) {                                                   # parabola through the minimum and neighbours
    x <- s[i + -1:1]; y <- v[i + -1:1]; f12 <- (y[2] - y[1]) / (x[2] - x[1])
    a <- ((y[3] - y[2]) / (x[3] - x[2]) - f12) / (x[3] - x[1])
    if (a > 0) { best <- (x[1] + x[2] - f12 / a) / 2; vmin <- y[1] + f12 * (best - x[1]) + a * (best - x[1]) * (best - x[2]) }
  }
  z <- (v - v[i]) / phi; crit <- qchisq(level, 1); inn <- z <= crit
  j <- range(which(inn)); r <- if (a > 0) sqrt(crit * phi / a)
  lo <- if (j[1] == 1) s[1] else if (j[1] == i && a > 0) best - r else
    s[j[1] - 1] + (z[j[1] - 1] - crit) / (z[j[1] - 1] - z[j[1]]) * (s[j[1]] - s[j[1] - 1])
  hi <- if (j[2] == n) s[n] else if (j[2] == i && a > 0) best + r else
    s[j[2]] + (crit - z[j[2]]) / (z[j[2] + 1] - z[j[2]]) * (s[j[2] + 1] - s[j[2]])
  list(best = best, vmin = vmin, lo = lo, hi = hi, edge_lo = inn[1], edge_hi = inn[n], disjoint = !all(inn[j[1]:j[2]]),
       sub_grid = j[1] == i && j[2] == i, min_out = s[i] < lo || s[i] > hi, z = z, inside = inn)
}

utils::globalVariables(c("ok", "N", "shift", "dev", "dev_nph", "lrt", "z", "inside"))
