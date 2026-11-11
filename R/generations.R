## Social generations from fitted cohort patterns (FCP), following Rosenberg & Miranda-Filho,
## JAMA Netw Open 2024 (doi:10.1001/jamanetworkopen.2024.15731): FCP at a reference age by birth
## cohort, mean FCP per generation, rate ratios between generations, sites combined, pooled ratios.

#' Covariance of an estimable function
#'
#' @inheritParams apc_ef
#' @return list with `x`, `est` (log scale) and covariance `V`.
#' @export
apc_ef_vcov <- function(M, comp = c("lac", "cac", "ftt", "fcp", "ld"), ref = NULL) {
  comp <- match.arg(comp)
  k <- apc_contrast(M$G, comp, apc_ref(M, ref))
  list(x = k$x, est = drop(k$X %*% M$B[k$i]), V = k$X %*% M$V[k$i, k$i] %*% t(k$X))
}

#' Fitted cohort pattern at any reference age
#'
#' FCP at `age` (log-linear interpolation between the adjacent age midpoints, e.g. 60 between 57.5
#' and 62.5 for 5-year groups), with its covariance on the rate scale (delta method).
#'
#' @param M fitted model ([apc_fit()]).
#' @param age reference age; `NULL` for the MATISSE default (central age group).
#' @param per rate multiplier.
#' @return list with cohort midpoints `x`, rates `f`, covariance `V`, `age`, and `status` of each cohort at `age`
#'   (`observed`, `forward` or `backward`; see [cohort_status()]).
#' @export
apc_fcp <- function(M, age = NULL, per = 1e5) {
  G <- M$G; r <- apc_ref(M)
  if (is.null(age)) age <- G$age[r["a"]]
  stopifnot(age >= min(G$age), age <= max(G$age))
  j <- findInterval(age, G$age, all.inside = TRUE); w <- (age - G$age[j]) / (G$age[j + 1] - G$age[j])
  k0 <- apc_contrast(G, "fcp", replace(r, "a", j)); k1 <- apc_contrast(G, "fcp", replace(r, "a", j + 1))
  X <- (1 - w) * k0$X + w * k1$X
  f <- per * exp(drop(X %*% M$B[k0$i]))
  list(x = k0$x, f = f, V = X %*% M$V[k0$i, k0$i] %*% t(X) * tcrossprod(f), age = age,
       status = cohort_status(M, age)$status)
}

#' Observed or projected: status of each cohort at a reference age
#'
#' The fitted cohort pattern at `age` ([apc_fcp()]) is interpolated between the age groups whose midpoints bracket
#' `age` (one group when `age` is a midpoint; for 5-year groups and age 50, the groups 45-49 and 50-54). A cohort
#' is `observed` at `age` if it has a Lexis cell in one of these age groups. Otherwise its rate at `age` is a
#' projection by the model along the cohort: `forward` for a cohort still younger than these age groups at the
#' end of the data (extrapolated from younger ages), `backward` for a cohort already older at the start of the
#' data (extrapolated from older ages). One rule for figures, tables and text: label projected cohorts, and
#' report generation means with and without them (`observed_only` in [fcp_generations()] and [fcp_contrast()]).
#'
#' @param M fitted model ([apc_fit()]) or comparative fit ([capricorn()]): a list with the design `G`.
#' @param age reference age; `NULL` for the MATISSE default (central age group).
#' @return data.table, one row per cohort: midpoint `x`, `status` (`observed`, `forward`, `backward`), the
#'   youngest and oldest age group midpoints in which the cohort is observed (`age_min`, `age_max`) and its
#'   number of Lexis cells (`n_cells`).
#' @examples
#' R <- rates_chunk(matisse_example(7), 5, 5)   # all cancers, ages 35-84, 1995-2019, 5 x 5 cells
#' cohort_status(apc_fit(rates_fill(R)), age = 60)
#' @export
cohort_status <- function(M, age = NULL) {
  G <- M$G
  if (is.null(age)) age <- G$age[floor((G$A + 1) / 2)]   # the central age group, as apc_ref()
  stopifnot(age >= min(G$age), age <= max(G$age))
  j <- findInterval(age, G$age, all.inside = TRUE); w <- (age - G$age[j]) / (G$age[j + 1] - G$age[j])
  used <- G$age[c(j, j + 1)[c(w < 1, w > 0)]]                # age groups that enter the interpolation at `age`
  tol <- 1e-6 * G$d
  has <- outer(G$coh, G$age, function(x, a) vapply(x + a, function(p) any(abs(p - G$per) < tol), TRUE))   # cohort x age cells
  obs <- apply(has[, match(used, G$age), drop = FALSE], 1, any)
  amin <- apply(has, 1, function(h) min(G$age[h])); amax <- apply(has, 1, function(h) max(G$age[h]))
  status <- ifelse(obs, "observed", ifelse(amax < min(used), "forward", "backward"))
  stopifnot(all(obs | amax < min(used) | amin > max(used)))
  data.table(x = G$coh, status = status, age_min = amin, age_max = amax, n_cells = rowSums(has))
}

#' Social generations and proxy parents (JAMA Netw Open 2024)
#'
#' @return data.table: `generation`, first and last birth year `y0`, `y1`.
#' @export
social_generations <- function() data.table(
  generation = c("Greatest", "Silent", "Boomers", "GenX", "Millennials",
                 "Proxy parents of Boomers", "Proxy parents of GenX", "Proxy parents of Millennials"),
  y0 = c(1908, 1928, 1946, 1965, 1981, 1917, 1936, 1952),
  y1 = c(1927, 1945, 1964, 1980, 1996, 1944, 1960, 1976))

#' Contrasts between successive generations and with proxy parents (JAMA Netw Open 2024)
#'
#' @return data.table: `contrast`, `from` (reference) and `to` generation.
#' @export
generation_contrasts <- function() data.table(
  contrast = c("Silent vs Greatest", "Boomers vs Silent", "GenX vs Boomers",
               "Boomers vs proxy parents", "GenX vs proxy parents", "Proxy parents: Millennials vs GenX"),
  from = c("Greatest", "Silent", "Boomers", "Proxy parents of Boomers", "Proxy parents of GenX", "Proxy parents of GenX"),
  to   = c("Silent", "Boomers", "GenX", "Boomers", "GenX", "Proxy parents of Millennials"))

#' Mean fitted cohort pattern by generation
#'
#' Average of the FCP over the cohorts (midpoints) in each generation, with its covariance.
#'
#' @param fcp result of [apc_fcp()] or [fcp_sum()].
#' @param gens generations ([social_generations()]).
#' @param min_cohorts fewest cohorts for a generation mean.
#' @param observed_only `TRUE` to average only the cohorts observed at the reference age (`fcp$status`, see
#'   [cohort_status()]); `FALSE` (default) to average all the cohorts of each generation, projected ones included.
#' @return list with `table` (`generation`, `n`, `mean`, `se`) and the covariance `V` of the means.
#' @export
fcp_generations <- function(fcp, gens = social_generations(), min_cohorts = 2, observed_only = FALSE) {
  keep <- if (observed_only) { if (is.null(fcp$status)) stop("fcp has no status: use apc_fcp() of this version"); fcp$status == "observed" }
          else rep(TRUE, length(fcp$x))
  W <- t(vapply(seq_len(nrow(gens)), function(g) { m <- fcp$x >= gens$y0[g] & fcp$x <= gens$y1[g] & keep
    if (sum(m) >= min_cohorts) m / sum(m) else rep(NA_real_, length(m)) }, numeric(length(fcp$x))))
  V <- W %*% fcp$V %*% t(W)
  list(table = data.table(generation = gens$generation, n = rowSums(W > 0), mean = drop(W %*% fcp$f), se = sqrt(diag(V))),
       V = V, W = W)
}

#' Rate ratios between generations
#'
#' Ratio of mean FCPs (`to` / `from`) with delta-method confidence intervals, accounting for the
#' covariance between the two means.
#'
#' @inheritParams fcp_generations
#' @param contrasts contrasts ([generation_contrasts()]).
#' @param alpha 1 - confidence level.
#' @return data.table: `contrast`, `rr`, `lo`, `hi`, `log_rr`, `se`.
#' @export
fcp_contrast <- function(fcp, contrasts = generation_contrasts(), gens = social_generations(), alpha = 0.05, min_cohorts = 2,
                         observed_only = FALSE) {
  g <- fcp_generations(fcp, gens, min_cohorts, observed_only); m <- g$table$mean; V <- g$V
  a <- match(contrasts$from, gens$generation); b <- match(contrasts$to, gens$generation)
  lr <- log(m[b] / m[a]); se <- sqrt(V[cbind(b, b)] / m[b]^2 + V[cbind(a, a)] / m[a]^2 - 2 * V[cbind(a, b)] / (m[a] * m[b]))
  q <- qnorm(1 - alpha / 2)
  data.table(contrast = contrasts$contrast, rr = exp(lr), lo = exp(lr - q * se), hi = exp(lr + q * se), log_rr = lr, se = se)
}

#' Sum of fitted cohort patterns (sites combined)
#'
#' @param fcps list of results of [apc_fcp()] on the same cohorts (independent strata).
#' @return FCP list: `x`, summed `f` and `V`, and `age` and `status` (of the first FCP; same Lexis diagram).
#' @export
fcp_sum <- function(fcps) {
  stopifnot(all(vapply(fcps, function(z) identical(z$x, fcps[[1]]$x), TRUE)))
  list(x = fcps[[1]]$x, f = Reduce(`+`, lapply(fcps, `[[`, "f")), V = Reduce(`+`, lapply(fcps, `[[`, "V")), age = fcps[[1]]$age,
       status = fcps[[1]]$status)
}

#' Mean local drift of early- and late-onset ages
#'
#' Average of the local drifts (log scale) over the age groups below `cut` (early onset) and from
#' `cut` (late onset), and their difference, with delta-method confidence intervals.
#'
#' @param M fitted model ([apc_fit()]).
#' @param cut age separating early and late onset (compared with age group midpoints).
#' @param alpha 1 - confidence level.
#' @return data.table: `onset` (early, late, early - late), `n` age groups, `est`, `se` (log scale),
#'   `drift`, `lo`, `hi` (% per year; for the difference, % per year of the early / late ratio).
#' @export
ld_onset <- function(M, cut = 50, alpha = 0.05) {
  v <- apc_ef_vcov(M, "ld"); e <- v$x < cut
  stopifnot(any(e), any(!e))
  W <- rbind(e / sum(e), (!e) / sum(!e)); W <- rbind(W, W[1, ] - W[2, ])
  est <- drop(W %*% v$est); se <- sqrt(diag(W %*% v$V %*% t(W))); q <- qnorm(1 - alpha / 2)
  data.table(onset = c("early", "late", "early - late"), n = c(sum(e), sum(!e), NA), est = est, se = se,
             drift = 100 * (exp(est) - 1), lo = 100 * (exp(est - q * se) - 1), hi = 100 * (exp(est + q * se) - 1))
}

#' Pooled rate ratio (random effects, DerSimonian-Laird)
#'
#' @param log_rr,se log rate ratios and standard errors (one per stratum, e.g. country).
#' @param alpha 1 - confidence level.
#' @return one-row data.table: `k`, pooled `rr`, `lo`, `hi`, `tau2`, `i2`, fixed-effect `rr_fixed`.
#' @export
pool_rr <- function(log_rr, se, alpha = 0.05) {
  ok <- is.finite(log_rr) & is.finite(se) & se > 0; y <- log_rr[ok]; v <- se[ok]^2; k <- length(y)
  if (!k) return(data.table(k = 0L, rr = NA_real_, lo = NA_real_, hi = NA_real_, tau2 = NA_real_, i2 = NA_real_, rr_fixed = NA_real_))
  w <- 1 / v; fe <- sum(w * y) / sum(w); Q <- sum(w * (y - fe)^2)
  tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
  ws <- 1 / (v + tau2); re <- sum(ws * y) / sum(ws); s <- sqrt(1 / sum(ws)); q <- qnorm(1 - alpha / 2)
  data.table(k = k, rr = exp(re), lo = exp(re - q * s), hi = exp(re + q * s), tau2 = tau2,
             i2 = if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) else 0, rr_fixed = exp(fe))
}
