## MATISSE @apc in R: the New Age-Period-Cohort Model (Rosenberg, MATISSE toolbox;
## Rosenberg & Anderson, Cancer Epidemiol Biomarkers Prev 2011, NCI APC Web Tool).
## apc_design mirrors GEODE.m, apc_fit APC.m, apc_ef the estimable functions of macaroons.m.

#' Age-period-cohort design matrix (MATISSE GEODE.m)
#'
#' Columns: intercept, linear age and cohort, quadratic age/period/cohort, then higher-order
#' deviations orthogonal to intercept, linear and quadratic terms.
#'
#' @param R rates object with square cells, at least 3 ages and 3 periods.
#' @return list with design `X` (cells ordered age fastest, as `c(R$Events)`), dimensions, midpoints,
#'   block pointers (`pa`, `pp`, `pc`) and edge-cell pointers (`inc_a`, `inc_p`, `inc_c`).
#' @export
apc_design <- function(R) {
  A <- nrow(R$Events); P <- ncol(R$Events); C <- A + P - 1; N <- A * P
  d <- diff(R$ages)[1]
  stopifnot(d == diff(R$periods)[1], A >= 3, P >= 3)
  age <- head(R$ages, -1) + d / 2; per <- head(R$periods, -1) + d / 2; coh <- seq(per[1] - age[A], per[P] - age[1], by = d)
  as <- rep(age, P); ps <- rep(per, each = A); cs <- ps - as       # cells, age fastest (as c(R$Events))
  a00 <- age[1] - d / 2; p00 <- per[1] - d / 2; c00 <- p00 - a00
  qa2 <- as^2 - (d * A + 2 * a00) * as + (a00 + d * A / 2)^2 - d^2 / 12 * (A - 1) * (A + 1)
  qp2 <- ps^2 - (d * P + 2 * p00) * ps + (p00 + d * P / 2)^2 - d^2 / 12 * (P - 1) * (P + 1)
  qc2 <- cs^2 - (d * (A + P) + 2 * (c00 - d * A)) * cs + (c00 - d * A)^2 + (A + P) * d * mean(cs) -
         d^2 / 6 * (2 * A^2 + 3 * A * P + 2 * P^2 - 1)
  dev <- function(x, lev, q) {                                         # deviations, dropping 1st and last 2 levels
    D <- outer(x, lev, `==`) * 1
    X0 <- cbind(1, x - mean(x), q)
    D0 <- D - X0 %*% solve(crossprod(X0), crossprod(X0, D))
    D0[, seq_len(max(0, length(lev) - 3)) + 1, drop = FALSE]
  }
  X <- cbind(1, as - mean(as), cs - mean(cs), qa2, qp2, qc2, dev(as, age, qa2), dev(ps, per, qp2), dev(cs, coh, qc2))
  na <- A - 3; np <- P - 3; nc <- C - 3
  list(X = X, A = A, P = P, C = C, d = d, age = age, per = per, coh = coh,
       pa = 6 + seq_len(na), pp = 6 + na + seq_len(np), pc = 6 + na + np + seq_len(nc),
       inc_a = ((P - 1) * A + 1):N, inc_p = A * seq_len(P), inc_c = c(A * seq_len(P), P * A - seq_len(A - 1)))
}

#' Fit the New APC Model (MATISSE APC.m)
#'
#' Weighted least squares on log rates (MATISSE default) or Poisson regression. With overdispersion,
#' the covariance is scaled by max(1, deviance / df).
#'
#' @param R rates object without zero events (see [rates_fill()]).
#' @param method "wls" or "poisson".
#' @param overdispersion scale the covariance by the dispersion.
#' @return list with the design `G`, coefficients `B`, covariance `V`, dispersion `s2`, deviance `DEV`
#'   and `net_drift` (% per year).
#' @examples
#' M <- apc_fit(rates_fill(matisse_example(16)))
#' M$net_drift
#' @export
apc_fit <- function(R, method = c("wls", "poisson"), overdispersion = TRUE) {
  method <- match.arg(method)
  G <- apc_design(R)
  Y <- c(R$Events); O <- c(R$Offset); X <- G$X
  stopifnot(all(Y > 0), all(O > 0))
  sc <- sqrt(colSums(X^2)); Xs <- sweep(X, 2, sc, "/")                  # unit-norm columns: well conditioned at 1 x 1
  if (method == "wls") {
    f <- lm.wfit(Xs, log(Y / O), Y); q <- f$qr; DEV <- sum(Y * f$residuals^2)
  } else {
    f <- suppressWarnings(glm.fit(Xs, Y, offset = log(O), family = poisson(), control = glm.control(maxit = 100))); stopifnot(f$converged)
    q <- qr(sqrt(f$weights) * Xs); DEV <- f$deviance
  }
  stopifnot(q$rank == ncol(X))
  B <- f$coefficients / sc
  V <- chol2inv(qr.R(q))[order(q$pivot), order(q$pivot)] / tcrossprod(sc)   # covariance via QR, not solve(X'WX)
  s2 <- if (overdispersion) max(1, DEV / (length(Y) - ncol(X))) else 1
  list(R = R, G = G, B = B, V = s2 * V, s2 = s2, DEV = DEV, method = method,
       net_drift = 100 * (exp(B[3]) - 1))                                # % per year
}

#' Reference age, period and cohort
#'
#' MATISSE defaults are the central age group and period, and the cohort through them.
#'
#' @param M fitted model (or `list(G = design)`).
#' @param ref `NULL` (defaults) or c(age, period, cohort) midpoints.
#' @return integer indices `a`, `p`, `c`.
#' @export
apc_ref <- function(M, ref = NULL) {
  G <- M$G
  if (is.null(ref)) { ia <- floor((G$A + 1) / 2); ip <- floor((G$P + 1) / 2); return(c(a = ia, p = ip, c = ip - ia + G$A)) }
  i <- c(a = match(ref[1], G$age), p = match(ref[2], G$per), c = match(ref[3], G$coh))
  if (anyNA(i)) stop("reference values must be age, period and cohort midpoints of the model")
  i
}

#' Estimable functions of the APC model (MATISSE macaroons.m)
#'
#' With delta-method confidence intervals:
#' * `lac` longitudinal age curve: rates by age in the reference cohort, period deviations removed;
#' * `cac` cross-sectional age curve: rates by age in the reference period, cohort deviations removed;
#' * `ftt` fitted temporal trends: rates by period at the reference age;
#' * `fcp` fitted cohort pattern: rates by cohort at the reference age;
#' * `ld` local drifts: annual % change by age (net drift plus the slope of the cohort deviations).
#'
#' @param M fitted model ([apc_fit()]).
#' @param comp estimable function.
#' @param ref reference values (see [apc_ref()]).
#' @param alpha 1 - confidence level.
#' @param per rate multiplier.
#' @return data.table: `ef`, `ref`, `x` (age, period or cohort midpoint), `value`, `lo`, `hi`.
#' @export
apc_ef <- function(M, comp = c("lac", "cac", "ftt", "fcp", "ld"), ref = NULL, alpha = 0.05, per = 1e5) {
  comp <- match.arg(comp)
  k  <- apc_contrast(M$G, comp, apc_ref(M, ref))
  ef <- drop(k$X %*% M$B[k$i])
  se <- sqrt(diag(k$X %*% M$V[k$i, k$i] %*% t(k$X)))
  q  <- qnorm(1 - alpha / 2); tr <- apc_scale(comp, per)
  data.table(ef = comp, ref = k$ref, x = k$x, value = tr(ef), lo = tr(ef - q * se), hi = tr(ef + q * se))
}

#' Several estimable functions, long format
#'
#' @inheritParams apc_ef
#' @param ... passed to [apc_ef()].
#' @export
apc_efs <- function(M, comp = c("lac", "cac", "ftt", "fcp", "ld"), ref = NULL, ...)
  rbindlist(lapply(comp, apc_ef, M = M, ref = ref, ...))

#' Scale of an estimable function: log scale to rate per `per` (basic EFs) or annual % change (ld)
#'
#' @param comp estimable function.
#' @param per rate multiplier.
#' @return a function.
#' @export
apc_scale <- function(comp, per = 1e5) if (comp == "ld") function(v) 100 * (exp(v) - 1) else function(v) per * exp(v)

#' Contrast of one estimable function
#'
#' EF = X %*% B\[i\] on the log scale, at reference cells `r`.
#'
#' @param G design ([apc_design()]).
#' @param comp estimable function.
#' @param r reference indices ([apc_ref()]).
#' @return list with `x`, `ref`, coefficient indices `i` and contrast matrix `X`.
#' @export
apc_contrast <- function(G, comp, r) {
  X <- G$X
  ab <- mean(G$age); pb <- mean(G$per); cb <- mean(G$coh)
  aref <- G$age[r["a"]]; pref <- G$per[r["p"]]; cref <- G$coh[r["c"]]
  rep_row <- function(x, n) matrix(x, n, length(x), byrow = TRUE)
  XA <- X[G$inc_a, c(4, G$pa), drop = FALSE]
  XAr <- X[G$inc_a[r["a"]], c(4, G$pa)]
  switch(comp,
    lac = list(x = G$age, ref = cref, i = c(1, 2, 3, 4, G$pa, 6, G$pc),
               X = cbind(1, G$age - ab, cref - cb, XA, rep_row(X[G$inc_c[r["c"]], c(6, G$pc)], G$A))),
    cac = list(x = G$age, ref = pref, i = c(1, 2, 3, 4, G$pa, 5, G$pp),
               X = cbind(1, G$age - ab, (pref - pb) - (G$age - ab), XA, rep_row(X[G$inc_p[r["p"]], c(5, G$pp)], G$A))),
    ftt = list(x = G$per, ref = aref, i = c(1, 2, 3, 5, G$pp, 4, G$pa),
               X = cbind(1, aref - ab, (G$per - pb) - (aref - ab), X[G$inc_p, c(5, G$pp), drop = FALSE], rep_row(XAr, G$P))),
    fcp = list(x = G$coh, ref = aref, i = c(1, 2, 3, 6, G$pc, 4, G$pa),
               X = cbind(1, aref - ab, G$coh - cb, X[G$inc_c, c(6, G$pc), drop = FALSE], rep_row(XAr, G$C))),
    ld  = {
      DP <- 12 / (G$d * (G$P - 1) * G$P * (G$P + 1)) * (seq_len(G$P) - (G$P + 1) / 2)   # least-squares slope over P cohorts
      K <- matrix(0, G$A, G$C)
      for (j in seq_len(G$A)) K[j, (G$A - j + 1):(G$A - j + G$P)] <- DP
      list(x = G$age, ref = NA_real_, i = c(6, G$pc, 3), X = cbind(K %*% X[G$inc_c, c(6, G$pc), drop = FALSE], 1))
    })
}
