## MATISSE @capricorn in R: hypothesis-based comparative age-period-cohort analysis
## (Rosenberg et al., BMC Med Res Methodol 2023, PMID 37853346). K rates objects over
## the same Lexis diagram are fitted jointly by Poisson regression. Proportional (PH)
## models make blocks of the New APC Model equal across strata:
##   PH-L  rates parallel along cohorts (diagonals): LAT, age and period deviations equal
##   PH-T  rates parallel over time within ages (rows): drift, period and cohort deviations equal
##   PH-X  cross-sectional age curves parallel (columns): CAT, age and cohort deviations equal
##   PH-A  whole Lexis diagrams parallel: all equal except the intercept
##   NPH   no constraint (separate APC models, fitted jointly)
## Mirrors CAPRICORN.m (fits, tofu expansion, global/homogeneity/composite tests, AICc) and EF.m.

#' Coefficient blocks of the APC design (GEODE pointers)
#'
#' 1 intercept, 2 LAT, 3 net drift, 4-6 quadratic age/period/cohort, 7-9 higher-order age/period/cohort deviations.
#'
#' @param G design ([apc_design()]).
#' @export
cap_blocks <- function(G) list(1, 2, 3, 4, 5, 6, G$pa, G$pp, G$pc)

#' Coefficient maps of a comparative model
#'
#' `E` (equal across strata) and `F` (free) maps, M x k: `X %*% E` are the shared columns. Models
#' `PHL`, `PHT`, `PHX`, `PHA`, `NPH`; single blocks (`LAT`, `NetDrift`, `QuadAge`, `QuadPer`, `QuadCoh`,
#' `HiOrdAge`, `HiOrdPer`, `HiOrdCoh`, `CAT`) make one block equal, for the homogeneity tests.
#'
#' @param G design ([apc_design()]).
#' @param model model or block name.
#' @export
cap_map <- function(G, model) {
  I <- diag(ncol(G$X)); Pt <- cap_blocks(G)
  sel <- function(b) I[, unlist(Pt[b]), drop = FALSE]
  one <- function(b) list(E = sel(b), F = sel(setdiff(1:9, b)))
  switch(model,
    NPH = list(E = I[, 0, drop = FALSE], F = I),
    PHL = list(E = sel(c(2, 4, 5, 7, 8)), F = sel(c(1, 3, 6, 9))),
    PHT = list(E = sel(c(3, 5, 6, 8, 9)), F = sel(c(1, 2, 4, 7))),
    PHX = list(E = cbind(I[, 2] - I[, 3], sel(c(4, 6, 7, 9))), F = cbind(I[, 1], I[, 2] + I[, 3], sel(c(5, 8)))),
    PHA = list(E = sel(2:9), F = sel(1)),
    LAT = one(2), NetDrift = one(3), QuadAge = one(4), QuadPer = one(5), QuadCoh = one(6),
    HiOrdAge = one(7), HiOrdPer = one(8), HiOrdCoh = one(9),
    CAT = list(E = cbind(I[, 2] - I[, 3]), F = cbind(I[, 1], I[, 2] + I[, 3], sel(4:9))))
}

#' Joint Poisson fit of a comparative model
#'
#' Coefficients and covariance expanded to one full APC set per stratum (MATISSE tofu).
#'
#' @param G design ([apc_design()]).
#' @param Y,O stacked events and offsets (one Lexis diagram per stratum, cells ordered age fastest).
#' @param K number of strata.
#' @param m coefficient map ([cap_map()]).
#' @param eta0 starting linear predictor (the NPH fit).
#' @export
cap_fit <- function(G, Y, O, K, m, eta0 = NULL) {
  Z  <- cbind(kronecker(rep(1, K), G$X %*% m$E), kronecker(diag(K), G$X %*% m$F))
  sc <- sqrt(colSums(Z^2)); Z <- sweep(Z, 2, sc, "/")                   # unit-norm columns: well conditioned at 1 x 1
  f  <- suppressWarnings(glm.fit(Z, Y, etastart = eta0, offset = log(O), family = poisson(),   # filled cells are non-integer
                                 control = glm.control(maxit = 100)))
  stopifnot(f$converged)
  q  <- qr(sqrt(f$weights) * Z); stopifnot(q$rank == ncol(Z))           # covariance via QR, not solve(X'WX)
  V  <- chol2inv(qr.R(q))[order(q$pivot), order(q$pivot)] / tcrossprod(sc)
  CM <- cbind(kronecker(rep(1, K), m$E), kronecker(diag(K), m$F))
  list(B = drop(CM %*% (f$coefficients / sc)), V = CM %*% V %*% t(CM), DEV = f$deviance, df = ncol(Z), eta = f$linear.predictors)
}

#' Comparative APC analysis (MATISSE CAPRICORN.m)
#'
#' Fits the PH-L, PH-T, PH-X, PH-A and NPH models jointly to K rates objects over the same Lexis
#' diagram; global, homogeneity and composite (Bonferroni) likelihood-ratio tests against NPH, and AICc.
#'
#' @param Rs named list of rates objects (zero events filled); the last stratum is the reference.
#' @param overdispersion `FALSE` as in MATISSE (Poisson tests); `TRUE` divides deviances by the NPH
#'   dispersion max(1, deviance / df) (quasi-Poisson tests, QAICc) and scales the EF covariances.
#' @return list with the design, fits, `tests` (`global`, `homogeneity`, `composite`), `aic` and `best` model.
#' @examples
#' \donttest{
#' Rs <- lapply(matisse_example(9:10), rates_fill)   # thyroid, male, NHB vs NHW
#' names(Rs) <- c("NHB", "NHW")
#' S <- capricorn(Rs)
#' S$aic
#' }
#' @export
capricorn <- function(Rs, overdispersion = FALSE) {
  K <- length(Rs); R1 <- Rs[[1]]
  stopifnot(K >= 2, !is.null(names(Rs)),
            all(vapply(Rs, function(R) identical(R$ages, R1$ages) && identical(R$periods, R1$periods), TRUE)))
  G <- apc_design(R1); N <- nrow(G$X)
  Y <- unlist(lapply(Rs, function(R) c(R$Events))); O <- unlist(lapply(Rs, function(R) c(R$Offset)))
  stopifnot(all(Y > 0), all(O > 0))
  mods <- c("PHL", "PHT", "PHX", "PHA", "NPH")
  comps <- c("LAT", "NetDrift", "QuadAge", "QuadPer", "QuadCoh", "HiOrdAge", "HiOrdPer", "HiOrdCoh", "CAT")
  nph  <- cap_fit(G, Y, O, K, cap_map(G, "NPH"))
  fits <- sapply(c(mods, comps), function(m) if (m == "NPH") nph else cap_fit(G, Y, O, K, cap_map(G, m), nph$eta), simplify = FALSE)
  phi <- if (overdispersion) max(1, fits$NPH$DEV / (N * K - fits$NPH$df)) else 1
  lrt <- function(m) { s <- (fits[[m]]$DEV - fits$NPH$DEV) / phi; d <- fits$NPH$df - fits[[m]]$df
    data.table(test = m, stat = s, df = d, p = if (d > 0) pchisq(s, d, lower.tail = FALSE) else NA_real_) }
  glob <- rbindlist(lapply(mods[1:4], lrt))
  hom  <- rbindlist(lapply(comps, lrt))
  ph <- setNames(hom$p, hom$test); pg <- setNames(glob$p, glob$test)
  bonf <- function(p, k) min(k * min(p, na.rm = TRUE), 1)
  comp <- data.table(test = c("PHL", "PHT", "PHX", "PHA", "Par-LAC", "Par-CAC", "Par-FTT", "Par-FCP"), p = c(
    bonf(c(ph[c("LAT", "QuadAge", "QuadPer", "HiOrdAge", "HiOrdPer")], pg["PHL"]), 6),
    bonf(c(ph[c("NetDrift", "QuadCoh", "HiOrdCoh", "QuadPer", "HiOrdPer")], pg["PHT"]), 6),
    bonf(c(ph[c("CAT", "QuadAge", "HiOrdAge", "QuadCoh", "HiOrdCoh")], pg["PHX"]), 6),
    bonf(c(ph[comps], pg["PHA"]), 10),
    bonf(ph[c("LAT", "QuadAge", "HiOrdAge")], 3), bonf(ph[c("CAT", "QuadAge", "HiOrdAge")], 3),
    bonf(ph[c("NetDrift", "QuadPer", "HiOrdPer")], 3), bonf(ph[c("NetDrift", "QuadCoh", "HiOrdCoh")], 3)))
  aic <- data.table(model = mods, df = vapply(fits[mods], `[[`, 0, "df"), DEV = vapply(fits[mods], `[[`, 0, "DEV"))
  aic[, `:=`(AICc = DEV / phi + 2 * (df + df * (df + 1) / (N * K - df - 1)))][, `:=`(delta = AICc - min(AICc), rank = frank(AICc))]
  list(G = G, K = K, labels = names(Rs), phi = phi, fits = fits[mods], tests = list(global = glob, homogeneity = hom, composite = comp),
       aic = aic[], best = aic[rank == 1, model])
}

#' Estimable functions of a comparative analysis (MATISSE capricorn/EF.m)
#'
#' Estimable functions by stratum, and contrasts of each stratum vs the last: rate ratios, or for
#' local drifts the ratio of annual changes as % per year.
#'
#' @param S result of [capricorn()].
#' @param comp estimable function (see [apc_ef()]).
#' @param model `NPH`, `PHL`, `PHT`, `PHX` or `PHA`.
#' @inheritParams apc_ef
#' @return data.table: `model`, `ef`, `stratum`, `type` (`ef` or `ratio`), `ref`, `x`, `value`, `lo`, `hi`.
#' @export
cap_ef <- function(S, comp = c("lac", "cac", "ftt", "fcp", "ld"), model = "NPH", ref = NULL, alpha = 0.05, per = 1e5) {
  comp <- match.arg(comp)
  G <- S$G; M <- ncol(G$X); K <- S$K; f <- S$fits[[model]]
  k <- apc_contrast(G, comp, apc_ref(list(G = G), ref))
  q <- qnorm(1 - alpha / 2); tr <- apc_scale(comp, per)
  ix <- function(g) (g - 1) * M + k$i
  one <- function(g, h = NULL) {
    D <- if (is.null(h)) cbind(k$X) else cbind(k$X, -k$X); j <- c(ix(g), if (!is.null(h)) ix(h))
    v <- drop(D %*% f$B[j]); se <- sqrt(pmax(0, diag(D %*% (S$phi * f$V[j, j]) %*% t(D))))
    t2 <- if (is.null(h)) tr else if (comp == "ld") apc_scale("ld") else exp
    data.table(model, ef = comp, stratum = if (is.null(h)) S$labels[g] else paste(S$labels[g], "vs", S$labels[h]),
               type = if (is.null(h)) "ef" else "ratio", ref = k$ref, x = k$x, value = t2(v), lo = t2(v - q * se), hi = t2(v + q * se))
  }
  rbind(rbindlist(lapply(seq_len(K), one)), rbindlist(lapply(seq_len(K - 1), one, h = K)))
}
