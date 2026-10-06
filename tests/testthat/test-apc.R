## The design spans the full APC model, the fits equal factor models, and each estimable function
## equals the fitted log rates along its reference line with the other deviations removed.

cells <- function(R) rates_long(R)[order(per_lo, age_lo)]           # column-major, as c(R$Events)

for (cell in c("1x1", "5x5")) test_that(paste("APC model and estimable functions,", cell), {
  R <- rates_fill(matisse_example(16))
  if (cell == "5x5") R <- rates_fill(rates_chunk(matisse_example(16), 5, 5))
  L <- cells(R)
  for (m in c("wls", "poisson")) {
    M <- apc_fit(R, m); G <- M$G; u <- drop(G$X %*% M$B)
    expect_equal(qr(G$X)$rank, G$A + G$P + G$C - 3)
    f <- if (m == "wls") fitted(lm(log(events / offset) ~ factor(age_mid) + factor(per_mid) + factor(coh_mid), weights = events, data = L))
         else predict(glm(events ~ factor(age_mid) + factor(per_mid) + factor(coh_mid) + offset(log(offset)), poisson, data = L)) - log(L$offset)
    expect_equal(unname(u), unname(f), tolerance = 1e-7)

    E <- apc_efs(M); r <- apc_ref(M); lu <- log(1e5) + u
    pd <- drop(G$X[, c(5, G$pp)] %*% M$B[c(5, G$pp)]); cd <- drop(G$X[, c(6, G$pc)] %*% M$B[c(6, G$pc)])
    chk <- function(e, on, xv, dev) { v <- E[ef == e]; max(abs(log(v$value[match(xv[on], v$x)]) - (lu - dev)[on])) }
    expect_lt(chk("lac", L$coh_mid == G$coh[r["c"]], L$age_mid, pd), 1e-8)
    expect_lt(chk("cac", L$per_mid == G$per[r["p"]], L$age_mid, cd), 1e-8)
    expect_lt(chk("ftt", L$age_mid == G$age[r["a"]], L$per_mid, cd), 1e-8)
    expect_lt(chk("fcp", L$age_mid == G$age[r["a"]], L$coh_mid, pd), 1e-8)

    sl <- sapply(seq_len(G$A), function(j) coef(lm(u[seq(j, G$A * G$P, G$A)] ~ G$per))[2])
    expect_equal(E[ef == "ld", value], unname(100 * (exp(sl) - 1)), tolerance = 1e-8)
    expect_equal(unname(M$net_drift), 100 * (exp(M$B[[3]]) - 1))
  }
})
