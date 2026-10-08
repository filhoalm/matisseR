M <- apc_fit(rates_fill(rates_chunk(matisse_example(16), 5, 5)))   # melanoma F NHW, 5 x 5

test_that("apc_fcp equals apc_ef at an age midpoint and interpolates log-linearly between midpoints", {
  G <- M$G; a <- G$age[4]
  f <- apc_fcp(M, a); e <- apc_ef(M, "fcp", ref = c(a, G$per[apc_ref(M)["p"]], G$coh[apc_ref(M)["c"]]))
  expect_equal(f$f, e$value, tolerance = 1e-10)
  expect_equal(sqrt(diag(f$V)) / f$f, (log(e$hi) - log(e$value)) / qnorm(0.975), tolerance = 1e-8)
  f60 <- apc_fcp(M, 60); f57 <- apc_fcp(M, 57.5); f62 <- apc_fcp(M, 62.5)
  expect_equal(log(f60$f), (log(f57$f) + log(f62$f)) / 2, tolerance = 1e-10)
  expect_equal(apc_ef_vcov(M, "fcp")$est, log(apc_fcp(M)$f / 1e5), tolerance = 1e-10)
})

test_that("generation means, contrasts and sums", {
  f <- apc_fcp(M, 60); g <- fcp_generations(f)
  m <- f$x >= 1946 & f$x <= 1964
  expect_equal(g$table[generation == "Boomers", mean], mean(f$f[m]))
  expect_equal(g$table[generation == "Boomers", se], sqrt(sum(f$V[m, m])) / sum(m))
  cc <- fcp_contrast(f)
  expect_equal(cc[contrast == "Boomers vs Silent", rr],
               g$table[generation == "Boomers", mean] / g$table[generation == "Silent", mean])
  s <- fcp_sum(list(f, f))
  expect_equal(s$f, 2 * f$f); expect_equal(fcp_contrast(s)$rr, cc$rr)
})

test_that("ld_onset averages local drifts below and above the cut", {
  ld <- apc_ef(M, "ld"); o <- ld_onset(M, 50); e <- ld$x < 50
  expect_equal(o$est[1], mean(log(1 + ld$value[e] / 100)), tolerance = 1e-10)
  expect_equal(o$est[2], mean(log(1 + ld$value[!e] / 100)), tolerance = 1e-10)
  expect_equal(o$est[3], o$est[1] - o$est[2])
  expect_equal(o$n[1:2], c(sum(e), sum(!e)))
  expect_true(all(o$se > 0)); expect_true(all(o$lo < o$drift & o$drift < o$hi))
})

test_that("pool_rr: fixed effect with homogeneous strata, DerSimonian-Laird otherwise", {
  p <- pool_rr(log(c(1.2, 1.2, 1.2)), c(0.1, 0.2, 0.1))
  expect_equal(p$rr, 1.2); expect_equal(p$tau2, 0); expect_equal(p$k, 3L)
  p <- pool_rr(log(c(0.8, 1.5)), c(0.05, 0.05))
  expect_gt(p$tau2, 0); expect_equal(p$rr_fixed, exp(mean(log(c(0.8, 1.5)))))
  expect_equal(pool_rr(numeric(0), numeric(0))$k, 0L)
})
