M <- apc_fit(rates_fill(rates_chunk(matisse_example(1), 5, 5)))    # synthetic melanoma-like, female, A, 5 x 5

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

test_that("cohort_status: observed at the reference age, or projected forward or backward", {
  R <- rates_chunk(matisse_example(7), 5, 5)
  M <- apc_fit(rates_fill(R))
  L <- data.table::as.data.table(rates_long(R))
  near <- function(a, b) abs(a - b) < 1e-9

  ## age 60 lies between the midpoints 57.5 and 62.5: a cohort is observed if it has a cell in either group
  s <- cohort_status(M, 60)
  ref <- L[, .(obs = any(near(age_mid, 57.5) | near(age_mid, 62.5)), amin = min(age_mid), amax = max(age_mid), n = .N), keyby = coh_mid]
  expect_equal(s$x, ref$coh_mid)
  expect_equal(s$status == "observed", ref$obs)
  expect_equal(s$age_min, ref$amin)
  expect_equal(s$age_max, ref$amax)
  expect_equal(s$n_cells, ref$n)
  expect_true(all(s[status == "forward", age_max] < 57.5))
  expect_true(all(s[status == "backward", age_min] > 62.5))
  expect_setequal(unique(s$status), c("observed", "forward", "backward"))

  ## at an age group midpoint only that group enters the interpolation
  s2 <- cohort_status(M, 62.5)
  expect_equal(s2$status == "observed", L[, any(near(age_mid, 62.5)), keyby = coh_mid]$V1)

  ## default reference age: the central age group, as apc_fcp()
  expect_equal(cohort_status(M)$status, apc_fcp(M)$status)
})

test_that("apc_fcp carries the status; observed_only restricts generation means and contrasts", {
  R <- rates_chunk(matisse_example(7), 5, 5)
  M <- apc_fit(rates_fill(R))
  f <- apc_fcp(M, 60)
  expect_identical(f$status, cohort_status(M, 60)$status)

  ## default unchanged: all cohorts of each generation, with or without a status
  f0 <- f[setdiff(names(f), "status")]
  expect_equal(fcp_generations(f)$table, fcp_generations(f0)$table)
  expect_equal(fcp_contrast(f), fcp_contrast(f0))
  expect_error(fcp_generations(f0, observed_only = TRUE))

  ## observed_only: the means use the observed cohorts only
  gens <- social_generations()
  g <- fcp_generations(f, observed_only = TRUE, min_cohorts = 1)
  for (i in seq_len(nrow(gens))) {
    k <- f$x >= gens$y0[i] & f$x <= gens$y1[i] & f$status == "observed"
    if (any(k)) expect_equal(g$table$mean[i], mean(f$f[k])) else expect_true(is.na(g$table$mean[i]))
  }
  ct <- fcp_contrast(f, observed_only = TRUE, min_cohorts = 1)
  cs <- generation_contrasts()
  m <- g$table$mean
  expect_equal(ct$rr, m[match(cs$to, gens$generation)] / m[match(cs$from, gens$generation)])

  ## sites combined keep the status
  expect_identical(fcp_sum(list(f, f))$status, f$status)
})

test_that("cohort_status works on a comparative fit (capricorn)", {
  R <- setNames(lapply(matisse_example(c(1, 2)), rates_chunk, 5, 5), c("A", "B"))
  S <- capricorn(lapply(R, rates_fill))
  expect_equal(cohort_status(S, 60), cohort_status(apc_fit(rates_fill(R[[1]])), 60))
})
