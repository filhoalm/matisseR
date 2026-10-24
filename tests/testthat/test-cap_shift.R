## cap_shift() on a synthetic reference (annual data, ages 30-84, 1950-2019, Poisson counts from an APC surface with
## cohort and period curvature): a target equal to the reference has lag 0 and a copy with years moved by +10 has lag
## +10, for each model; a lagged target with cohort curvature 0.6 x the reference's is recovered by CohScale; profile
## deviances are those of the matching factor models; interval rule; input validation.

set.seed(1)
ref <- CJ(age = 7:17, year = 1950:2019)
ref[, `:=`(py = 1e6, cases = rpois(.N, 1e6 * exp(-12 + 0.09 * (5 * age - 2.5) + 0.15 * sin((year - 1950) / 7) +
                                                 0.6 * exp(-((year + 0.5 - (5 * age - 2.5) - 1925) / 15)^2))))]
tgt <- copy(ref)[, year := year + 10]
models <- c("CohScale", "CohCurv", "PHL", "PHT")

test_that("self-lag is 0 and a copy shifted by +10 has lag +10, for each model", {
  for (m in models) {
    s0 <- cap_shift(ref, ref, -3:3, 7:17, model = m)
    s1 <- cap_shift(tgt, ref, 7:13, 7:17, model = m)
    for (x in list(list(s = s0, L = 0), list(s = s1, L = 10))) {
      sm <- x$s$summary
      expect_lt(abs(sm$shift - x$L), 0.5)
      expect_true(sm$lo <= x$L && x$L <= sm$hi && sm$lo <= sm$shift && sm$shift <= sm$hi)
      expect_lt(abs(x$s$profile[shift == x$L, lrt]), 1e-6)                    # identical strata: model = NPH
      expect_true(all(x$s$profile[shift != x$L, lrt] > 1))
      expect_false(sm$edge_lo || sm$edge_hi || sm$min_out)
      expect_equal(sm$kappa, if (m == "CohScale") 1 else NA_real_, tolerance = 1e-3)
    }
    s1d <- cap_shift(tgt, ref, 7:13, 7:17, model = m, stat = "dev")$summary
    expect_lt(abs(s1d$shift - 10), 0.5)
    expect_equal(s1$window, c(1967, 2026))                                    # reference at every shift
    expect_equal(s0$window, c(1957, 2016))
    expect_equal(s1$summary$periods, 12)
  }
})

test_that("CohScale recovers a lag of 10 when the target cohort curvature is 0.6 x the reference's", {
  t6 <- CJ(age = 7:17, year = 1960:2029)[, `:=`(py = 5e5, c = year + 0.5 - (5 * age - 2.5))]
  t6[, cases := rpois(.N, 5e5 * exp(-12.5 + 0.1 * (5 * age - 2.5) + 0.1 * cos((year - 1950) / 5) +
                                      0.6 * 0.6 * exp(-((c - 10 - 1925) / 15)^2)))]
  s <- cap_shift(t6, ref, 5:15, 7:17)$summary
  expect_lt(abs(s$shift - 10), 0.3); expect_lt(abs(s$kappa - 0.6), 0.05)
  expect_true(s$lo <= 10 && 10 <= s$hi); expect_lt(s$misfit, 3)
  s1 <- cap_shift(t6, ref, 5:15, 7:17, model = "CohCurv")$summary         # amplitude difference: biased, flagged
  expect_gt(s1$shift, 10.5); expect_gt(s1$misfit, 20)
})

test_that("profile deviances equal factor models with shared, scaled and stratum-specific effects", {
  t2 <- copy(ref)[, year := year + 4][, cases := rpois(.N, cases * exp(0.1 * cos((year - 1950) / 5) +
                                                                          0.005 * (year - 1985)))]
  fm <- list(NPH = events ~ s:(fa + fp + fc), CohCurv = events ~ s:(fa + fp) + fc, PHL = events ~ fa + fp + s:fc,
             PHT = events ~ s:fa + fp + fc, CohScale = events ~ s:(fa + fp) + fc:w)
  pr <- list()
  for (m in models) {
    x <- cap_shift(t2, ref, 0:8, 7:17, model = m); P <- pr[[m]] <- x$profile; i <- which.min(P$lrt)
    L <- rbindlist(lapply(names(x$Rs), function(s) rates_long(x$Rs[[s]])[order(per_lo, age_lo)][, s := s]))
    L[, `:=`(fa = factor(age_mid), fp = factor(per_mid), fc = factor(coh_mid), s = factor(s, names(x$Rs)),
             w = ifelse(s == "target", P$kappa[i], 1))]
    dev <- function(f) suppressWarnings(glm(f, poisson, data = L, offset = log(offset)))$deviance
    expect_equal(c(P$dev[i], P$dev_nph[i]), c(dev(fm[[m]]), dev(fm$NPH)), tolerance = 1e-7)
    expect_equal(x$phi, max(1, P$dev_nph[i] / (2 * length(x$Rs$target$Events) - P$df_nph[i])))
    expect_true(x$summary$misfit >= 0 && x$summary$misfit <= P$lrt[i] / (P$df_nph[i] - P$df[i]))
  }
  expect_true(all(pr$CohScale$dev <= pr$CohCurv$dev + 1e-6))                 # kappa free: CohScale <= CohCurv
  S <- capricorn(x$Rs)
  expect_equal(P$dev[i], S$fits$PHT$DEV, tolerance = 1e-7)
})

test_that("interval rule: parabola next to the minimum, linear further out, phi, edges", {
  s <- -5:10
  r <- sqrt(qchisq(0.95, 1) / 50)
  a <- matisseR:::shift_ci(s, 50 * (s - 2.1)^2 + 100, 1)                    # sharp: inside the grid step
  expect_equal(c(a$best, a$vmin, a$lo, a$hi), c(2.1, 100, 2.1 - r, 2.1 + r))
  expect_true(a$sub_grid); expect_false(a$min_out)
  a <- matisseR:::shift_ci(s, 50 * (s - 2.4)^2, 1)                          # the grid minimum is not forced inside
  expect_equal(c(a$best, a$lo, a$hi), c(2.4, 2.4 - r, 2.4 + r)); expect_true(a$min_out)
  b <- matisseR:::shift_ci(s, 0.5 * (s - 2.3)^2, 1); w <- sqrt(qchisq(0.95, 1) / 0.5)
  expect_equal(b$best, 2.3); expect_lt(max(abs(c(b$lo, b$hi) - (2.3 + c(-w, w)))), 0.1); expect_false(b$sub_grid)
  b4 <- matisseR:::shift_ci(s, 0.5 * (s - 2.3)^2, 4)
  expect_lt(max(abs(c(b4$lo, b4$hi) - (2.3 + 2 * c(-w, w)))), 0.25)
  e <- matisseR:::shift_ci(s, 10 * s, 1)                                    # minimum at the grid edge
  expect_equal(c(e$best, e$lo), c(-5, -5)); expect_true(e$edge_lo); expect_false(e$edge_hi)
  g <- matisseR:::shift_ci(s, c(0, 10, 0.5, 10, rep(20, 12)), 1)             # two minima
  expect_true(g$disjoint); expect_equal(g$lo, -5)
  v <- matisseR:::shift_ci(9:14, c(50, 0, 3, 2, 3, 30), 1)                  # grid points judged on the observed minimum
  expect_equal(v$inside, c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)); expect_equal(v$hi, 13 + (qchisq(0.95, 1) - 3) / 27)
})

test_that("input validation", {
  expect_error(cap_shift(ref[, !"py"], ref, -3:3, 7:17), "lacks column")
  expect_error(cap_shift(as.matrix(ref), ref, -3:3, 7:17), "data.frame")
  expect_error(cap_shift(ref, ref, c(0, 0.5, 1), 7:17), "whole numbers")
  expect_error(cap_shift(ref, ref, 0:1, 7:17), "at least 3")
  expect_error(cap_shift(ref, ref, -3:3, c(7, 9, 11)), "consecutive")
  expect_error(cap_shift(ref, ref, -3:3, 7:17, model = "PHX"))
  expect_error(cap_shift(ref, ref, -3:3, 7:17, level = 1.5), "level")
  expect_error(cap_shift(ref[age != 9], ref, -3:3, 7:17), "lacks age group\\(s\\) 9")
  expect_error(cap_shift(copy(ref)[, age := 5 * age], ref, -3:3, 7:17), "lacks age group")   # ages coded in years
  expect_error(cap_shift(copy(ref)[, age := as.character(age)], ref, -3:3, 7:17), "numeric")
  expect_equal(cap_shift(copy(ref)[, ages := 1], ref, -3:3, 7:17, model = "CohCurv")$summary,
               cap_shift(ref, ref, -3:3, 7:17, model = "CohCurv")$summary)  # a column 'ages' does not mask the argument
  expect_error(cap_shift(ref, ref, -28:28, 7:17), "at least 15")            # 1978-1991: 14 years
  expect_error(cap_shift(ref[!(age == 9 & year == 1985)], ref, -3:3, 7:17, years = 1980:1994), "at least 15")
  expect_equal(cap_shift(ref, ref, -3:3, 7:17, years = 1990:2010)$window, c(1991, 2010))
})
