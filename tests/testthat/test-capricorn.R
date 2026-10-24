## NPH equals separate APC fits; each PH model equals a factor model with shared or stratum-specific
## age, period and cohort effects; stratum estimable functions and ratios match the single-model ones.

test_that("CAPRICORN, synthetic populations C (zero cells) vs A (1 x 1)", {
  Rs <- lapply(matisse_example(c(3, 1)), rates_fill); names(Rs) <- c("C", "A")
  S <- capricorn(Rs); G <- S$G
  sep <- lapply(Rs, apc_fit, method = "poisson")
  expect_equal(unname(S$fits$NPH$B), unname(unlist(lapply(sep, `[[`, "B"))), tolerance = 1e-6)
  expect_equal(S$fits$NPH$DEV, sum(sapply(sep, `[[`, "DEV")), tolerance = 1e-8)

  L <- rbindlist(lapply(names(Rs), function(s) rates_long(Rs[[s]])[order(per_lo, age_lo)][, s := s]))
  L[, `:=`(fa = factor(age_mid), fp = factor(per_mid), fc = factor(coh_mid), s = factor(s, names(Rs)))]
  dev <- function(fm) suppressWarnings(glm(fm, poisson, data = L, offset = log(offset)))$deviance
  ref <- c(PHL = dev(events ~ fa + fp + s:fc), PHT = dev(events ~ s:fa + fp + fc), PHX = dev(events ~ fa + s:fp + fc),
           PHA = dev(events ~ s + fa + fp + fc), NPH = dev(events ~ s:(fa + fp + fc)))
  expect_equal(S$aic[match(names(ref), model), DEV], unname(ref), tolerance = 1e-7)

  for (cp in c("lac", "cac", "ftt", "fcp", "ld")) {
    e <- cap_ef(S, cp)
    for (g in 1:2) expect_equal(e[type == "ef" & stratum == names(Rs)[g], value], apc_ef(sep[[g]], cp)$value, tolerance = 1e-6)
    a <- e[stratum == "C", value]; b <- e[stratum == "A", value]
    expect_equal(e[type == "ratio", value], if (cp == "ld") 100 * ((1 + a / 100) / (1 + b / 100) - 1) else a / b, tolerance = 1e-8)
  }
  expect_setequal(S$tests$global$test, c("PHL", "PHT", "PHX", "PHA"))
  expect_true(all(S$tests$composite$p >= 0 & S$tests$composite$p <= 1))
})

test_that("AICc picks the model that holds by construction: B vs A PH-X, C vs A PH-A, D vs A PH-L (5 x 5)", {
  R5 <- lapply(matisse_example(1:4), function(R) rates_fill(rates_chunk(R, 5, 5))); names(R5) <- LETTERS[1:4]
  for (x in list(c("B", "PHX"), c("C", "PHA"), c("D", "PHL"))) expect_equal(capricorn(R5[c(x[1], "A")])$best, x[2])
})

test_that("CAPRICORN with four strata, synthetic populations B, C, D and A (5 x 5)", {
  Rs <- lapply(matisse_example(c(2, 3, 4, 1)), function(R) rates_fill(rates_chunk(R, 5, 5))); names(Rs) <- c("B", "C", "D", "A")
  S <- capricorn(Rs); S1 <- capricorn(Rs, overdispersion = TRUE)
  expect_equal(S$K, 4); expect_equal(nrow(S$aic), 5)
  expect_equal(nrow(cap_ef(S, "fcp")), 4 * S$G$C + 3 * S$G$C)
  expect_gte(S1$phi, 1)
  expect_equal(S1$tests$global$stat, S$tests$global$stat / S1$phi)
})
