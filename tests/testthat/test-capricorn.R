## NPH equals separate APC fits; each PH model equals a factor model with shared or stratum-specific
## age, period and cohort effects; stratum estimable functions and ratios match the single-model ones.

test_that("CAPRICORN, thyroid NHB vs NHW (1 x 1)", {
  Rs <- lapply(matisse_example(9:10), rates_fill); names(Rs) <- c("NHB", "NHW")
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
    a <- e[stratum == "NHB", value]; b <- e[stratum == "NHW", value]
    expect_equal(e[type == "ratio", value], if (cp == "ld") 100 * ((1 + a / 100) / (1 + b / 100) - 1) else a / b, tolerance = 1e-8)
  }
  expect_setequal(S$tests$global$test, c("PHL", "PHT", "PHX", "PHA"))
  expect_true(all(S$tests$composite$p >= 0 & S$tests$composite$p <= 1))
})

test_that("CAPRICORN with four strata, thyroid API/HIS/NHB/NHW (5 x 5)", {
  Rs <- lapply(matisse_example(7:10), function(R) rates_fill(rates_chunk(R, 5, 5))); names(Rs) <- c("API", "HIS", "NHB", "NHW")
  S <- capricorn(Rs); S1 <- capricorn(Rs, overdispersion = TRUE)
  expect_equal(S$K, 4); expect_equal(nrow(S$aic), 5)
  expect_equal(nrow(cap_ef(S, "fcp")), 4 * S$G$C + 3 * S$G$C)
  expect_gte(S1$phi, 1)
  expect_equal(S1$tests$global$stat, S$tests$global$stat / S1$phi)
})
