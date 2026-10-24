## Fitted cohort share on synthetic diagrams of population A: examples 1, 6 and 7 (melanoma-like, lung-like and all
## cancers, female) share their person-years, and all cancers contains the melanoma and lung cases (5 x 5).

ex <- lapply(matisse_example(c(1, 6, 7)), rates_chunk, age_block = 5, per_block = 5)
Ra <- ex[[3]]
Rc <- ex[[1]]

test_that("direct estimator on all / all is exactly 1, with and without zero cells", {
  for (R in list(Ra, matisse_example(8))) {                         # 8: thyroid-like, male, C, 1 x 1, mostly zero cells
    s <- suppressWarnings(fcp_share(R, R, 60, "direct"))
    expect_lt(max(abs(s$share$share - 1)), 1e-10)
    expect_lt(max(abs(s$generations[!is.na(mean), mean] - 1)), 1e-10)
    expect_true(all(s$share$lo < 1 & s$share$hi > 1))
  }
  expect_equal(suppressWarnings(fcp_share(matisse_example(8), matisse_example(8), 60, "direct"))$check$n_filled,
               matisse_example()$zero_cells[8])
})

test_that("shares of a subset lie in (0, 1); direct = FCP of rates_share(), ratio = FCP ratio", {
  s <- fcp_share(Rc, Ra, 60); d <- s$share[estimator == "direct"]; r <- s$share[estimator == "ratio"]
  expect_true(all(s$share$share > 0 & s$share$share < 1 & s$share$lo > 0))
  expect_equal(d$share, apc_fcp(apc_fit(rates_share(Rc, Ra)), 60, per = 1)$f, tolerance = 1e-10)
  fc <- apc_fcp(apc_fit(rates_fill(Rc)), 60); fa <- apc_fcp(apc_fit(rates_fill(Ra)), 60)
  expect_equal(r$share, fc$f / fa$f, tolerance = 1e-10)
  expect_equal(r$se_log, sqrt(diag(fc$V) / fc$f^2 + diag(fa$V) / fa$f^2), tolerance = 1e-10)
  expect_equal(r$hi / r$share, exp(qnorm(0.975) * r$se_log), tolerance = 1e-10)
  expect_true(all(d$share >= r$lo & d$share <= r$hi & r$share >= d$lo & r$share <= d$hi))
  expect_equal(s$generations[estimator == "direct", mean], fcp_generations(s$fcp$direct)$table$mean)
  expect_equal(s$check$frac_half, rep(mean(Rc$Events / Ra$Events > 0.5), 2))
  expect_equal(s$check$max_cell, c(max(exp(s$models$direct$G$X %*% s$models$direct$B)), NA))
  expect_equal(s$check$max_fcp_all, c(max(d$share), max(r$share)))
  expect_true(all(s$check$max_fcp_observed <= s$check$max_fcp_all & s$check$max_fcp_all < 0.95))
})

test_that("the warning ignores cohorts projected to ref_age", {
  a <- rep(head(Ra$ages, -1) + 2.5, ncol(Ra$Events)); co <- rep(head(Ra$periods, -1) + 2.5, each = nrow(Ra$Events)) - a
  Rx <- Ra; Rx$Events[] <- Ra$Events * 0.4 * exp(0.02 * (a - 60)) * ifelse(co > 1961, 2.5, 1)   # observed at 55-64: <= 1961
  expect_no_warning(s <- fcp_share(Rx, Ra, 60, "direct"))
  expect_equal(s$check$max_cell, max(Rx$Events / Ra$Events), tolerance = 1e-10)
  expect_equal(c(s$check$max_fcp_observed, s$check$max_fcp_all), c(0.4, 1), tolerance = 1e-10)
  expect_warning(fcp_share(Rx, Ra, 60, "direct", warn_above = 0.8), "above 0.8")
})

test_that("warning above 0.95; rates_share() checks the subset", {
  R97 <- Ra; R97$Events <- 0.97 * Ra$Events
  expect_warning(s <- fcp_share(R97, Ra, 60, "direct"), "above 0.95")
  expect_equal(s$share$share, rep(0.97, nrow(s$share)), tolerance = 1e-10)
  expect_error(rates_share(Ra, Rc))
})

test_that("capricorn() takes cause deaths per all-cause death directly", {
  Rs <- lapply(list(melanoma = ex[[1]], lung = ex[[2]]), rates_share, R_all = Ra)
  S <- capricorn(Rs); e <- cap_ef(S, "fcp", per = 1)
  for (g in names(Rs)) expect_equal(e[stratum == g, value], apc_ef(apc_fit(Rs[[g]], "poisson"), "fcp", per = 1)$value, tolerance = 1e-6)
  expect_equal(e[type == "ratio", value], e[stratum == "melanoma", value] / e[stratum == "lung", value], tolerance = 1e-8)
})
