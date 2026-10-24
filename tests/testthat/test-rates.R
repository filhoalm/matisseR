test_that("csv2rates reads every example; all cancers contains its sites", {
  idx <- matisse_example()
  expect_equal(nrow(idx), 8)
  for (k in idx$number) {
    R <- matisse_example(k)
    expect_equal(dim(R$Events), c(idx$ages[k], idx$periods[k]))
    expect_equal(R$ages[1], 35); expect_equal(R$periods[1], 1991)
    expect_true(all(R$Offset > 0)); expect_false(anyNA(R$Events))
    expect_equal(c(sum(R$Events), sum(R$Events == 0)), c(idx$cases[k], idx$zero_cells[k]))
  }
  ex <- matisse_example(c(1, 5, 6, 7))
  for (R in ex[1:3]) expect_identical(R$Offset, ex[[4]]$Offset)
  expect_true(all(ex[[4]]$Events >= ex[[1]]$Events + ex[[2]]$Events + ex[[3]]$Events))
})

test_that("rates2csv and csv2rates round-trip", {
  R <- matisse_example(1); f <- tempfile(fileext = ".csv")
  rates2csv(R, f, title = "melanoma F A")
  R2 <- csv2rates(f)
  expect_equal(unname(R2$Events), unname(R$Events)); expect_equal(unname(R2$Offset), unname(R$Offset))
  expect_equal(R2$ages, R$ages); expect_equal(R2$periods, R$periods); expect_equal(R2$fullname, "melanoma F A")
})

test_that("chunk keeps totals and the most recent periods (LTF)", {
  R <- matisse_example(1); R5 <- rates_chunk(R, 5, 5)
  expect_equal(dim(R5$Events), c(10, 5))
  expect_equal(tail(R5$periods, 1), tail(R$periods, 1))
  keep <- R$periods[-length(R$periods)] >= R5$periods[1]
  expect_equal(sum(R5$Events), sum(R$Events[, keep]))
})

test_that("fill replaces zeros by a value >= 0.5", {
  R <- rates_fill(matisse_example(8))
  expect_true(all(R$Events > 0)); expect_gte(R$fill_value, 0.5)
  expect_equal(R$n_filled, matisse_example()$zero_cells[8])
})

test_that("pchip interpolates the knots and preserves monotonicity", {
  x <- c(20, 25, 30, 35, 40, 45); y <- c(0, 3, 10, 30, 31, 80); xi <- seq(20, 45, 0.5)
  expect_equal(pchip(x, y, x), y)
  expect_true(all(diff(pchip(x, y, xi)) >= -1e-12))
})

test_that("i5 keeps period totals (both methods) and group totals (linear)", {
  R5 <- rates_chunk(matisse_example(1), age_block = 5)             # 5-year ages, single years
  for (m in c("linear", "cubic")) {
    R1 <- rates_i5(R5, m)
    expect_equal(dim(R1$Events), c(50, ncol(R5$Events)))
    expect_equal(colSums(R1$Events), colSums(R5$Events)); expect_equal(colSums(R1$Offset), colSums(R5$Offset))
  }
  R1 <- rates_i5(R5, "linear"); g <- rep(seq_len(nrow(R5$Events)), each = 5)
  expect_equal(unname(rowsum(R1$Events, g)), unname(R5$Events))
})
