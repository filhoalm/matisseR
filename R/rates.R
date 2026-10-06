## MATISSE `rates` objects (@rates class): a Lexis diagram as Events and Offset matrices
## (A age groups x P periods) with age and period cutpoints. Methods mirror MATISSE:
## csv2rates.m, CHUNK.m, FILL.m, I5.m, UNI.m, RATES2CSV.m and LEXIS.m.

#' Rates object from long data
#'
#' Builds a MATISSE-style rates object from one row per age group and year.
#'
#' @param dt data.table with `age` (5-year group index 1-18: 1 = 0-4, ..., 18 = 85+), `year`, `cases`, `py`.
#' @param age_groups age group indices to keep (consecutive).
#' @param fullname,event_label labels.
#' @return list with `Events`, `Offset` (age x year matrices), `ages` and `periods` (cutpoints, length A + 1 and P + 1),
#'   `fullname`, `event_label`.
#' @export
rates_make <- function(dt, age_groups, fullname = "", event_label = "Cases") {
  dt <- dt[age %in% age_groups, .(cases = sum(cases), py = sum(py)), by = .(age, year)]
  yrs <- sort(unique(dt$year))
  stopifnot(all(diff(yrs) == 1))
  dt <- merge(CJ(age = age_groups, year = yrs), dt, by = c("age", "year"), all.x = TRUE)
  E  <- matrix(dt$cases, length(age_groups), length(yrs), byrow = TRUE, dimnames = list(age_groups, yrs))
  O  <- matrix(dt$py, length(age_groups), length(yrs), byrow = TRUE, dimnames = list(age_groups, yrs))
  list(Events = E, Offset = O, ages = c(5 * (age_groups - 1), 5 * max(age_groups)), periods = c(yrs, max(yrs) + 1),
       fullname = fullname, event_label = event_label)
}

#' Read a rates object from a CSV file (MATISSE csv2rates.m)
#'
#' Reads the NCI Age Period Cohort Web Tool format used by MATISSE: header lines (`Title:`,
#' `Description:`, `Start Year:`, `Start Age:`, `Interval (Years):` or `Age Bin (Years):` and
#' `Period Bin (Years):`), then one row per age group with paired events and offset columns per period.
#'
#' @param file path to the CSV file.
#' @return a rates object (see [rates_make()]).
#' @export
csv2rates <- function(file) {
  L <- readLines(file, warn = FALSE)
  hd <- c(title = "Title:", description = "Description:", year = "Start Year:", age = "Start Age:",
          interval = "Interval (Years):", age_bin = "Age Bin (Years):", per_bin = "Period Bin (Years):")
  h <- lapply(hd, function(k) { x <- grep(k, head(L, 8), fixed = TRUE, value = TRUE)[1]
    if (is.na(x)) NULL else trimws(gsub('"', "", sub(",+\\s*$", "", substring(x, regexpr(k, x, fixed = TRUE) + nchar(k))))) })
  n_hd <- sum(vapply(hd, function(k) any(grepl(k, head(L, 8), fixed = TRUE)), TRUE))
  N <- as.matrix(fread(text = L[-seq_len(n_hd)], header = FALSE))
  N <- N[, colSums(!is.na(N)) > 0, drop = FALSE]
  E <- N[, c(TRUE, FALSE), drop = FALSE]; O <- N[, c(FALSE, TRUE), drop = FALSE]
  da <- as.numeric(if (is.null(h$age_bin)) h$interval else h$age_bin)
  dp <- as.numeric(if (is.null(h$per_bin)) h$interval else h$per_bin)
  ages <- floor(as.numeric(h$age)) + da * (0:nrow(E)); periods <- floor(as.numeric(h$year)) + dp * (0:ncol(E))
  dimnames(E) <- dimnames(O) <- list(head(ages, -1), head(periods, -1))
  list(Events = E, Offset = O, ages = ages, periods = periods, fullname = if (is.null(h$title)) "Rates" else h$title,
       event_label = "Cases", description = h$description)
}

#' MATISSE example datasets
#'
#' The 16 example Lexis diagrams shipped with the MATISSE toolbox (`ratesdata.m`): SEER incidence,
#' single years of age 35-84 by single calendar years from 1992 (27 or 29 years).
#'
#' @param i example number(s), 1-16; `NULL` returns the index table.
#' @return a rates object, a list of rates objects, or the index (data.table).
#' @examples
#' matisse_example()
#' R <- matisse_example(16)   # melanoma, female, non-Hispanic white
#' @export
matisse_example <- function(i = NULL) {
  idx <- fread(system.file("extdata", "examples.csv", package = "matisseR"))
  if (is.null(i)) return(idx)
  R <- lapply(i, function(k) csv2rates(system.file("extdata", idx[number == k, file], package = "matisseR")))
  if (length(i) == 1) R[[1]] else R
}

#' Aggregate a rates object into larger cells (MATISSE CHUNK.m)
#'
#' Sums successive periods (and/or ages) into blocks; leftover periods are dropped from the start
#' ('LTF', MATISSE default) so the most recent years are kept.
#'
#' @param R rates object.
#' @param age_block,per_block number of age groups / periods per block.
#' @return rates object.
#' @export
rates_chunk <- function(R, age_block = 1, per_block = 1) {
  blocks <- function(M, cut, k, margin) {
    n <- if (margin == 1) nrow(M) else ncol(M)
    keep <- (n %% k + 1):n
    g <- rep(seq_len(n %/% k), each = k)
    S <- if (margin == 1) rowsum(M[keep, , drop = FALSE], g, na.rm = TRUE)
         else t(rowsum(t(M[, keep, drop = FALSE]), g, na.rm = TRUE))
    list(M = S, cut = cut[keep[1]] + (0:(n %/% k)) * k * diff(cut)[1])
  }
  for (x in list(list(k = age_block, m = 1, cut = "ages"), list(k = per_block, m = 2, cut = "periods"))) {
    if (x$k == 1) next
    e <- blocks(R$Events, R[[x$cut]], x$k, x$m)
    o <- blocks(R$Offset, R[[x$cut]], x$k, x$m)
    R$Events <- e$M; R$Offset <- o$M; R[[x$cut]] <- e$cut
  }
  dimnames(R$Events) <- dimnames(R$Offset) <- list(head(R$ages, -1), head(R$periods, -1))
  R
}

#' Fill zero or missing events (MATISSE FILL.m, zero_fill = 'adp')
#'
#' Replaces zero/missing event counts by MATISSE's adaptive tiny value (>= 0.5), needed before
#' fitting models on log rates.
#'
#' @param R rates object.
#' @return rates object with `n_filled` and `fill_value`.
#' @export
rates_fill <- function(R) {
  e <- R$Events
  U <- sort(unique(e[!is.na(e) & e > 0]))
  if (!length(U)) U <- 0.1
  C <- min(length(U), 5)
  delta <- diff(R$ages)[1]
  zfv <- if (C < 2) 0.5 else {
    L <- 12 / (delta * (C - 1) * C * (C + 1)) * (seq_len(C) - (C + 1) / 2)
    max(0.5 * (U[1] - sum(L * U[seq_len(C)])) / delta, 0.5)
  }
  stopifnot(!anyNA(R$Offset))
  R$n_filled <- sum(is.na(e) | e == 0)
  R$fill_value <- zfv
  R$Events[is.na(e) | e == 0] <- zfv
  R
}

#' Long format: one row per Lexis cell
#'
#' @param R rates object.
#' @param per rate multiplier (per 100,000 by default).
#' @return data.table with age, period and cohort (midpoint c = p - a) of each cell, events, offset and rate.
#' @export
rates_long <- function(R, per = 1e5) {
  da <- diff(R$ages)[1]; dp <- diff(R$periods)[1]
  dt <- CJ(i = seq_len(nrow(R$Events)), j = seq_len(ncol(R$Events)))
  dt[, `:=`(age_lo = R$ages[i], age_mid = R$ages[i] + da / 2, per_lo = R$periods[j], per_mid = R$periods[j] + dp / 2,
            events = R$Events[cbind(i, j)], offset = R$Offset[cbind(i, j)])]
  dt[, `:=`(coh_mid = per_mid - age_mid, rate = per * events / offset)][, c("i", "j") := NULL][]
}

#' MATLAB pchip interpolation
#'
#' Shape-preserving piecewise cubic Hermite interpolation with MATLAB's slopes (pchipslopes).
#'
#' @param x,y knots (x increasing).
#' @param xi points to interpolate.
#' @return interpolated values.
#' @export
pchip <- function(x, y, xi) {
  n <- length(x); h <- diff(x); del <- diff(y) / h
  d <- numeric(n)
  k <- which(sign(del[-(n - 1)]) * sign(del[-1]) > 0)
  if (length(k)) {
    hs <- h[k] + h[k + 1]
    w1 <- (h[k] + hs) / (3 * hs); w2 <- (hs + h[k + 1]) / (3 * hs)
    dmax <- pmax(abs(del[k]), abs(del[k + 1])); dmin <- pmin(abs(del[k]), abs(del[k + 1]))
    d[k + 1] <- dmin / (w1 * (del[k] / dmax) + w2 * (del[k + 1] / dmax))
  }
  end_slope <- function(h1, h2, d1, d2) {
    s <- ((2 * h1 + h2) * d1 - h1 * d2) / (h1 + h2)
    if (sign(s) != sign(d1)) 0 else if (sign(d1) != sign(d2) && abs(s) > abs(3 * d1)) 3 * d1 else s
  }
  d[1] <- end_slope(h[1], h[2], del[1], del[2])
  d[n] <- end_slope(h[n - 1], h[n - 2], del[n - 1], del[n - 2])
  j <- pmin(findInterval(xi, x, all.inside = TRUE), n - 1)
  s <- xi - x[j]
  c3 <- (3 * del[j] - 2 * d[j] - d[j + 1]) / h[j]
  b3 <- (d[j] - 2 * del[j] + d[j + 1]) / h[j]^2
  y[j] + s * (d[j] + s * (c3 + s * b3))
}

#' Single-year ages from grouped ages (MATISSE I5.m)
#'
#' Interpolates cumulative events and offsets within each period. 'linear' keeps the group totals
#' (flat rates within a group); 'cubic' (pchip) gives smooth single-year rates and keeps the period totals.
#'
#' @param R rates object with grouped ages.
#' @param method "linear" or "cubic".
#' @return rates object with single-year ages.
#' @export
rates_i5 <- function(R, method = c("linear", "cubic")) {
  method <- match.arg(method)
  a <- R$ages; sya <- seq(a[1], tail(a, 1))
  one <- function(y, normalise) {
    cy <- c(0, cumsum(y))
    v  <- if (method == "linear") pmax(0, approx(a, cy, sya)$y) else pchip(a, cy, sya)
    d  <- pmax(0, diff(v))
    if (normalise) sum(y) * d / sum(d) else d
  }
  R$Events <- apply(R$Events, 2, one, normalise = method == "cubic")
  R$Offset <- apply(R$Offset, 2, one, normalise = TRUE)
  R$ages   <- sya
  dimnames(R$Events) <- dimnames(R$Offset) <- list(head(sya, -1), head(R$periods, -1))
  R
}

#' Summary statistics of a rates object (cf. MATISSE uni)
#'
#' @param R rates object.
#' @return one-row data.table.
#' @export
rates_uni <- function(R) {
  A <- nrow(R$Events); P <- ncol(R$Events)
  data.table(object = R$fullname, ages = sprintf("%g-%g", R$ages[1], tail(R$ages, 1) - 1),
             periods = sprintf("%g-%g", R$periods[1], tail(R$periods, 1) - 1), cell = sprintf("%gx%g", diff(R$ages)[1], diff(R$periods)[1]),
             A = A, P = P, C = A + P - 1, events = round(sum(R$Events)), person_years = round(sum(R$Offset)),
             crude_rate = round(1e5 * sum(R$Events) / sum(R$Offset), 1), zero_cells_filled = if (is.null(R$n_filled)) 0L else R$n_filled)
}

#' Write a rates object as CSV (inverse of MATISSE csv2rates.m)
#'
#' NCI APC Web Tool format: header lines, then paired events/offset columns per period, one row per age group.
#'
#' @param R rates object with equal age and period intervals.
#' @param file output path.
#' @param title,description header text.
#' @export
rates2csv <- function(R, file, title = R$fullname, description = "") {
  d <- diff(R$ages)[1]
  stopifnot(d == diff(R$periods)[1])
  M <- matrix(0, nrow(R$Events), 2 * ncol(R$Events))
  M[, c(TRUE, FALSE)] <- R$Events
  M[, c(FALSE, TRUE)] <- R$Offset
  writeLines(c(paste0("Title: ", title), paste0("\"Description: ", description, "\""),
               paste0("Start Year: ", R$periods[1]), paste0("Start Age: ", R$ages[1]), paste0("Interval (Years): ", d),
               apply(M, 1, function(x) paste(sprintf("%f", x), collapse = ", "))), file)
}

#' Lexis diagram plot (MATISSE LEXIS.m)
#'
#' Rates over age x period beside the same cells over age x birth cohort; the period panel is padded
#' to the cohort span so both panels share one aspect ratio. Requires ggplot2, patchwork and scales.
#'
#' @param R rates object with square cells.
#' @param cols low and high colours of the (log) scale.
#' @param legend legend title.
#' @return a patchwork object.
#' @export
lexis_plot <- function(R, cols = c("#e3eefc", "#0d366b"), legend = "Rate per 100,000") {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE), requireNamespace("patchwork", quietly = TRUE))
  d <- diff(R$ages)[1]
  stopifnot(d == diff(R$periods)[1])
  dt <- rates_long(R)
  yc <- range(dt$coh_mid) + c(-d, d) / 2
  yp <- mean(range(dt$per_mid)) + c(-1, 1) * diff(yc) / 2
  per_lab <- if (d == 1) head(R$periods, -1) else sprintf("%d-%02d", head(R$periods, -1), (tail(R$periods, -1) - 1) %% 100)
  base <- list(
    ggplot2::geom_tile(width = d, height = d, colour = "white", linewidth = 0.3),
    ggplot2::scale_fill_gradient(low = cols[1], high = cols[2], trans = "log10", labels = scales::label_comma(), name = legend),
    ggplot2::scale_x_continuous(breaks = seq(R$ages[1], tail(R$ages, 1), 10), expand = c(0, 0)),
    ggplot2::coord_fixed(), ggplot2::theme_minimal(base_size = 10),
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom", legend.key.width = ggplot2::unit(1.6, "cm")))
  p1 <- ggplot2::ggplot(dt, ggplot2::aes(age_mid, per_mid, fill = rate)) + base +
    ggplot2::scale_y_continuous(limits = yp, breaks = unique(dt$per_mid), labels = per_lab, expand = c(0, 0)) +
    ggplot2::labs(x = "Cross-sectional age", y = "Cross-sectional period")
  p2 <- ggplot2::ggplot(dt, ggplot2::aes(age_mid, coh_mid, fill = rate)) + base +
    ggplot2::scale_y_continuous(limits = yc, breaks = seq(ceiling(yc[1] / 10) * 10, yc[2], 10), expand = c(0, 0)) +
    ggplot2::labs(x = "Longitudinal age", y = "Birth cohort (midpoint)")
  patchwork::wrap_plots(p1, p2, nrow = 1) + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}
