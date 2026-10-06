## matisseR outputs for the MATISSE examples, in the layout written by validate_matisse.m (MATLAB),
## on MATISSE's scale (log rates per 100,000; local drifts and contrasts as log differences).
## Usage (package root): Rscript validation/validate_r.R   -> validation/r/*.csv ; then Rscript validation/compare.R

library(data.table); library(matisseR)
out <- file.path("validation", "r"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
comps <- c("lac", "cac", "ftt", "fcp", "ld")
mlog <- function(ef, v) suppressWarnings(fifelse(rep_len(ef == "ld", length(v)), log1p(v / 100), log(v)))   # MATISSE: lot + X * B (lot = log(1e5); 0 for ld)

## APC (WLS, overdispersion) on the 16 examples at 1 x 1, zero cells filled: estimable functions and key parameters
efs <- list(); key <- list()
for (i in 1:16) {
  M <- apc_fit(rates_fill(matisse_example(i)))
  e <- apc_efs(M)[, k := seq_len(.N), by = ef]
  efs[[i]] <- e[, .(example = i, comp = match(ef, comps), k, value = mlog(ef, value))]
  key[[i]] <- data.table(example = i, k = 1:7, value = c(M$B[1:3], M$B[2] - M$B[3], M$B[4:6]))   # Intercept, LAT, NetDrift, CAT, THETAa/p/c
}
fwrite(rbindlist(efs), file.path(out, "apc_ef.csv"), col.names = FALSE)
fwrite(rbindlist(key), file.path(out, "apc_key.csv"), col.names = FALSE)

## CAPRICORN: thyroid NHB vs NHW and ER+ breast NHB vs NHW (1 x 1), thyroid API/HIS/NHB/NHW (5 x 5)
cmp <- list(list(ex = 9:10, d = 1), list(ex = 1:2, d = 1), list(ex = 7:10, d = 5))
hom_order <- c("CAT", "HiOrdCoh", "QuadCoh", "NetDrift", "HiOrdPer", "QuadPer", "LAT", "HiOrdAge", "QuadAge")   # MATISSE row order
tab <- list(global = list(), homogeneity = list(), composite = list(), aic = list(), ef = list())
for (j in seq_along(cmp)) {
  Rs <- lapply(matisse_example(cmp[[j]]$ex), function(R) rates_fill(if (cmp[[j]]$d > 1) rates_chunk(R, cmp[[j]]$d, cmp[[j]]$d) else R))
  names(Rs) <- paste0("ex", cmp[[j]]$ex)
  S <- capricorn(Rs); NK <- nrow(S$G$X) * S$K
  tab$global[[j]] <- S$tests$global[, .(comparison = j, k = .I, stat, df, p)]
  tab$homogeneity[[j]] <- S$tests$homogeneity[match(hom_order, test), .(comparison = j, k = .I, stat, df, p)]
  tab$composite[[j]] <- S$tests$composite[, .(comparison = j, k = .I, p)]
  tab$aic[[j]] <- S$aic[, .(comparison = j, k = .I, df, bcdf = df + df * (df + 1) / (NK - df - 1), DEV, AICc, delta, rank, model = .I)]
  tab$ef[[j]] <- rbindlist(lapply(comps, function(cp) cap_ef(S, cp)[type == "ratio"][
    , .(comparison = j, comp = match(cp, comps), g = match(sub(" vs .*", "", stratum), names(Rs)), value = mlog(cp, value)), by = stratum][
    , k := seq_len(.N), by = g][, .(comparison, comp, g, k, value)]))
}
for (t in names(tab)) fwrite(rbindlist(tab[[t]]), file.path(out, paste0("cap_", t, ".csv")), col.names = FALSE)
cat("wrote", length(list.files(out)), "files to", out, "\n")
