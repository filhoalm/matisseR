## Compare matisseR (validation/r) with MATISSE (validation/matisse) outputs.
## Usage (package root): Rscript validation/compare.R

library(data.table)
keys <- list(apc_ef = c("example", "comp", "k"), apc_key = c("example", "k"),
             cap_global = c("comparison", "k"), cap_homogeneity = c("comparison", "k"), cap_composite = c("comparison", "k"),
             cap_aic = c("comparison", "k"), cap_ef = c("comparison", "comp", "g", "k"))
vals <- list(apc_ef = "value", apc_key = "value", cap_global = c("stat", "df", "p"), cap_homogeneity = c("stat", "df", "p"),
             cap_composite = "p", cap_aic = c("df", "bcdf", "DEV", "AICc", "delta", "rank", "model"), cap_ef = "value")
rd <- function(dir, t) { f <- file.path("validation", dir, paste0(t, ".csv"))
  if (!file.exists(f)) return(NULL); setnames(fread(f, header = FALSE), c(keys[[t]], vals[[t]])) }
res <- rbindlist(lapply(names(keys), function(t) {
  r <- rd("r", t); m <- rd("matisse", t)
  if (is.null(r) || is.null(m)) return(data.table(table = t, quantity = NA_character_, n = NA, max_abs = NA, max_rel = NA))
  x <- merge(r, m, by = keys[[t]], suffixes = c("_r", "_m"))
  rbindlist(lapply(vals[[t]], function(v) { a <- x[[paste0(v, "_r")]]; b <- x[[paste0(v, "_m")]]
    data.table(table = t, quantity = v, n = sprintf("%d/%d", nrow(x), max(nrow(r), nrow(m))),
               max_abs = max(abs(a - b), na.rm = TRUE), max_rel = max(abs(a - b) / pmax(abs(b), 1e-12), na.rm = TRUE)) }))
}))
res[, status := fifelse(is.na(max_rel), "missing", fifelse(max_rel < 1e-6 | max_abs < 1e-8, "agree", "DIFFER"))]
print(res, digits = 3)
