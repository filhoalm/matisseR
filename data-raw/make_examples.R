## Synthetic example Lexis diagrams for matisseR. They replace the SEER example datasets of MATISSE, which may not be
## redistributed under the SEER Research Data Use Agreement. Counts are Poisson draws from known age-period-cohort
## surfaces (log rate = function of age + function of period + function of cohort), single years of age 35-84 by single
## calendar years 1991-2019, with synthetic person-years. The shapes are set by hand to resemble common cancer
## patterns; they are not fitted to any registry data.
## Usage (package root): Rscript data-raw/make_examples.R   -> inst/extdata/syn*.csv and inst/extdata/examples.csv

library(data.table)
pkgload::load_all(quiet = TRUE)                                   # csv2rates (check) and apc_fit (true net drift)
set.seed(20261024)
ages <- 35:84; years <- 1991:2019; A <- length(ages); P <- length(years)
a <- rep(ages + 0.5, P); p <- rep(years + 0.5, each = A); co <- p - a   # cell midpoints, age fastest (as c(Events))

## person-years: survival by age, a baby boom around the 1962 cohort, slow growth over time
py <- function(size, sex) round(size * exp(-exp((a - c(F = 88, M = 83)[[sex]]) / 9)) *
                                  (1 + 0.3 * exp(-((co - 1962) / 9)^2)) * (1 + 0.004 * (p - 2005)))
## log rate per person-year at age 60 in 2005 (r60 per 100,000), with a linear and quadratic age curve
agec <- function(r60, b1, b2) log(r60 * 1e-5) + b1 * (a - 60) + b2 * (a - 60)^2
mel_coh <- function(c0) 1 / (1 + exp(-(co - c0) / 9)) - 0.5      # melanoma-like rise over cohorts, centred on c0
mel <- agec(30, 0.035, -3e-4) + mel_coh(1945) + 0.04 * sin((p - 1991) / 3)

S <- list(
  list(site = "melanoma", sex = "F", pop = "A", size = 40000, lr = mel,
       design = "main example: rising cohort pattern"),
  list(site = "melanoma", sex = "F", pop = "B", size = 30000, lr = mel + log(0.8) + 0.012 * (p - 2005),
       design = "A with a different period trend (vs A: PH-X holds)"),
  list(site = "melanoma", sex = "F", pop = "C", size = 6000, lr = mel + log(1.3),
       design = "A x 1.3 in a small population, zero cells (vs A: PH-A holds)"),
  list(site = "melanoma", sex = "F", pop = "D", size = 25000, lr = mel + log(0.9) + 0.015 * (co - 1945),
       design = "A with a steeper cohort trend (vs A: PH-L holds)"),
  list(site = "breast", sex = "F", pop = "A", size = 40000,
       lr = agec(250, 0.03, -1.2e-3) + 0.006 * (co - 1945) + 0.03 * cos((p - 1991) / 4),
       design = "population A: steady cohort increase"),
  list(site = "lung", sex = "F", pop = "A", size = 40000,
       lr = agec(60, 0.07, -8e-4) + 0.9 * exp(-((co - 1948) / 14)^2) - 0.5,
       design = "population A: cohort peak around 1948"),
  list(site = "all cancers", sex = "F", pop = "A", size = 40000, lr = NULL,
       design = "population A: examples 1, 5 and 6 plus other cancers (same person-years)"),
  list(site = "thyroid", sex = "M", pop = "C", size = 6000, lr = agec(10, 0.01, -5e-4) + 0.035 * (p - 2005),
       design = "small population, mostly zero cells, rising period trend"))
other <- agec(450, 0.05, -4e-4) + 0.003 * (p - 2005)               # other cancers of population A (in 'all cancers')

ext <- file.path("inst", "extdata")
unlink(Sys.glob(file.path(ext, "*.csv")))
write_rates <- function(E, O, file, title) {                       # NCI APC Web Tool format, as MATISSE's files
  M <- matrix(0L, A, 2 * P); M[, c(TRUE, FALSE)] <- as.integer(E); M[, c(FALSE, TRUE)] <- as.integer(O)
  pad <- strrep(",", 2 * P - 1)
  con <- file(file, "wb"); on.exit(close(con))
  writeLines(c(paste0(c(paste("Title:", title), "Description: synthetic data (matisseR data-raw/make_examples.R)",
                        paste("Start Year:", years[1]), paste("Start Age:", ages[1]), "Interval (Years): 1"), pad),
               apply(M, 1, paste, collapse = ",")), con)
}
true_drift <- function(mu, O) {                                    # noise-free surface: the APC fit is exact
  R <- list(Events = matrix(mu, A), Offset = matrix(O, A), ages = c(ages, max(ages) + 1), periods = c(years, max(years) + 1))
  apc_fit(R)$net_drift
}

E <- list(); idx <- list()
for (k in seq_along(S)) {
  s <- S[[k]]; O <- py(s$size, s$sex)
  if (is.null(s$lr)) {                                             # all cancers = sum of the site counts + other
    mu <- NA; E[[k]] <- E[[1]] + E[[5]] + E[[6]] + rpois(A * P, O * exp(other))
  } else {
    mu <- O * exp(s$lr); E[[k]] <- rpois(A * P, mu)
  }
  file <- sprintf("syn%02d_%s_%s_%s.csv", k, gsub(" ", "_", s$site), s$sex, s$pop)
  title <- sprintf("Synthetic %s %s %s", s$site, s$sex, s$pop)
  write_rates(E[[k]], O, file.path(ext, file), title)
  R <- csv2rates(file.path(ext, file))
  stopifnot(identical(as.numeric(R$Events), as.numeric(E[[k]])), identical(as.numeric(R$Offset), as.numeric(O)),
            R$fullname == title)
  idx[[k]] <- data.table(number = k, file, site = s$site, sex = s$sex, population = s$pop, ages = A, periods = P,
                         cases = sum(E[[k]]), zero_cells = sum(E[[k]] == 0),
                         net_drift = if (anyNA(mu)) NA_real_ else round(true_drift(mu, O), 6), design = s$design)
}
idx <- rbindlist(idx)
fwrite(idx, file.path(ext, "examples.csv"))
print(idx[, !"design"])
