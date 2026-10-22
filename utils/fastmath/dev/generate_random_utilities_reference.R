## R part of the reference data for `fastmath.random.utilities-test` (Random Audit WP 7.3, 9.6).
##
## Tool: R with gamlss.dist (`dNBI`, `pNBI`, `qNBI`, `dNBII`, `pNBII`, `qNBII`), jsonlite for the output.
## The numbers are written as strings with 17 significant digits ("%.17g"), so that Python reads back exactly the doubles of R.
## Precision: double precision (gamlss uses `lgamma` differences, the relative error of the pmf is ~1e-11 for mu = 1e4).
## Domain: (mu, sigma) in (1,1), (100,.01), (.01,10), (1e4,1e-3), (5,2e-4), (3,.5) for both families; sigma >= 1e-4 only
## (below that gamlss itself switches to the Poisson law, which is not the exact negative binomial). The quantiles are
## computed at fixed levels (not .5: for mu = 100, sigma = .01 the exact cdf(99) = .5 + 2e-19, a tie in double precision). No random numbers are used: the output is deterministic.
##
## It is called by `generate_random_utilities_reference.py` (which merges the output into the EDN file):
##     Rscript generate_random_utilities_reference.R /tmp/utilities_gamlss.json
suppressMessages({library(jsonlite); library(gamlss.dist)})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "/tmp/utilities_gamlss.json"

fmt <- function(v) vapply(v, function(z) sprintf("%.17g", z), "")

PS <- c(1e-10, 1e-3, 0.1, 0.3, 0.9, 0.999, 1 - 1e-10)
LEVELS <- c(1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)

mk <- function(fam, mu, sigma) {
  d <- get(paste0("d", fam)); p <- get(paste0("p", fam)); q <- get(paste0("q", fam))
  qs <- sapply(LEVELS, function(l) q(l, mu, sigma))
  ks <- unique(c(qs - 1, qs, qs + 1, 0:4))
  ks <- sort(ks[ks >= 0])
  pmf <- sapply(ks, function(k) d(k, mu, sigma))
  lpmf <- sapply(ks, function(k) d(k, mu, sigma, log = TRUE))
  cdf <- sapply(ks, function(k) p(k, mu, sigma))
  sf <- sapply(ks, function(k) p(k, mu, sigma, lower.tail = FALSE))
  icdf <- sapply(PS, function(pp) q(pp, mu, sigma))
  var <- if (fam == "NBI") mu + sigma * mu^2 else mu * (1 + sigma)
  list(mu = fmt(mu), sigma = fmt(sigma), ks = fmt(ks), pmf = fmt(pmf), lpmf = fmt(lpmf), cdf = fmt(cdf), sf = fmt(sf),
       ps = fmt(PS), icdf = fmt(icdf), mean = fmt(mu), var = fmt(var))
}

sets <- list(c(1, 1), c(100, 0.01), c(0.01, 10), c(1e4, 1e-3), c(5, 2e-4), c(3, 0.5))
res <- list(
  nbi = lapply(sets, function(s) mk("NBI", s[1], s[2])),
  nbii = lapply(sets, function(s) mk("NBII", s[1], s[2])),
  version = list(R = R.version.string, gamlss.dist = as.character(packageVersion("gamlss.dist")))
)
write_json(res, out, auto_unbox = TRUE, digits = NA)
