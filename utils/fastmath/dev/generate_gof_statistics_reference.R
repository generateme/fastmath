# Reference values for fastmath.stats/ks-statistic, cvm-statistic, ad-statistic.
# Tools: R (stats::ks.test, goftest::cvm.test, goftest::ad.test, fitdistrplus::mgedist).
# KS, CvM and AD come from stats/goftest. The AD variants (ADR, ADL, AD2R, AD2L, AD2) come from the
# objective of fitdistrplus::mgedist, captured at the given parameters with a custom optimizer;
# fitdistrplus divides every objective by n, so the reference is objective * n (the scale of A2).
# Datasets are stored in the output, so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_gof_statistics_reference.R
# Output: test/resources/stats/gof_statistics_reference.edn
suppressMessages({library(goftest); library(fitdistrplus)})

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

captured <- NULL
capture <- function(fn, par, fix.arg, obs, pdistnam, hessian, method, ...) {
  captured <<- fn(par, fix.arg, obs, pdistnam)
  list(par = par, value = captured, hessian = NULL, counts = c(NA, NA), convergence = 0, message = NULL)
}
objective <- function(x, distr, start, gof) {
  mgedist(x, distr, gof = gof, start = start, custom.optim = capture)
  captured * length(x)
}

# fm: fastmath key and parameter map (as EDN text); distr: R distribution name; args: R parameters
case <- function(x, fm, distr, args) {
  p <- paste0("p", distr)
  ks <- suppressWarnings(do.call(ks.test, c(list(x, p), args)))$statistic
  cvm <- do.call(cvm.test, c(list(x, null = p), args))$statistic
  ad <- do.call(ad.test, c(list(x, null = p), args))$statistic
  v <- sapply(c("ADR", "ADL", "AD2R", "AD2L", "AD2"), function(g) objective(x, distr, args, g))
  # fitdistrplus' own CvM/AD objectives, same scale check as a guard on the capture
  stopifnot(abs(objective(x, distr, args, "AD") - ad) < 1e-9 * max(1, abs(ad)))
  stopifnot(abs(objective(x, distr, args, "CvM") - cvm) < 1e-9 * max(1, abs(cvm)))
  paste0("{:distribution ", fm$key, " :params ", fm$params, "\n :data ", vec(x),
         "\n :ks ", fmt(ks), " :cvm ", fmt(cvm), " :ad ", fmt(ad),
         "\n :adr ", fmt(v["ADR"]), " :adl ", fmt(v["ADL"]), " :ad2r ", fmt(v["AD2R"]),
         " :ad2l ", fmt(v["AD2L"]), " :ad2 ", fmt(v["AD2"]), "}")
}

set.seed(20261029)
out <- c(
  case(rnorm(25), list(key = ":normal", params = "{:mu 0.1 :sd 1.2}"), "norm", list(mean = 0.1, sd = 1.2)),
  case(rexp(30), list(key = ":exponential", params = paste0("{:mean ", fmt(1 / 1.3), "}")), "exp", list(rate = 1.3)),
  case(rgamma(40, 2, 1), list(key = ":gamma", params = paste0("{:shape 1.8 :scale ", fmt(1 / 1.1), "}")),
       "gamma", list(shape = 1.8, rate = 1.1)),
  # ties: data rounded to one decimal place
  case(round(rexp(30), 1) + 0.1, list(key = ":exponential", params = "{:mean 1.0}"), "exp", list(rate = 1)))

dir.create("test/resources/stats", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(out, collapse = "\n"), "]"), "test/resources/stats/gof_statistics_reference.edn")
