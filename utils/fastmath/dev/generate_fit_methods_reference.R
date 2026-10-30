# Reference values for fastmath.fitdistr/fit with the methods other than :mle.
# Tools: R fitdistrplus (mgedist: CvM KS AD ADR ADL AD2R AD2L AD2; msedist; qmedist; mmedist), lmomco::lmoms (sample L-moments).
# :lmme has no ready R function for the two-parameter parametrisation used here, so its reference is an independent
# solution of the L-moment equations: the theoretical L-moments are integrate(q(u) * P*_{r-1}(u), 0, 1) (shifted Legendre
# polynomials) and the system is solved with nleqslv.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_methods_reference.R
# Output: test/resources/fitdistr/fit_methods_reference.edn
suppressMessages({library(fitdistrplus); library(lmomco); library(nleqslv)})

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

# qmedist matches as many quantiles as there are parameters (two here)
probs <- c(0.25, 0.75)

mweibull <- function(order, shape, scale) scale^order * gamma(1 + order / shape)

# shifted Legendre polynomials P*_0 .. P*_2
legendre <- list(function(u) 1, function(u) 2 * u - 1, function(u) 6 * u^2 - 6 * u + 1)

lmme <- function(x, qdist, start) {
  l <- lmoms(x)$lambdas[1:2]
  theo <- function(par) sapply(1:2, function(r)
    integrate(function(u) do.call(qdist, c(list(u), as.list(par))) * legendre[[r]](u), 0, 1,
              rel.tol = 1e-12, subdivisions = 1000)$value)
  nleqslv(start, function(par) theo(par) - l, control = list(xtol = 1e-12, ftol = 1e-12))$x
}

case <- function(key, x, distr, qdist, mapper, mme.args = list(), start = NULL) {
  fit <- function(f, ...) f(x, distr, ...)
  est <- list()
  for (g in c("CvM", "KS", "AD", "ADR", "ADL", "AD2R", "AD2L", "AD2"))
    est[[tolower(g)]] <- mgedist(x, distr, gof = g, start = start)$estimate
  est[["mps"]] <- msedist(x, distr, start = start)$estimate
  est[["qme"]] <- qmedist(x, distr, probs = probs, start = start)$estimate
  est[["mme"]] <- do.call(mmedist, c(list(x, distr), mme.args))$estimate
  mle <- fitdist(x, distr)$estimate
  est[["lmme"]] <- setNames(lmme(x, qdist, mle), names(mle))
  entries <- sapply(names(est), function(m) paste0(" :", m, " ", mapper(est[[m]])))
  paste0("{:distribution ", key, " :probs ", vec(probs), "\n :data ", vec(x), "\n :estimates {\n",
         paste(entries, collapse = "\n"), "}}")
}

set.seed(20261030)
out <- c(
  case(":normal", rnorm(150, 1.5, 2), "norm", "qnorm",
       function(e) sprintf("{:mu %.17g :sd %.17g}", e["mean"], e["sd"])),
  case(":gamma", rgamma(200, 2.5, 1 / 1.7), "gamma", "qgamma",
       function(e) sprintf("{:shape %.17g :scale %.17g}", e["shape"], 1 / e["rate"])),
  case(":weibull", rweibull(200, 1.8, 3), "weibull", "qweibull",
       function(e) sprintf("{:alpha %.17g :beta %.17g}", e["shape"], e["scale"]),
       mme.args = list(order = 1:2, memp = function(x, order) mean(x^order), start = list(shape = 1.5, scale = 2))))

dir.create("test/resources/fitdistr", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(out, collapse = "\n"), "]"), "test/resources/fitdistr/fit_methods_reference.edn")
