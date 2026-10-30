# Reference values for the standard errors of fastmath.fitdistr/fit (maximum likelihood).
# Tools: R fitdistrplus::fitdist, standard errors from sqrt(diag(vcov)) (inverse of the numerical Hessian of optim).
# The densities are written in the fastmath parameterization, so the estimates and the errors are the ones of fastmath.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_se_reference.R
# Output: test/resources/fitdistr/fit_se_reference.edn
suppressMessages(library(fitdistrplus))

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

dgamma_scale <- function(x, shape, scale, log = FALSE) dgamma(x, shape = shape, scale = scale, log = log)
pgamma_scale <- function(q, shape, scale, ...) pgamma(q, shape = shape, scale = scale, ...)
dexp_mean <- function(x, mean, log = FALSE) dexp(x, rate = 1 / mean, log = log)
pexp_mean <- function(q, mean, ...) pexp(q, rate = 1 / mean, ...)

# names: fastmath names of the parameters in the order of the R estimates
case <- function(key, x, distr, start, names) {
  f <- if (is.null(start)) fitdist(x, distr) else fitdist(x, distr, start = start)
  se <- sqrt(diag(vcov(f)))
  est <- paste0("{", paste(paste0(":", names, " ", fmt(f$estimate)), collapse = " "), "}")
  ses <- paste0("{", paste(paste0(":", names, " ", fmt(se)), collapse = " "), "}")
  paste0("{:distribution ", key, " :params ", est, " :se ", ses, "\n :data ", vec(x), "}")
}

set.seed(20261030)
out <- c(
  case(":normal", rnorm(200, 1.5, 2), "norm", NULL, c("mu", "sd")),
  case(":gamma", rgamma(300, 2.5, 1 / 1.7), "gamma_scale", list(shape = 2, scale = 1.5), c("shape", "scale")),
  case(":weibull", rweibull(250, 1.8, 3), "weibull", NULL, c("alpha", "beta")),
  case(":beta", rbeta(200, 2, 5), "beta", NULL, c("alpha", "beta")),
  case(":log-normal", rlnorm(200, 0.5, 0.8), "lnorm", NULL, c("scale", "shape")),
  case(":exponential", rexp(150, 0.7), "exp_mean", list(mean = 1.5), c("mean")),
  case(":poisson", rpois(300, 4.2), "pois", NULL, c("p")),
  case(":cauchy", rcauchy(300, -1, 2), "cauchy", list(location = -1, scale = 2), c("median", "scale")))

dir.create("test/resources/fitdistr", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(out, collapse = "\n"), "]"), "test/resources/fitdistr/fit_se_reference.edn")
