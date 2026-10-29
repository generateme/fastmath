# Reference values for fastmath.random.fitting/fit with the :mle method.
# Tools: R fitdistrplus::fitdist(method = "mle"), fitdistrplus::gofstat (KS, CvM, AD statistics).
# Parameters are mapped to the fastmath names; AICc is computed here from the log likelihood.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_mle_reference.R
# Output: test/resources/random/fit_mle_reference.edn
suppressMessages(library(fitdistrplus))

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

# key: fastmath key; map: function from the R estimates to the fastmath parameter map text
case <- function(key, x, distr, start, map) {
  f <- if (is.null(start)) fitdist(x, distr) else fitdist(x, distr, start = start)
  g <- gofstat(f)
  n <- length(x); k <- length(f$estimate)
  aicc <- f$aic + 2 * k * (k + 1) / (n - k - 1)
  paste0("{:distribution ", key, " :params ", map(f$estimate), "\n :data ", vec(x),
         "\n :loglik ", fmt(f$loglik), " :aic ", fmt(f$aic), " :aicc ", fmt(aicc), " :bic ", fmt(f$bic),
         " :k ", k, " :n ", n,
         "\n :ks ", fmt(g$ks), " :cvm ", fmt(g$cvm), " :ad ", fmt(g$ad), "}")
}

set.seed(20261029)
out <- c(
  case(":normal", rnorm(200, 1.5, 2), "norm", NULL,
       function(e) sprintf("{:mu %.17g :sd %.17g}", e["mean"], e["sd"])),
  case(":exponential", rexp(150, 0.7), "exp", NULL,
       function(e) sprintf("{:mean %.17g}", 1 / e["rate"])),
  case(":gamma", rgamma(300, 2.5, 1 / 1.7), "gamma", NULL,
       function(e) sprintf("{:shape %.17g :scale %.17g}", e["shape"], 1 / e["rate"])),
  case(":weibull", rweibull(250, 1.8, 3), "weibull", NULL,
       function(e) sprintf("{:alpha %.17g :beta %.17g}", e["shape"], e["scale"])),
  case(":beta", rbeta(200, 2, 5), "beta", NULL,
       function(e) sprintf("{:alpha %.17g :beta %.17g}", e["shape1"], e["shape2"])),
  case(":log-normal", rlnorm(200, 0.5, 0.8), "lnorm", NULL,
       function(e) sprintf("{:scale %.17g :shape %.17g}", e["meanlog"], e["sdlog"])),
  case(":cauchy", rcauchy(300, -1, 2), "cauchy", list(location = -1, scale = 2),
       function(e) sprintf("{:median %.17g :scale %.17g}", e["location"], e["scale"])))

dir.create("test/resources/random", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(out, collapse = "\n"), "]"), "test/resources/random/fit_mle_reference.edn")
