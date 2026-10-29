# Reference values for fastmath.random.fitting/fit of the families of the metadata batch D which R supports.
# Tools: R fitdistrplus::fitdist (method "mle") with the probability functions of stats (noncentral chisq, f, t, beta),
# evd (gev, gpd), gnorm (generalized normal), gamlss.dist (GG, ZIP, ZIP2, ZAP, ZIBI, ZABI, ZINBI, ZANBI, ZIBB, ZABB, ZAGA,
# ZAIG), GeneralizedHyperbolic (generalized hyperbolic and normal inverse Gaussian, optim of the log likelihood).
# Parameters that R's version fits differently are fixed: the location of the generalized Pareto distribution at 0 and the
# integer binomial denominators `bd`: `:fixed` is part of the output.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_batch_d_reference.R
# Output: test/resources/random/fit_batch_d_reference.edn
# gamlss.dist is loaded last: the packages share some names
suppressMessages({library(fitdistrplus); library(evd); library(gnorm); library(GeneralizedHyperbolic); library(gamlss.dist)})

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

case <- function(key, x, fixed, params, loglik) {
  paste0("{:distribution ", key, " :method :mle :fixed ", fixed, " :params ", params,
         "\n :data ", vec(x), "\n :loglik ", fmt(loglik), "}")
}

# bounded fits (L-BFGS-B): probabilities in (0, 1), scales and shapes positive
eps <- 1e-6
fit <- function(x, distr, start, fix.arg = NULL, discrete = FALSE, lower = -Inf, upper = Inf) {
  f <- fitdist(x, distr, start = start, fix.arg = fix.arg, discrete = discrete, lower = lower, upper = upper)
  list(e = f$estimate, loglik = f$loglik)
}

set.seed(20261101)
out <- list()
add <- function(s) out[[length(out) + 1]] <<- s

x <- evd::rgev(400, 1, 2, 0.2); f <- fit(x, "gev", list(loc = 0.8, scale = 1.8, shape = 0.1), lower = c(-Inf, eps, -Inf))
add(case(":generalized-extreme-value", x, "{}", sprintf("{:mu %.17g :sigma %.17g :xi %.17g}", f$e["loc"], f$e["scale"], f$e["shape"]), f$loglik))

x <- evd::rgpd(400, 0, 2, 0.2); f <- fit(x, "gpd", list(scale = 1.5, shape = 0.1), list(loc = 0), lower = c(eps, -Inf))
add(case(":generalized-pareto", x, "{:mu 0.0}", sprintf("{:sigma %.17g :xi %.17g}", f$e["scale"], f$e["shape"]), f$loglik))

x <- gnorm::rgnorm(400, 1, 2, 3); f <- fit(x, "gnorm", list(mu = 0.5, alpha = 1.5, beta = 2.5), lower = c(-Inf, eps, eps))
add(case(":generalized-normal", x, "{}", sprintf("{:mu %.17g :alpha %.17g :beta %.17g}", f$e["mu"], f$e["alpha"], f$e["beta"]), f$loglik))

x <- gamlss.dist::rGG(400, mu = 2, sigma = 0.5, nu = 1.5); f <- fit(x, "GG", list(mu = 1.5, sigma = 0.6, nu = 1), lower = c(eps, eps, -Inf))
add(case(":generalized-gamma", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rchisq(400, 4, ncp = 2); f <- fit(x, "chisq", list(df = 3, ncp = 1), lower = c(eps, 0))
add(case(":chi-squared-noncentral", x, "{}", sprintf("{:nu %.17g :lambda %.17g}", f$e["df"], f$e["ncp"]), f$loglik))

x <- rf(300, 5, 12, ncp = 2); f <- fit(x, "f", list(df1 = 4, df2 = 10, ncp = 1), lower = c(eps, eps, 0))
add(case(":f-noncentral", x, "{}", sprintf("{:df1 %.17g :df2 %.17g :ncp %.17g}", f$e["df1"], f$e["df2"], f$e["ncp"]), f$loglik))

x <- rt(300, 6, ncp = 1); f <- fit(x, "t", list(df = 5, ncp = 0.5), lower = c(eps, -Inf))
add(case(":t-noncentral", x, "{}", sprintf("{:df %.17g :ncp %.17g}", f$e["df"], f$e["ncp"]), f$loglik))

x <- rbeta(300, 2, 4, ncp = 1.5); f <- fit(x, "beta", list(shape1 = 1.5, shape2 = 3, ncp = 1), lower = c(eps, eps, 0))
add(case(":beta-noncentral", x, "{}", sprintf("{:alpha %.17g :beta %.17g :ncp %.17g}", f$e["shape1"], f$e["shape2"], f$e["ncp"]), f$loglik))

# the zero-inflated and zero-adjusted families of gamlss.dist
x <- rZIP(500, mu = 3, sigma = 0.2); f <- fit(x, "ZIP", list(mu = 2, sigma = 0.1), discrete = TRUE, lower = c(eps, eps), upper = c(Inf, 1 - eps))
add(case(":zero-inflated-poisson", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rZIP2(500, mu = 3, sigma = 0.2); f <- fit(x, "ZIP2", list(mu = 2, sigma = 0.1), discrete = TRUE, lower = c(eps, eps), upper = c(Inf, 1 - eps))
add(case(":zero-inflated-poisson2", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rZAP(500, mu = 3, sigma = 0.2); f <- fit(x, "ZAP", list(mu = 2, sigma = 0.1), discrete = TRUE, lower = c(eps, eps), upper = c(Inf, 1 - eps))
add(case(":zero-adjusted-poisson", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rZIBI(500, mu = 0.4, sigma = 0.2, bd = 10); f <- fit(x, "ZIBI", list(mu = 0.5, sigma = 0.1), list(bd = 10), TRUE, lower = c(eps, eps), upper = c(1 - eps, 1 - eps))
add(case(":zero-inflated-binomial", x, "{:bd 10}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rZABI(500, mu = 0.4, sigma = 0.2, bd = 10); f <- fit(x, "ZABI", list(mu = 0.5, sigma = 0.1), list(bd = 10), TRUE, lower = c(eps, eps), upper = c(1 - eps, 1 - eps))
add(case(":zero-adjusted-binomial", x, "{:bd 10}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rZINBI(500, mu = 3, sigma = 0.5, nu = 0.2); f <- fit(x, "ZINBI", list(mu = 2, sigma = 0.8, nu = 0.1), discrete = TRUE, lower = c(eps, eps, eps), upper = c(Inf, Inf, 1 - eps))
add(case(":zero-inflated-negative-binomial", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rZANBI(500, mu = 3, sigma = 0.5, nu = 0.2); f <- fit(x, "ZANBI", list(mu = 2, sigma = 0.8, nu = 0.1), discrete = TRUE, lower = c(eps, eps, eps), upper = c(Inf, Inf, 1 - eps))
add(case(":zero-adjusted-negative-binomial", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rZIBB(500, mu = 0.4, sigma = 0.5, nu = 0.2, bd = 10); f <- fit(x, "ZIBB", list(mu = 0.5, sigma = 0.8, nu = 0.1), list(bd = 10), TRUE, lower = c(eps, eps, eps), upper = c(1 - eps, Inf, 1 - eps))
add(case(":zero-inflated-beta-binomial", x, "{:bd 10}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rZABB(500, mu = 0.4, sigma = 0.5, nu = 0.2, bd = 10); f <- fit(x, "ZABB", list(mu = 0.5, sigma = 0.8, nu = 0.1), list(bd = 10), TRUE, lower = c(eps, eps, eps), upper = c(1 - eps, Inf, 1 - eps))
add(case(":zero-adjusted-beta-binomial", x, "{:bd 10}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rZAGA(500, mu = 2, sigma = 0.5, nu = 0.2); f <- fit(x, "ZAGA", list(mu = 1.5, sigma = 0.6, nu = 0.1), lower = c(eps, eps, eps), upper = c(Inf, Inf, 1 - eps))
add(case(":zero-adjusted-gamma", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- rZAIG(500, mu = 2, sigma = 0.5, nu = 0.2); f <- fit(x, "ZAIG", list(mu = 1.5, sigma = 0.6, nu = 0.1), lower = c(eps, eps, eps), upper = c(Inf, Inf, 1 - eps))
add(case(":zero-adjusted-inverse-gaussian", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

# the generalized hyperbolic family: an optimized log likelihood (these fits are slow in R as well)
x <- GeneralizedHyperbolic::rghyp(300, param = c(0, 1, 2, 0.5, 1))
nll <- function(p) { v <- try(-sum(log(GeneralizedHyperbolic::dghyp(x, param = c(p[1], exp(p[2]), exp(p[3]), p[4], p[5])))), silent = TRUE)
                     if (inherits(v, "try-error") || !is.finite(v)) 1e10 else v }
o <- optim(c(0, 0, log(2), 0.3, 0.5), nll, control = list(maxit = 2000, reltol = 1e-12))
add(case(":generalized-hyperbolic", x, "{}", sprintf("{:mu %.17g :delta %.17g :alpha %.17g :beta %.17g :lambda %.17g}",
                                                    o$par[1], exp(o$par[2]), exp(o$par[3]), o$par[4], o$par[5]), -o$value))

dir.create("test/resources/random", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(unlist(out), collapse = "\n"), "]"), "test/resources/random/fit_batch_d_reference.edn")
cat("cases:", length(out), "\n")
