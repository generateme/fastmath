# Reference values for fastmath.fitdistr/fit of the discrete families of the metadata batch B.
# Tools: R fitdistrplus::fitdist (method "mle") with the probability functions of stats (pois, geom, binom, nbinom),
# extraDistr (lgser, bbinom), VGAM (zipf) and gamlss.dist (NBI, NBII, BB). The discrete uniform distribution has no
# numerical fit: its estimate is the sample range.
# Parameters R does not fit, the integer ones (binomial trials, Pascal successes, Zipf elements, beta-binomial trials),
# are fixed: `:fixed` is part of the output; fastmath fits `:trials` and `:r` without fixing them (see the tests).
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_batch_b_reference.R
# Output: test/resources/fitdistr/fit_batch_b_reference.edn
# gamlss.dist is loaded last: the packages share some names
suppressMessages({library(fitdistrplus); library(extraDistr); library(VGAM); library(gamlss.dist)})

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

case <- function(key, x, fixed, params, loglik) {
  paste0("{:distribution ", key, " :method :mle :fixed ", fixed, " :params ", params,
         "\n :data ", vec(x), "\n :loglik ", fmt(loglik), "}")
}

fit <- function(x, distr, start = NULL, fix.arg = NULL) {
  f <- fitdist(x, distr, start = start, fix.arg = fix.arg, discrete = TRUE)
  list(e = f$estimate, loglik = f$loglik)
}

set.seed(20261030)
out <- list()

x <- rpois(300, 3); f <- fit(x, "pois")
out[[1]] <- case(":poisson", x, "{}", sprintf("{:p %.17g}", f$e["lambda"]), f$loglik)

x <- rgeom(300, 0.3); f <- fit(x, "geom")
out[[2]] <- case(":geometric", x, "{}", sprintf("{:p %.17g}", f$e["prob"]), f$loglik)

x <- rbinom(300, 25, 0.3); f <- fit(x, "binom", list(prob = 0.5), list(size = 25))
out[[3]] <- case(":binomial", x, "{:trials 25}", sprintf("{:p %.17g}", f$e["prob"]), f$loglik)

x <- rnbinom(400, size = 3.5, prob = 0.4); f <- fit(x, "nbinom")
out[[4]] <- case(":negative-binomial", x, "{}",
                 sprintf("{:r %.17g :p %.17g}", f$e["size"], f$e["size"] / (f$e["size"] + f$e["mu"])), f$loglik)

x <- rnbinom(300, size = 4, prob = 0.4); f <- fit(x, "nbinom", list(prob = 0.5), list(size = 4))
out[[5]] <- case(":pascal", x, "{:r 4}", sprintf("{:p %.17g}", f$e["prob"]), f$loglik)

x <- rNBI(400, mu = 3, sigma = 0.5); f <- fit(x, "NBI", list(mu = 2, sigma = 1))
out[[6]] <- case(":nbi", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik)

x <- rNBII(400, mu = 3, sigma = 1.5); f <- fit(x, "NBII", list(mu = 2, sigma = 1))
out[[7]] <- case(":nbii", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik)

x <- rlgser(300, 0.7); f <- fit(x, "lgser", list(theta = 0.5))
out[[8]] <- case(":logarithmic", x, "{}", sprintf("{:p %.17g}", f$e["theta"]), f$loglik)

# fitdist fails for zipf (the optimizer stops with the error code 100): one parameter, so optimize the log likelihood
x <- rzipf(300, N = 20, shape = 1.5)
zipf.loglik <- function(s) sum(-s * log(x)) - length(x) * log(sum((1:20)^(-s)))
o <- optimize(zipf.loglik, c(0.01, 10), maximum = TRUE, tol = 1e-12)
out[[9]] <- case(":zipf", x, "{:number-of-elements 20}", sprintf("{:exponent %.17g}", o$maximum), o$objective)

x <- rbbinom(300, 12, 2, 3); f <- fit(x, "bbinom", list(alpha = 1, beta = 1), list(size = 12))
out[[10]] <- case(":beta-binomial", x, "{:n 12}", sprintf("{:alpha %.17g :beta %.17g}", f$e["alpha"], f$e["beta"]), f$loglik)

x <- rBB(300, mu = 0.4, sigma = 0.5, bd = 12); f <- fit(x, "BB", list(mu = 0.5, sigma = 1), list(bd = 12))
out[[11]] <- case(":bb", x, "{:bd 12}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mu"], f$e["sigma"]), f$loglik)

x <- sample(3:12, 200, TRUE)
out[[12]] <- case(":uniform-int", x, "{}", sprintf("{:lower %d :upper %d}", min(x), max(x)), -length(x) * log(max(x) - min(x) + 1))

dir.create("test/resources/fitdistr", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(unlist(out), collapse = "\n"), "]"), "test/resources/fitdistr/fit_batch_b_reference.edn")
