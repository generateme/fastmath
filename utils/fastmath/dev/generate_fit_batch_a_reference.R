# Reference values for fastmath.random.fitting/fit of the families of the metadata batch A that R supports.
# Tools: R fitdistrplus::fitdist (method "mle") with the density, cdf and quantile functions of stats (chisq, f, t),
# statmod (invgauss), actuar (invgamma), extraDistr (rayleigh, hnorm), VGAM (naka); fitdistrplus::mmedist for the beta
# method of moments. The uniform distribution has no numerical fit: its estimate is the sample range.
# Families with a location or an integer parameter that R's version does not have are fitted with that parameter
# fixed (rayleigh and half-normal at 0, erlang with the integer shape): `:fixed` is part of the output.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_batch_a_reference.R
# Output: test/resources/random/fit_batch_a_reference.edn
# actuar is loaded last: extraDistr has an inverse gamma with other parameter names
suppressMessages({library(fitdistrplus); library(extraDistr); library(VGAM); library(statmod); library(actuar)})

# VGAM masks the Rayleigh functions of extraDistr
drayleigh <- function(x, sigma, log = FALSE) extraDistr::drayleigh(x, sigma, log)
prayleigh <- function(q, sigma, lower.tail = TRUE, log.p = FALSE) extraDistr::prayleigh(q, sigma, lower.tail, log.p)
qrayleigh <- function(p, sigma, lower.tail = TRUE, log.p = FALSE) extraDistr::qrayleigh(p, sigma, lower.tail, log.p)

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

case <- function(key, x, fixed, params, loglik, method = ":mle") {
  paste0("{:distribution ", key, " :method ", method, " :fixed ", fixed, " :params ", params,
         "\n :data ", vec(x), "\n :loglik ", if (is.na(loglik)) "nil" else fmt(loglik), "}")
}

fit <- function(x, distr, start = NULL, fix.arg = NULL, ...) {
  f <- fitdist(x, distr, start = start, fix.arg = fix.arg, ...)
  list(e = f$estimate, loglik = f$loglik)
}

set.seed(20261030)
out <- list()

x <- runif(100, -1, 3)
out[[1]] <- case(":uniform-real", x, "{}", sprintf("{:lower %.17g :upper %.17g}", min(x), max(x)),
                 -length(x) * log(max(x) - min(x)))

x <- rchisq(200, 4); f <- fit(x, "chisq", list(df = 3))
out[[2]] <- case(":chi-squared", x, "{}", sprintf("{:degrees-of-freedom %.17g}", f$e["df"]), f$loglik)

x <- rf(300, 5, 12); f <- fit(x, "f", list(df1 = 3, df2 = 8))
out[[3]] <- case(":f", x, "{}", sprintf("{:numerator-degrees-of-freedom %.17g :denominator-degrees-of-freedom %.17g}",
                                          f$e["df1"], f$e["df2"]), f$loglik)

x <- rt(300, 6); f <- fit(x, "t", list(df = 10))
out[[4]] <- case(":t", x, "{}", sprintf("{:degrees-of-freedom %.17g}", f$e["df"]), f$loglik)

x <- rinvgauss(200, mean = 2, shape = 3); f <- fit(x, "invgauss", list(mean = 1, shape = 1))
out[[5]] <- case(":inverse-gaussian", x, "{}", sprintf("{:mu %.17g :lambda %.17g}", f$e["mean"], f$e["shape"]), f$loglik)

x <- actuar::rinvgamma(300, shape = 3, scale = 2); f <- fit(x, "invgamma", list(shape = 2, scale = 1))
out[[6]] <- case(":inverse-gamma", x, "{}", sprintf("{:alpha %.17g :beta %.17g}", f$e["shape"], f$e["scale"]), f$loglik)

x <- extraDistr::rrayleigh(300, sigma = 2); f <- fit(x, "rayleigh", list(sigma = 1))
out[[7]] <- case(":rayleigh", x, "{:a 0.0}", sprintf("{:beta %.17g}", f$e["sigma"]), f$loglik)

x <- rhnorm(300, sigma = 1.5); f <- fit(x, "hnorm", list(sigma = 1))
out[[8]] <- case(":half-normal", x, "{:mu 0.0}", sprintf("{:sigma %.17g}", f$e["sigma"]), f$loglik)

x <- rgamma(200, shape = 3, rate = 0.8); f <- fit(x, "gamma", list(rate = 1), fix.arg = list(shape = 3))
out[[9]] <- case(":erlang", x, "{:k 3}", sprintf("{:lambda %.17g}", f$e["rate"]), f$loglik)

x <- rnaka(300, scale = 2, shape = 1.5); f <- fit(x, "naka", list(scale = 1, shape = 1))
out[[10]] <- case(":nakagami", x, "{}", sprintf("{:mu %.17g :omega %.17g}", f$e["shape"], f$e["scale"]), f$loglik)

x <- rbeta(200, 2, 5); e <- mmedist(x, "beta")$estimate
out[[11]] <- case(":beta", x, "{}", sprintf("{:alpha %.17g :beta %.17g}", e["shape1"], e["shape2"]), NA, ":mme")

dir.create("test/resources/random", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(unlist(out), collapse = "\n"), "]"), "test/resources/random/fit_batch_a_reference.edn")
