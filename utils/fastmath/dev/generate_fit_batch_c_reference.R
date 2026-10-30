# Reference values for fastmath.fitdistr/fit of the heavy-tailed, location-scale and bounded families of batch C.
# Tools: R fitdistrplus::fitdist (method "mle") with the probability functions of stats (logis), extraDistr (laplace,
# frechet, triang, hcauchy, betapr), actuar (pareto1, llogis), VGAM (gumbel, bisa, foldnorm, levy), gamlss.dist (exGAUS),
# SuppDists (Johnson), circular::mle.vonmises; the one-parameter or independent versions (symmetrical beta, half-logistic,
# hyperbolic secant) are an optimized log likelihood written here from the density.
# A parameter that R's version has no counterpart for (the location of the half-Cauchy, Levy and fatigue life
# distributions, the minimum of the Pareto distribution) is fixed: `:fixed` is part of the output.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_batch_c_reference.R
# Output: test/resources/fitdistr/fit_batch_c_reference.edn
# The packages share names: the later ones win, actuar and gamlss.dist are loaded last
suppressMessages({library(fitdistrplus); library(extraDistr); library(VGAM); library(SuppDists); library(circular)
                  library(actuar); library(gamlss.dist)})

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

case <- function(key, x, fixed, params, loglik, method = ":mle") {
  paste0("{:distribution ", key, " :method ", method, " :fixed ", fixed, " :params ", params,
         "\n :data ", vec(x), "\n :loglik ", fmt(loglik), "}")
}

fit <- function(x, distr, start, fix.arg = NULL) {
  f <- fitdist(x, distr, start = start, fix.arg = fix.arg)
  list(e = f$estimate, loglik = f$loglik)
}

# wrappers of the functions of other packages with the arguments fitdist needs
dlaplace <- extraDistr::dlaplace; plaplace <- extraDistr::plaplace; qlaplace <- extraDistr::qlaplace
dpareto1 <- actuar::dpareto1; ppareto1 <- actuar::ppareto1; qpareto1 <- actuar::qpareto1
dllogis <- actuar::dllogis; pllogis <- actuar::pllogis; qllogis <- actuar::qllogis
dfrechet <- extraDistr::dfrechet; pfrechet <- extraDistr::pfrechet; qfrechet <- extraDistr::qfrechet
dtriang <- extraDistr::dtriang; ptriang <- extraDistr::ptriang; qtriang <- extraDistr::qtriang
dhcauchy <- extraDistr::dhcauchy; phcauchy <- extraDistr::phcauchy; qhcauchy <- extraDistr::qhcauchy
dbetapr <- extraDistr::dbetapr; pbetapr <- extraDistr::pbetapr; qbetapr <- extraDistr::qbetapr
dbisa <- VGAM::dbisa; pbisa <- VGAM::pbisa; qbisa <- VGAM::qbisa
dfoldnorm <- VGAM::dfoldnorm; pfoldnorm <- VGAM::pfoldnorm; qfoldnorm <- VGAM::qfoldnorm
dlevy <- VGAM::dlevy; plevy <- VGAM::plevy; qlevy <- VGAM::qlevy
dgumbel <- VGAM::dgumbel; pgumbel <- VGAM::pgumbel; qgumbel <- VGAM::qgumbel
# the Gumbel distribution of VGAM is the one of the maxima
stopifnot(abs(pgumbel(0.5, 1, 2) - exp(-exp(-(0.5 - 1) / 2))) < 1e-12)

set.seed(20261031)
out <- list()
add <- function(s) out[[length(out) + 1]] <<- s

x <- extraDistr::rlaplace(300, 1, 2); f <- fit(x, "laplace", list(mu = 0, sigma = 1))
add(case(":laplace", x, "{}", sprintf("{:mu %.17g :beta %.17g}", f$e["mu"], f$e["sigma"]), f$loglik))

x <- rlogis(300, 1, 2); f <- fit(x, "logis", list(location = 0, scale = 1))
add(case(":logistic", x, "{}", sprintf("{:mu %.17g :s %.17g}", f$e["location"], f$e["scale"]), f$loglik))

x <- VGAM::rgumbel(300, 1, 2); f <- fit(x, "gumbel", list(location = 0, scale = 1))
add(case(":gumbel", x, "{}", sprintf("{:mu %.17g :beta %.17g}", f$e["location"], f$e["scale"]), f$loglik))

x <- actuar::rpareto1(300, shape = 3, min = 1.5); f <- fit(x, "pareto1", list(shape = 2), list(min = min(x)))
add(case(":pareto", x, "{}", sprintf("{:scale %.17g :shape %.17g}", min(x), f$e["shape"]), f$loglik))

x <- extraDistr::rfrechet(300, lambda = 3, sigma = 2, mu = 0.5); f <- fit(x, "frechet", list(lambda = 2.5, sigma = 1.8, mu = 0.3))
add(case(":frechet", x, "{}", sprintf("{:alpha %.17g :beta %.17g :delta %.17g}", f$e["lambda"], f$e["sigma"], f$e["mu"]), f$loglik))

x <- actuar::rllogis(300, shape = 3, scale = 2); f <- fit(x, "llogis", list(shape = 2, scale = 1.5))
add(case(":log-logistic", x, "{}", sprintf("{:alpha %.17g :beta %.17g}", f$e["shape"], f$e["scale"]), f$loglik))

x <- VGAM::rbisa(300, scale = 2, shape = 0.7); f <- fit(x, "bisa", list(scale = 1.5, shape = 1))
add(case(":fatigue-life", x, "{:alpha 0.0}", sprintf("{:beta %.17g :gamma %.17g}", f$e["scale"], f$e["shape"]), f$loglik))

x <- extraDistr::rhcauchy(300, sigma = 2); f <- fit(x, "hcauchy", list(sigma = 1))
add(case(":half-cauchy", x, "{:mu 0.0}", sprintf("{:scale %.17g}", f$e["sigma"]), f$loglik))

x <- VGAM::rfoldnorm(300, mean = 1.5, sd = 1); f <- fit(x, "foldnorm", list(mean = 1, sd = 1.2), list(a1 = 1, a2 = 1))
add(case(":folded-normal", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", f$e["mean"], f$e["sd"]), f$loglik))

x <- gamlss.dist::rexGAUS(300, mu = 1, sigma = 0.5, nu = 1); f <- fit(x, "exGAUS", list(mu = 0.5, sigma = 0.3, nu = 0.5))
add(case(":exgaus", x, "{}", sprintf("{:mu %.17g :sigma %.17g :nu %.17g}", f$e["mu"], f$e["sigma"], f$e["nu"]), f$loglik))

x <- extraDistr::rbetapr(300, 3, 5, 2); f <- fit(x, "betapr", list(shape1 = 2, shape2 = 4, scale = 1.5))
add(case(":pearson-6", x, "{}", sprintf("{:alpha1 %.17g :alpha2 %.17g :beta %.17g}", f$e["shape1"], f$e["shape2"], f$e["scale"]), f$loglik))

x <- VGAM::rlevy(300, location = 0, scale = 2); f <- fit(x, "levy", list(scale = 1), list(location = 0))
add(case(":levy", x, "{:mu 0.0}", sprintf("{:c %.17g}", f$e["scale"]), f$loglik))

# the maximum likelihood estimate of the circular mean and the concentration
# fastmath's distribution lives on [mu - pi, mu + pi] and is not wrapped: the angles are moved into the window around the
# circular mean
th <- circular::rvonmises(300, mu = circular::circular(0.5), kappa = 2)
e <- circular::mle.vonmises(th)
mu <- as.numeric(e$mu)
x <- mu + atan2(sin(as.numeric(th) - mu), cos(as.numeric(th) - mu))
ll <- sum(circular::dvonmises(circular::circular(x), circular::circular(mu), e$kappa, log = TRUE))
add(case(":von-mises", x, "{}", sprintf("{:mu %.17g :kappa %.17g}", mu, e$kappa), ll))

# one parameter or independent versions: the log likelihood of the density is maximized
x <- rbeta(300, 3, 3)
o <- optimize(function(a) sum(dbeta(x, a, a, log = TRUE)), c(0.05, 50), maximum = TRUE, tol = 1e-12)
add(case(":beta-symmetrical", x, "{}", sprintf("{:alpha %.17g}", o$maximum), o$objective))

x <- abs(rlogis(300, 0, 2))
hl <- function(s) sum(log(2 / s) - x / s - 2 * log1p(exp(-x / s)))
o <- optimize(hl, c(0.05, 50), maximum = TRUE, tol = 1e-12)
add(case(":half-logistic", x, "{}", sprintf("{:scale %.17g}", o$maximum), o$objective))

u <- runif(300); x <- 1 + (2 * 2 / pi) * log(tan(pi * u / 2))
hs <- function(p) { s <- exp(p[2]); z <- pi * (x - p[1]) / (2 * s); sum(-log(2 * s) - log(cosh(z))) }
o <- optim(c(0, 0), function(p) -hs(p), method = "BFGS", control = list(reltol = 1e-14))
add(case(":hyperbolic-secant", x, "{}", sprintf("{:mu %.17g :sigma %.17g}", o$par[1], exp(o$par[2])), -o$value))

dir.create("test/resources/fitdistr", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(unlist(out), collapse = "\n"), "]"), "test/resources/fitdistr/fit_batch_c_reference.edn")
cat("cases:", length(out), "\n")
