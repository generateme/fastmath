# Generator of R reference values for the regression audit
# (test/fastmath/ml/regression_test.clj and test/fastmath/ml/regression/contrast_test.clj).
#
# Writes test/resources/regression/audit_reference.edn
# Run from the repository root:
#   Rscript utils/fastmath/dev/regression_audit_ref_gen.R
#
# Needs R packages: stats, MASS. All data are fixed vectors, except the
# negative binomial samples (set.seed(42), embedded in the output).

suppressMessages(library(MASS))
options(digits = 17, warn = -1)

# ---- minimal EDN writer ---------------------------------------------------

num <- function(v) {
  if (is.na(v) && !is.nan(v)) "nil"
  else if (is.nan(v)) "##NaN"
  else if (v == Inf) "##Inf"
  else if (v == -Inf) "##-Inf"
  else sprintf("%.17g", v)
}

edn <- function(x) {
  if (is.null(x)) return("nil")
  if (is.character(x)) {
    s <- sprintf("\"%s\"", x)
    return(if (length(x) == 1) s else paste0("[", paste(s, collapse = " "), "]"))
  }
  if (is.logical(x)) {
    s <- ifelse(x, "true", "false")
    return(if (length(x) == 1) s else paste0("[", paste(s, collapse = " "), "]"))
  }
  if (is.matrix(x)) {
    return(paste0("[", paste(sapply(seq_len(nrow(x)), function(i) edn(I(as.vector(x[i, ])))),
                             collapse = " "), "]"))
  }
  if (is.list(x)) {
    if (!is.null(names(x))) {
      keys <- gsub("[._]", "-", names(x))
      return(paste0("{", paste(sprintf(":%s %s", keys, sapply(x, edn)), collapse = " "), "}"))
    }
    return(paste0("[", paste(sapply(x, edn), collapse = " "), "]"))
  }
  s <- sapply(as.vector(x), num)
  # I() forces a vector for length one
  if (length(x) == 1 && !inherits(x, "AsIs")) s else paste0("[", paste(s, collapse = " "), "]")
}

ctl <- glm.control(epsilon = 1e-14, maxit = 100)

# a vector, also when it has one element
V <- function(x) I(as.vector(x))

# ---- model summaries -------------------------------------------------------

lm_info <- function(l, influence = TRUE) {
  s <- summary(l)
  out <- list(
    intercept = if ("(Intercept)" %in% names(coef(l))) unname(coef(l)["(Intercept)"]) else 0,
    beta = V(unname(coef(l)[names(coef(l)) != "(Intercept)"])),
    stderr = V(unname(s$coefficients[, 2])),
    sigma = s$sigma,
    df_residual = l$df.residual,
    r_squared = s$r.squared,
    adjusted_r_squared = s$adj.r.squared,
    f_statistic = if (is.null(s$fstatistic)) NaN else unname(s$fstatistic[1]),
    f_p_value = if (is.null(s$fstatistic)) NaN else
      unname(pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)),
    log_likelihood = as.numeric(logLik(l)),
    aic = AIC(l),
    bic = BIC(l)
  )
  if (influence) {
    im <- influence.measures(l)$infmat
    out$hat <- V(unname(im[, "hat"]))
    out$studentized <- V(unname(rstudent(l)))
    out$cooks_distance <- V(unname(im[, "cook.d"]))
    out$dffits <- V(unname(im[, "dffit"]))
    out$covratio <- V(unname(im[, "cov.r"]))
    out$dfbetas <- unname(t(dfbetas(l)))
  }
  out
}

glm_info <- function(f) {
  s <- summary(f)
  list(
    intercept = unname(coef(f)[1]),
    beta = V(unname(coef(f)[-1])),
    stderr = V(unname(s$coefficients[, 2])),
    deviance = deviance(f),
    null_deviance = f$null.deviance,
    df_residual = f$df.residual,
    df_null = f$df.null,
    dispersion = s$dispersion,
    aic = f$aic,
    log_likelihood = as.numeric(logLik(f)),
    # BIC of R counts the rows with a zero weight, here only the positive weights (as for lm)
    bic = -2 * as.numeric(logLik(f)) + attr(logLik(f), "df") * log(sum(f$prior.weights > 0)),
    hat = V(unname(hatvalues(f))),
    cooks_distance = V(unname(cooks.distance(f)))
  )
}

out <- list()
x <- 1:10

# ---- 1. zero weights ----------------------------------------------------------
yg <- c(1.3, 2.1, 3.8, 3.9, 5.4, 5.7, 7.9, 8.6, 9.2, 10.9)
w  <- c(1, 2, 0, 1, 1.5, 1, 0, 2, 1, 1)
pb <- c(0.1, 0.17, 0.5, 0.33, 0.5, 0.5, 0.5, 0.75, 0.9, 0.89)
nb <- c(10, 12, 0, 15, 10, 8, 0, 12, 10, 9)
yc <- c(2, 3, 6, 5, 9, 12, 0, 15, 19, 22)
wc <- c(1, 1, 1, 0, 2, 1, 1, 0, 1, 1)
zw <- list(
  data = list(x = x, y = yg, w = w, pb = pb, nb = nb, yc = yc, wc = wc),
  lm = lm_info(lm(yg ~ x, weights = w)),
  glm_gaussian = glm_info(glm(yg ~ x, weights = w, control = ctl)),
  glm_binomial = glm_info(glm(pb ~ x, weights = nb, family = binomial, control = ctl)),
  glm_poisson = glm_info(glm(yc ~ x, weights = wc, family = poisson, control = ctl))
)
# the gaussian AIC of glm is infinite in R for zero weights (log(0)): use the lm value
zw$glm_gaussian$aic <- AIC(lm(yg ~ x, weights = w))
zw$glm_gaussian$bic <- BIC(lm(yg ~ x, weights = w))
zw$glm_gaussian$log_likelihood <- as.numeric(logLik(lm(yg ~ x, weights = w)))
out$zero_weights <- zw

# ---- 2. intercept only ----------------------------------------------------------
io <- function(l) {
  p <- predict(l, data.frame(x = 1), se.fit = TRUE, interval = "confidence")
  pp <- predict(l, data.frame(x = 1), interval = "prediction", weights = 1)
  list(info = lm_info(l, FALSE),
       predict = list(fit = unname(p$fit[1, "fit"]), stderr = p$se.fit,
                      confidence_interval = unname(p$fit[1, c("lwr", "upr")]),
                      prediction_interval = unname(pp[1, c("lwr", "upr")])))
}
out$intercept_only <- list(
  data = list(y = yg, w = c(1, 2, 1, 1, 1.5, 1, 3, 2, 1, 1)),
  ols = io(lm(yg ~ 1)),
  wls = io(lm(yg ~ 1, weights = c(1, 2, 1, 1, 1.5, 1, 3, 2, 1, 1))),
  zero_weights = io(lm(yg ~ 1, weights = w))
)

# ---- 3. binomial pairs with offset -------------------------------------------------
xp <- 1:8; s <- c(1, 2, 3, 5, 6, 8, 9, 9); n <- rep(10, 8)
off <- seq(0, 0.7, length.out = 8)
out$pairs_offset <- list(
  data = list(x = xp, successes = s, failures = n - s, offset = off),
  with_offset = glm_info(glm(cbind(s, n - s) ~ xp, family = binomial, offset = off, control = ctl)),
  no_offset = glm_info(glm(cbind(s, n - s) ~ xp, family = binomial, control = ctl)),
  probit_offset = glm_info(glm(cbind(s, n - s) ~ xp, family = binomial(link = "probit"),
                               offset = off, control = ctl))
)

# binomial pairs with a zero number of trials: weight 0, R sets the proportion to 0
ztx <- 1:8; zts <- c(0, 2, 3, 5, 6, 8, 9, 9); ztf <- c(0, 8, 7, 5, 4, 2, 1, 1)
out$zero_trials <- list(
  data = list(x = ztx, successes = zts, failures = ztf),
  fit = glm_info(glm(cbind(zts, ztf) ~ ztx, family = binomial, control = ctl))
)

# ---- 4. prediction intervals ----------------------------------------------------------
xn <- c(3.5, 12)
offn <- c(2, 11)
pi_info <- function(l, nd) {
  p <- predict(l, nd, se.fit = TRUE, interval = "confidence")
  pp <- predict(l, nd, interval = "prediction", weights = rep(1, nrow(nd)))
  list(fit = V(unname(p$fit[, "fit"])), stderr = V(unname(p$se.fit)),
       confidence_interval = lapply(seq_len(nrow(nd)), function(i) unname(p$fit[i, c("lwr", "upr")])),
       prediction_interval = lapply(seq_len(nrow(nd)), function(i) unname(pp[i, c("lwr", "upr")])))
}
offx <- seq(0, 9)
out$prediction_interval <- list(
  data = list(x = x, y = yg, w = c(1, 2, 1, 1, 1.5, 1, 3, 2, 1, 1), offset = offx, new_x = xn, new_offset = offn),
  ols = pi_info(lm(yg ~ x), data.frame(x = xn)),
  wls = pi_info(lm(yg ~ x, weights = c(1, 2, 1, 1, 1.5, 1, 3, 2, 1, 1)), data.frame(x = xn)),
  no_intercept = pi_info(lm(yg ~ x - 1), data.frame(x = xn)),
  offset = pi_info(lm(yg ~ x + offset(offx)), data.frame(x = xn, offx = offn))
)

# ---- 5. negative binomial: theta and its standard error -------------------------------------
set.seed(42)
gen <- function(n, p, theta, off = FALSE) {
  X <- matrix(rnorm(n * p, sd = 0.5), n, p); b <- c(1.0, 0.5, -0.3)[1:(p + 1)]
  o <- if (off) log(runif(n, 1, 5)) else rep(0, n)
  mu <- exp(b[1] + X %*% b[-1] + o); y <- rnbinom(n, mu = mu, size = theta)
  list(X = X, y = y, o = o)
}
nbs <- list()
for (cfg in list(c(60, 2, 1.5, 0), c(200, 2, 3, 0), c(300, 1, 0.7, 1), c(80, 2, 10, 0))) {
  d <- gen(cfg[1], cfg[2], cfg[3], cfg[4] == 1)
  f <- if (cfg[4] == 1) glm.nb(d$y ~ d$X + offset(d$o)) else glm.nb(d$y ~ d$X)
  nbs[[length(nbs) + 1]] <- list(
    offset = cfg[4] == 1, y = d$y, x = d$X, o = d$o,
    intercept = unname(coef(f)[1]), beta = V(unname(coef(f)[-1])),
    stderr = V(unname(sqrt(diag(vcov(f))))), theta = f$theta, theta_stderr = f$SE.theta,
    log_likelihood = as.numeric(logLik(f)), aic = AIC(f), deviance = deviance(f),
    null_deviance = f$null.deviance)
}
out$nbinomial <- nbs

# ---- 6. links and families ------------------------------------------------------------------
xl <- 1:8
yl <- c(1.2, 2.1, 3.5, 3.9, 5.2, 5.8, 7.9, 8.1)
ycl <- c(2, 3, 6, 7, 8, 9, 10, 12)
sl <- c(1, 2, 3, 5, 6, 8, 9, 9)
lk <- function(f) {
  s <- summary(f)
  list(intercept = unname(coef(f)[1]), beta = V(unname(coef(f)[-1])), stderr = V(unname(s$coefficients[, 2])),
       deviance = deviance(f), null_deviance = f$null.deviance, aic = f$aic, dispersion = s$dispersion)
}
out$links <- list(
  data = list(x = xl, y = yl, yc = ycl, successes = sl, failures = 10 - sl),
  poisson_sqrt = lk(glm(ycl ~ xl, family = poisson(link = "sqrt"), control = ctl)),
  poisson_identity = lk(glm(ycl ~ xl, family = poisson(link = "identity"), control = ctl)),
  binomial_cauchit = lk(glm(cbind(sl, 10 - sl) ~ xl, family = binomial(link = "cauchit"), control = ctl)),
  quasibinomial = lk(glm(cbind(sl, 10 - sl) ~ xl, family = quasibinomial, control = ctl)),
  gaussian_log = lk(glm(yl ~ xl, family = gaussian(link = "log"), control = ctl)),
  inverse_gaussian_inverse = lk(glm(yl ~ xl, family = inverse.gaussian(link = "inverse"), control = ctl)),
  gamma_inverse = lk(glm(yl ~ xl, family = Gamma(link = "inverse"), control = ctl)),
  poisson_power0 = lk(glm(ycl ~ xl, family = poisson(link = "log"), control = ctl)),
  poisson_power05 = lk(glm(ycl ~ xl, family = poisson(link = power(0.5)), control = ctl)),
  # non-integer numbers of successes: R warns and rounds half to even in the AIC
  binomial_noninteger = lk(glm(c(0.2, 0.3, 0.4, 0.55, 0.6, 0.8, 0.9, 0.95) ~ xl, family = binomial,
                               weights = c(3, 5, 2, 4, 6, 3, 5, 2), control = ctl))
)
# dispersion estimators of a poisson model: Pearson and mean deviance
fp <- glm(ycl ~ xl, family = poisson, control = ctl)
out$dispersion <- list(
  pearson = sum(residuals(fp, "pearson")^2) / fp$df.residual,
  mean_deviance = deviance(fp) / fp$df.residual,
  stderr_pearson = V(unname(sqrt(diag(vcov(fp)) * sum(residuals(fp, "pearson")^2) / fp$df.residual))),
  stderr_mean_deviance = V(unname(sqrt(diag(vcov(fp)) * deviance(fp) / fp$df.residual)))
)

# ---- 7. transformers -------------------------------------------------------------------------------
xt <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10.5)
pp <- poly(xt, 3)
newt <- c(0.5, 2.5, 7, 12)
library(splines)
xl_b <- 0; xr_b <- 1; nseg <- 4; degb <- 3; dx <- (xr_b - xl_b) / nseg
knots <- seq(xl_b - degb * dx, xr_b + degb * dx, by = dx)
xb <- c(0, 0.1, 0.3, 0.5, 0.77, 0.999, 1)
out$transformers <- list(
  poly = list(x = xt, degree = 3, new_x = newt,
              train = pp[, 1:3], new = predict(pp, newt)),
  bspline = list(xl = xl_b, xr = xr_b, nseg = nseg, degree = degb, x = xb,
                 basis = splineDesign(knots, xb, ord = degb + 1, outer.ok = TRUE)),
  bspline_degrees = lapply(1:6, function(deg) {
    kn <- seq(xl_b - deg * dx, xr_b + deg * dx, by = dx)
    list(degree = deg, basis = splineDesign(kn, xb, ord = deg + 1, outer.ok = TRUE))
  })
)

# ---- 8. contrasts -----------------------------------------------------------------------------------
cs <- list()
for (k in 2:6) {
  cs[[length(cs) + 1]] <- list(n = k, treatment = contr.treatment(k), sum = contr.sum(k),
                               helmert = contr.helmert(k), SAS = contr.SAS(k),
                               poly = if (k >= 2) contr.poly(k) else NULL)
}
out$contrasts <- cs

dir.create("test/resources/regression", showWarnings = FALSE, recursive = TRUE)
con <- file("test/resources/regression/audit_reference.edn", "w")
writeLines(c(";; generated by utils/fastmath/dev/regression_audit_ref_gen.R, do not edit",
             paste0(";; R ", R.version.string, ", MASS ", packageVersion("MASS")),
             edn(out)), con)
close(con)
cat("done\n")
