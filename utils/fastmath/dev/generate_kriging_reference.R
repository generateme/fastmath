#!/usr/bin/env Rscript
# Generator of reference values for `fastmath.interpolation.kriging` from the R package gstat.
#
# Tool: R 4.5 with gstat 2.1-6 and sp 2.2-1 (`krige`, `vgm`, `variogramLine`); values are printed with
# 17 significant digits (round trip doubles).
#
# Domain:
#   * :cases - ordinary kriging (formula z ~ 1) and universal kriging with a linear drift (z ~ x + y, in 1-d
#     z ~ x) with the explicit semivariogram models Sph (spherical), Exp (exponential), Gau (gaussian), each with
#     a nugget, psill and range (`vgm(psill, model, range, nugget)`, the same parametrisation as
#     `fastmath.kernel.variogram/{spherical,exponential,gaussian}` with {:nugget :psill :range}), on a 2-d data set
#     (25 random nodes) and on a 1-d data set (15 nodes, points on the line y = 0); 30 and 25 random queries
#     inside the data region (never at a node). Additionally the gaussian ordinary case with the measurement error
#     `Err` (gstat measurement error VARIANCE, added to the diagonal of the covariance: `fastmath`'s `:error` is the
#     standard deviation, so Err = error^2).
#   * :error-sequence - ordinary kriging with a measurement error per node (a vector `error`): gstat has only a
#     scalar `Err`, so the system is solved here with R's `solve` and the gstat semivariogram `variogramLine`
#     (gamma matrix minus the diagonal of the error variances, the unbiasedness row/column, the same system as
#     `krige`); the script asserts that this solution equals `krige` for a constant error before it is written.
# Why: the audit measured the kriging values (3.6e-15 for the same variogram) but the former test suite had no
# reference data; `:error` was documented as a nugget although it is a standard deviation.
#
# Produces `test/resources/interpolation/kriging_reference.edn` (a map with the keys :cases and :error-sequence).
# The data sets use `set.seed` and `runif`/`rnorm` of the default Mersenne-Twister generator, stable between R versions.
#
# Run from anywhere:
#     Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_kriging_reference.R

suppressPackageStartupMessages({
  library(gstat)
  library(sp)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
out_file <- file.path(root, "test", "resources", "interpolation", "kriging_reference.edn")

fmt <- function(x) {
  s <- sprintf("%.17g", x)
  ifelse(grepl("[.eEni]", s), s, paste0(s, ".0"))
}
vec <- function(x) paste0("[", paste(fmt(x), collapse = " "), "]")
mat <- function(m) paste0("[", paste(apply(m, 1, vec), collapse = " "), "]")

# ---------------------------------------------------------------- data
set.seed(21)
n2 <- 25
x2 <- runif(n2, 0, 10); y2 <- runif(n2, 0, 10)
z2 <- sin(x2 / 3) + cos(y2 / 4) + 0.05 * x2 + rnorm(n2, 0, 0.1)
q2 <- cbind(runif(30, 1, 9), runif(30, 1, 9))

set.seed(22)
n1 <- 15
x1 <- sort(runif(n1, 0, 10))
z1 <- sin(x1 / 2) + 0.1 * x1 + rnorm(n1, 0, 0.1)
q1 <- sort(runif(25, 0.5, 9.5))

datasets <- list(
  d2 = list(dim = 2, coords = cbind(x2, y2), z = z2, queries = q2),
  d1 = list(dim = 1, coords = cbind(x1, 0 * x1), z = z1, queries = cbind(q1, 0 * q1)))

models <- list(
  spherical = list(gstat = "Sph", nugget = 0.2, psill = 1.2, range = 6.0),
  exponential = list(gstat = "Exp", nugget = 0.05, psill = 0.8, range = 2.5),
  gaussian = list(gstat = "Gau", nugget = 0.1, psill = 1.5, range = 4.0))

krige_values <- function(d, drift, model, err = 0.0) {
  df <- data.frame(x = d$coords[, 1], y = d$coords[, 2], z = d$z)
  coordinates(df) <- ~x + y
  nd <- data.frame(x = d$queries[, 1], y = d$queries[, 2])
  coordinates(nd) <- ~x + y
  vm <- vgm(psill = model$psill, model = model$gstat, range = model$range, nugget = model$nugget, Err = err)
  frm <- if (drift == "ordinary") z ~ 1 else if (d$dim == 1) z ~ x else z ~ x + y
  krige(frm, df, nd, model = vm, nmax = Inf, debug.level = 0)$var1.pred
}

case_lines <- c()
add_case <- function(name, d, drift, model_name, err) {
  model <- models[[model_name]]
  vals <- krige_values(d, drift, model, err)
  xs <- if (d$dim == 1) vec(d$coords[, 1]) else mat(d$coords)
  qs <- if (d$dim == 1) vec(d$queries[, 1]) else mat(d$queries)
  case_lines <<- c(case_lines, sprintf(
    "  {:name :%s :dim %d :model :%s :gstat \"%s\" :nugget %s :psill %s :range %s :drift :%s :err %s :xs %s :ys %s :queries %s :values %s}",
    name, d$dim, model_name, model$gstat, fmt(model$nugget), fmt(model$psill), fmt(model$range), drift,
    fmt(err), xs, vec(d$z), qs, vec(vals)))
}

for (dn in names(datasets)) {
  for (mn in names(models)) {
    for (drift in c("ordinary", "linear")) {
      add_case(paste(dn, mn, drift, sep = "-"), datasets[[dn]], drift, mn, 0.0)
    }
  }
}
add_case("d2-gaussian-ordinary-err", datasets$d2, "ordinary", "gaussian", 0.04)
add_case("d1-spherical-ordinary-err", datasets$d1, "ordinary", "spherical", 0.09)

# ---------------------------------------------------------------- error per node
manual_values <- function(d, model, errs) {
  n <- nrow(d$coords)
  gam <- function(h) variogramLine(vgm(psill = model$psill, model = model$gstat, range = model$range,
                                       nugget = model$nugget), dist_vector = h)$gamma
  D <- as.matrix(dist(d$coords))
  G <- matrix(gam(as.vector(D)), n, n)
  diag(G) <- 0.0
  G <- G - diag(errs^2)
  A <- rbind(cbind(G, 1), c(rep(1, n), 0))
  w <- solve(A, c(d$z, 0))
  sapply(seq_len(nrow(d$queries)), function(i) {
    dq <- sqrt(rowSums((d$coords - matrix(d$queries[i, ], n, 2, byrow = TRUE))^2))
    sum(w[1:n] * gam(dq)) + w[n + 1]
  })
}

# sanity: the manual system equals gstat for a constant error (Err = error^2)
chk <- manual_values(datasets$d2, models$gaussian, rep(0.2, n2))
stopifnot(max(abs(chk - krige_values(datasets$d2, "ordinary", models$gaussian, 0.04))) < 1e-9)

set.seed(23)
errs2 <- runif(n2, 0.05, 0.4)
errs1 <- runif(n1, 0.02, 0.3)
err_lines <- c(
  sprintf("  {:name :d2-gaussian-ordinary :model :gaussian :nugget %s :psill %s :range %s :errors %s :xs %s :ys %s :queries %s :values %s}",
          fmt(models$gaussian$nugget), fmt(models$gaussian$psill), fmt(models$gaussian$range), vec(errs2),
          mat(datasets$d2$coords), vec(datasets$d2$z), mat(datasets$d2$queries),
          vec(manual_values(datasets$d2, models$gaussian, errs2))),
  sprintf("  {:name :d1-spherical-ordinary :model :spherical :nugget %s :psill %s :range %s :errors %s :xs %s :ys %s :queries %s :values %s}",
          fmt(models$spherical$nugget), fmt(models$spherical$psill), fmt(models$spherical$range), vec(errs1),
          vec(datasets$d1$coords[, 1]), vec(datasets$d1$z), vec(datasets$d1$queries[, 1]),
          vec(manual_values(datasets$d1, models$spherical, errs1))))

lines <- c("{",
           " :cases [", case_lines, " ]",
           " :error-sequence [", err_lines, " ]",
           "}")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_file)
cat("gstat", as.character(packageVersion("gstat")), "wrote", length(case_lines), "cases and",
    length(err_lines), "error sequence cases to", out_file, "\n")
