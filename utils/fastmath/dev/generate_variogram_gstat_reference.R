#!/usr/bin/env Rscript
# Generator of reference values for `fastmath.kernel.variogram` from the R package gstat.
#
# Tool: R 4.5 with gstat 2.1-6 and sp 2.2-1 (`variogramLine`, `variogram`); values are printed with
# 17 significant digits (round trip doubles).
#
# Domain:
#   * :models - the semivariogram models which gstat has with the same parametrisation as fastmath:
#     Sph (spherical), Exp (exponential), Gau (gaussian), Pen (pentaspherical), Cir (circular),
#     Lin (linear), Wav (hole = normalised sinc); two parameter sets, 17 lags per set (a tiny lag, fractions
#     of the range, the range itself, just beyond, far beyond). Lag zero is not stored (gstat gives 0).
#   * :empirical - the classical (Matheron) estimator of `variogram` on 6 deterministic data sets (2-d
#     random, 2-d clustered, 3-d, jittered grid, 1-d line, 2-d with an outlier), with explicit `boundaries`.
#     `fastmath.kernel.variogram/empirical` bins the pairs with the distance within the cutoff on
#     (0, dmax] with `bins` equal parts (dmax is the largest distance within the cutoff, left-open
#     right-closed), so the boundaries here are those edges. Why: the default bins of gstat (cutoff / 15)
#     are different, the comparison is exact only for the same edges.
#
# Produces `test/resources/kernel/variogram_gstat_reference.edn` (a map with the keys :models and :empirical).
# The data sets use `set.seed` and `runif`/`rnorm` of the default Mersenne-Twister generator, stable between R versions.
#
# Run from anywhere:
#     Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_variogram_gstat_reference.R

suppressPackageStartupMessages({
  library(gstat)
  library(sp)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
out_file <- file.path(root, "test", "resources", "kernel", "variogram_gstat_reference.edn")

fmt <- function(x) {
  s <- sprintf("%.17g", x)
  ifelse(grepl("[.eEni]", s), s, paste0(s, ".0"))
}
vec <- function(x) paste0("[", paste(fmt(x), collapse = " "), "]")
ivec <- function(x) paste0("[", paste(sprintf("%d", as.integer(x)), collapse = " "), "]")

# ---------------------------------------------------------------- models
param_sets <- list(c(nugget = 0.0, psill = 1.0, range = 1.0),
                   c(nugget = 0.3, psill = 2.5, range = 7.5))
rel_lags <- c(1e-9, 1e-3, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99, 1.0, 1.0000001, 1.2, 1.5, 3.0, 10.0)
models <- c(spherical = "Sph", exponential = "Exp", gaussian = "Gau", pentaspherical = "Pen",
            circular = "Cir", linear = "Lin", hole = "Wav")

model_lines <- c()
for (name in names(models)) {
  for (p in param_sets) {
    lags <- rel_lags * p[["range"]]
    vgm_model <- vgm(psill = p[["psill"]], model = models[[name]], range = p[["range"]], nugget = p[["nugget"]])
    g <- variogramLine(vgm_model, dist_vector = lags)$gamma
    model_lines <- c(model_lines, sprintf(
      "  {:model :%s :gstat \"%s\" :nugget %s :psill %s :range %s :lag %s :gamma %s}",
      name, models[[name]], fmt(p[["nugget"]]), fmt(p[["psill"]]), fmt(p[["range"]]), vec(lags), vec(g)))
  }
}

# ---------------------------------------------------------------- empirical
datasets <- list()
set.seed(11)
n <- 40
x <- runif(n, 0, 10); y <- runif(n, 0, 10)
datasets$d1 <- list(coords = cbind(x, y), z = sin(x / 2) + rnorm(n, 0, 0.3), bins = 5, frac = 0.5)

set.seed(12)
cx <- rep(c(2, 8, 5), each = 15); cy <- rep(c(2, 3, 8), each = 15)
x <- cx + rnorm(45, 0, 1); y <- cy + rnorm(45, 0, 1)
datasets$d2 <- list(coords = cbind(x, y), z = 0.3 * x - 0.2 * y + rnorm(45, 0, 0.5), bins = 8, frac = 0.6)

set.seed(13)
x <- runif(35, 0, 5); y <- runif(35, 0, 5); w <- runif(35, 0, 5)
datasets$d3 <- list(coords = cbind(x, y, w), z = cos(x) + sin(y) * w / 5 + rnorm(35, 0, 0.2), bins = 7, frac = 0.7)

set.seed(14)
gx <- rep(seq(0, 12, length.out = 8), times = 8) + runif(64, -0.3, 0.3)
gy <- rep(seq(0, 12, length.out = 8), each = 8) + runif(64, -0.3, 0.3)
datasets$d4 <- list(coords = cbind(gx, gy), z = sin(gx / 3) * cos(gy / 4) + rnorm(64, 0, 0.1), bins = 10, frac = 0.45)

set.seed(15)
x <- sort(runif(60, 0, 100))
datasets$d5 <- list(coords = cbind(x, 0 * x), z = cumsum(rnorm(60, 0, 0.5)), bins = 12, frac = 0.5)

set.seed(16)
x <- runif(30, 0, 8); y <- runif(30, 0, 8)
z <- x / 4 + rnorm(30, 0, 0.2); z[7] <- z[7] + 6
datasets$d6 <- list(coords = cbind(x, y), z = z, bins = 6, frac = 0.8)

emp_lines <- c()
for (name in names(datasets)) {
  d <- datasets[[name]]
  m <- nrow(d$coords)
  dm <- as.matrix(dist(d$coords))
  all_d <- dm[upper.tri(dm)]
  cutoff <- d$frac * max(all_d)
  kept <- all_d[all_d > 0 & all_d <= cutoff]
  dmax <- max(kept)
  edges <- seq(0, dmax, length.out = d$bins + 1)
  edges[length(edges)] <- dmax * (1 + 1e-12)
  df <- as.data.frame(d$coords)
  names(df) <- c("x", "y", "w")[seq_len(ncol(d$coords))]
  df$z <- d$z
  coordinates(df) <- as.formula(paste("~", paste(names(df)[seq_len(ncol(d$coords))], collapse = "+")))
  v <- variogram(z ~ 1, df, boundaries = edges, cutoff = max(edges))
  stopifnot(sum(v$np) == length(kept))
  emp_lines <- c(emp_lines, sprintf(
    "  {:name :%s :xs [%s] :ys %s :cutoff %s :bins %d :n %s :h %s :gamma %s}",
    name,
    paste(apply(d$coords, 1, function(r) vec(r)), collapse = " "),
    vec(d$z), fmt(cutoff), d$bins, ivec(v$np), vec(v$dist), vec(v$gamma)))
}

lines <- c("{",
           " :models [", model_lines, " ]",
           " :empirical [", emp_lines, " ]",
           "}")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_file)
cat("gstat", as.character(packageVersion("gstat")), "wrote", length(model_lines), "model cases and",
    length(emp_lines), "data sets to", out_file, "\n")
