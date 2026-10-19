#!/usr/bin/env Rscript
# Generator of reference values for the vector kernels of `fastmath.kernel` from the R package kernlab.
#
# Tool: R 4.5 with kernlab 0.9.33 (`rbfdot`, `laplacedot`, `polydot`, `tanhdot`, `vanilladot`, `anovadot`,
# `splinedot`, `besseldot`); values are printed with 17 significant digits (round trip doubles).
#
# Domain: 8 deterministic pairs of points per kernel (`set.seed(20260514)`, `runif` of the default
# Mersenne-Twister generator, stable between R versions), a pair of identical points included. The
# mappings of the parameters:
#   rbfdot(sigma)                  exp(-sigma |x - y|^2)         gaussian :sigma 1/sqrt(2 sigma)
#   laplacedot(sigma)              exp(-sigma |x - y|)           laplacian :sigma 1/sigma
#   polydot(degree, scale, offset) (scale <x,y> + offset)^degree polynomial :p :alpha :c
#   tanhdot(scale, offset)         tanh(scale <x,y> + offset)    hyperbolic-tangent :alpha :c
#   anovadot(sigma, degree 1)      sum_i exp(-sigma (x_i - y_i)^2)  anova :sigma :k 1 :d 1
#   splinedot()                    product of the spline expression (non-negative input)  spline
#   besseldot(sigma, order, degree) (J_order(b) b^-order / lim)^degree, b = sigma |x - y|  bessel2
#                                  kernlab switches to the limit 1 at b < 1e-4 and fastmath at b < 1e-5, so the
#                                  pairs are chosen with b >= 0.02 (identical points: both give 1).
# Why: pin the parametrisation of the kernels which come from kernlab (spline, anova, bessel2) and of the
# kernlab kernels with a mapping (rbf, laplace, poly, tanh). Note that the `anovadot` degree acts on the sum,
# fastmath's `:d` acts on every term, hence only `degree = 1` is stored.
#
# Produces `test/resources/kernel/kernel_vector_kernlab_reference.edn`, a map `{:cases [...]}`; a case is
# `{:kernel :params :oracle :x [[...]] :y [[...]] :expected [...]}` (the format of
# `kernel_vector_reference.edn`).
#
# Run from anywhere:
#     Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_kernel_vector_kernlab_reference.R

suppressPackageStartupMessages({
  library(kernlab)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- normalizePath(file.path(dirname(script), "..", "..", ".."))
out_file <- file.path(root, "test", "resources", "kernel", "kernel_vector_kernlab_reference.edn")

fmt <- function(x) {
  s <- sprintf("%.17g", x)
  ifelse(grepl("[.eEni]", s), s, paste0(s, ".0"))
}
vec <- function(x) paste0("[", paste(fmt(x), collapse = " "), "]")
mat <- function(m) paste0("[", paste(apply(m, 1, vec), collapse = " "), "]")

set.seed(20260514)
pairs <- function(d, lo = -1.5, hi = 1.5) {
  x <- matrix(runif(7 * d, lo, hi), ncol = d)
  y <- matrix(runif(7 * d, lo, hi), ncol = d)
  list(x = rbind(x, x[1, ]), y = rbind(y, x[1, ]))
}

cases <- character(0)
add <- function(kernel, params, p, kern, oracle) {
  expected <- sapply(seq_len(nrow(p$x)), function(i) kern(p$x[i, ], p$y[i, ]))
  cases <<- c(cases, sprintf(" {:kernel :%s :params %s :oracle \"%s\" :x %s :y %s :expected %s}",
                             kernel, params, oracle, mat(p$x), mat(p$y), vec(expected)))
}

p <- pairs(3)
for (sigma in c(0.3, 1.2)) {
  add("gaussian", sprintf("{:sigma %s}", fmt(1 / sqrt(2 * sigma))), p, rbfdot(sigma = sigma), "kernlab rbfdot")
  add("laplacian", sprintf("{:sigma %s}", fmt(1 / sigma)), p, laplacedot(sigma = sigma), "kernlab laplacedot")
}
add("linear", "{}", p, vanilladot(), "kernlab vanilladot")
for (cfg in list(c(2, 1, 0), c(3, 0.5, 1), c(2, 2, 0.5))) {
  add("polynomial", sprintf("{:p %s :alpha %s :c %s}", fmt(cfg[1]), fmt(cfg[2]), fmt(cfg[3])), p,
      polydot(degree = cfg[1], scale = cfg[2], offset = cfg[3]), "kernlab polydot")
}
q <- pairs(3, -0.8, 0.8)
for (cfg in list(c(1, 0), c(0.5, 0.3))) {
  add("hyperbolic-tangent", sprintf("{:alpha %s :c %s}", fmt(cfg[1]), fmt(cfg[2])), q,
      tanhdot(scale = cfg[1], offset = cfg[2]), "kernlab tanhdot")
}
a <- pairs(3)
for (sigma in c(0.5, 1.3)) {
  add("anova", sprintf("{:sigma %s :k 1.0 :d 1.0}", fmt(sigma)), a, anovadot(sigma = sigma, degree = 1), "kernlab anovadot")
}
s <- pairs(3, 0, 2)
add("spline", "{}", s, splinedot(), "kernlab splinedot")
b <- pairs(2, -1.5, 1.5)
for (cfg in list(c(1, 1, 1), c(0.7, 1, 2), c(1.5, 2, 1), c(2, 0.5, 1))) {
  add("bessel2", sprintf("{:sigma %s :order %s :degree %s}", fmt(cfg[1]), fmt(cfg[2]), fmt(cfg[3])), b,
      besseldot(sigma = cfg[1], order = cfg[2], degree = cfg[3]), "kernlab besseldot")
}

writeLines(c("{:cases [", cases, "]}"), out_file)
cat("kernlab", as.character(packageVersion("kernlab")), "wrote", length(cases), "cases to", out_file, "\n")
