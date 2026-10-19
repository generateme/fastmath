#!/usr/bin/env Rscript
# Reference values of the Wendland functions from the R package fields for `fastmath.kernel.rbf/wendland`.
#
# Tool: R with fields (`Wendland(d, theta = 1, dimension, k)`, normalised to 1 at d = 0, support radius theta = 1);
# every double is printed with 17 significant digits (round trip).
#
# Domain: dimension 1..8 (fields uses l = floor(dimension / 2) + k + 1, Fasshauer Definition 4.4), k = 1..4
# (fields supports k 2, 3, 4 and k = 1 - the closed forms only), 8 distances: 0, a tiny one, interior points,
# the last point before the support radius, the radius and a point beyond.
#
# Output: an EDN vector on stdout, one map `{:s :k :r [...] :v [...]}` per line. It is embedded verbatim under
# the key :fields-wendland by `generate_rbf_kernel_reference.py`, which owns the file
# `test/resources/kernel/rbf_kernel_reference.edn`.
#
# Run (normally through the python generator, but it works alone):
#     Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_rbf_kernel_reference.R

suppressPackageStartupMessages(library(fields))

fmt <- function(x) {
  s <- sprintf("%.17g", x)
  if (!grepl("[.eE]", s)) s <- paste0(s, ".0")
  s
}

rs <- c(0, 0.001, 0.1, 0.25, 0.5, 0.9, 0.999, 1, 1.5)

cat("[\n")
rows <- character(0)
for (s in 1:8) {
  for (k in 1:4) {
    v <- Wendland(rs, theta = 1, dimension = s, k = k)
    rows <- c(rows, sprintf(" {:s %d :k %d :r [%s] :v [%s]}", s, k,
                            paste(sapply(rs, fmt), collapse = " "),
                            paste(sapply(v, fmt), collapse = " ")))
  }
}
cat(paste(rows, collapse = "\n"), "\n", sep = "")
cat("]\n")
