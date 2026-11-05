# Reference values of the regularized incomplete beta function I_x(a, b) for large shapes (fastmath.special,
# upstream group 2): one large shape (the incomplete gamma limit region, both orientations), moderate-to-large shapes
# where the continued fraction is kept, and two large shapes (min(a, b) >= 1e7, where fastmath returns NaN; kept for
# documentation). x = mean + z * sd of Beta(a, b), z in {-5, -1, 0, 1, 5}.
#
# Tool: R pbeta (TOMS 708), checked against mpmath quadrature of the density at (1e10, 1e5) and (1e20, 1e10)
# (agreement 1e-15 and 7e-12). Numbers are written with 17 significant digits (the default 15 changes the inputs).
#
#   Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_special_gamma_beta_large_reference.R
#
# Output: test/resources/special/gamma_beta_large_reference.edn  {:beta-large {:x [..] :a [..] :b [..] :ref [..]}}

f <- function(v) sprintf("%.17g", v)
xs <- c(); as <- c(); bs <- c(); refs <- c()
add <- function(a, b) {
  m <- a / (a + b); s <- sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  for (z in c(-5, -1, 0, 1, 5)) {
    x <- m + z * s
    if (x > 0 && x < 1) {
      i <- pbeta(x, a, b)
      if (!is.na(i) && i >= 1e-300) {
        xs <<- c(xs, x); as <<- c(as, a); bs <<- c(bs, b); refs <<- c(refs, i)
      }
    }
  }
}
for (la in c(8, 10, 12, 15, 20, 100, 300)) for (b in c(0.5, 1, 10, 1e3, 1e5, 1e6)) {
  a <- 10^la
  if (b <= 1e-5 * a) { add(a, b); add(b, a) }
}
for (la in c(3, 4, 5, 6)) for (lr in c(0, 1, 2)) add(10^la, 10^(la + lr))
for (la in c(7, 8, 10)) add(10^la, 10^(la + 1))
dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
out <- file.path(dir, "..", "..", "..", "test", "resources", "special", "gamma_beta_large_reference.edn")
v <- function(name, x) paste0(":", name, " [", paste(f(x), collapse = " "), "]")
writeLines(paste0("{:beta-large {", v("x", xs), "\n  ", v("a", as), "\n  ", v("b", bs), "\n  ", v("ref", refs), "}}"), out)
