# R oracle of the time-series reference of `fastmath.random` (Random Audit, work package 9.5).
#
# Called by generate_random_series_reference.py (do not run by hand):
#   Rscript generate_random_series_reference.R <signal.txt> <out.txt>
# <signal.txt>: the innovations e_0..e_{N-1}, one number per line (written by the python script).
# <out.txt>:    the echo of the signal (to check the text round trip) and, per model, the output of
#               stats::arima.sim with explicit innovations, start innovations and burn-in:
#                 ns = 2 (p + q), n = N - ns
#                 arima.sim(model = list(ar = phi, ma = theta), n = n, innov = e[(ns + 1):N], n.start = ns, start.innov = e[1:ns])
#               R convention: ar = +phi (Y_t = phi Y_{t-1} + e_t), ma = +theta (Y_t = e_t + theta e_{t-1}), the same as fastmath.
#               The first R output is the process value at signal index ns, as the first output of the library functions.
# Line format of <out.txt>:  "E v1 v2 ..." for the echo, "M name ar|... ; ma|... ; values" for models (coefficients with %.17g).

args <- commandArgs(trailingOnly = TRUE)
e <- scan(args[1], quiet = TRUE)
N <- length(e)

models <- list(
  list("ma-1", numeric(0), c(0.5)),
  list("ma-2", numeric(0), c(0.5, -0.3)),
  list("ma-3", numeric(0), c(0.4, -0.3, 0.2)),
  list("ma-4", numeric(0), c(0.5, -0.3, 0.2, 0.1)),
  list("ar-1", c(0.6), numeric(0)),
  list("ar-2", c(0.5, -0.3), numeric(0)),
  list("ar-3", c(0.2, 0.3, 0.1), numeric(0)),
  list("ar-4", c(0.3, 0.2, 0.1, 0.05), numeric(0)),
  list("arma-11", c(0.7), c(0.4)),
  list("arma-22", c(0.5, -0.2), c(0.3, 0.1)),
  list("arma-33", c(0.2, 0.3, 0.1), c(0.4, -0.3, 0.2)),
  list("arma-44", c(0.3, 0.2, 0.1, 0.05), c(0.5, -0.3, 0.2, 0.1)),
  list("arma-21", c(0.5, -0.3), c(0.4)),
  list("arma-12", c(0.5), c(0.4, 0.2))
)

fmt <- function(x) paste(sprintf("%.17g", x), collapse = " ")

con <- file(args[2], "w")
writeLines(paste("E", fmt(e)), con)
writeLines(paste("V", R.version.string), con)
for (m in models) {
  name <- m[[1]]; ar <- m[[2]]; ma <- m[[3]]
  ns <- 2 * (length(ar) + length(ma))
  n <- N - ns
  model <- list()
  if (length(ar) > 0) model$ar <- ar
  if (length(ma) > 0) model$ma <- ma
  x <- as.numeric(arima.sim(model = model, n = n, innov = e[(ns + 1):N], n.start = ns, start.innov = e[1:ns]))
  stopifnot(length(x) == n)
  writeLines(paste("M", name, fmt(ar), ";", fmt(ma), ";", fmt(x)), con)
}
close(con)
