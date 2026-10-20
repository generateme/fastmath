# Generator of the R part of test/resources/kernel/bandwidth_reference.edn  (step 1 of 2).
#
# Regeneration (byte-identical; run in this order, from anywhere):
#   Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bandwidth_reference.R
#   cd /home/ts/penv && uv run python \
#       /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bandwidth_reference.py
# Step 2 (the Python generator) reads the file written here (datasets + R outputs), adds the independent
# Python values and the cross-validation battery, and rewrites the whole file.  Run alone, this script
# produces a valid file with the R part only.
#
# Tools: R 4.5.2 `stats` (bw.nrd, bw.nrd0, bw.SJ), RNGkind("Mersenne-Twister", "Inversion", "Rejection"),
#        set.seed(20261020) once, datasets drawn in the fixed order below.  Precision: double; numbers are
#        written with the shortest decimal string that reads back to the same double.
#
# Domain: 32 datasets (n = 1 ... 1500): normal samples (n 5 ... 1000), n = 1, 2, 3, constants (sd = 0), IQR = 0
# with sd > 0, integer data and heavy ties, three distinct values, sd 1000, offset 1e6 (also with sd 1000),
# tiny sd 1e-6, lognormal (skewed), bimodal, trimodal mixture, t(3), exponential, uniform.  n > 500
# exercises the binned pair-count path of bw.SJ (n > nb/2).
#
# For every dataset: sample mean, sd, IQR (type 7), the `bw.SJ` scale min(sd, IQR/1.349) and
# hmax = 1.144 scale n^-0.2, and the outputs of
#   bw.nrd(x)                            1.06 min(sd, IQR/1.34) n^-0.2 (type 7 quantiles; no fallback)
#   bw.nrd0(x)                           0.9 lo n^-0.2 with the chain lo = min(sd, IQR/1.34), else sd, else |x[1]|, else 1
#   bw.SJ(x, method = "ste")             defaults: nb 1000, lower 0.1 hmax, upper hmax, tol 0.1 lower
#   bw.SJ(x, method = "dpi")
#   bw.SJ(x, method = "ste", tol = 1e-10 hmax)   the same equation with a tight root-finding tolerance
# An outcome is {:value v} or, when R signals an error or returns NA/NaN/Inf, {:error true :message "..."} (the
# message of the R condition, "NA" for a missing value).  `:warning` lists warnings (none were observed).
#
# Why: KB-07 (type 6 vs type 7 IQR), KB-08 (fallback chain), D6/D7 of the Kernel Density audit, WP 2.1 and 2.6.

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20261020)

OUT <- file.path("/home/ts/clojure/fastmath/test/resources/kernel/bandwidth_reference.edn")

num <- function(x) {
  x <- as.numeric(x)
  if (is.na(x)) return("##NaN")
  if (is.infinite(x)) return(if (x > 0) "##Inf" else "##-Inf")
  for (d in 1:17) {
    s <- sprintf(paste0("%.", d, "g"), x)
    if (as.numeric(s) == x) break
  }
  if (!grepl("[.eE]", s)) s <- paste0(s, ".0")
  s
}
vec <- function(x) paste0("[", paste(vapply(x, num, ""), collapse = " "), "]")
str1 <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")

outcome <- function(expr) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "bwerr"))
  s <- if (inherits(res, "bwerr")) {
    paste0("{:error true :message ", str1(res$msg), "}")
  } else if (is.na(res) || !is.finite(res)) {
    paste0("{:error true :message ", str1(if (is.na(res)) "NA" else "Inf"), "}")
  } else paste0("{:value ", num(res), "}")
  if (length(warns)) s <- sub("}$", paste0(" :warning ", str1(paste(unique(warns), collapse = "; ")), "}"), s)
  s
}

# ---------------------------------------------------------------- datasets (fixed order of draws)
ds <- list()
add <- function(name, x) ds[[length(ds) + 1L]] <<- list(name = name, x = as.numeric(x))
add("n1", 5)
add("n2", c(1, 4))
add("n2-same", c(2, 2))
add("n3", c(-1, 0, 1))
add("n5-norm", rnorm(5))
add("n10-norm", rnorm(10))
add("n20-norm", rnorm(20))
add("n33-norm", rnorm(33))
add("n100-norm", rnorm(100))
add("n1000-norm", rnorm(1000))
add("n501-norm", rnorm(501))
add("const20", rep(5, 20))
add("const0", rep(0, 20))
x <- numeric(100); x[1:10] <- rnorm(10, 0, 3)
add("iqr0-zeros", x)
add("iqr0-outliers", c(rep(7, 40), 1, 2, 30))
add("ties-int100", round(rnorm(100, 0, 2)))
add("poisson200", rpois(200, 3))
add("integers1000", sample(1:50, 1000, replace = TRUE))
add("three-values60", rep(c(0, 1, 5), each = 20))
add("sd1000-n100", 1000 * rnorm(100))
add("offset1e6-n50", 1e6 + rnorm(50))
add("offset1e6-sd1000-n200", 1e6 + 1000 * rnorm(200))
add("negative-mean-small-sd", rnorm(80, -50, 0.1))
add("tiny-sd-1e-6", 1e-6 * rnorm(40))
add("lognorm100", rlnorm(100))
add("lognorm500", rlnorm(500))
add("lognorm1500", rlnorm(1500))
u <- runif(200); add("bimodal200", ifelse(u < 0.5, rnorm(200, -3, 0.5), rnorm(200, 3, 0.5)))
u <- runif(500); add("trimix500", ifelse(u < 0.5, rnorm(500), ifelse(u < 0.8, rnorm(500, 5, 0.1), -4 + 0.3 * (rexp(500) - rexp(500)))))
add("t3-n300", rt(300, 3))
add("expo100", rexp(100))
add("unif50", runif(50))

# ---------------------------------------------------------------- results
lines <- character(0)
for (d in ds) {
  x <- d$x
  n <- length(x)
  sdx <- if (n > 1) sd(x) else NA_real_
  iqr7 <- if (n > 1) IQR(x) else NA_real_
  scale <- min(sdx, iqr7 / 1.349)
  hmax <- 1.144 * scale * n^(-1/5)
  tight <- if (is.finite(hmax) && hmax > 0) 1e-10 * hmax else NA_real_
  rec <- paste0("{:name ", str1(d$name), " :n ", n,
                " :mean ", num(mean(x)), " :sd ", num(sdx), " :iqr7 ", num(iqr7),
                " :sj-scale ", num(scale), " :hmax ", num(hmax),
                "\n   :nrd ", outcome(bw.nrd(x)),
                "\n   :nrd0 ", outcome(bw.nrd0(x)),
                "\n   :sj-ste ", outcome(bw.SJ(x, method = "ste")),
                "\n   :sj-dpi ", outcome(bw.SJ(x, method = "dpi")),
                "\n   :sj-ste-tight ", outcome(bw.SJ(x, method = "ste", tol = tight)),
                "\n   :data ", vec(x), "}")
  lines <- c(lines, rec)
}

hdr <- c(
  ";; bandwidth_reference.edn - generated by utils/fastmath/dev/generate_bandwidth_reference.R then .py (do not edit)",
  paste0(";; [R] ", R.version.string, ", stats::bw.nrd, bw.nrd0, bw.SJ; 32 datasets, set.seed(20261020), Mersenne-Twister/Inversion/Rejection."),
  ";; [R] case = {:name :n :mean :sd :iqr7 :sj-scale :hmax :nrd :nrd0 :sj-ste :sj-dpi :sj-ste-tight :data}; outcome = {:value v} or {:error true :message s}."
)
writeLines(c(hdr, "{", " :cases [", paste0("  ", lines), " ]", "}"), OUT)
cat("wrote", OUT, file.info(OUT)$size, "bytes,", length(ds), "datasets\n")
