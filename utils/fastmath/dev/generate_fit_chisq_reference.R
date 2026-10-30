# Reference values for the chi-square statistic of fastmath.fitdistr/fit for discrete distributions.
# Tools: R fitdistrplus::fitdist and fitdistrplus::gofstat (automatic cells, and the cells given by chisqbreaks).
# Output per case: the data, the fitted parameters (fastmath names), the automatic breaks of gofstat with its
# chisq and chisqdf, and the explicit breaks with its chisq and chisqdf.
# Datasets are stored in the output so the Clojure tests do not depend on the R random generator.
# Run: Rscript utils/fastmath/dev/generate_fit_chisq_reference.R
# Output: test/resources/fitdistr/fit_chisq_reference.edn
suppressMessages(library(fitdistrplus))

fmt <- function(v) sprintf("%.17g", v)
vec <- function(v) paste0("[", paste(fmt(v), collapse = " "), "]")

case <- function(key, x, distr, start, params, breaks) {
  f <- if (is.null(start)) fitdist(x, distr, discrete = TRUE) else fitdist(x, distr, start = start, discrete = TRUE)
  a <- gofstat(f)
  b <- gofstat(f, chisqbreaks = breaks)
  paste0("{:distribution ", key, " :params ", params(f$estimate), "\n :data ", vec(x),
         "\n :auto-breaks ", vec(a$chisqbreaks), " :auto-chisq ", fmt(a$chisq), " :auto-df ", a$chisqdf,
         "\n :breaks ", vec(breaks), " :chisq ", fmt(b$chisq), " :df ", b$chisqdf, "}")
}

set.seed(20261102)
out <- c(
  case(":poisson", rpois(400, 4.5), "pois", NULL,
       function(e) sprintf("{:p %.17g}", e["lambda"]), c(2, 4, 6, 8)),
  case(":geometric", rgeom(500, 0.25), "geom", NULL,
       function(e) sprintf("{:p %.17g}", e["prob"]), c(1, 3, 6, 10, 15)),
  case(":negative-binomial", rnbinom(600, size = 3.5, prob = 0.4), "nbinom", NULL,
       function(e) sprintf("{:r %.17g :p %.17g}", e["size"], e["size"] / (e["size"] + e["mu"])), c(1, 3, 5, 8, 12)),
  case(":poisson", rpois(2000, 12), "pois", NULL,
       function(e) sprintf("{:p %.17g}", e["lambda"]), c(8, 11, 14, 17)))

dir.create("test/resources/fitdistr", showWarnings = FALSE, recursive = TRUE)
writeLines(c("[", paste(out, collapse = "\n"), "]"), "test/resources/fitdistr/fit_chisq_reference.edn")
