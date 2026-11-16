# Reference data for the hierarchical clustering tests in test/fastmath/ml/clustering_test.clj
#
# Tools: R 4.5.2, stats::hclust, stats::cutree, stats::cophenetic, cluster 2.1.8.2 (diana).
# Run: Rscript utils/fastmath/dev/generate_hclust_reference.R (from the fastmath directory)
# Output: test/resources/clustering/hclust_reference.edn (regeneration is byte-identical:
# fixed data, %.17g printing).
#
# Why: fastmath.ml.clustering/hclust ports R hclust (:generic algorithm), cutree and diana
# (mapped as as.hclust), so R is the primary oracle for merges, heights, orders, labels,
# cophenetic distances and the cophenetic correlation.
# - :centroid and :median follow SciPy (squared distances): the reference is R on d^2, with
#   square roots of heights and cophenetic distances (equal to SciPy, checked in the SciPy
#   generator);
# - :weighted is R "mcquitty", :ward is "ward.D2", :ward-d is "ward.D";
# - tie sets (collinear, square, dup) check the R tie-breaking of the generic algorithm;
# - cut by height: at every height and between heights, only for monotone trees;
# - diana: merge rows are canonicalised (single observation first, two observations by
#   position, two clusters by step), as fastmath does; partitions for every k.
# The datasets are formulas repeated in generate_hclust_reference.py and .jl.

suppressMessages(library(cluster))

out <- file.path("test", "resources", "clustering", "hclust_reference.edn")

# ---- data (same formulas in the Python and Julia generators)
x8 <- matrix(c(-0.626, 0.576, 0.184, -0.305, -0.836, 1.512, 1.595, 0.39,
               0.33, -0.621, -0.82, -2.215, 0.487, 1.125, 0.738, -0.045), ncol = 2, byrow = TRUE)
x30 <- t(sapply(1:30, function(i) sapply(1:3, function(j) round(sin(i * 1.3 + j * 2.1) * 3 + cos(i * j * 0.7), 3))))
collinear <- matrix(0:5, ncol = 1)
square <- matrix(c(0, 0, 0, 1, 1, 0, 1, 1), ncol = 2, byrow = TRUE)
dup <- rbind(x8[1:5, ], x8[1, ], x8[3, ])
triangle <- matrix(c(0, 0, 1, 0, 0.51, 0.85), ncol = 2, byrow = TRUE)
pair <- matrix(c(0, 0, 3, 4), ncol = 2, byrow = TRUE)

stopifnot(!anyDuplicated(dist(x8)), !anyDuplicated(dist(x30)), !anyDuplicated(dist(x8, "manhattan")))

# ---- EDN printing
num <- function(x) ifelse(is.na(x), "##NaN", sprintf("%.17g", x))
int <- function(x) sprintf("%d", as.integer(x))
vec <- function(xs, f = num) paste0("[", paste(f(xs), collapse = " "), "]")
mat <- function(m, f = num) paste0("[", paste(apply(m, 1, function(r) vec(r, f)), collapse = " "), "]")
kv <- function(...) {
  args <- list(...)
  paste0("{", paste(paste0(":", names(args), " ", unlist(args)), collapse = "\n "), "}")
}

methods <- c(single = "single", complete = "complete", average = "average", weighted = "mcquitty",
             centroid = "centroid", median = "median", ward = "ward.D2", "ward-d" = "ward.D")
squared <- c("centroid", "median")

canonical <- function(m) t(apply(m, 1, function(r) {
  a <- r[1]; b <- r[2]
  if (a < 0 && b < 0) { if (a > b) c(a, b) else c(b, a) }
  else if (a < 0) c(a, b) else if (b < 0) c(b, a) else c(min(a, b), max(a, b))
}))

cuts <- function(tree, n) {
  ks <- t(sapply(1:n, function(k) cutree(tree, k)))
  res <- list(k = mat(ks - 1, int))
  if (!is.unsorted(tree$height)) {
    hs <- sort(unique(c(0, tree$height, (head(tree$height, -1) + tail(tree$height, -1)) / 2, max(tree$height) + 1)))
    res$h <- vec(hs)
    res[["h-labels"]] <- mat(t(sapply(hs, function(h) cutree(tree, h = h))) - 1, int)
  }
  res
}

agglomerative <- function(d, n, method) {
  sq <- method %in% squared
  tree <- hclust(if (sq) d^2 else d, methods[[method]])
  coph <- as.matrix(cophenetic(tree))
  if (sq) {
    tree$height <- sqrt(tree$height)
    coph <- sqrt(coph)
  }
  cc <- if (n >= 3) cor(as.vector(as.dist(coph)), as.vector(d)) else NA
  do.call(kv, c(list(merge = mat(tree$merge, int), height = vec(tree$height),
                     order = vec(tree$order - 1, int), cophenetic = mat(coph),
                     "cophenetic-correlation" = num(cc)),
                cuts(tree, n)))
}

divisive <- function(d, n) {
  dv <- diana(d)
  tree <- as.hclust(dv)
  kv(merge = mat(canonical(tree$merge), int), height = vec(tree$height),
     order = vec(dv$order - 1, int), dc = num(dv$dc),
     cophenetic = mat(as.matrix(cophenetic(tree))),
     "cophenetic-correlation" = num(cor(as.vector(cophenetic(tree)), as.vector(d))),
     k = mat(t(sapply(1:n, function(k) cutree(tree, k))) - 1, int))
}

dataset <- function(x, metric = "euclidean", precomputed = FALSE, diana = TRUE) {
  d <- dist(x, metric)
  n <- attr(d, "Size")
  res <- lapply(names(methods), function(m) agglomerative(d, n, m))
  names(res) <- names(methods)
  entries <- c(list(data = mat(if (precomputed) as.matrix(d) else x),
                    distance = if (precomputed) ":precomputed" else paste0(":", metric)),
               res)
  if (diana) entries$diana <- divisive(d, n)
  do.call(kv, entries)
}

sets <- list(x8 = dataset(x8),
             x30 = dataset(x30),
             collinear = dataset(collinear),
             square = dataset(square),
             dup = dataset(dup),
             manhattan = dataset(x8, "manhattan", precomputed = TRUE),
             triangle = dataset(triangle),
             pair = dataset(pair))

writeLines(c(";; generated by utils/fastmath/dev/generate_hclust_reference.R, do not edit",
             do.call(kv, sets)), out)
