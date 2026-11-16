# Reference values of fastmath.kernel-test 2d KDE (test/fastmath/kernel/density_2d_test.clj): MASS::kde2d, default bandwidth.
# Run: Rscript utils/fastmath/dev/kde2d_reference.R
library(MASS)
x <- c(1.2, 2.5, 3.1, 0.4, 2.2, 4.0, 3.3, 1.9, 2.8, 0.9)
y <- c(5.0, 3.2, 4.4, 6.1, 2.7, 3.9, 5.5, 4.8, 3.0, 5.9)
k <- kde2d(x, y, n = 3, lims = c(1, 3, 3, 5))
cat(sprintf("%.15g", bandwidth.nrd(x)/4), sprintf("%.15g", bandwidth.nrd(y)/4), "\n")
cat(sprintf("%.15g", as.vector(k$z)), "\n")
