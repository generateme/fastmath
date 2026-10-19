(ns fastmath.kernel.vector
  "Vector kernels: functions of two vectors (or two numbers) returning a double.

  A vector kernel measures similarity (or, for the conditionally positive definite ones, a negated dissimilarity) of two points. Every constructor takes an optional map of parameters, returns a function `(fn [x y] ...)` and accepts one dimensional input as plain numbers as well as vectors. Constructors which use a distance take a `:distance` key, a function of two points returning a double (default: [[fastmath.vector/dist]], euclidean).

  Positive definiteness is stated for every kernel: positive definite (PD, every Gram matrix is positive semi-definite), conditionally positive definite (CPD, the Gram matrix is positive semi-definite on vectors with zero sum) or neither. Only PD kernels are valid covariance functions of a Gaussian process; CPD ones need a constant (or polynomial) term.

  Kernels: [[linear]], [[polynomial]], [[gaussian]], [[exponential]], [[laplacian]], [[anova]], [[hyperbolic-tangent]], [[rational-quadratic]], [[multiquadratic]], [[inverse-multiquadratic]], [[geometric]], [[wave]], [[periodic]], [[power]], [[log]], [[spline]], [[b-spline]], [[bessel]], [[bessel2]], [[cauchy]], [[chi-square]], [[chi-square2]], [[histogram]], [[generalized-histogram]], [[generalized-t-student]], [[dirichlet]], [[hellinger]], [[pearson]], [[hyperbolic-secant]], [[matern]]; [[rbf->kernel]] converts a radial basis function into a vector kernel.

  Literature: Genton, Classes of Kernels for Machine Learning: A Statistics Perspective, JMLR 2 (2001), http://www.jmlr.org/papers/volume2/genton01a/genton01a.pdf; http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/ ."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special :as special]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- invalid-option
  "Creates the `ex-info` for an invalid kernel option."
  [kernel option value expected]
  (ex-info (str "Kernel " kernel ": invalid " option ", expected " expected ", got " (pr-str value))
           {:type :invalid-input
            :kernel kernel
            :option option
            :value value
            :expected expected}))

(defn linear
  "Creates the linear kernel, the dot product of two vectors.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning `<x, y>` as a double. The kernel is positive definite. The values are not normalised: K(x, x) is the squared norm of `x`.

  See also [[polynomial]], [[hyperbolic-tangent]]."
  ([] (fn ^double [x y] (v/dot x y)))
  ([_] (fn ^double [x y] (v/dot x y))))

(defn polynomial
  "Creates the polynomial kernel, a power of the shifted and scaled dot product: `(c + alpha <x, y>)^p`.

  Parameters:

  - `:alpha` - dot product multiplier (default: `1.0`)
  - `:c` - shift (default: `0.0`)
  - `:p` - exponent (default: `2.0`)

  Returns a function `(fn [x y])` returning a double. Without arguments the kernel is `<x, y>^2`. The kernel is positive definite for a positive integer `:p`, `alpha >= 0` and `c >= 0`; a non-integer `:p` gives `NaN` for a negative base and no guarantee of positive definiteness. Corresponds to the polynomial kernel of scikit-learn with `gamma = alpha`, `coef0 = c` and `degree = p`.

  See also [[linear]], [[hyperbolic-tangent]]."
  ([] (fn ^double [x y] (m/sq (v/dot x y))))
  ([{:keys [^double alpha ^double c ^double p]
     :or {alpha 1.0 c 0.0 p 2.0}}]
   (fn ^double [x y] (m/pow (m/+ c (m/* alpha (v/dot x y))) p))))

(defn gaussian
  "Creates the Gaussian (squared exponential, RBF) kernel `exp(-d^2 / (2 sigma^2))`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - length scale (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance. Equal to the scikit-learn RBF kernel with `length_scale = sigma` (the `gamma` of `rbf_kernel` is `1 / (2 sigma^2)`).

  See also [[exponential]], [[laplacian]], [[cauchy]], [[matern]]."
  ([] (gaussian nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv (m/sq (double (distance x y)))))))))

(defn exponential
  "Creates the exponential kernel `exp(-d / (2 sigma^2))`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`); note that `sigma` enters squared
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance. It equals [[laplacian]] with `sigma_L = 2 sigma^2` and the Matern kernel of order 1 (the sklearn Matern with `nu = 0.5`) with a length scale `2 sigma^2`. Genton writes the kernel as `exp(-d / theta)`; use [[laplacian]] with `:sigma theta` for that form.

  See also [[laplacian]], [[gaussian]], [[matern]]."
  ([] (exponential nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-s2halfinv (m// (m/* -2.0 sigma sigma))]
     (fn ^double [x y] (m/exp (m/* -s2halfinv ^double (distance x y)))))))

(defn laplacian
  "Creates the Laplacian kernel `exp(-d / sigma)`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance and for the manhattan distance. The default distance is euclidean, while the scikit-learn `laplacian_kernel` uses the manhattan (L1) distance with `gamma = 1 / sigma`: pass `:distance v/dist-abs` to get it. With the euclidean distance the kernel equals the Matern kernel of order 1 with the length scale `sigma`.

  See also [[exponential]], [[gaussian]], [[matern]]."
  ([] (laplacian nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (let [-sinv (m// -1.0 sigma)]
     (fn ^double [x y] (m/exp (m/* -sinv ^double (distance x y)))))))

(defn anova
  "Creates the ANOVA kernel, a sum over the components of a Gaussian term: `sum_i exp(-sigma (x_i^k - y_i^k)^2)^d`.

  Parameters:

  - `:sigma` - multiplier of the squared difference (default: `1.0`)
  - `:k` - power applied to every component before the difference is taken (default: `1.0`)
  - `:d` - power applied to every term of the sum (default: `1.0`)

  Returns a function `(fn [x y])` returning a double; K(x, x) is the dimension of the vector. The power `:d` acts on every term, not on the sum as in the anovadot kernel of the R package kernlab: both agree for `:d 1.0`. The power `:k` of a negative component is `NaN` unless `:k` is an integer. The sum of Gaussian terms is positive definite for positive `sigma`. This is not the elementary symmetric polynomial form of the ANOVA decomposition.

  See also [[gaussian]]."
  ([] (anova nil))
  ([{:keys [^double sigma ^double k ^double d]
     :or {sigma 1.0 k 1.0 d 1.0}}]
   (let [powk (fn ^double [^double v] (m/pow v k))
         powd (fn ^double [^double v] (m/pow v d))
         -s (m/- sigma)]
     (fn ^double [x y]
       (let [xk (v/fmap x powk)
             yk (v/fmap y powk)]
         (-> (v/sub xk yk)
             (v/sq)
             (v/mult -s)
             (v/exp)
             (v/fmap powd)
             (v/sum)))))))

(defn hyperbolic-tangent
  "Creates the hyperbolic tangent (sigmoid) kernel `tanh(c + alpha <x, y>)`.

  Parameters:

  - `:alpha` - dot product multiplier (default: `1.0`)
  - `:c` - shift (default: `0.0`)

  Returns a function `(fn [x y])` returning a double in `(-1, 1)`. The kernel is not positive definite in general (Gram matrices with a minimum eigenvalue near -1 were found), it is only conditionally positive definite for some parameters. Corresponds to the sigmoid kernel of scikit-learn with `gamma = alpha` and `coef0 = c`.

  See also [[polynomial]], [[linear]]."
  ([] (hyperbolic-tangent nil))
  ([{:keys [^double alpha ^double c]
     :or {alpha 1.0 c 0.0}}]
   (fn ^double [x y] (m/tanh (m/+ c (m/* alpha (v/dot x y)))))))

(defn rational-quadratic
  "Creates the rational quadratic kernel `1 - d^2 / (d^2 + c)`, which equals `c / (d^2 + c)`.

  Parameters:

  - `:c` - shift (default: `1.0`); it is a squared distance, compare [[multiquadratic]] where `c` is a distance
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance. Equal to [[cauchy]] with `sigma = sqrt(c)` and to the scikit-learn RationalQuadratic kernel with `alpha = 1` and `length_scale = sqrt(c / 2)`.

  See also [[cauchy]], [[multiquadratic]], [[inverse-multiquadratic]]."
  ([] (rational-quadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                      (m// c (m/+ dist c))))))

(defn multiquadratic
  "Creates the multiquadratic kernel `sqrt(d^2 + c^2)`.

  Parameters:

  - `:c` - shift (default: `1.0`); it is a distance, compare [[rational-quadratic]] where `c` is a squared distance
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double; K(x, x) = `abs(c)`. The kernel increases with the distance and is not positive definite: it is conditionally negative definite, so its negation is conditionally positive definite of order 1.

  See also [[inverse-multiquadratic]], [[rational-quadratic]]."
  ([] (multiquadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (let [c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                        (m/sqrt (m/+ dist c2)))))))

(defn inverse-multiquadratic
  "Creates the inverse multiquadratic kernel `1 / sqrt(d^2 + c^2)`.

  Parameters:

  - `:c` - shift (default: `1.0`); it is a distance, compare [[rational-quadratic]] where `c` is a squared distance
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1/abs(c)]`; K(x, x) = `1 / abs(c)`. The kernel is positive definite for a euclidean distance.

  See also [[multiquadratic]], [[rational-quadratic]]."
  ([] (inverse-multiquadratic nil))
  ([{:keys [^double c distance]
     :or {c 1.0 distance v/dist}}]
   (let [c2 (m/* c c)]
     (fn ^double [x y] (let [dist (m/sq (double (distance x y)))]
                        (m// (m/sqrt (m/+ dist c2))))))))

;; http://perso.lcpc.fr/tarel.jean-philippe/publis/jpt-icann05b.pdf

(defn- geometric-internal
  "Unnormalised volume recursion of the geometric kernel for `n >= 2`, `ds = d / (2r)`.

  The seed `phi_2 = (acos(s) - s sqrt(1 - s^2)) / 2` is half of the value printed in Tarel et al. (eq. 4)."
  ^double [^long n ^double ds]
  (loop [k (long 3)
         p1 (m/* 0.5 (m/- (m/acos ds) (m/* ds (m/sqrt (m/- 1.0 (m/* ds ds))))))
         p2 (m/- 1.0 ds)]
    (if (m/> k n)
      p1
      (let [k- (m/- k 1.0)
            curr (m/- (m/* (m// k- k) p2)
                      (m/* (m// 1.0 k) ds
                           (m/pow (m/- 1.0 (m/* ds ds)) (m/* 0.5 k-))))]
        (recur (m/inc k) curr p1)))))

(defn geometric
  "Creates the geometric compactly supported kernel: the volume of the intersection of two n-dimensional balls of radius `r`, divided by the volume of one ball.

  The kernel is a function of the distance `d` of the points, it decreases from `1` at `d = 0` to `0` at `d = 2 r` and is zero beyond. Dimensions `1`, `2` and `3` have names in the registry: triangular, circular and spherical.

  Parameters:

  - `:n` - dimension of the balls (required, integer, at least `1`)
  - `:r` - radius of the balls (default: `1.0`, must be positive); the support radius of the kernel is `2 r`
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `[0, 1]`; K(x, x) = 1 for every `n`, in particular the circular kernel (`n = 2`) is normalised by `2 / pi`. The kernel is positive definite in the euclidean space of dimension `n` and lower. For `n >= 4` the value is computed with a recursion in `n`, so the cost of one evaluation grows linearly with `n`; results which are slightly negative by rounding near the edge of the support are returned as `0.0`.

  Throws `ex-info` when `:n` is missing or smaller than `1`, or when `:r` is not positive.

  See also [[rational-quadratic]], [[pearson]]."
  ([] (geometric nil))
  ([{:keys [n ^double r distance]
     :or {r 1.0 distance v/dist}}]
   (when (nil? n) (throw (invalid-option "geometric" :n n "an integer >= 1 (required)")))
   (when-not (and (number? n) (m/>= (double n) 1.0) (m/== (double n) (m/floor (double n))))
     (throw (invalid-option "geometric" :n n "an integer >= 1")))
   (when-not (m/> r 0.0) (throw (invalid-option "geometric" :r r "a positive number")))
   (let [n (long n)
         r2 (m/* 2.0 r)]
     (case n
       1 (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (m/>= dist r2) 0.0
                                (m/- 1.0 (m// dist r2)))))
       2 (let [c (m// 2.0 m/PI)]
           (fn ^double [x y] (let [^double dist (distance x y)]
                              (if (m/>= dist r2) 0.0
                                  (let [ds (m// dist r2)]
                                    (m/* c (m/- (m/acos ds)
                                                (m/* ds (m/sqrt (m/- 1.0 (m/* ds ds)))))))))))
       3 (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (m/>= dist r2) 0.0
                                (let [ds (m// dist r2)]
                                  (m/+ (m/- 1.0 (m/* 1.5 ds)) (m/* 0.5 ds ds ds))))))
       (let [z (geometric-internal n 0.0)]
         (fn ^double [x y] (let [^double dist (distance x y)]
                            (if (m/>= dist r2) 0.0
                                (m/max 0.0 (m// (geometric-internal n (m// dist r2)) z))))))))))

(defn wave
  "Creates the wave (sinc) kernel `(sigma / d) sin(d / sigma)`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `[-0.22, 1]`; K(x, x) = 1. The kernel is positive definite in the euclidean space of dimension 3 and lower only; in higher dimensions Gram matrices with negative eigenvalues exist.

  See also [[periodic]], [[bessel2]]."
  ([] (wave nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (fn ^double [x y] (let [^double dist (distance x y)]
                      (if (m/zero? dist) 1.0
                          (m/* (m// sigma dist) (m/sin (m// dist sigma))))))))

(defn periodic
  "Creates the periodic kernel `exp(-2 sin^2(pi d / periodicity) / sigma^2)`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:periodicity` - period (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1 and K = 1 at every distance which is a multiple of the period. The kernel is positive definite for one dimensional input; for vectors of higher dimension the euclidean distance gives a kernel which is not positive definite (Gram matrices with negative eigenvalues exist). Equal to the scikit-learn ExpSineSquared kernel with `length_scale = sigma`.

  See also [[gaussian]], [[wave]]."
  ([] (periodic nil))
  ([{:keys [^double sigma ^double periodicity distance]
     :or {sigma 1.0 periodicity 1.0 distance v/dist}}]
   (let [p (m// m/PI periodicity)
         s2 (m/* sigma sigma)]
     (fn ^double [x y] (let [^double dist (distance x y)]
                        (m/exp (m// (m/* -2.0 (m/sq (m/sin (m/* p dist)))) s2)))))))

(defn power
  "Creates the power kernel, the negative power of the distance `-d^p`.

  Parameters:

  - `:p` - exponent (default: `2.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double; K(x, x) = `-0.0` for a positive `p`. The kernel is conditionally positive definite for `0 < p <= 2` and is not positive definite.

  See also [[log]], [[generalized-t-student]]."
  ([] (power nil))
  ([{:keys [^double p distance]
     :or {p 2.0 distance v/dist}}]
   (fn ^double [x y] (m/- (m/pow (double (distance x y)) p)))))

(defn log
  "Creates the logarithmic kernel `-log(1 + d^p)`, where `d` is the distance of the points.

  Parameters:

  - `:p` - exponent (default: `2.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double; K(x, x) = `-0.0`. The kernel is conditionally positive definite for `0 < p <= 2` and is not positive definite.

  See also [[power]], [[generalized-t-student]]."
  ([] (log nil))
  ([{:keys [^double p distance]
     :or {p 2.0 distance v/dist}}]
   (fn ^double [x y] (m/- (m/log (m/inc (m/pow (double (distance x y)) p)))))))

(defn spline
  "Creates the spline kernel, a product over the components of `1 + x y + x y min(x, y) - (x + y) min(x, y)^2 / 2 + min(x, y)^3 / 3`.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning a double. The components are expected to be non-negative. The kernel is positive definite on non-negative input and equals the splinedot kernel of the R package kernlab.

  See also [[b-spline]]."
  ([] (spline nil))
  ([_] (fn ^double [x y]
         (double
          (reduce m/* 1.0 (map (fn [^double xi ^double yi]
                                 (let [xiyi (m/* xi yi)
                                       m (min xi yi)
                                       m2 (m/* m m)]
                                   (inc (m/+ xiyi
                                             (m/* xiyi m)
                                             (m/* -0.5 (m/+ xi yi) m2)
                                             (m/* m/THIRD m2 m)))))
                               (if (number? x) [x] x)
                               (if (number? y) [y] y)))))))

;; https://www.researchgate.net/publication/2538959_Inverse_B-spline_interpolation

(defn- b-spline-beta
  "Centered cardinal B-spline of degree `n` at `x >= 0`, supported on `abs(x) < (n + 1) / 2`.

  Returns exactly `0.0` at and beyond the edge of the support, where the alternating sum is zero only in exact arithmetic."
  ^double [^long n ^double x]
  (let [n!inv (m/inv-factorial n)
        n+ (m/long-inc n)
        nx (m/+ x (m/* 0.5 n+))]
    (if (m/>= nx n+)
      0.0
      (loop [k (long 0)
             neg 1.0
             curr 0.0]
        (if (or (m/> k n+)
                (m/<= nx k))
          (m/* n!inv curr)
          (recur (m/inc k)
                 (m/* -1.0 neg)
                 (m/+ curr (m/* neg (m/combinations n+ k)
                                (m/fpow (m/- nx k) n)))))))))

(defn b-spline
  "Creates the B-spline kernel, a product over the components of the centered cardinal B-spline of degree `2 n + 1` evaluated at the difference of the components.

  Parameters:

  - `:n` - half of the degree minus one, a non-negative integer (default: `2`); the degree of the spline is `2 n + 1`

  Returns a function `(fn [x y])` returning a double. The kernel is supported on the components with `abs(x_i - y_i) < n + 1`; it is exactly `0.0` at and beyond the edge of the support in every component. The kernel is positive definite (the Fourier transform of the odd degree cardinal B-spline is a non-negative power of sinc). It is not normalised: K(x, x) is the central value of the spline for one component, for example `2/3` for `n = 1`.

  Throws `ex-info` when `:n` is negative.

  See also [[spline]]."
  ([] (b-spline nil))
  ([{:keys [^long n]
     :or {n 2.0}}]
   (when (m/neg? n) (throw (invalid-option "b-spline" :n n "a non-negative integer")))
   (let [nn (m/long-inc (m/long-mult 2 n))]
     (fn ^double [x y]
       (-> (v/sub x y)
           (v/abs)
           (v/fmap (fn [^double x] (b-spline-beta nn x)))
           (v/prod))))))

(defn bessel
  "Creates the Bessel kernel `J_{v+1}(sigma d) d^(n (v+1))`, where `J` is the Bessel function of the first kind and `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:n` - exponent factor (default: `2.0`)
  - `:v` - order of the Bessel function minus one (default: `-1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double. The kernel is the one of http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/ and is not normalised: for `v > -1` the kernel is `0` at `d = 0`, K(x, x) = 0, only for the default `v = -1` it is `1`. The kernel is not positive definite in general (11 of 12 settings tried had Gram matrices with negative eigenvalues).

  See also [[bessel2]]."
  ([] (bessel nil))
  ([{:keys [^double sigma ^double n ^double v distance]
     :or {sigma 1.0 n 2.0 v -1.0 distance v/dist}}]
   (fn ^double [x y] (let [^double dist (distance x y)
                          v+ (m/inc v)]
                      (/ (special/bessel-J v+ (m/* sigma dist))
                         (m/pow dist (m/- (m/* n v+))))))))

;; R kernlab
(defn bessel2
  "Creates the Bessel kernel of the R package kernlab: `(J_order(b) b^(-order) / lim)^degree`, where `b = sigma d`, `d` is the distance of the points and `lim = 1 / (Gamma(order + 1) 2^order)` is the limit of the numerator at zero.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:degree` - exponent (default: `1.0`)
  - `:order` - order of the Bessel function of the first kind (default: `0.0`); the default of kernlab is `1`
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double; K(x, x) = 1. For `b < 1e-5` the limit `1.0` is returned (kernlab switches at `1e-4`, which gives differences of about `1e-9`). The kernel is positive definite in the euclidean space of dimension up to `2 order + 2` (for an integer `degree`); with the default order this is the plane, so Gram matrices of points in 3 or 5 dimensions are not positive semi-definite.

  See also [[bessel]], [[wave]]."
  ([] (bessel2 nil))
  ([{:keys [^double sigma ^double order ^double degree distance]
     :or {sigma 1.0 order 0.0 degree 1.0 distance v/dist}}]
   (fn ^double [x y]
     (let [lim (m// (m/* (special/gamma (m/inc order)) (m/pow 2.0 order)))
           bkt (m/* sigma ^double (distance x y))]
       (if (m/< bkt 1.0e-5)
         1.0
         (m/pow (m// (m/* (special/bessel-J order bkt)
                          (m/pow bkt (m/- order))) lim) degree))))))

(defn cauchy
  "Creates the Cauchy kernel `1 / (1 + (d / sigma)^2)`, where `d` is the distance of the points.

  Parameters:

  - `:sigma` - scale (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance. Equal to [[rational-quadratic]] with `c = sigma^2`.

  See also [[rational-quadratic]], [[generalized-t-student]]."
  ([] (cauchy nil))
  ([{:keys [^double sigma distance]
     :or {sigma 1.0 distance v/dist}}]
   (fn ^double [x y] (m// (m/inc (m/sq (m// ^double (distance x y) sigma)))))))

(defn chi-square
  "Creates the chi-square kernel `sum_i 2 x_i y_i / (x_i + y_i)` for non-negative vectors, for example histograms.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning a double. A component with `x_i + y_i = 0` contributes `0` (not `NaN`), as in the additive chi-square kernel of scikit-learn. With `A = additive_chi2_kernel`, which is the negative `sum (x_i - y_i)^2 / (x_i + y_i)`, the kernel equals `sum (x_i + y_i) / 2 + A / 2`. The kernel is positive definite on non-negative input.

  See also [[chi-square2]], [[histogram]], [[hellinger]]."
  ([] (chi-square nil))
  ([_] (fn ^double [x y]
         (double
          (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                 (let [s (m/+ xi yi)]
                                   (if (m/zero? s)
                                     0.0
                                     (m// (m/* 2.0 xi yi) s))))
                               (if (number? x) [x] x)
                               (if (number? y) [y] y)))))))

(defn chi-square2
  "Creates the second chi-square kernel `1 - sum_i (x_i - y_i)^2 / (0.5 (x_i + y_i))` for non-negative vectors.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning a double. A component with `x_i + y_i = 0` contributes `0` (not `NaN`), as in the additive chi-square kernel of scikit-learn. With `A = additive_chi2_kernel`, which is the negative `sum (x_i - y_i)^2 / (x_i + y_i)`, the kernel equals `1 + 2 A`. The kernel is conditionally positive definite and is not positive definite.

  See also [[chi-square]], [[histogram]]."
  ([] (chi-square2 nil))
  ([_] (fn ^double [x y]
         (- 1.0 ^double (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                               (let [s (m/+ xi yi)]
                                                 (if (m/zero? s)
                                                   0.0
                                                   (m// (m/sq (m/- xi yi))
                                                        (m/* 0.5 s)))))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y)))))))

(defn histogram
  "Creates the histogram intersection kernel `sum_i min(x_i, y_i)`.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning a double. The kernel is positive definite on non-negative input.

  See also [[generalized-histogram]], [[chi-square]]."
  ([] (histogram nil))
  ([_] (fn ^double [x y] (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                               (min xi yi))
                                             (if (number? x) [x] x)
                                             (if (number? y) [y] y))))))

(defn generalized-histogram
  "Creates the generalized histogram intersection kernel `sum_i min(abs(x_i)^p, abs(y_i)^p)`.

  Parameters:

  - `:p` - exponent (default: `2.0`)

  Returns a function `(fn [x y])` returning a double. The kernel is positive definite for a positive `p`.

  See also [[histogram]]."
  ([] (generalized-histogram nil))
  ([{:keys [^double p]
     :or {p 2.0}}]
   (fn ^double [x y] (reduce m/+ 0.0 (map (fn [^double xi ^double yi]
                                           (min (m/pow (m/abs xi) p)
                                                (m/pow (m/abs yi) p)))
                                         (if (number? x) [x] x)
                                         (if (number? y) [y] y))))))

(defn generalized-t-student
  "Creates the generalized t-student kernel `1 / (1 + d^p)`, where `d` is the distance of the points.

  Parameters:

  - `:p` - exponent (default: `2.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance and `0 < p <= 2`. With `p = 2` it equals [[cauchy]] with `sigma = 1`.

  See also [[cauchy]], [[power]], [[log]]."
  ([] (generalized-t-student nil))
  ([{:keys [^double p distance]
     :or {p 2.0 distance v/dist}}]
   (fn ^double [x y] (m// (m/inc (m/pow (double (distance x y)) p))))))

(defn- dirichlet-dim
  ^double [^double nhalf ^double delta]
  (if (m/zero? delta)
    nhalf
    (let [den (m/* 2.0 (m/sin (m/* 0.5 delta)))]
      (if (m/zero? den)
        nhalf
        (m// (m/sin (m/* nhalf delta)) den)))))

(defn dirichlet
  "Creates the Dirichlet kernel, a product over the components of `sin((n + 1/2) delta) / (2 sin(delta / 2))`, where `delta = x_i - y_i`.

  Parameters:

  - `:n` - order (default: `1.0`)

  Returns a function `(fn [x y])` returning a double. The kernel is not normalised: for every component the value at `delta = 0` is `n + 1/2`. For an integer `n` the kernel is positive definite (a half of the Dirichlet kernel is the sum of cosines with non-negative coefficients), for a non-integer `n` positive definiteness is not guaranteed.

  See also [[periodic]]."
  ([] (dirichlet nil))
  ([{:keys [^double n]
     :or {n 1.0}}]
   (let [nhalf (m/+ 0.5 n)]
     (fn ^double [x y]
       (let [diff (v/sub x y)]
         (if (number? diff)
           (dirichlet-dim nhalf diff)
           (reduce m/* 1.0 (map (partial dirichlet-dim nhalf) diff))))))))

(defn hellinger
  "Creates the Hellinger kernel `sum_i sqrt(x_i y_i)` for non-negative vectors.

  Parameters:

  - `params` (optional) - ignored, accepted for a uniform signature.

  Returns a function `(fn [x y])` returning a double. A negative product `x_i y_i` contributes `0`. The kernel is the dot product of the square roots of the vectors, so it is positive definite on non-negative input.

  See also [[chi-square]], [[histogram]]."
  ([] (hellinger nil))
  ([_] (fn ^double [x y] (-> (v/emult x y) v/safe-sqrt v/sum))))

(defn pearson
  "Creates the Pearson VII kernel (PUK) `1 / (1 + (2 d sqrt(2^(1/omega) - 1) / sigma)^2)^omega`, where `d` is the distance of the points.

  The scale is chosen so that `sigma` is the full width at half maximum of the peak for every `omega`: the kernel is `1/2` at `d = sigma / 2`. This is the form of Ustun et al. (2006) and of the PUK kernel in Weka.

  Parameters:

  - `:sigma` - width of the peak (default: `1.0`)
  - `:omega` - tailing factor (default: `1.0`); `omega = 1` gives the Lorentzian, [[cauchy]] with `sigma / 2` and a large `omega` approaches [[gaussian]]
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance.

  See also [[cauchy]], [[gaussian]]."
  ([] (pearson nil))
  ([{:keys [^double sigma ^double omega distance]
     :or {sigma 1.0 omega 1.0 distance v/dist}}]
   (let [c (m// (m/* 4.0 (m/dec (m/pow 2.0 (m// omega))))
                (m/sq sigma))]
     (fn ^double [x y]
       (m// (m/pow (m/inc (m/* c (m/sq (double (distance x y))))) omega))))))

(defn hyperbolic-secant
  "Creates the hyperbolic secant kernel `sech(a d)`, where `d` is the distance of the points.

  Parameters:

  - `:a` - scaling factor (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for one dimensional input.

  See also [[gaussian]], [[cauchy]]."
  ([] (hyperbolic-secant nil))
  ([{:keys [^double a distance]
     :or {a 1.0 distance v/dist}}]
   (fn ^double [x y]
     (m/sech (m/* a ^double (distance x y))))))

(defn matern
  "Creates the Matern kernel of the half-integer smoothness `nu = order / 2`, `2^(1-nu) / Gamma(nu) v^nu K_nu(v)` with `v = 2 sqrt(nu) d / theta`, where `K` is the modified Bessel function of the second kind and `d` is the distance of the points.

  This is the parametrisation of Handcock and Stein: `theta = sqrt(2) ell`, where `ell` is the length scale of the Matern kernel of scikit-learn (`Matern(length_scale = theta / sqrt(2), nu = order / 2)`, whose argument is `sqrt(2 nu) d / ell`). The registry presets are `:matern-12` (`order` 1, equal to [[exponential]] in the form `exp(-sqrt(2) d / theta)`), `:matern-32` (`order` 3) and `:matern-52` (`order` 5).

  Parameters:

  - `:order` - twice the smoothness, an odd positive integer (default: `1`); `5` gives the Matern 5/2 kernel
  - `:theta` - scale (default: `1.0`)
  - `:distance` - distance function (default: euclidean)

  Returns a function `(fn [x y])` returning a double in `(0, 1]`; K(x, x) = 1. The kernel is positive definite for a euclidean distance and an odd `:order`. At a distance so large that the intermediate power overflows (above `v = 700`) the result is `0.0`.

  Throws `ex-info` when `:order` is not an odd positive integer: the even orders (integer smoothness) are not supported by the half-integer Bessel function.

  See also [[exponential]], [[laplacian]], [[gaussian]]."
  ([] (matern nil))
  ([{:keys [order ^double theta distance]
     :or {order 1 theta 1.0 distance v/dist}}]
   (when-not (and (number? order)
                  (m/== (double order) (m/floor (double order)))
                  (m/pos? (double order))
                  (m/< (double order) 1.0e9)
                  (m/odd? (long order)))
     (throw (invalid-option "matern" :order order "an odd positive integer")))
   (let [order (long order)
         mu (m// order 2.0)
         gf (m// (m/* (m/pow 2.0 (m/- 1.0 mu))) (special/gamma mu))
         s (m// (m/* 2.0 (m/sqrt mu)) theta)]
     (fn ^double [x y]
       (let [v (m/* s (double (distance x y)))]
         (if (m/< v 1.0e-16)
           1.0
           (let [res (m/* gf (m/pow v mu) (special/bessel-K-half-odd order v))]
             ;; v^mu overflows before the Bessel function underflows (0 * Inf): the limit is 0
             (if (and (m/nan? res) (m/> v 700.0)) 0.0 res))))))))

(defn rbf->kernel
  "Converts a radial basis function into a vector kernel using a distance function.

  Parameters:

  - `rbf-kernel` - a function of one double (the distance), for example created by `fastmath.kernel/rbf`
  - `distance` (optional) - distance function of two points (default: euclidean)

  Returns a function `(fn [x y])` which applies `rbf-kernel` to the distance of `x` and `y`.

  See also [[gaussian]], [[matern]]."
  ([rbf-kernel] (rbf->kernel rbf-kernel v/dist))
  ([rbf-kernel distance]
   (fn [x y] (rbf-kernel (distance x y)))))
