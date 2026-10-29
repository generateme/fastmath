;; Radial Basis Functions
(ns fastmath.kernel.rbf
  "Radial basis function (RBF) families.

  Every family is a constructor which takes a map of parameters and returns a function `r -> double` of the distance `r` (a negative `r` is treated as `abs(r)` unless stated otherwise). Constructors are called with a map, with `nil` or without arguments; parameters which are not optional have no default and a missing (or non-numeric, or non-finite) one throws `ex-info` with `:type :invalid-input`, `:family` (the family keyword) and `:property` (the missing key).

  Every family accepts `:scale` (default `1.0`), which divides the distance: the kernel is evaluated at `rho = r / scale`. For the compactly supported families (`truncated-power`, `wendland`, `gneiting`, `wu`) `scale` is the support radius, for the others it is a length scale.

  Families by kind:

  - positive definite in every dimension: [[gaussian]], [[matern]] (odd order), [[generalized-multiquadratic]] with negative `beta` (inverse multiquadric),
  - compactly supported: [[truncated-power]], [[wendland]], [[gneiting]] (oscillating), [[wu]],
  - positive definite up to a dimension: [[gaussians-laguerre]], [[poisson]], [[whittaker]],
  - conditionally positive definite (use them with polynomial terms): [[linear]], [[generalized-multiquadratic]], [[radial-powers]], [[thin-plate-splines]], [[shifted-surface-splines]].

  The registry `fastmath.kernel/rbf` creates these functions by keyword and adds presets.

  Literature: [1] Meshfree Approximation Methods with MATLAB (Fasshauer, 2007) and MATH 590 lecture notes, chapter 11; [2] Wendland, Scattered Data Approximation (2005); [3] Luh, The Shape Parameter in the Shifted Surface Spline (arXiv:1011.0124); [4] Gneiting, Compactly supported correlation functions (2002); [5] Wu, Compactly supported positive definite radial functions (1995)."
  (:require [fastmath.core :as m]
            [fastmath.calculus.quadrature :as calc]
            [fastmath.calculus.finite :as finite]
            [fastmath.special :as special]
            [fastmath.polynomials :as poly]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Literature
;; http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol4/Choosing_Basis_Functions_and_Shape_Parameters.pdf
;; https://www.math.ucdavis.edu/~saito/data/jim/buhmann-actanumerica.pdf
;; https://www.math.utah.edu/vigre/pubs/pdfs/flw08.pdf
;; https://num.math.uni-goettingen.de/schaback/teaching/sc.pdf

;; parameter helpers

(defn- invalid-input!
  "Throws `ex-info` for an invalid parameter `property` of `family`."
  [family property message value]
  (throw (ex-info (str "rbf " (name family) ", " (name property) ": " message)
                  {:type :invalid-input
                   :family family
                   :property property
                   :value value})))

(defn- required-double
  "Returns the finite number under `property` in `params` as a double or throws `ex-info`."
  ^double [family params property]
  (let [v (get params property)]
    (cond
      (not (number? v)) (invalid-input! family property "required parameter is missing or is not a number" v)
      (not (m/valid-double? (double v))) (invalid-input! family property "must be finite" v)
      :else (double v))))

(defn- required-long
  "Returns the integer valued number under `property` in `params` as a long or throws `ex-info`."
  ^long [family params property]
  (let [v (required-double family params property)]
    (if (m/== v (m/floor v))
      (long v)
      (invalid-input! family property "must be an integer" v))))

(defn- scale-of
  "Returns the `:scale` of `params` (default 1.0, also for an explicit nil)."
  ^double [params]
  (let [s (get params :scale)]
    (if (nil? s) 1.0 (double s))))

;; simple kernels

(defn linear
  "Linear kernel, the distance itself: `abs(r / scale)`.

  Parameters:

  - `params` (map, optional):
    - `:scale` (double, default `1.0`): length scale, the distance is divided by it.

  Returns a function `r -> double`. The value at `r = 0` is `0.0`, the function is symmetric in `r` and grows without bound.

  This is not a positive definite function: `abs(r)` is conditionally negative definite of order 1, so interpolation with it needs a constant polynomial term (and the sign is positive here, while [[radial-powers]] with `beta = 1` gives `-abs(r)`, the conditionally positive definite form).

  See also [[radial-powers]], [[thin-plate-splines]], [[gaussian]]."
  ([] (linear nil))
  ([params]
   (let [scale (scale-of params)]
     (fn ^double [^double x] (m/abs (m// x scale))))))

(defn gaussian
  "Gaussian kernel: `exp(-(r / scale)^2)`.

  Parameters:

  - `params` (map, optional):
    - `:scale` (double, default `1.0`): length scale. The value at `r = scale` is `exp(-1)`; the standard deviation of the matching Gaussian density is `scale / sqrt(2)`.

  Returns a function `r -> double`. The value at `r = 0` is `1.0`; the function is symmetric, positive, decreasing in `abs(r)` and not compactly supported (it underflows to `0.0` for `abs(r / scale)` above about 27).

  Strictly positive definite in every dimension.

  See also [[matern]], [[wendland]], [[generalized-multiquadratic]]."
  ([] (gaussian nil))
  ([params]
   (let [scale (scale-of params)]
     (fn ^double [^double x] (m/exp (m/- (m/sq (m// x scale))))))))

;; RBF kernels
;; [1] https://www.math.unipd.it/~demarchi/RBF/LectureNotes.pdf
;; [2] Meshfree Approximation Methods with MATLAB, Gregory E. Fasshauer

;; [1] p.32

(defn truncated-power
  "Truncated power kernel: `(1 - r / scale)^k` for `r` below `scale` and `0` from there on.

  Parameters:

  - `params` (map):
    - `:k` (double, required): the exponent. Closed forms are used for `1/3`, `0.5`, `1`, `2` and `3`.
    - `:scale` (double, default `1.0`): the support radius.

  Returns a function `r -> double`. The value at `r = 0` is `1.0`. The kernel is symmetric, decreasing, vanishes at `abs(r) = scale` (for `k > 0`; for `k = 0` the value at `abs(r) = scale` is `1.0` and the jump is just beyond) and is `0.0` beyond it. A missing `:k` throws `ex-info`.

  Positive definite in R^d only for `d <= 2k - 1` (Askey: the exponent must be at least `(d + 1) / 2`): `k = 1` in R^1, `k = 2` in R^1 to R^3, `k = 3` in R^1 to R^5. The exponents `1/3` and `0.5` are not positive definite in any dimension.

  See also [[wendland]], [[gneiting]], [[wu]]."
  ([] (truncated-power nil))
  ([params]
   (let [k (required-double :truncated-power params :k)
         scale (scale-of params)]
     (cond
       (m/== k m/THIRD) (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0)
                                                                             (m/cbrt (m/- 1.0 ex)) 0.0)))
       (m/== k 0.5) (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0) (m/sqrt (m/- 1.0 ex)) 0.0)))
       (m/== k 1.0) (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0) (m/- 1.0 ex) 0.0)))
       (m/== k 2.0) (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0) (m/sq (m/- 1.0 ex)) 0.0)))
       (m/== k 3.0) (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0) (m/cb (m/- 1.0 ex)) 0.0)))
       :else (fn ^double [^double x] (let [ex (m/abs (m// x scale))] (if (m/<= ex 1.0) (m/pow (m/- 1.0 ex) k) 0.0)))))))

;; [1] p.33

(defn gaussians-laguerre
  "Gaussians-Laguerre kernel: `exp(-rho^2) * L_n^(d/2)(rho^2)` with `rho = r / scale`, the generalized Laguerre polynomial of degree `n` and order `d/2`.

  Parameters:

  - `params` (map):
    - `:dimension` (double, required, positive): the dimension `d`, it sets the order `d/2` of the polynomial.
    - `:degree` (long, required, non-negative): the degree `n`. A whole double such as `2.0` is accepted.
    - `:scale` (double, default `1.0`): length scale.

  Returns a function `r -> double`. The kernel is not normalised: the value at `r = 0` is `L_n^(d/2)(0)`, the binomial coefficient `C(n + d/2, n)` (for example `1.5` for `:dimension 1.0`, `:degree 1`; `1.0` for degree 0). It is symmetric and decays like a Gaussian times a polynomial, so it is not compactly supported and oscillates for a positive degree. A missing `:dimension` or `:degree` throws `ex-info`.

  Positive definite in R^s for every `s <= dimension` (and not beyond it).

  See also [[gaussian]], [[poisson]]."
  ([] (gaussians-laguerre nil))
  ([params]
   (let [dimension (required-double :gaussians-laguerre params :dimension)
         degree (required-long :gaussians-laguerre params :degree)
         scale (scale-of params)
         hdimension (m/* 0.5 dimension)]
     (when-not (m/pos? dimension)
       (invalid-input! :gaussians-laguerre :dimension "must be positive" dimension))
     (when (m/neg? degree)
       (invalid-input! :gaussians-laguerre :degree "must be non-negative" degree))
     (fn ^double [^double x]
       (let [ex (m/sq (m// x scale))
             e (m/exp (m/- ex))]
         ;; the polynomial overflows where the exponential underflows: avoid Inf * 0
         (if (m/zero? e) 0.0 (m/* e (poly/eval-laguerre-L degree hdimension ex))))))))

;; [1] p.33

(defn poisson
  "Poisson (Bessel) kernel: `2^p * Gamma(p + 1) * J_p(rho) / rho^p` with `p = d/2 - 1`, `rho = r / scale` and the Bessel function of the first kind `J_p`.

  Parameters:

  - `params` (map):
    - `:d` (double, required, at least `1`): the dimension, it sets the Bessel order `d/2 - 1`. Closed forms are used for `1` (cosine), `2` (`J_0`), `3` (`sin(rho)/rho`) and `4` (`2 J_1(rho)/rho`).
    - `:scale` (double, default `1.0`): length scale.

  Returns a function `r -> double`. The function is normalised: the value at `r = 0` is `1.0` for every `d` (the limit of `J_p(rho)/rho^p`, divided by its value `1/(2^p Gamma(p + 1))`). It is symmetric, oscillating, and decays only like `rho^(-(d - 1)/2)`, so it is not compactly supported. A missing `:d` throws `ex-info`.

  Positive definite in R^s for every `s <= d`.

  See also [[gaussians-laguerre]], [[matern]]."
  ([] (poisson nil))
  ([params]
   (let [d (required-double :poisson params :d)
         scale (scale-of params)]
     (when (m/< d 1.0)
       (invalid-input! :poisson :d "must be at least 1" d))
     (cond
       (m/== d 1.0) (fn ^double [^double x] (m/cos (m// x scale)))
       (m/== d 2.0) (fn ^double [^double x] (special/bessel-J0 (m/abs (m// x scale))))
       (m/== d 3.0) (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                           (if (m/< ex 1.0e-13) 1.0 (m// (m/sin ex) ex))))
       (m/== d 4.0) (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                           (if (m/< ex 1.0e-13) 1.0 (m// (m/* 2.0 (special/bessel-J1 ex)) ex))))
       :else (let [p (m/dec (m/* 0.5 d))
                   c (m/exp (m/+ (m/* p m/LN2) (special/log-gamma (m/inc p))))
                   c2 (m// -0.5 d)]
               (fn ^double [^double x]
                 (let [ex (m/abs (m// x scale))]
                   (if (m/< ex 1.0e-5)
                     (m/inc (m/* c2 ex ex)) ;; series, the error is below 1e-20
                     (m/* c (m// (special/bessel-J p ex) (m/pow ex p)))))))))))

;; [1] p. 34

;; note: matern52, from formula it's 1/3 of the example from paper

(defn matern
  "Matern kernel for an odd `order`: `2^(1 - v) / Gamma(v) * rho^v * K_v(rho)` with `v = order / 2`, `rho = r / scale` and the modified Bessel function of the second kind `K_v`.

  Parameters:

  - `params` (map):
    - `:order` (long, required, odd and positive): twice the smoothness `v`. `1` gives `exp(-rho)`, `3` gives `(1 + rho) exp(-rho)` and `5` gives `(1 + rho + rho^2/3) exp(-rho)`. A whole double such as `3.0` is accepted.
    - `:scale` (double, default `1.0`): length scale, applied without the `sqrt(order)` factor of other conventions: the length scale `ell` of the usual form (for example in scikit-learn) corresponds to `scale = ell / sqrt(order)`.

  Returns a function `r -> double`. The value at `r = 0` is `1.0`; it is symmetric, positive, decreasing and not compactly supported (it is `0.0` where it underflows). A missing `:order` throws `ex-info`; an even, zero, negative or non-integer `:order` throws `ex-info` with `:property :order`.

  Strictly positive definite in every dimension.

  See also [[gaussian]], [[wendland]]."
  ([] (matern nil))
  ([params]
   (let [order (required-long :matern params :order)
         scale (scale-of params)]
     (when-not (and (m/pos? order) (m/odd? order))
       (invalid-input! :matern :order "must be a positive odd integer" order))
     (let [v (m/* 0.5 order)
           rg (m/exp (m/- (m/* (m/- 1.0 v) m/LN2) (special/log-gamma v)))]
       (fn ^double [^double x]
         (let [ex (m/abs (m// x scale))]
           (if (m/< ex 1.0e-20)
             1.0
             (let [k (special/bessel-K-half-odd order ex)]
               ;; K underflows to 0 where ex^v overflows: avoid 0 * Inf
               (if (m/zero? k) 0.0 (m/* rg (m/pow ex v) k))))))))))

;; [1] p. 43

(defn- infer-negate
  "Sign factor for the families with `negate?`. `kind` is `:gmq`, `:radial-powers` or `:tps`."
  ^double [kind ^double beta negate?]
  (cond
    (or (nil? negate?) (= :auto negate?))
    (case kind
      ;; (-1)^ceil(beta) for beta > 0, positive (positive definite) otherwise
      :gmq (if (and (m/pos? beta) (m/odd? (long (m/ceil beta)))) -1.0 1.0)
      ;; (-1)^ceil(beta/2)
      :radial-powers (if (m/even? (long (m/ceil (m/* 0.5 beta)))) 1.0 -1.0)
      ;; (-1)^(ceil(beta) + 1)
      :tps (if (m/even? (long (m/ceil beta))) -1.0 1.0))
    negate? -1.0
    :else 1.0))

(defn generalized-multiquadratic
  "Generalized multiquadratic kernel: `sign * (1 + rho^2)^beta` with `rho = r / scale`.

  Parameters:

  - `params` (map):
    - `:beta` (double, required): the exponent. `0.5` is the multiquadric, `-0.5` the inverse multiquadric, `-1` the inverse quadric. Closed forms are used for `1/3`, `0.5`, `1`, `2`, `3` and the negatives of them.
    - `:negate?` (`:auto`, `true` or `false`, default `:auto`): the sign. `true` multiplies by `-1`, `false` by `+1`. `:auto` (also for `nil`) chooses the sign of the conditionally positive definite form: `(-1)^ceil(beta)` for `beta > 0` (negative for `beta` in `(0, 1]`, positive for `(1, 2]`, negative for `(2, 3]`, ...) and `+1` for `beta <= 0`.
    - `:scale` (double, default `1.0`): length scale (shape parameter is `1 / scale`).

  Returns a function `r -> double`. The value at `r = 0` is the sign. The function is symmetric, smooth and grows like `abs(r)^(2 beta)` for a positive `beta` (it overflows to infinity for a huge `r`), and decays for a negative `beta`. A missing `:beta` throws `ex-info`.

  Positive definiteness: strictly positive definite in every dimension for `beta < 0` (sign `+1`). For `beta > 0` conditionally positive definite of order `ceil(beta)` (needs polynomial terms of degree `ceil(beta) - 1`) with the `:auto` sign, checked numerically for `beta` up to 3.5; beyond that the sign rule follows the literature but is not verified. A non-negative integer `beta` gives a polynomial in `r` which is not a usable basis function on its own.

  See also [[radial-powers]], [[thin-plate-splines]], [[shifted-surface-splines]]."
  ([] (generalized-multiquadratic nil))
  ([params]
   (let [beta (required-double :generalized-multiquadratic params :beta)
         scale (scale-of params)
         fact (infer-negate :gmq beta (get params :negate?))]
     (cond
       (m/== beta m/THIRD) (fn ^double [^double x] (m/* fact (m/cbrt (m/inc (m/sq (m// x scale))))))
       (m/== beta (m/- m/THIRD)) (fn ^double [^double x] (m/* fact (m/cbrt (m// (m/inc (m/sq (m// x scale)))))))
       (m/== beta 0.5) (fn ^double [^double x] (m/* fact (m/sqrt (m/inc (m/sq (m// x scale))))))
       (m/== beta -0.5) (fn ^double [^double x] (m/* fact (m/sqrt (m// (m/inc (m/sq (m// x scale)))))))
       (m/== beta 1.0) (fn ^double [^double x] (m/* fact (m/inc (m/sq (m// x scale)))))
       (m/== beta -1.0) (fn ^double [^double x] (m/* fact (m// (m/inc (m/sq (m// x scale))))))
       (m/== beta 2.0) (fn ^double [^double x] (m/* fact (m/sq (m/inc (m/sq (m// x scale))))))
       (m/== beta -2.0) (fn ^double [^double x] (m/* fact (m/sq (m// (m/inc (m/sq (m// x scale)))))))
       (m/== beta 3.0) (fn ^double [^double x] (m/* fact (m/cb (m/inc (m/sq (m// x scale))))))
       (m/== beta -3.0) (fn ^double [^double x] (m/* fact (m/cb (m// (m/inc (m/sq (m// x scale)))))))
       :else (fn ^double [^double x] (m/* fact (m/pow (m/inc (m/sq (m// x scale))) beta)))))))

;; [1] p. 43

(defn radial-powers
  "Radial powers kernel: `sign * rho^beta` with `rho = abs(r / scale)`.

  Parameters:

  - `params` (map):
    - `:beta` (double, required): the exponent. Closed forms are used for `1` and `3`.
    - `:negate?` (`:auto`, `true` or `false`, default `:auto`): the sign. `true` multiplies by `-1`, `false` by `+1`. `:auto` (also for `nil`) uses `(-1)^ceil(beta / 2)`: negative for `beta` in `(0, 2]`, positive for `(2, 4]`, negative for `(4, 6]`, ...
    - `:scale` (double, default `1.0`): length scale (it only scales the kernel by `scale^(-beta)`).

  Returns a function `r -> double`. The value at `r = 0` is `0.0` for a positive `beta` (`1.0` for `beta = 0`, infinity for a negative `beta`). The function is symmetric. A missing `:beta` throws `ex-info`.

  Conditionally positive definite of order `ceil(beta / 2)` for a positive `beta` which is not an even integer (needs polynomial terms of degree `ceil(beta / 2) - 1`); an even integer `beta` is a polynomial in the coordinates and is not usable. For a negative `beta` the `:auto` sign is just the same formula (`+1` for `beta` in `(-2, 0)`, `-1` for `(-4, -2]`) and is not meant for it: the positive definite inverse powers (`beta` above `-d` in R^d, the Riesz kernels) need `:negate? false`.

  See also [[linear]], [[thin-plate-splines]], [[generalized-multiquadratic]]."
  ([] (radial-powers nil))
  ([params]
   (let [beta (required-double :radial-powers params :beta)
         scale (scale-of params)
         fact (infer-negate :radial-powers beta (get params :negate?))]
     (cond
       (m/== beta 1.0) (fn ^double [^double x] (m/* fact (m/abs (m// x scale))))
       (m/== beta 3.0) (fn ^double [^double x] (m/* fact (m/cb (m/abs (m// x scale)))))
       :else (fn ^double [^double x] (m/* fact (m/pow (m/abs (m// x scale)) beta)))))))

(defn thin-plate-splines
  "Thin plate splines kernel: `sign * rho^(2 beta) * log(rho)` with `rho = abs(r / scale)`.

  Parameters:

  - `params` (map):
    - `:beta` (double, required, positive): the exponent. `1` is the classical thin plate spline `rho^2 log(rho)`.
    - `:negate?` (`:auto`, `true` or `false`, default `:auto`): the sign. `true` multiplies by `-1`, `false` by `+1`. `:auto` (also for `nil`) uses `(-1)^(ceil(beta) + 1)`: positive for `beta = 1`, negative for `2`, positive for `3`, negative for `4`.
    - `:scale` (double, default `1.0`): length scale. For an integer `beta` a change of the scale adds only the term `rho^(2 beta) * log(scale)`, a polynomial which interpolation with polynomial terms of the needed degree does not see.

  Returns a function `r -> double`. The value at `r = 0` is `0.0` (the limit). The function is symmetric and negative for `rho < 1` when the sign is positive. A missing `:beta` throws `ex-info`; a non-positive `:beta` throws `ex-info` with `:property :beta`.

  Conditionally positive definite of order `beta + 1` for an integer `beta` (needs polynomial terms of degree `beta`); the sign was checked numerically for `beta` 1, 2 and 3 (and 4 as an observation). For a non-integer `beta` the rule is uncertain: random point sets gave signs which depend on the dimension, so the function is not a valid conditionally positive definite function in all dimensions and the `:auto` sign has no guarantee.

  See also [[radial-powers]], [[generalized-multiquadratic]]."
  ([] (thin-plate-splines nil))
  ([params]
   (let [beta (required-double :thin-plate-splines params :beta)
         scale (scale-of params)
         beta2 (m/* 2.0 beta)
         fact (infer-negate :tps beta (get params :negate?))]
     (when-not (m/pos? beta)
       (invalid-input! :thin-plate-splines :beta "must be positive" beta))
     (fn ^double [^double x]
       (let [ex (m/abs (m// x scale))]
         (if (m/zero? ex) 0.0 (m/* fact (m/pow ex beta2) (m/log ex))))))))

;; [2] p.44

(defn- whittaker-integration [^double alpha ^double k ^double beta]
  (let [k- (m/dec k)]
    (if (m/zero? alpha)
      (fn ^double [^double r]
        (let [f (if (m/zero? r)
                  (fn ^double [^double t] (m/exp (m/- (m/* t beta))))
                  (fn ^double [^double t] (m/* (m/pow (m/max 0.0 (m/- 1.0 (m/* r t))) k-)
                                            (m/exp (m/- (m/* t beta))))))]
          (calc/gk-quadrature f 0.0 ##Inf {:max-iters 1000 :abs 1.0e-6 :rel 1.0e-6})))
      (fn ^double [^double r]
        (let [f (if (m/zero? r)
                  (fn ^double [^double t] (m/* (m/pow t alpha)
                                            (m/exp (m/- (m/* t beta)))))
                  (fn ^double [^double t] (m/* (m/pow (m/max 0.0 (m/- 1.0 (m/* r t))) k-)
                                            (m/pow t alpha)
                                            (m/exp (m/- (m/* t beta))))))]
          (calc/gk-quadrature f 0.0 ##Inf {:max-iters 1000 :abs 1.0e-6 :rel 1.0e-6}))))))

(defn- whittaker-direct
  "Closed forms of the Whittaker kernel, `kind` 0 to 3 for alpha, k = (0, 2), (1, 2), (0, 3), (1, 3); `ex` is rho."
  ^double [^long kind ^double beta ^double ex]
  (let [e (m/exp (m/- (m// beta ex)))]
    (case kind
      0 (m// (m/+ (m/- beta ex) (m/* ex e)) beta)
      1 (let [ex2 (m/* 2.0 ex)]
          (m// (m/+ (m/- beta ex2) (m/* (m/+ beta ex2) e)) beta))
      2 (let [b (m/sq beta)
              ex2 (m/* 2.0 (m/sq ex))]
          (m// (m/+ (m/- b (m/* 2.0 beta ex)) ex2 (m/- (m/* ex2 e))) b))
      3 (let [b (m/sq beta)
              ex2 (m/* 6.0 (m/sq ex))
              bex (m/* 2.0 beta ex)]
          (m// (m/+ (m/- b (m/* 2.0 bex)) ex2 (m/- (m/* (m/+ bex ex2) e))) b)))))

(defn- whittaker-tail
  "Series in `u = beta / rho` of the closed forms (they cancel catastrophically for a small `u`)."
  ^double [^long kind ^double u]
  (let [n0 (case kind 0 2 1 3 2 3 3 4)
        shift (if (m/< kind 2) 1 2)
        fact0 (case kind 0 2.0 1 6.0 2 6.0 3 24.0)]
    ;; value = sum over n >= n0 of c(n) u^(n - shift) / n!
    (loop [n n0
           t (m// (m/pow u (m/- n0 shift)) fact0)
           acc 0.0]
      (let [sgn (if (m/even? n) 1.0 -1.0)
            c (case kind
                0 sgn
                1 (m/* -1.0 sgn (m/- n 2))
                2 (m/* -2.0 sgn)
                3 (m/* 2.0 sgn (m/- n 3)))
            term (m/* c t)
            acc (m/+ acc term)]
        (if (or (m/> n 60) (and (m/> n (m/inc n0)) (m/<= (m/abs term) (m/* 1.0e-18 (m/abs acc)))))
          acc
          (recur (m/inc n) (m// (m/* t u) (m/inc n)) acc))))))

(defn- whittaker-closed [^long kind ^double beta ^double scale]
  (fn ^double [^double x]
    (let [ex (m/abs (m// x scale))
          u (m// beta ex)]
      (if (m/< u 1.5)
        (whittaker-tail kind u)
        (whittaker-direct kind beta ex)))))

(defn whittaker
  "Whittaker kernel family: a positive mixture of truncated powers, `W(rho) / W(0)` with `W(rho) = integral over t >= 0 of (1 - rho t)_+^(k - 1) * t^alpha * exp(-beta t) dt` and `rho = r / scale`.

  Parameters:

  - `params` (map):
    - `:alpha` (double, required): the exponent of the power weight (`alpha > -1`).
    - `:k` (double, required, positive): the order of the truncated powers (the exponent is `k - 1`).
    - `:beta` (double, required, positive): the rate of the exponential weight.
    - `:scale` (double, default `1.0`): length scale.

  Closed forms are used for `alpha` equal to 0 or 1 together with `k` equal to 2 or 3 (whatever the numeric type of the arguments: `0`, `0.0` and `1` and `1.0` are the same). Every other combination uses numerical integration (adaptive Gauss-Kronrod, tolerance 1e-6), which is much slower and less accurate.

  Returns a function `r -> double`. The function is normalised: the value at `r = 0` is `1.0`. It is symmetric, positive and decreasing, and not compactly supported: for large `rho` it decays algebraically, like `rho^(-(alpha + 1))`, and underflows to `0.0` only for an enormous `rho`. A missing `:alpha`, `:k` or `:beta` throws `ex-info`; a non-positive `:beta` or `:k` or an `:alpha` not above `-1` throws `ex-info`.

  As a positive mixture of truncated powers with exponent `k - 1` (see [[truncated-power]]) it is positive definite in R^d at least for `d <= 2k - 3` (Askey, a sufficient condition: `d = 1` for `k = 2`).

  See also [[truncated-power]], [[gaussian]]."
  ([] (whittaker nil))
  ([params]
   (let [alpha (required-double :whittaker params :alpha)
         k (required-double :whittaker params :k)
         beta (required-double :whittaker params :beta)
         scale (scale-of params)]
     (when-not (m/> alpha -1.0)
       (invalid-input! :whittaker :alpha "must be greater than -1" alpha))
     (when-not (m/pos? k)
       (invalid-input! :whittaker :k "must be positive" k))
     (when-not (m/pos? beta)
       (invalid-input! :whittaker :beta "must be positive" beta))
     (cond
       (and (m/== alpha 0.0) (m/== k 2.0)) (whittaker-closed 0 beta scale)
       (and (m/== alpha 1.0) (m/== k 2.0)) (whittaker-closed 1 beta scale)
       (and (m/== alpha 0.0) (m/== k 3.0)) (whittaker-closed 2 beta scale)
       (and (m/== alpha 1.0) (m/== k 3.0)) (whittaker-closed 3 beta scale)
       :else
       (let [w (whittaker-integration alpha k beta)
             w0 (double (w 0.0))]
         (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                  (m// ^double (w ex) w0))))))))

;; [2] p.131

(defn shifted-surface-splines
  "Shifted surface splines kernel family (Dyn, Levin and Rippa; in the form of Luh, equations (1) and (2) with `beta` in place of lambda).

  With `t = beta - s/2` and `rho = r / scale` the kernel is, for an odd `s`, `(-1)^ceil(t) * (1 + rho^2)^t` and, for an even `s`, `(-1)^(1 + t) * (1 + rho^2)^t * log(sqrt(1 + rho^2))`.

  Parameters:

  - `params` (map):
    - `:s` (long, required, positive): the space dimension. A whole double such as `3.0` is accepted.
    - `:beta` (double, required): the exponent, must be greater than `s / 2` (`t > 0`). For an even `s` it must be an integer (so that `t` is a positive integer), otherwise the sign is not defined.
    - `:scale` (double, default `1.0`): length scale; it is the shape parameter `c` of Luh (his function is `(r^2 + c^2)^t`, the same as here up to the constant factor `scale^(2 t)`).

  Returns a function `r -> double`. The value at `r = 0` is the sign for an odd `s` and `0.0` for an even `s`. The function is symmetric and smooth (including at 0, the logarithm is of `1 + rho^2`), and grows without bound (it overflows to infinity for a huge `r`). For example `s = 3`, `beta = 2.5`, `r = 1` gives `-2.0` (`t = 1`, sign `-1`). A missing `:s` or `:beta`, a non-integer `:s`, `:s` below 1, `:beta` not above `s / 2` or a non-integer `:beta` for an even `s` throws `ex-info`.

  Conditionally positive definite of order `ceil(t)` for an odd `s` (the generalized multiquadric with exponent `t`, see [[generalized-multiquadratic]]) and of order `t + 1` for an even `s`; use polynomial terms of that order minus one degree.

  See also [[generalized-multiquadratic]], [[thin-plate-splines]]."
  ([] (shifted-surface-splines nil))
  ([params]
   (let [s (required-long :shifted-surface-splines params :s)
         beta (required-double :shifted-surface-splines params :beta)
         scale (scale-of params)
         t (m/- beta (m/* 0.5 s))
         ct (long (m/ceil t))]
     (when (m/< s 1)
       (invalid-input! :shifted-surface-splines :s "must be positive" s))
     (when-not (m/pos? t)
       (invalid-input! :shifted-surface-splines :beta "must be greater than half of s" beta))
     (if (m/odd? s)
       (let [fact (if (m/even? ct) 1.0 -1.0)]
         (fn ^double [^double x]
           (m/* fact (m/pow (m/inc (m/sq (m// x scale))) t))))
       (let [fact (if (m/even? ct) -1.0 1.0)]
         (when-not (m/== t (m/floor t))
           (invalid-input! :shifted-surface-splines :beta "must be an integer for an even s" beta))
         (fn ^double [^double x]
           (let [r2 (m/sq (m// x scale))]
             (m/* fact 0.5 (m/pow (m/inc r2) t) (m/log1p r2)))))))))

;; [2] p.88
;; normalized to have 1.0 at 0.0

(defn wendland
  "Wendland kernel family: `I^k (1 - rho)_+^l` normalised to `1` at `0`, with `l = floor(s / 2) + k + 1`, `rho = r / scale` and `I f(rho) = integral from rho to 1 of t f(t) dt` (Wendland and Fasshauer: the function of smoothness `2k` which is positive definite in R^s).

  Parameters:

  - `params` (map):
    - `:s` (double, required, at least `1`): the space dimension. Only `floor(s / 2)` matters (`s = 3` and `s = 2` give the same function).
    - `:k` (long, required, from `0` to `4`): the smoothness, the function is `2k` times continuously differentiable. A whole double such as `2.0` is accepted.
    - `:scale` (double, default `1.0`): the support radius.

  Returns a function `r -> double`. On its support the function is a polynomial of degree `floor(s / 2) + 3k + 1`; for example `k = 1`, `s = 3` is `(1 - rho)^4 (4 rho + 1)` (value `0.1875` at `rho = 0.5`). The value at `r = 0` is `1.0`; it is symmetric, decreasing, equal to `0.0` at `abs(r) = scale` and beyond. A missing `:s` or `:k`, `:s` below 1, a negative or non-integer `:k` or `:k` above 4 throws `ex-info`.

  Strictly positive definite in R^d for every `d <= s`.

  See also [[wu]], [[gneiting]], [[truncated-power]], [[matern]]."
  ([] (wendland nil))
  ([params]
   (let [s (required-double :wendland params :s)
         k (required-long :wendland params :k)
         scale (scale-of params)]
     (when (m/< s 1.0)
       (invalid-input! :wendland :s "must be at least 1" s))
     (when-not (m/<= 0 k 4)
       (invalid-input! :wendland :k "must be an integer from 0 to 4" k))
     (let [l (m/+ (m/floor (m/* s 0.5)) (double k) 1.0)
           p (m/+ l k)]
       (case (int k)
         0 (fn ^double [^double x] (let [r (m/abs (m// x scale))]
                                    (if (m/< r 1.0) (m/pow (m/- 1.0 r) l) 0.0)))
         1 (fn ^double [^double x] (let [r (m/abs (m// x scale))]
                                    (if (m/< r 1.0)
                                      (m/* (m/pow (m/- 1.0 r) p) (m/inc (m/* p r))) 0.0)))
         2 (let [f2 (poly/mevalpoly l 3.0 4.0 1.0)
                 f1 (poly/mevalpoly l 6.0 3.0)]
             (fn ^double [^double x] (let [r (m/abs (m// x scale))]
                                      (if (m/< r 1.0)
                                        (m// (m/* (m/pow (m/- 1.0 r) p)
                                              (poly/mevalpoly r 3.0 f1 f2)) 3.0) 0.0))))
         3 (let [f3 (poly/mevalpoly l 15.0 23.0 9.0 1.0)
                 f2 (poly/mevalpoly l 45.0 36.0 6.0)
                 f1 (poly/mevalpoly l 45.0 15.0)]
             (fn ^double [^double x] (let [r (m/abs (m// x scale))]
                                      (if (m/< r 1.0)
                                        (m// (m/* (m/pow (m/- 1.0 r) p)
                                              (poly/mevalpoly r 15.0 f1 f2 f3)) 15.0) 0.0))))
         4 (let [f4 (poly/mevalpoly l 105.0 176.0 86.0 16.0 1.0)
                 f3 (m/* 5.0 (poly/mevalpoly l 84.0 85.0 24.0 2.0))
                 f2 (m/* 45.0 (poly/mevalpoly l 14.0 8.0 1.0))
                 f1 (m/* 105.0 (poly/mevalpoly l 4.0 1.0))]
             (fn ^double [^double x] (let [r (m/abs (m// x scale))]
                                      (if (m/< r 1.0)
                                        (m// (m/* (m/pow (m/- 1.0 r) p)
                                              (poly/mevalpoly r 105.0 f1 f2 f3 f4)) 105.0) 0.0)))))))))

;; [2] p. 91

(defn gneiting
  "Gneiting kernel family: the oscillating compactly supported function `(1 - rho)_+^l * (1 + l rho - (l + 1)(l + s + 2) / s * rho^2)` with `rho = r / scale` (the turning bands operator applied to a Wendland function).

  Parameters:

  - `params` (map):
    - `:s` (double, required, positive): the space dimension.
    - `:l` (double, required): the exponent.
    - `:scale` (double, default `1.0`): the support radius.

  Returns a function `r -> double`. The value at `r = 0` is `1.0`; it is symmetric, takes negative values (it oscillates) and is `0.0` at `abs(r) = scale` and beyond. A missing `:s` or `:l` or a non-positive `:s` throws `ex-info` (the formula divides by `s`).

  Positive definite in R^s if and only if `l >= (s + 5) / 2`.

  See also [[wendland]], [[wu]], [[truncated-power]]."
  ([] (gneiting nil))
  ([params]
   (let [s (required-double :gneiting params :s)
         l (required-double :gneiting params :l)
         scale (scale-of params)]
     (when-not (m/pos? s)
       (invalid-input! :gneiting :s "must be positive" s))
     (let [f (m/- (m// (m/* (m/inc l) (m/+ l s 2.0)) s))]
       (fn ^double [^double x]
         (let [ex (m/abs (m// x scale))]
           (if (m/< ex 1.0)
             (m/* (m/pow (m/- 1.0 ex) l)
                (poly/mevalpoly ex 1.0 l f))
             0.0)))))))

;;

;; https://www.researchgate.net/publication/226871599_Compactly_supported_positive_definite_radial_function
;; https://clojurians.zulipchat.com/#narrow/stream/151924-data-science/topic/.E2.9C.94.20Calculus.20help

(defn- wu-integration [^double l]
  (fn ^double [^double r]
    (let [rr (m/* 2.0 r)]
      (calc/gk-quadrature (fn [^double t] (m/* (m/pow (m/- 1.0 (m/* t t)) l)
                                            (m/pow (m/- 1.0 (m/sq (m/- rr t))) l)))
                          (m/dec rr) 1.0))))

(defn- wu-differentiation [f ^long k]
  (if (m/zero? k)
    f
    (let [Dk (finite/derivative f) ;; derivative
          fDk (fn ^double [^double r ^double h] (let [hr (m/+ r h)]
                                                 (m// ^double (Dk hr) (m/- hr)))) ;; D operator
          eDk (finite/extrapolate fDk {:power 2.0})] ;; extrapolation
      (recur (fn ^double [^double r]
               (let [r (m/abs r)]
                 (if (m/< r 1.0e-3)
                   (eDk r) ;; extrapolate at 0 for small `r`
                   (fDk r 0.0))))
             (m/dec k)))))

(defn- wu-closed-form
  "Closed form of the Wu function for the integer `l` (0 to 3) and `k` (0 to `l`) or nil."
  [^long l ^long k ^double scale]
  (case (m/+ (m/* 10 l) k)
    0 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                               (if (m/< ex 1.0) (m/- 1.0 ex) 0.0)))
    10 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m/* (m/cb (m/- 1.0 ex))
                                     (poly/mevalpoly ex 1.0 3.0 1.0))
                                  0.0)))
    11 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m/* 0.5 (m/sq (m/- 1.0 ex))
                                     (poly/mevalpoly ex 2.0 1.0))
                                  0.0)))
    20 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m/* (m/fpow (m/- 1.0 ex) 5)
                                     (poly/mevalpoly ex 1.0 5.0 9.0 5.0 1.0))
                                  0.0)))
    21 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m/* 0.25 (m/fpow (m/- 1.0 ex) 4)
                                     (poly/mevalpoly ex 4.0 16.0 12.0 3.0))
                                  0.0)))
    22 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m/* 0.125 (m/cb (m/- 1.0 ex))
                                     (poly/mevalpoly ex 8.0 9.0 3.0))
                                  0.0)))
    30 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m// (m/* (m/fpow (m/- 1.0 ex) 7)
                                        (poly/mevalpoly ex 5.0 35.0 101.0 147.0 101.0 35.0 5.0)) 5.0)
                                  0.0)))
    31 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m// (m/* (m/fpow (m/- 1.0 ex) 6)
                                        (poly/mevalpoly ex 6.0 36.0 82.0 72.0 30.0 5.0)) 6.0)
                                  0.0)))
    32 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m// (m/* (m/fpow (m/- 1.0 ex) 5)
                                        (poly/mevalpoly ex 8.0 40.0 48.0 25.0 5.0)) 8.0)
                                  0.0)))
    33 (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                (if (m/< ex 1.0)
                                  (m// (m/* (m/fpow (m/- 1.0 ex) 4)
                                        (poly/mevalpoly ex 16.0 29.0 20.0 5.0)) 16.0)
                                  0.0)))
    nil))

(defn wu
  "Wu kernel family: the compactly supported functions `psi_{l,k} = D^k psi_l` normalised to `1` at `0`, where `psi_l(rho)` is the convolution of `(1 - t^2)_+^l` with itself evaluated at `2 rho` (the integral over t of `(1 - t^2)_+^l * (1 - (2 rho - t)^2)_+^l`), `D f(rho) = -f'(rho) / rho` and `rho = r / scale`.

  Parameters:

  - `params` (map):
    - `:l` (double, required, non-negative): the exponent of the base function.
    - `:k` (long, required, non-negative): the number of `D` steps. A whole double such as `2.0` is accepted (the presets do that).
    - `:scale` (double, default `1.0`): the support radius.

  Closed form polynomials are used for an integer `l` from 0 to 3 and `k` from 0 to `l` (whatever the numeric type of the arguments: `2`, `2.0` are the same). Their values equal the exact rational ones to about 1e-14. For other combinations (a larger or non-integer `l`, or `k` above `l`) the base function is integrated and differentiated numerically (Gauss-Kronrod and finite differences): this is accurate for `k = 0`, but unstable for `k >= 2` (errors grow quickly with `k`, up to garbage for `k >= 3`), so prefer the closed forms.

  Returns a function `r -> double`. The value at `r = 0` is `1.0`; it is symmetric, equal to `0.0` at `abs(r) = scale` and beyond. A missing `:l` or `:k`, a negative `:l` or `:k` or a non-integer `:k` throws `ex-info`.

  Positive definite in R^d for `d <= 2k + 1` (the closed forms).

  See also [[wendland]], [[gneiting]], [[truncated-power]]."
  ([] (wu nil))
  ([params]
   (let [l (required-double :wu params :l)
         k (required-long :wu params :k)
         scale (scale-of params)]
     (when (m/neg? l)
       (invalid-input! :wu :l "must be non-negative" l))
     (when (m/neg? k)
       (invalid-input! :wu :k "must be non-negative" k))
     (or (when (and (m/== l (m/floor l)) (m/<= l 3.0) (m/<= k l))
           (wu-closed-form (long l) k scale))
         (let [integr (wu-integration l)
               phi (wu-differentiation integr k)
               fact (m// ^double (phi 0.0))]
           (fn ^double [^double x] (let [ex (m/abs (m// x scale))]
                                    (if (m/<= ex 1.0) (m/* fact ^double (phi ex)) 0.0))))))))

