(ns fastmath.kernel.density
  "One dimensional kernel density estimation: kernel functions, bandwidth selectors and the density estimators.

  A kernel density estimate of the sample `x_1 ... x_n` is `f(x) = 1/(n h) * sum K((x - x_i) / h)` for a kernel `K` (a symmetric probability density) and a bandwidth `h > 0`.

  Kernels (standard forms, the compact ones have support `[-1, 1]`; `gaussian` has standard deviation `1`, which is not the convention of `density` in R, where the bandwidth is the standard deviation of the kernel): [[uniform]], [[triangular]], [[epanechnikov]], [[quartic]], [[triweight]], [[tricube]], [[gaussian]], [[cosine]], [[logistic]], [[sigmoid]], [[silverman]], [[laplace]], [[wigner]], [[cauchy]]. Their constants (roughness, second moment, canonical bandwidth, efficiency, radius) are in [[kde-data]].

  Estimators: [[kernel-density]] (the function), [[kernel-density+]] (the function with the fitting details), [[kernel-density-ci]] (with pointwise confidence intervals).

  Bandwidth selection: [[bandwidth]] with the rules of thumb `:nrd`, `:nrd0`, `:nrd-adjust`, the plug-in selectors `:sj`, `:sj-dpi` (as `bw.SJ` in R) and the cross-validation selectors `:lscv`, `:lcv`, `:rlcv`.

  Evaluation is exact: every data point within `1.01 * h * radius` of `x` contributes, no point outside does. Option `:binned?` replaces the data by a linear binning on a regular grid for large samples."
  (:require [fastmath.core :as m]
            [fastmath.optimization.acm :as acm])
  (:import [org.apache.commons.math3.distribution NormalDistribution]
           [clojure.lang IFn$DD]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- invalid-input!
  "Throws `ex-info` for an invalid `property` of the density `function`."
  [function property message value expected]
  (throw (ex-info (str "kernel density " (name function) ", " (name property) ": " message)
                  {:type :invalid-input
                   :function function
                   :property property
                   :value value
                   :expected expected})))

;; kernels

(defn uniform
  "Uniform kernel, `0.5` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) 0.5 0.0))

(def ^{:const true :private true :tag 'double} gaussian-factor (m// (m/sqrt m/TWO_PI)))

(defn gaussian
  "Gaussian kernel, the standard normal density.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (m/* gaussian-factor (m/exp (m/* -0.5 x x))))

(defn triangular
  "Triangular kernel, `1 - abs(x)` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0) (m/- 1.0 absx) 0.0)))

(defn epanechnikov
  "Epanechnikov kernel, `0.75 * (1 - x^2)` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0) (m/* 0.75 (m/- 1.0 absx) (m/+ 1.0 absx)) 0.0)))

(defn quartic
  "Quartic (biweight) kernel, `15/16 * (1 - x^2)^2` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0)
      (m/* 0.9375 (m/sq (m/* (m/- 1.0 absx) (m/+ 1.0 absx))))
      0.0)))

(defn triweight
  "Triweight kernel, `35/32 * (1 - x^2)^3` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0)
      (let [v (m/* (m/- 1.0 absx) (m/+ 1.0 absx))]
        (m/* 1.09375 v v v))
      0.0)))

(defn tricube
  "Tricube kernel, `70/81 * (1 - abs(x)^3)^3` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0)
      (let [v (m/- 1.0 (m/* absx absx absx))]
        (m/* 0.8641975308641975 v v v))
      0.0)))

(defn cosine
  "Cosine kernel, `pi/4 * cos(pi/2 * x)` on `[-1, 1]`, `0.0` outside.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (if (m/<= (m/abs x) 1.0) (m/* m/QUARTER_PI (m/cos (m/* m/HALF_PI x))) 0.0))

(defn logistic
  "Logistic kernel, the density of the standard logistic distribution `1 / (2 + e^x + e^-x)`.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`; `0.0` for very large `abs(x)`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (m// (m/+ 2.0 (m/exp x) (m/exp (m/- x)))))

(defn sigmoid
  "Sigmoid kernel, `2/pi / (e^x + e^-x)`.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`; `0.0` for very large `abs(x)`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (m/* m/TWO_INV_PI (m// (m/+ (m/exp x) (m/exp (m/- x))))))

(defn silverman
  "Silverman kernel, `0.5 * e^-u * sin(u + pi/4)` for `u = abs(x) / sqrt(2)`.

  The kernel takes negative values (a higher order kernel), so a density estimate built with it can be negative.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`; `0.0` at `+-Inf`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [xx (m// (m/abs x) m/SQRT2)]
    (if (m/inf? xx)
      0.0
      (m/* 0.5 (m/exp (m/- xx)) (m/sin (m/+ xx m/QUARTER_PI))))))

(defn laplace
  "Laplace kernel, `0.5 * e^-abs(x)`, experimental.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (m/* 0.5 (m/exp (m/- (m/abs x)))))

(defn wigner
  "Wigner (semicircle) kernel, `2/pi * sqrt(1 - x^2)` on `[-1, 1]`, `0.0` outside, experimental.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (m/<= absx 1.0) (m/* m/TWO_INV_PI (m/sqrt (m/* (m/- 1.0 absx) (m/+ 1.0 absx)))) 0.0)))

(defn cauchy
  "Cauchy kernel, `1 / (pi/2 * (1 + 4 x^2))`, experimental.

  The density has no finite variance, so the canonical bandwidth is not defined, the exact evaluation sums every data point and the `:nrd-adjust` bandwidth is not available.

  Parameters: `x` - a point.

  Returns the value of the kernel at `x`.

  See also [[kde-data]], [[kernel-density]]."
  ^double [^double x]
  (m// (m/* m/HALF_PI (m/inc (m/sq (m/* 2.0 x))))))

;; Canonical kernels for density estimation
;; Marron, Nolan

;; factor used to scale nrd kernel size
(def ^{:private true :const true :tag 'double} GAUSS-FACT (m// (m/pow (m/* 0.5 (m// m/SQRTPI)) 0.2)))

(defn- calc-additional-consts
  [{:keys [^double k2 ^double x2k radius] :as d}]
  (let [delta0 (m/* (m/pow k2 0.2)
                    (m/pow x2k -0.4))]
    (assoc d :delta0 delta0
           :nrd-scale (m/* delta0 GAUSS-FACT)
           :efficiency (m/* k2 (m/sqrt x2k))
           :extent (get d :extent radius))))

(def ^{:doc "Constants of the predefined kernels.

  Keys - kernel keywords, values - maps with:

  * `:kernel` - the kernel function
  * `:k2` - integral of `K(x)^2` (the roughness)
  * `:x2k` - integral of `x^2 K(x)` (the second moment); `0.0` for `:silverman` and `##Inf` for `:cauchy`
  * `:delta0` - the canonical bandwidth `k2^0.2 * x2k^-0.4`; `##Inf` for `:silverman` and `0.0` for `:cauchy`, which means undefined (the bandwidth adjusted to the kernel is not available)
  * `:nrd-scale` - `delta0` relative to the gaussian kernel, the factor of the `:nrd-adjust` bandwidth (`1.0` for `:gaussian`); undefined as `:delta0`
  * `:efficiency` - the efficiency constant `k2 * sqrt(x2k)`
  * `:radius` - the distance from the data point beyond which the kernel is not used in the evaluation: `1.0` for the compact kernels, a distance with `abs(K) < 1.0e-10` for the other ones, `##Inf` for `:cauchy` (every data point is summed)
  * `:extent` - the distance beyond which the density is treated as zero when the support of the estimate (`:mn`, `:mx`) is inferred; the same as `:radius` except `:cauchy` (`1000.0`, the mass cut off is about `4.0e-4`)

  A custom kernel function has none of these constants; `:radius` (default `##Inf`) and `:extent` (default `5.0`) can be passed as options to [[kernel-density+]] and [[bandwidth]].

  See also [[kernel-density]], [[bandwidth]]."}
  kde-data
  (let [in {:uniform {:k2 0.5 :x2k m/THIRD :kernel uniform :radius 1.0}
            :triangular {:k2 m/TWO_THIRD :x2k m/SIXTH :kernel triangular :radius 1.0}
            :epanechnikov {:k2 0.6 :x2k 0.2 :kernel epanechnikov :radius 1.0}
            :quartic {:k2 (m// 5.0 7.0) :x2k (m// 7.0) :kernel quartic :radius 1.0}
            :triweight {:k2 (m// 350.0 429.0) :x2k (m// 9.0) :kernel triweight :radius 1.0}
            :tricube {:k2 (m// 175.0 247.0) :x2k (m// 35.0 243.0) :kernel tricube :radius 1.0}
            :gaussian {:k2 (m/* 0.5 (m// m/SQRTPI)) :x2k 1.0 :kernel gaussian :radius 6.65}
            :cosine {:k2 (m/* 0.0625 m/PI m/PI) :x2k (m/- 1.0 (m// 8.0 (m/* m/PI m/PI)))
                     :kernel cosine :radius 1.0}
            :logistic {:k2 m/SIXTH :x2k (m// (m/* m/PI m/PI) 3.0) :kernel logistic :radius 23.1}
            :sigmoid {:k2 (m// 2.0 (m/* m/PI m/PI)) :x2k (m// (m/* m/PI m/PI) 4.0) :kernel sigmoid :radius 22.6}
            :silverman {:k2 (m/* 0.0625 3.0 m/SQRT2) :x2k 0.0 :kernel silverman :radius 29.7}
            :wigner {:k2 (m// 16.0 (m/* 3.0 m/PI m/PI)) :x2k 0.25 :kernel wigner :radius 1.0}
            :cauchy {:k2 m/M_1_PI :x2k ##Inf :kernel cauchy :radius ##Inf :extent 1000.0}
            :laplace {:k2 0.25 :x2k 2.0 :kernel laplace :radius 22.4}}]
    (reduce (fn [m k] (update m k calc-additional-consts)) in (keys in))))

;; validation

(def ^:private kernel-keywords (vec (sort (keys kde-data))))

(def ^:private bandwidth-keywords [:nrd :nrd0 :nrd-adjust :sj :sj-dpi :rlcv :lcv :lscv])

(def ^{:private true :const true :tag 'long} MAX-BINNED-NODES 16777216)

(defn- finite-positive?
  [v]
  (and (number? v) (m/valid-double? (double v)) (m/pos? (double v))))

(defn- check-binned!
  "Returns `binned?` validated: `nil`, `false`, `true` or a finite number greater than zero."
  [function binned?]
  (cond
    (or (nil? binned?) (false? binned?)) false
    (true? binned?) true
    (finite-positive? binned?) binned?
    :else (invalid-input! function :binned? "must be nil, false, true or a finite number greater than 0" binned?
                          "nil, false, true or a number > 0")))

(defn- sorted-copy
  "Returns a sorted copy of finite `data` as `double[]` (the input is not modified)."
  ^doubles [function data]
  (let [^doubles a (try (m/seq->double-array data)
                        (catch Exception _
                          (invalid-input! function :data "must be a sequence of numbers" data "a sequence of finite numbers")))]
    (when (or (nil? a) (m/zero? (alength a)))
      (invalid-input! function :data "must not be empty" data "a non-empty sequence of finite numbers"))
    (let [^doubles c (java.util.Arrays/copyOf a (alength a))]
      (dotimes [i (alength c)]
        (when-not (m/valid-double? (Array/aget c i))
          (invalid-input! function :data (str "contains a non finite value at index " i) (Array/aget c i)
                          "finite numbers")))
      (java.util.Arrays/sort c)
      c)))

(defn- ->prim-kernel
  ^IFn$DD [f]
  (if (instance? IFn$DD f)
    f
    (reify IFn$DD (invokePrim [_ x] (double (f x))))))

(defn- ->kdata
  "Validates the input and returns the preprocessed map.

  Keys: `:data` sorted `double[]`, `:len`, `:kernel-name`, `:kernel` (`IFn$DD`), `:radius`, `:extent`, `:binned?`."
  [function data kernel binned? {:keys [radius extent]}]
  (let [known (when (keyword? kernel) (get kde-data kernel))]
    (when-not (or known (fn? kernel))
      (invalid-input! function :kernel (str "unknown kernel, expected a function or one of " kernel-keywords)
                      kernel "a kernel function or one of the keywords of kde-data"))
    (let [data (sorted-copy function data)
          radius (cond
                   known (:radius known)
                   (nil? radius) ##Inf
                   (and (number? radius) (not (m/nan? (double radius))) (m/pos? (double radius))) (double radius)
                   :else (invalid-input! function :radius "must be a number greater than 0 or ##Inf" radius "> 0"))
          extent (cond
                   known (:extent known)
                   (nil? extent) 5.0
                   (finite-positive? extent) (double extent)
                   :else (invalid-input! function :extent "must be a finite number greater than 0" extent "> 0"))]
      {:data data
       :len (alength data)
       :kernel-name kernel
       :kernel (->prim-kernel (if known (:kernel known) kernel))
       :radius (double radius)
       :extent (double extent)
       :binned? (check-binned! function binned?)})))

(defn- check-bandwidth!
  "Returns `h` as a double if it is finite and greater than 0, throws otherwise."
  ^double [function h]
  (if (finite-positive? h)
    (double h)
    (invalid-input! function :bandwidth "must be finite and greater than 0" h "a finite number > 0")))

;; evaluation

(defn- lower-bound
  "Index of the first element of sorted `a` not less than `v`."
  ^long [^doubles a ^double v]
  (loop [lo 0
         hi (alength a)]
    (if (m/< lo hi)
      (let [mid (m/unsigned-bit-shift-right (m/long-add lo hi) 1)]
        (if (m/< (Array/aget a mid) v)
          (recur (m/inc mid) hi)
          (recur lo mid)))
      lo)))

(defn- upper-bound
  "Index of the first element of sorted `a` greater than `v`."
  ^long [^doubles a ^double v]
  (loop [lo 0
         hi (alength a)]
    (if (m/< lo hi)
      (let [mid (m/unsigned-bit-shift-right (m/long-add lo hi) 1)]
        (if (m/<= (Array/aget a mid) v)
          (recur (m/inc mid) hi)
          (recur lo mid)))
      lo)))

(defmacro ^:private window-sum
  "Sum of the kernel over the data points within `span` of `x`, the point of index `skip` is left out."
  [data kernel h span x skip]
  `(let [^"[D" data# ~data
         ^IFn$DD kernel# ~kernel
         h# (double ~h)
         span# (double ~span)
         x# (double ~x)
         skip# (long ~skip)
         to# (upper-bound data# (m/+ x# span#))]
     (double
      (loop [i# (lower-bound data# (m/- x# span#))
             s# 0.0]
        (if (m/< i# to#)
          (recur (m/inc i#)
                 (if (m/== i# skip#)
                   s#
                   (m/+ s# (.invokePrim kernel# (m// (m/- x# (Array/aget data# i#)) h#)))))
          s#)))))

(defn- binned-nodes
  "Linear binning of sorted `data` on the regular grid `g0 + j * delta`, returns the node weights."
  ^doubles [^doubles data ^double g0 ^double delta ^long nodes]
  (let [counts (double-array nodes)
        last-node (m/long-dec nodes)]
    (dotimes [i (alength data)]
      (let [t (m// (m/- (Array/aget data i) g0) delta)
            j (m/long-min last-node (long (m/floor t)))
            frac (m/- t (double j))]
        (Array/add counts j (m/- 1.0 frac))
        (when (m/< j last-node)
          (Array/add counts (m/inc j) frac))))
    counts))

(defn- kde-
  "Builds the density estimator for the preprocessed `kdata` and bandwidth `h`.

  Returns a map with `:kde`, `:factor` (`1/(n h)`), `:h`, `:mn` and `:mx`."
  [{:keys [^doubles data ^long len ^IFn$DD kernel ^double radius ^double extent binned?]} ^double h]
  (let [span (m/* 1.01 h radius)
        factor (m// (m/* (double len) h))
        mn (Array/aget data 0)
        mx (Array/aget data (m/dec len))
        kde-fn (if binned?
                 (let [divisor (if (number? binned?) (double binned?) 5.0)
                       delta (m// h divisor)
                       nodes-d (m/inc (m/ceil (m// (m/- mx mn) delta)))]
                   (when (m/> nodes-d (double MAX-BINNED-NODES))
                     (invalid-input! :kernel-density :binned?
                                     (str "the grid would have " (long nodes-d) " nodes, more than " MAX-BINNED-NODES
                                          ", use a larger bandwidth or divisor, or no binning")
                                     binned? (str "at most " MAX-BINNED-NODES " grid nodes")))
                   (let [^doubles counts (binned-nodes data mn delta (long nodes-d))
                         last-node (double (m/long-dec (alength counts)))]
                     (fn ^double [^double x]
                       (if (m/nan? x)
                         ##NaN
                         (let [lo (long (m/max 0.0 (m/ceil (m// (m/- (m/- x span) mn) delta))))
                               hi (long (m/min last-node (m/floor (m// (m/- (m/+ x span) mn) delta))))]
                           (m/* factor
                                (double
                                 (loop [j lo
                                        s 0.0]
                                   (if (m/<= j hi)
                                     (let [c (Array/aget counts j)]
                                       (recur (m/inc j)
                                              (if (m/pos? c)
                                                (m/+ s (m/* c (.invokePrim kernel (m// (m/- x (m/+ mn (m/* (double j) delta))) h))))
                                                s)))
                                     s)))))))))
                 (fn ^double [^double x]
                   (if (m/nan? x)
                     ##NaN
                     (m/* factor (window-sum data kernel h span x -1)))))
        reach (m/* 1.01 h extent)]
    {:kde kde-fn
     :factor factor
     :h h
     :mn (m/- mn reach)
     :mx (m/+ mx reach)}))

;; bandwidth estimation

(defn- mean-sd
  "Mean and the sample standard deviation of `data`, two pass."
  [^doubles data]
  (let [n (alength data)
        mean (m// (double (loop [i 0 s 0.0] (if (m/< i n) (recur (m/inc i) (m/+ s (Array/aget data i))) s))) (double n))
        ss (double
            (loop [i 0 s 0.0]
              (if (m/< i n)
                (recur (m/inc i) (m/+ s (m/sq (m/- (Array/aget data i) mean))))
                s)))]
    [mean (m/sqrt (m// ss (double (m/long-dec n))))]))

(defn- quantile-7
  "Quantile `p` of sorted `data`, linear interpolation of the order statistics (type 7, the default of R)."
  ^double [^doubles data ^double p]
  (let [pos (m/* p (m/dec (alength data)))
        lo (long (m/floor pos))
        hi (m/min (m/inc lo) (m/dec (alength data)))
        frac (m/- pos lo)
        a (Array/aget data lo)
        b (Array/aget data hi)]
    (m/+ a (m/* frac (m/- b a)))))

(defn- iqr-7
  ^double [^doubles data]
  (m/- (quantile-7 data 0.75) (quantile-7 data 0.25)))

(defn- check-sample-size!
  [function ^long len what]
  (when (m/< len 2)
    (invalid-input! function :data (str what " needs at least 2 data points") len "at least 2 data points")))

;; nrd0 0.9
;; nrd 1.06
(defn- nrd
  "Rule of thumb bandwidth `scale * s * n^-0.2` with `s = min(sd, IQR/1.34)` and the fallback chain of `bw.nrd0` in R:
  `sd`, then `abs(x_1)`, then `1.0` when the spread is zero."
  ^double [{:keys [^doubles data ^long len]} ^double scale]
  (check-sample-size! :bandwidth len "the rule of thumb")
  (let [[_ ^double sd] (mean-sd data)
        lo (m/min sd (m// (iqr-7 data) 1.34))
        s (double
           (cond
             (m/pos? lo) lo
             (m/pos? sd) sd
             (m/pos? (m/abs (Array/aget data 0))) (m/abs (Array/aget data 0))
             :else 1.0))]
    (m/* scale s (m/pow len -0.2))))

(defn- kernel-scale
  "Canonical bandwidth factor of the kernel relative to gaussian or throws when it is not defined."
  ^double [{:keys [kernel-name]} selector]
  (let [nrd-scale (:nrd-scale (when (keyword? kernel-name) (get kde-data kernel-name)))]
    (if (and nrd-scale (finite-positive? nrd-scale))
      (double nrd-scale)
      (invalid-input! :bandwidth :bandwidth
                      (str selector " needs a kernel with a canonical bandwidth, not available for "
                           (if (keyword? kernel-name) kernel-name "a custom kernel function"))
                      kernel-name "a predefined kernel other than :silverman and :cauchy"))))

;; solve-the-equation and direct plug-in selectors (Sheather, Jones 1991), port of bw.SJ from R

(def ^{:private true :const true :tag 'double} SJ-DELMAX 1000.0)

(def ^{:private true :const true :tag 'long} SJ-MAX-BINS 10000)

(defn- check-bins!
  "Returns the number of bins of `:sj` and `:sj-dpi` as a long: `nil` is the default `1000`, otherwise an integer valued number in `[2, 10000]` (the cost of the pair counts is quadratic in the number of bins)."
  ^long [bins]
  (cond
    (nil? bins) 1000
    (and (number? bins)
         (let [d (double bins)]
           (and (m/valid-double? d) (m/== d (m/floor d)) (m/<= 2.0 d (double SJ-MAX-BINS)))))
    (long (double bins))
    :else (invalid-input! :bandwidth :bins (str "must be an integer in [2, " SJ-MAX-BINS "]") bins
                          (str "integer in [2, " SJ-MAX-BINS "]"))))

(defn- sj-pair-counts
  "Binned pair counts of the sorted `data`: returns `[d counts]` with the bin width `d` and the `bins` counts of pairs
  by the distance in bins (the `bw_pair_cnts` of R: exact pairs for `n <= bins / 2`, pairs of bin counts otherwise)."
  [^doubles data ^long bins]
  (let [n (alength data)
        mn (Array/aget data 0)
        mx (Array/aget data (m/dec n))
        d (m// (m/* (m/- mx mn) 1.01) (double bins))]
    (when-not (m/pos? d)
      (invalid-input! :bandwidth :data "data are constant in bandwidth calculation" mn "non-constant data"))
    (let [counts (double-array bins)]
      (if (m/> n (m// (double bins) 2.0))
        (let [idx (long-array n)
              _ (dotimes [i n]
                  (let [x (Array/aget data i)]
                    (aset idx i (long (m/* (m/signum x) (m/floor (m// (m/abs x) d)))))))
              mn-idx (reduce m/min (seq idx))
              w (double-array bins)]
          (dotimes [i n]
            (let [b (m/- (aget idx i) (long mn-idx))]
              (when (m/< b bins)
                (Array/add w b 1.0))))
          (dotimes [ii bins]
            (let [wi (Array/aget w ii)]
              (Array/add counts 0 (m/* wi (m/dec wi)))
              (dotimes [jj ii]
                (Array/add counts (m/- ii jj) (m/* wi (Array/aget w jj))))))
          (Array/mult counts 0 0.5))
        (let [idx (long-array n)]
          (dotimes [i n]
            (aset idx i (long (m// (Array/aget data i) d))))
          (dotimes [i n]
            (dotimes [j i]
              (Array/inc counts (int (m/abs (m/- (aget idx i) (aget idx j)))))))))
      [d counts])))

(defn- sj-phi-fn
  "Returns the function `(fn [h order])`: the estimates of the functionals `phi4` (`order` 4) and `phi6` (`order` 6)
  of the gaussian kernel at the bandwidth `h` from the pair counts."
  [^doubles counts ^double n ^double d]
  (fn ^double [^double h ^long order]
    (let [len (alength counts)
          total (double
                 (loop [i 0
                        s 0.0]
                   (if (m/< i len)
                     (let [delta (m/sq (m// (m/* (double i) d) h))]
                       (if (m/>= delta SJ-DELMAX)
                         s
                         (recur (m/inc i)
                                (m/+ s (m/* (Array/aget counts i)
                                            (m/* (m/exp (m/* -0.5 delta))
                                                 (if (m/== order 4)
                                                   (m/+ (m/- (m/* delta delta) (m/* 6.0 delta)) 3.0)
                                                   (m/+ (m/- (m/* delta delta delta) (m/* 15.0 delta delta))
                                                        (m/- (m/* 45.0 delta) 15.0)))))))))
                     s)))
          sum (if (m/== order 4)
                (m/+ (m/* 2.0 total) (m/* n 3.0))
                (m/- (m/* 2.0 total) (m/* n 15.0)))]
      (m// sum (m/* n (m/- n 1.0) (m/fpow h (m/long-inc order)) m/SQRT2PI)))))

(defn- zeroin
  "Brent root finding `R_zeroin2` of R for `f` on `[ax, bx]` where `f` has values of opposite signs."
  ^double [f ^double ax ^double bx ^double tol]
  (let [fa0 (double (f ax))
        fb0 (double (f bx))]
  (cond
    (m/zero? fa0) ax
    (m/zero? fb0) bx
    :else
    (loop [a ax
           b bx
           c ax
           fa fa0
           fb fb0
           fc fa0
           iter 0]
      (let [prev-step (m/- b a)
            swap? (m/< (m/abs fc) (m/abs fb))
            a1 (if swap? b a)
            b1 (if swap? c b)
            c1 (if swap? b c)
            fa1 (if swap? fb fa)
            fb1 (if swap? fc fb)
            fc1 (if swap? fb fc)
            tol-act (m/+ (m/* 2.0 2.220446049250313e-16 (m/abs b1)) (m/* 0.5 tol))
            new-step (m/* 0.5 (m/- c1 b1))]
        (cond
          (or (m/<= (m/abs new-step) tol-act) (m/zero? fb1)) b1
          (m/> iter 1000) (invalid-input! :bandwidth :data "the root finding of the Sheather-Jones equation did not converge" iter
                                          "a converging root finding")
          :else
          (let [interpolate? (and (m/>= (m/abs prev-step) tol-act) (m/> (m/abs fa1) (m/abs fb1)))
                cb (m/- c1 b1)
                linear? (m/== a1 c1)
                t1 (if linear? (m// fb1 fa1) (m// fb1 fc1))
                t2 (m// fb1 fa1)
                qq (m// fa1 fc1)
                p0 (if linear?
                     (m/* cb t1)
                     (m/* t2 (m/- (m/* cb qq (m/- qq t1)) (m/* (m/- b1 a1) (m/- t1 1.0)))))
                q0 (if linear?
                     (m/- 1.0 t1)
                     (m/* (m/- qq 1.0) (m/- t1 1.0) (m/- t2 1.0)))
                p (if (m/> p0 0.0) p0 (m/- p0))
                q (if (m/> p0 0.0) (m/- q0) q0)
                step (if (and interpolate?
                              (m/< p (m/- (m/* 0.75 cb q) (m/* 0.5 (m/abs (m/* tol-act q)))))
                              (m/< p (m/abs (m/* 0.5 prev-step q))))
                       (m// p q)
                       new-step)
                step (if (m/< (m/abs step) tol-act)
                       (if (m/> step 0.0) tol-act (m/- tol-act))
                       step)
                b2 (m/+ b1 step)
                fb2 (double (f b2))
                same-sign? (or (and (m/> fb2 0.0) (m/> fc1 0.0)) (and (m/< fb2 0.0) (m/< fc1 0.0)))]
            (if same-sign?
              (recur b1 b2 b1 fb1 fb2 fb1 (m/inc iter))
              (recur b1 b2 c1 fb1 fb2 fc1 (m/inc iter))))))))))

(defn- sj
  "Sheather-Jones bandwidth for the gaussian kernel, `method` is `:ste` (solve-the-equation) or `:dpi`."
  ^double [{:keys [^doubles data ^long len]} method bins]
  (check-sample-size! :bandwidth len "the Sheather-Jones selector")
  (let [[d counts] (sj-pair-counts data bins)
        d (double d)
        ^doubles counts counts
        n (double len)
        [_ ^double sd] (mean-sd data)
        scale (m/min sd (m// (iqr-7 data) 1.349))
        a (m/* 1.24 scale (m/pow n (m// -1.0 7.0)))
        b (m/* 1.23 scale (m/pow n (m// -1.0 9.0)))
        c1 (m// 1.0 (m/* 2.0 m/SQRTPI n))
        phi (sj-phi-fn counts n d)
        phi4 (fn ^double [^double h] (double (phi h 4)))
        td (m/- (double (phi b 6)))
        sparse! (fn [what] (invalid-input! :bandwidth :data (str "sample is too sparse to find " what) len
                                          "a sample with a positive spread"))]
    (when-not (and (m/valid-double? td) (m/pos? td))
      (sparse! "TD"))
    (if (= method :dpi)
      (m/pow (m// c1 (double (phi4 (m/pow (m// 2.394 (m/* n td)) (m// 1.0 7.0))))) 0.2)
      (let [hmax (m/* 1.144 scale (m/pow n -0.2))
            alph2 (m/* 1.357 (m/pow (m// (double (phi4 a)) td) (m// 1.0 7.0)))
            _ (when-not (m/valid-double? alph2) (sparse! "alph2"))
            fsd (fn ^double [^double h]
                  (m/- (m/pow (m// c1 (double (phi4 (m/* alph2 (m/pow h (m// 5.0 7.0)))))) 0.2) h))
            [lower upper] (loop [lower (m/* 0.1 hmax)
                                 upper hmax
                                 itry 1]
                            (let [p (m/* (double (fsd lower)) (double (fsd upper)))]
                              (cond
                                (m/nan? p) (invalid-input! :bandwidth :data "the Sheather-Jones equation is not finite" len
                                                           "a sample with a positive spread")
                                (m/<= p 0.0) [lower upper]
                                (m/> itry 99) (invalid-input! :bandwidth :data "no solution in the specified range of bandwidths" len
                                                              "a sample with a positive spread")
                                (m/odd? itry) (recur lower (m/* upper 1.2) (m/inc itry))
                                :else (recur (m// lower 1.2) upper (m/inc itry)))))]
        (zeroin fsd (double lower) (double upper) (m/* 0.1 (double lower)))))))

;; cross validation

(def ^{:private true :const true :tag 'double} CV-GRID-DENSITY 40.0)
(def ^{:private true :const true :tag 'long} CV-MAX-INTEGRAL-STEPS 4194304)
(def ^{:private true :const true :tag 'double} CV-MAX-EXTENT 100.0)

(defn- loo-densities
  "Leave-one-out density at every data point: the point itself is left out and the sum is rescaled by `1/((n-1) h)`."
  ^doubles [{:keys [^doubles data ^long len ^IFn$DD kernel ^double radius]} ^double h]
  (let [span (m/* 1.01 h radius)
        factor (m// (m/* (m/dec len) h))
        out (double-array len)]
    (dotimes [i len]
      (Array/set out i (m/* factor (window-sum data kernel h span (Array/aget data i) i))))
    out))

(defn- integral-over-grid
  "Trapezoid integral of `g(f(t))` for the density `f` of `kdata` and bandwidth `h` over the range of the data extended
  by the extent of the kernel; the step is `h / 40` or larger when there would be too many steps."
  ^double [{:keys [^doubles data ^long len ^double extent] :as kdata} ^double h ^IFn$DD g]
  (let [^IFn$DD f (:kde (kde- (assoc kdata :binned? false) h))
        reach (m/* h (m/min extent CV-MAX-EXTENT))
        lo (m/- (Array/aget data 0) reach)
        hi (m/+ (Array/aget data (m/dec len)) reach)
        steps (m/long-min CV-MAX-INTEGRAL-STEPS (long (m/ceil (m/* (m// (m/- hi lo) h) CV-GRID-DENSITY))))
        dt (m// (m/- hi lo) steps)]
    (m/* dt (double
             (loop [i 1
                    s (m/* 0.5 (m/+ (.invokePrim g (.invokePrim f lo)) (.invokePrim g (.invokePrim f hi))))]
               (if (m/< i steps)
                 (recur (m/inc i) (m/+ s (.invokePrim g (.invokePrim f (m/+ lo (m/* (double i) dt))))))
                 s))))))

(defn- rlcv-a
  "The threshold `a_n` of the robust likelihood cross validation for the dimension `1`."
  ^double [^long len]
  (m// (m/* m/INV_SQRT_2 (m/sqrt (m/ln len))) (double len)))

(defn- mean-of
  ^double [^doubles v]
  (let [n (alength v)]
    (m// (double (loop [i 0 s 0.0] (if (m/< i n) (recur (m/inc i) (m/+ s (Array/aget v i))) s))) (double n))))

(defn- exact-sum
  "`kdata` with the sum over every data point for the predefined kernels of infinite support: the logarithm of a tiny
  density needs the tail of the kernel."
  [{:keys [kernel-name ^double radius] :as kdata}]
  (if (and (keyword? kernel-name) (m/> radius 1.0))
    (assoc kdata :radius ##Inf)
    kdata))

(defn- cv-criterion
  "Returns `[value sign]` of the cross validation criterion at `h`: the natural value of the criterion
  (to minimise for `:lscv`, to maximise for `:lcv` and `:rlcv`) in the units of the standardised data."
  [kdata criterion rlcv-form ^double h]
  (let [^doubles loo (loo-densities kdata h)]
    (case criterion
      :lscv (m/- (integral-over-grid kdata h (reify IFn$DD (invokePrim [_ v] (m/* v v))))
                 (m/* 2.0 (mean-of loo)))
      :lcv (let [^doubles loo (loo-densities (exact-sum kdata) h)
                 logs (double-array (alength loo))]
             (dotimes [i (alength loo)]
               (Array/set logs i (m/ln (Array/aget loo i))))
             (mean-of logs))
      :rlcv (let [a (rlcv-a (:len kdata))
                  la- (m/dec (m/ln a))
                  half-a (m/* 0.5 a)
                  continuous? (= rlcv-form :continuous)
                  gs (double-array (alength loo))
                  _ (dotimes [i (alength loo)]
                      (let [x (Array/aget loo i)]
                        (Array/set gs i (if (m/>= x a) (m/ln x) (m/+ la- (m// x a))))))
                  b (integral-over-grid kdata h
                                        (reify IFn$DD
                                          (invokePrim [_ v]
                                            (if (m/>= v a)
                                              (if continuous? (m/- v half-a) v)
                                              (m// (m/* v v) (m/* 2.0 a))))))]
              (m/- (mean-of gs) b)))))

(defn- check-bounds!
  [bounds]
  (if (and (sequential? bounds) (m/== (count bounds) 2)
           (every? finite-positive? bounds) (m/< (double (first bounds)) (double (second bounds))))
    [(double (first bounds)) (double (second bounds))]
    (invalid-input! :bandwidth :bounds "must be a pair of finite numbers 0 < lo < hi, factors of the rule of thumb" bounds
                    "[lo hi]")))

(defn- cv
  "Cross validation bandwidth.

  The data are standardised (`z = (x - mean) / sd`), the criterion is evaluated on a log-spaced grid of bandwidths
  over `bounds * s` (`s` is the rule of thumb adjusted to the kernel), the best grid point is refined by Brent's method
  in the neighbouring grid interval. Returns the info map with the bandwidth in the units of `data`."
  [kdata criterion {:keys [bounds grid rlcv-form]}]
  (check-sample-size! :bandwidth (:len kdata) "cross validation")
  (let [[lo-f hi-f] (check-bounds! (or bounds [0.01 5.0]))
        lo-f (double lo-f)
        hi-f (double hi-f)
        rlcv-form (or rlcv-form :printed)
        _ (when-not (#{:printed :continuous} rlcv-form)
            (invalid-input! :bandwidth :rlcv-form "must be :printed or :continuous" rlcv-form ":printed or :continuous"))
        ^doubles data (:data kdata)
        len (long (:len kdata))
        [^double mean ^double sd] (mean-sd data)
        _ (when-not (m/pos? sd)
            (invalid-input! :bandwidth :data "cross validation needs data which are not constant" mean "non-constant data"))
        z (let [z (double-array len)]
            (dotimes [i len]
              (Array/set z i (m// (m/- (Array/aget data i) mean) sd)))
            z)
        zk (assoc kdata :data z :binned? false)
        scale (double
               (let [raw (nrd zk 1.06)
                     nrd-scale (:nrd-scale (when (keyword? (:kernel-name kdata)) (get kde-data (:kernel-name kdata))))]
                 (if (and nrd-scale (finite-positive? nrd-scale)) (m/* raw (double nrd-scale)) raw)))
        g (long (if (nil? grid) (if (m/> len 5000) 30 60) grid))
        _ (when (m/< g 3)
            (invalid-input! :bandwidth :grid "must be an integer of at least 3" grid ">= 3"))
        h-lo (m/* scale lo-f)
        h-hi (m/* scale hi-f)
        log-ratio (m/ln (m// h-hi h-lo))
        hs (mapv (fn [^long k] (m/* h-lo (m/exp (m/* log-ratio (m// (double k) (double (m/long-dec g)))))))
                 (range g))
        sign (if (= criterion :lscv) -1.0 1.0)
        score (fn ^double [^double h]
                (let [v (double (cv-criterion zk criterion rlcv-form h))]
                  (if (m/valid-double? v) (m/* sign v) ##-Inf)))
        scores (mapv score hs)
        best (long (reduce (fn [^long bi ^long i]
                             (if (m/> (double (scores i)) (double (scores bi))) i bi))
                           0 (range 1 g)))
        best-score (double (scores best))]
    (if (m/== best-score ##-Inf)
      {:h (m/* sd h-lo) :value ##NaN :at-bound? true :status :degenerate}
      (let [[h-best score-best]
            (if (or (m/zero? best) (m/== best (m/dec g)))
              [(hs best) best-score]
              (let [[^double x ^double y] (acm/brent score {:bounds [(hs (m/dec best)) (hs (m/inc best))]
                                                            :initial (hs best)
                                                            :goal :maximize})]
                (if (m/> y best-score) [x y] [(hs best) best-score])))
            h-best (double h-best)]
        {:h (m/* sd h-best)
         :value (m/* sign (double score-best))
         :at-bound? (boolean (or (m/< (m/abs (m/dec (m// h-best h-lo))) 1.0e-6)
                                 (m/< (m/abs (m/dec (m// h-best h-hi))) 1.0e-6)))
         :status :ok}))))

(defn- infer-h
  "Returns the info map `{:h :value :at-bound? :status}` for the bandwidth specification `h`."
  [function kdata h opts]
  (let [direct (fn [^double v] {:h v :value ##NaN :at-bound? false :status :closed-form})]
    (cond
      (nil? h) (direct (nrd kdata 1.06))
      (number? h) (direct (check-bandwidth! function h))
      (keyword? h)
      (let [info (case h
                   :nrd (direct (nrd kdata 1.06))
                   :nrd0 (direct (nrd kdata 0.9))
                   :nrd-adjust (direct (m/* (nrd kdata 1.06) (kernel-scale kdata :nrd-adjust)))
                   :sj (direct (m/* (sj kdata :ste (check-bins! (:bins opts))) (kernel-scale kdata :sj)))
                   :sj-dpi (direct (m/* (sj kdata :dpi (check-bins! (:bins opts))) (kernel-scale kdata :sj-dpi)))
                   :rlcv (cv kdata :rlcv opts)
                   :lcv (cv kdata :lcv opts)
                   :lscv (cv kdata :lscv opts)
                   (invalid-input! function :bandwidth (str "unknown bandwidth selector, expected a number or one of " bandwidth-keywords)
                                   h "a number or a selector keyword"))]
        (check-bandwidth! function (:h info))
        (when (and (= :degenerate (:status info)) (not (:info? opts)))
          (invalid-input! function :bandwidth
                          "the criterion is not finite for any bandwidth of the search interval, use :info? to get the bound, or another selector or :bounds"
                          h "a selector with a finite criterion"))
        info)
      :else (invalid-input! function :bandwidth (str "must be a number or one of " bandwidth-keywords) h
                            "a number or a selector keyword"))))

(defn bandwidth
  "Estimates the bandwidth `h` of the kernel density estimate of `data`.

  The bandwidth controls how much each data point is spread out by the kernel: a smaller value follows the data more closely (a noisy estimate), a larger value gives a smoother and more biased estimate.

  Parameters:

  - `kernel` - a keyword naming a predefined kernel of [[kde-data]], or a kernel function (a function of a double returning a double; the selectors which need the canonical bandwidth of the kernel, `:nrd-adjust`, `:sj` and `:sj-dpi`, do not support it).
  - `data` - a sequence of numbers (finite, at least 2 for the selectors), not modified.
  - `h` - the number to use, or the keyword of a selector (`nil` means `:nrd`):
      - `:nrd` - the normal reference rule `1.06 * s * n^-0.2` with `s = min(sd, IQR / 1.34)`; the IQR uses the quantile of type 7 as R. When `s` is `0` the fallback is `sd`, then `abs(x_1)`, then `1.0` (as `bw.nrd0` in R).
      - `:nrd0` - Silverman's rule of thumb `0.9 * s * n^-0.2`, the same `s` (`bw.nrd0` in R).
      - `:nrd-adjust` - `:nrd` multiplied by the canonical bandwidth factor of the kernel (`:nrd-scale` in [[kde-data]]), so the kernels are comparable; throws for `:silverman`, `:cauchy` and custom kernels.
      - `:sj` - the Sheather-Jones solve-the-equation plug-in selector, the gaussian kernel result of `bw.SJ(method = 'ste')` in R (binned pair counts, the root found by `uniroot` with the same tolerance) multiplied by the canonical bandwidth factor of the kernel. It is the recommended selector for multimodal or non gaussian data. Throws when the sample is constant or has a zero interquartile range with too sparse data, as R.
      - `:sj-dpi` - the direct plug-in variant, `bw.SJ(method = 'dpi')`.
      - `:lscv` - least squares cross validation: minimises the integral of `f^2` minus twice the mean of the leave-one-out density.
      - `:lcv` - likelihood cross validation: maximises the mean log of the leave-one-out density (a data point with a non positive leave-one-out density makes the criterion `-Inf`; it is `-Inf` for every bandwidth of a compact kernel when a data point is isolated, then the selector throws `ex-info`, or with `:info?` returns the lower bound of the search interval with status `:degenerate`).
      - `:rlcv` - robust likelihood cross validation of Wu (*Journal of Business and Economic Statistics* 37(4), 2019): the likelihood criterion with the logarithm replaced by a linear function below `a_n = sqrt(ln n) / (sqrt(2) n)` and the bias term `b*(h)` the integral of `f` for `f >= a_n` and of `f^2 / (2 a_n)` below it.
  - `opts` (optional map):
      - `:info?` - when `true`, return the map described below instead of the number. Default: `false`.
      - `:bounds` - the cross validation search interval, a pair `[lo hi]` of factors of the rule of thumb `s`. Default: `[0.01 5.0]`.
      - `:grid` - the number of log spaced points of the first search stage, at least 3. Default: `60`, `30` for more than 5000 points.
      - `:bins` - the number of bins of `:sj` and `:sj-dpi`, an integer in `[2, 10000]`. Default: `1000`.
      - `:rlcv-form` - the bias term of `:rlcv`: `:printed` (the formula as printed in Wu 2019) or `:continuous`. The printed `b*` jumps by `a_n / 2` at `f = a_n`, the continuous form uses `f - a_n / 2` there, which is the one that follows from the definition of the paper. Default: `:printed`.
      - `:radius`, `:extent` - for a custom kernel function: the distance beyond which the kernel is ignored in the sum (default `##Inf`: every point) and beyond which the density is treated as zero (default `5.0`).

  The cross validation selectors standardise the data, so the result is equivariant to a shift and a scale of the data, search the criterion on a deterministic grid, refine the best grid point with Brent's method and cost `O(n^2)` evaluations of the kernel for each bandwidth tried (the `:binned?` option is ignored). Measured for the gaussian kernel and `n` points (default search, one core, including the JIT warm-up): `n = 1000` takes 0.5 s (`:lscv`), 0.6 s (`:lcv`, `:rlcv`); `n = 5000` takes 3.6 s, 8.3 s and 4.7 s; `n = 10000` takes 8 s, 25 s and 16 s. A kernel which sums every point (`:cauchy`, a custom kernel without `:radius`) is slower: `n = 1000` takes 10 s for `:lscv` against 0.14 s for `:epanechnikov`. The integral terms use a grid with a step of at most `h / 40` over the range of the data extended by `min(extent, 100) h`, so for `:cauchy` the tail beyond `100 h` is left out. `:sj` and `:sj-dpi` need `O(bins^2)` operations (below 0.1 s). When the best point is a bound of the search interval the bandwidth may be outside of it: see `:at-bound?`.

  Returns the bandwidth as a double. With `:info?` returns a map with `:h`, `:value` (the value of the criterion at `h` in standardised units, `##NaN` for the closed form selectors), `:at-bound?` (`true` when `h` is within `1.0e-6` relative of an end of the search interval) and `:status` (`:closed-form`, `:ok` or `:degenerate`).

  Throws `ex-info` (`:type :invalid-input`) for invalid data, kernel, selector or options.

  See also [[kernel-density]], [[kernel-density+]], [[kernel-density-ci]], [[kde-data]]."
  ([kernel data h] (bandwidth kernel data h nil))
  ([kernel data h {:keys [info?] :as opts}]
   (let [kdata (->kdata :bandwidth data kernel nil opts)
         info (infer-h :bandwidth kdata h opts)]
     (if info? info (double (:h info))))))

(defn kernel-density+
  "Builds a 1d kernel density estimator (KDE) for `data` and returns it together with the fitting details.

  A kernel density estimate is a smoothed histogram: it places a scaled copy of `kernel` at every data point and sums the contributions, weighted by a bandwidth `h` that controls how much each point spreads out.

  Parameters:

  - `kernel` - a keyword naming one of the predefined kernels of [[kde-data]], or a custom kernel function accepting and returning a double (same shape as [[gaussian]], [[epanechnikov]], etc.).
  - `data` - a sequence of numbers (finite, not empty), not modified.
  - `params` (optional map):
      - `:bandwidth` - the bandwidth `h`, a number greater than 0 or a keyword accepted by [[bandwidth]] (`:nrd`, `:nrd0`, `:nrd-adjust`, `:sj`, `:sj-dpi`, `:rlcv`, `:lcv`, `:lscv`). Default: `:nrd`.
      - `:bandwidth-options` - a map of the options of [[bandwidth]] (`:bounds`, `:grid`, `:bins`, `:rlcv-form`).
      - `:binned?` - `false` (default) for the exact evaluation. When `true` or a number, the data are replaced by a linear binning on a regular grid with the spacing `h / divisor` (the divisor is `5` for `true`, the number otherwise) and the density is a sum over the grid nodes; the error is of the second order in the spacing, and the cost does not depend on the number of data points. A grid of more than `2^24` nodes throws.
      - `:radius`, `:extent` - for a custom kernel function: the distance beyond which the kernel is ignored (default `##Inf`: the sum is over every data point) and beyond which the density is treated as zero (default `5.0`). Ignored for the predefined kernels (see [[kde-data]]).

  The evaluation is exact: the density at `x` sums the kernel over every data point within `1.01 * h * radius` of `x`, so it is not zero in a gap of the data wider than the kernel radius when the kernel is not compact.

  Returns a map:

  - `:kde` - the density function, a function of a double returning the estimate (`##NaN` at `##NaN`); it is negative at some points for `:silverman`.
  - `:factor` - the normalizing factor `1/(n*h)`, `n` is the number of data points (also with `:binned?`).
  - `:h` - the bandwidth used.
  - `:mn` and `:mx` - the support of the estimator: the data range extended by `1.01 * h * extent`; useful as integration bounds.

  Throws `ex-info` (`:type :invalid-input`) for empty data, a non finite value in the data, an unknown kernel or bandwidth selector, a bandwidth which is not finite and positive, or invalid `:binned?`.

  See also [[kernel-density]] (returns only the function), [[kernel-density-ci]] (adds confidence intervals), [[bandwidth]], [[kde-data]]."
  ([kernel data] (kernel-density+ kernel data nil))
  ([kernel data {:keys [bandwidth binned? bandwidth-options] :as params}]
   (let [kdata (->kdata :kernel-density data kernel binned? params)
         h (:h (infer-h :kernel-density kdata bandwidth (merge bandwidth-options (select-keys params [:radius :extent]))))]
     (kde- kdata h))))

(defn kernel-density
  "Returns a 1d kernel density estimation (KDE) function for `data`.

  A kernel density estimate is a smoothed histogram: it places a scaled copy of `kernel` at every data point and sums the contributions, weighted by a bandwidth that controls how much each point spreads out.

  Parameters:

  - `kernel` - a keyword naming one of the predefined kernels of [[kde-data]] (`:uniform`, `:triangular`, `:epanechnikov`, `:quartic`, `:triweight`, `:tricube`, `:gaussian`, `:cosine`, `:logistic`, `:sigmoid`, `:silverman`, `:wigner`, `:cauchy`, `:laplace`), or a custom kernel function accepting and returning a double.
  - `data` - a sequence of numbers (finite, not empty).
  - `params` (optional map): `:bandwidth`, `:bandwidth-options`, `:binned?`, `:radius`, `:extent`, with the meaning and the defaults of [[kernel-density+]].

  Returns a function of a single double returning the estimated density at that point as a double.

  Throws `ex-info` (`:type :invalid-input`) for invalid input, see [[kernel-density+]].

  See also [[kernel-density+]] (the same function plus the bandwidth, the factor and the support), [[kernel-density-ci]] (density with confidence intervals), [[bandwidth]], [[kde-data]]."
  ([kernel data] (kernel-density kernel data nil))
  ([kernel data params] (:kde (kernel-density+ kernel data params))))

(defn kernel-density-ci
  "Returns a 1d KDE function that also reports pointwise asymptotic confidence intervals.

  Uses the asymptotic normal approximation of the estimate (section 6.1.5 of http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html): the variance of the density estimate `f(x)` is approximated by `k2 * f(x) / (n*h)`, where `k2` is the roughness of the kernel (integral of the squared kernel), and the bounds are `f(x)` plus and minus the standard normal quantile of `1 - alpha/2` times the standard deviation. The interval is for the expected value of the estimate (the bias is ignored); the variance is computed from `max(f(x), 0)`.

  Parameters:

  - `kernel` - a keyword naming one of the predefined kernels of [[kde-data]]; a custom kernel function is not accepted, since its roughness `k2` is not known.
  - `data` - a sequence of numbers (finite, not empty).
  - `params` (optional map), forwarded to [[kernel-density+]], plus:
      - `:alpha` - the interval has the coverage `1 - alpha`, a number in `(0, 1)`. Default: `0.05`.

  Returns a function of a single double `x` returning a 3-element vector `[fx lower upper]`: the estimate and the confidence bounds. The lower bound is not clipped at zero.

  Throws `ex-info` (`:type :invalid-input`) for a kernel which is not in [[kde-data]], `alpha` outside of `(0, 1)` and the invalid input of [[kernel-density+]].

  See also [[kernel-density]], [[kernel-density+]], [[kde-data]]."
  ([kernel data] (kernel-density-ci kernel data nil))
  ([kernel data {:keys [alpha] :as params}]
   (when-not (and (keyword? kernel) (contains? kde-data kernel))
     (invalid-input! :kernel-density-ci :kernel (str "no roughness is known for the kernel, expected one of " kernel-keywords)
                     kernel "a keyword of kde-data"))
   (let [alpha (if (nil? alpha)
                 0.05
                 (if (and (number? alpha) (m/< 0.0 (double alpha) 1.0))
                   (double alpha)
                   (invalid-input! :kernel-density-ci :alpha "must be a number in (0, 1)" alpha "(0, 1)")))
         alpha (double alpha)
         kded (double (:k2 (kde-data kernel)))
         ^NormalDistribution local-normal (NormalDistribution.)
         za (.inverseCumulativeProbability local-normal (m/- 1.0 (m/* 0.5 alpha)))
         {:keys [kde ^double factor]} (kernel-density+ kernel data params)]
     (fn [^double x]
       (let [^double fx (kde x)
             band (m/* za (m/sqrt (m/* factor kded (m/max fx 0.0))))]
         [fx (m/- fx band) (m/+ fx band)])))))
