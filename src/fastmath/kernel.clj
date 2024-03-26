(ns fastmath.kernel
  "Various kernel functions.

  * RBF (double -> double functions)
  * vector kernels (vector x vector -> double function; may be positive definite, conditional positive definite, positive semi-definite, mercer)
  * density estimation
  * some kernel operations"
  (:require [fastmath.core :as m]
            [fastmath.kernel.rbf :as rbf]
            [fastmath.kernel.vector :as vk])
  (:import [org.apache.commons.math3.distribution NormalDistribution]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defmacro ^:private emit
  ([kind nm f] `(emit ~kind ~nm ~f nil))
  ([kind nm f params]
   (if params
     `(defmethod ~kind ~nm
        ([_#] (~f ~params))
        ([_# params#] (~f (merge params# ~params))))
     `(defmethod ~kind ~nm
        ([_#] (~f))
        ([_# params#] (~f params#))))))


(defmulti rbf
  "RBF kernel creator. RBF is double->double function.

  Refer `fastmath.kernel.rbf` namespace for details."
  (fn [k & _] k))

(defmacro ^:private emit-rbf
  ([nm f] `(emit ~'rbf ~nm ~f))
  ([nm f params] `(emit ~'rbf ~nm ~f ~params)))

(emit-rbf :linear rbf/linear)
(emit-rbf :gaussian rbf/gaussian)

(emit-rbf :truncated-power rbf/truncated-power)
(emit-rbf :truncated-power-1 rbf/truncated-power {:k 1.0})
(emit-rbf :truncated-power-2 rbf/truncated-power {:k 2.0})
(emit-rbf :truncated-power-3 rbf/truncated-power {:k [3.0]})
(emit-rbf :truncated-power-half rbf/truncated-power {:k 0.5})
(emit-rbf :truncated-power-third rbf/truncated-power {:k m/THIRD})

(emit-rbf :gaussians-laguerre rbf/gaussians-laguerre)
(emit-rbf :gaussians-laguerre-11 rbf/gaussians-laguerre {:dimension 1.0 :degree 1.0})
(emit-rbf :gaussians-laguerre-12 rbf/gaussians-laguerre {:dimension 1.0 :degree 2.0})
(emit-rbf :gaussians-laguerre-21 rbf/gaussians-laguerre {:dimension 2.0 :degree 1.0})
(emit-rbf :gaussians-laguerre-22 rbf/gaussians-laguerre {:dimension 2.0 :degree 2.0})

(emit-rbf :poisson rbf/poisson)
(emit-rbf :poisson-2 rbf/poisson {:d 2.0})
(emit-rbf :poisson-3 rbf/poisson {:d 3.0})
(emit-rbf :poisson-4 rbf/poisson {:d 4.0})

(emit-rbf :matern rbf/matern)
(emit-rbf :matern-c0 rbf/matern {:order 1.0})
(emit-rbf :matern-c2 rbf/matern {:order 3.0})
(emit-rbf :matern-c4 rbf/matern {:order 5.0})

(emit-rbf :generalized-multiquadratic rbf/generalized-multiquadratic)
(emit-rbf :multiquadratic rbf/generalized-multiquadratic {:beta 0.5 :negate? false})
(emit-rbf :inverse-multiquadratic rbf/generalized-multiquadratic {:beta -0.5 :negate? false})

(emit-rbf :radial-powers rbf/radial-powers)
(emit-rbf :radial-powers-3 rbf/radial-powers {:beta 3.0 :negate? false})

(emit-rbf :thin-plate-splines rbf/thin-plate-splines)
(emit-rbf :thin-plate rbf/thin-plate-splines {:beta 1.0 :negate? false})

(emit-rbf :shifted-surface-splines rbf/shifted-surface-splines)

(emit-rbf :wendland rbf/wendland)
(emit-rbf :gneiting rbf/gneiting)

(emit-rbf :wu rbf/wu)
(emit-rbf :wu-10 rbf/wu {:l 1.0 :k 0.0})
(emit-rbf :wu-11 rbf/wu {:l 1.0 :k 1.0})
(emit-rbf :wu-20 rbf/wu {:l 2.0 :k 0.0})
(emit-rbf :wu-21 rbf/wu {:l 2.0 :k 1.0})
(emit-rbf :wu-22 rbf/wu {:l 2.0 :k 2.0})
(emit-rbf :wu-30 rbf/wu {:l 3.0 :k 0.0})
(emit-rbf :wu-31 rbf/wu {:l 3.0 :k 1.0})
(emit-rbf :wu-32 rbf/wu {:l 3.0 :k 2.0})
(emit-rbf :wu-33 rbf/wu {:l 3.0 :k 3.0})

(emit-rbf :whittaker rbf/whittaker)

;;

(defmulti kernel
  "Create vector kernel.

  Vector kernel returns a number for two vectors (or numbers)

  Kernels can be Mercer, positive definite, conditional positive definite, positive semi-definite or other."
  (fn [k & _] k))

(defmacro ^:private emit-kernel
  ([nm f] `(emit ~'kernel ~nm ~f))
  ([nm f params] `(emit ~'kernel ~nm ~f ~params)))

(emit-kernel :linear vk/linear)
(emit-kernel :polynomial vk/polynomial)
(emit-kernel :gaussian vk/gaussian)
(emit-kernel :exponential vk/exponential)
(emit-kernel :laplacian vk/laplacian)
(emit-kernel :anova vk/anova)

(emit-kernel :hyperbolic-tangent vk/hyperbolic-tangent)
(emit-kernel :hyperbolic-secant vk/hyperbolic-secant)

(emit-kernel :rational-quadratic vk/rational-quadratic)
(emit-kernel :multiquadratic vk/multiquadratic)
(emit-kernel :inverse-multiquadratic vk/inverse-multiquadratic)

(emit-kernel :triangular vk/geometric {:n 1})
(emit-kernel :circular vk/geometric {:n 2})
(emit-kernel :spherical vk/geometric {:n 3})
(emit-kernel :geometric vk/geometric)

(emit-kernel :wave vk/wave)
(emit-kernel :periodic vk/periodic)

(emit-kernel :power vk/power)
(emit-kernel :log vk/log)

(emit-kernel :spline vk/spline)
(emit-kernel :b-spline vk/b-spline)

(emit-kernel :bessel vk/bessel)
(emit-kernel :bessel2 vk/bessel2)

(emit-kernel :cauchy vk/cauchy)

(emit-kernel :chi-square vk/chi-square)
(emit-kernel :chi-square2 vk/chi-square2)

(emit-kernel :histogram vk/histogram)
(emit-kernel :generalized-histogram vk/generalized-histogram)
(emit-kernel :generalized-t-student vk/generalized-t-student)

(emit-kernel :dirichlet vk/dirichlet)
(emit-kernel :hellinger vk/hellinger)
(emit-kernel :pearson vk/pearson)

(emit-kernel :matern vk/matern)
(emit-kernel :matern-12 vk/matern {:order 1})
(emit-kernel :matern-32 vk/matern {:order 3})
(emit-kernel :matern-52 vk/matern {:order 5})

(def kernels-list ^{:doc "List of available vector kernels."} (sort (keys (methods kernel))))

(defn exp
  "Kernel wraper. exp of kernel `k` with optional scaling value `t`."
  ([k] (exp k 1.0))
  ([k ^double t]
   (fn [x y] (m/exp (* t ^double (k x y))))))

(defn scale
  "Kernel wrapper. Scale kernel result."
  [k ^double scale]
  (fn [x y] (* scale ^double (k x y))))

(defn mult
  "Kernel wrapper. Multiply two or more kernels."
  ([k1] k1)
  ([k1 k2] (fn [x y] (* ^double (k1 x y) ^double (k2 x y))))
  ([k1 k2 k3] (fn [x y] (* ^double (k1 x y) ^double (k2 x y) ^double (k3 x y))))
  ([k1 k2 k3 & r]
   (let [k (mult k1 k2 k3)]
     (if-not (seq r) k
             (apply mult k r)))))

(defn wmean
  "Kernel wrapper. (Weighted) mean of kernel results."
  ([kernels] (wmean kernels (repeat (count kernels) 1.0)))
  ([kernels weights]
   (fn [x y] (wmean (map (fn [k] (k x y)) kernels) weights))))


;;;;;;;;; density

(def ^{:const true :private true :tag 'double} gaussian-factor (/ (m/sqrt m/TWO_PI)))

(defn uniform-density-kernel ^double [^double x] (if (<= (m/abs x) 1.0) 0.5 0.0))
(defn gaussian-density-kernel ^double [^double x] (* gaussian-factor (m/exp (* -0.5 x x))))
(defn triangular-density-kernel ^double [^double x] (let [absx (m/abs x)]
                                                   (if (<= absx 1.0) (- 1.0 absx) 0.0)))
(defn epanechnikov-density-kernel ^double [^double x] (if (<= (m/abs x) 1.0) (* 0.75 (- 1.0 (* x x))) 0.0))
(defn quartic-density-kernel ^double [^double x] (if (<= (m/abs x) 1.0) (* 0.9375 (m/sq (- 1.0 (* x x)))) 0.0))
(defn triweight-density-kernel
  ^double [^double x]
  (if (<= (m/abs x) 1.0)
    (let [v (- 1.0 (* x x))]
      (* 1.09375 v v v)) 0.0))

(defn tricube-density-kernel
  ^double [^double x]
  (let [absx (m/abs x)]
    (if (<= absx 1.0)
      (let [v (- 1.0 (* absx absx absx))]
        (* 0.8875 v v v)) 0.0)))

(defn cosine-density-kernel ^double [^double x] (if (<= (m/abs x) 1.0) (* m/QUARTER_PI (m/cos (* m/HALF_PI x))) 0.0))
(defn logistic-density-kernel ^double [^double x] (/ (+ 2.0 (m/exp x) (m/exp (- x)))))
(defn sigmoid-density-kernel ^double [^double x] (/ (/ 2.0 (+ (m/exp x) (m/exp (- x)))) m/PI))
(defn silverman-density-kernel
  ^double [^double x]
  (let [xx (/ (m/abs x) m/SQRT2)]
    (* 0.5 (m/exp (- xx)) (m/sin (+ xx m/QUARTER_PI)))))

;; experimental
(defn laplace-density-kernel ^double [^double x] (* 0.5 (m/exp (- (m/abs x)))))
(defn wigner-density-kernel ^double [^double x] (if (<= (m/abs x) 1.0) (/ (* 2.0 (m/sqrt (- 1.0 (* x x)))) m/PI) 0.0))
(defn cauchy-density-kernel ^double [^double x] (/ (* m/HALF_PI (inc (m/sq (/ x 0.5))))))

;;

(defn- nrd
  ^double [data]
  (let [adata (m/seq->double-array data)
        sd (smile.math.MathEx/sd adata)
        iqr (- (smile.math.MathEx/q3 adata)
               (smile.math.MathEx/q1 adata))
        res (double (cond
                      (and (pos? sd) (pos? iqr)) (min sd (/ iqr 1.34))
                      (pos? sd) sd
                      (pos? iqr) (/ iqr 1.34)
                      :else 1.0))]
    (* 1.06 res (m/pow (alength ^doubles adata) -0.2))))

(defn- kde
  "Return kernel density estimation function"
  ([data k] (kde data k nil))
  ([data k h]
   (let [data (let [a (m/seq->double-array data)] (java.util.Arrays/sort a) a)
         last-idx (dec (alength data))
         h (double (or h (nrd data)))
         hrev (/ h)
         span (* 6.0 h)
         factor (/ (* (alength data) h))
         mn (aget data 0)
         mx (aget data last-idx)]
     [(fn [^double x]
        (let [start (java.util.Arrays/binarySearch data (- x span))
              start (long (if (neg? start) (dec (- start)) start))
              end (java.util.Arrays/binarySearch data (+ x span))
              end (min last-idx (long (if (neg? end) (dec (- end)) end)))]
          (loop [i start
                 sum 0.0]
            (if (<= i end)
              (recur (inc i) (+ sum ^double (k (* hrev (- x (aget data i))))))
              (* factor sum)))))
      factor h (- mn span) (+ mx span)])))

(defonce ^:private kde-integral
  {:uniform 0.5
   :triangular m/TWO_THIRD
   :epanechnikov 0.6
   :quartic (/ 5.0 7.0)
   :triweight (/ 350.0 429.0)
   :tricube (/ 175.0 247.0)
   :gaussian (* 0.5 (/ m/SQRTPI))
   :cosine (* 0.0625 m/PI m/PI)
   :logistic m/SIXTH
   :sigmoid (/ 2.0 (* m/PI m/PI))
   :silverman (* 0.0625 3.0 m/SQRT2)
   :wigner (/ 16.0 (* 3 m/PI m/PI))
   :cauchy m/M_1_PI
   :laplace 0.25})

(defmulti kernel-density
  "Create kernel density estimator.

  Parameters:

  * kernel name, see [[kernel-density-list]].
  * sequence of data values
  * optional: bandwidth (by default, bandwidth is estimated using nrd method)"
  (fn [k & _] k))

(defmacro ^:private make-kernel-density-fns
  [lst]
  `(do ~@(for [v (eval lst)
               :let [n (symbol (str (name v) "-density-kernel"))]]
           `(defmethod kernel-density ~v
              ([k# vs#] (first (kde vs# ~n)))
              ([k# vs# h#] (first (kde vs# ~n h#)))
              ([k# vs# h# all?#] (let [kded# (kde vs# ~n h#)]
                                   (if all?# kded# (first kded#))))))))

(make-kernel-density-fns (keys kde-integral))

(defmethod kernel-density :default [_ & r] (apply kernel-density :gaussian r))

(def kernel-density-list ^{:doc "List of available density kernels."} (sort (keys (methods kernel-density))))

(defn kernel-density-ci
  "Create function which returns confidence intervals for given kde method.

  Check 6.1.5 http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html

  Parameters:

  * `method` - kernel name
  * `data` - sequence of data values
  * `bandwidth`
  * `alpha` - confidence level parameter

  Returns three values: density, lower confidence, upper confidence"
  ([method data] (kernel-density-ci method data nil))
  ([method data bandwidth] (kernel-density-ci method data bandwidth 0.05))
  ([method data bandwidth ^double alpha]
   (if (contains? kde-integral method)
     (let [^NormalDistribution local-normal (NormalDistribution.)
           za (.inverseCumulativeProbability local-normal (- 1.0 (* 0.5 (or alpha 0.05))))
           [kde-f ^double factor] (kernel-density method data bandwidth true)]
       (fn [^double x]
         (let [^double fx (kde-f x)
               band (* za (m/sqrt (* factor ^double (kde-integral method) fx)))]
           [fx (- fx band) (+ fx band)])))
     (let [kde-f (kernel-density method data bandwidth)]
       (fn [x] (let [fx (kde-f x)]
                [fx fx fx]))))))

(m/unuse-primitive-operators)
