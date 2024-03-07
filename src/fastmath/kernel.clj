(ns fastmath.kernel
  "Various kernel functions.

  * RBF (double -> double functions)
  * vector kernels (vector x vector -> double function; may be positive definite, conditional positive definite, positive semi-definite, mercer)
  * density estimation
  * some kernel operations"
  (:require [fastmath.core :as m]
            [fastmath.distance :as d]
            [fastmath.vector :as v]
            [fastmath.kernel.rbf :as rbf])
  (:import [smile.math.kernel MercerKernel]
           [smile.stat.distribution KernelDensity]
           [clojure.lang IFn]
           [org.apache.commons.math3.distribution NormalDistribution]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defmulti rbf
  "RBF kernel creator. RBF is double->double function.

  Parameters:

  All kernels accept `scale` parameter (as last parameter).

  Following kernels also accept `beta`: `:multiquadratic`, `:inverse-multiquadratic`, `:truncated-power`, `:radial-powers` and `:thin-plate`."
  (fn [k & _] k))

(defmacro ^:private emit-rbf
  ([nm f] `(emit-rbf ~nm ~f []))
  ([nm f argnames] `(emit-rbf ~nm ~f ~argnames nil))
  ([nm f argnames params]
   (let [argnames (conj argnames 'shape)
         orargnames '{shape 1.0}]
     `(defmethod rbf ~nm
        ([_#] (rbf ~nm nil))
        ([_# {:keys ~argnames :or ~orargnames}] (~f ~@params ~@argnames))))))

(emit-rbf :linear rbf/linear)
(emit-rbf :gaussian rbf/gaussian)

(emit-rbf :truncated-power rbf/truncated-power [k])
(emit-rbf :truncated-power-1 rbf/truncated-power [] [1.0])
(emit-rbf :truncated-power-2 rbf/truncated-power [] [2.0])
(emit-rbf :truncated-power-3 rbf/truncated-power [] [3.0])
(emit-rbf :truncated-power-half rbf/truncated-power [] [0.5])
(emit-rbf :truncated-power-third rbf/truncated-power [] [m/THIRD])

(emit-rbf :gaussians-laguerre rbf/gaussians-laguerre [dimension degree])
(emit-rbf :gaussians-laguerre-11 rbf/gaussians-laguerre [] [1.0 1.0])
(emit-rbf :gaussians-laguerre-12 rbf/gaussians-laguerre [] [1.0 2.0])
(emit-rbf :gaussians-laguerre-21 rbf/gaussians-laguerre [] [2.0 1.0])
(emit-rbf :gaussians-laguerre-22 rbf/gaussians-laguerre [] [2.0 1.0])

(emit-rbf :poisson rbf/poisson [d])
(emit-rbf :poisson-2 rbf/poisson [] [2.0])
(emit-rbf :poisson-3 rbf/poisson [] [3.0])
(emit-rbf :poisson-4 rbf/poisson [] [4.0])

(emit-rbf :matern rbf/matern [order])
(emit-rbf :matern-c0 rbf/matern [] [1.0])
(emit-rbf :matern-c2 rbf/matern [] [3.0])
(emit-rbf :matern-c4 rbf/matern [] [5.0])

(emit-rbf :generalized-multiquadratic rbf/generalized-multiquadratic [beta negate?])
(emit-rbf :multiquadratic rbf/generalized-multiquadratic [] [0.5 false])
(emit-rbf :inverse-multiquadratic rbf/generalized-multiquadratic [] [-0.5 false])

(emit-rbf :radial-powers rbf/radial-powers [beta negate?])
(emit-rbf :radial-powers-3 rbf/radial-powers [] [3.0 false])

(emit-rbf :thin-plate-splines rbf/thin-plate-splines [beta negate?])
(emit-rbf :thin-plate rbf/thin-plate-splines [] [1.0 false])

(emit-rbf :shifted-surface-splines rbf/shifted-surface-splines [s beta])

(emit-rbf :wendland rbf/wendland [s k])
(emit-rbf :gneiting rbf/gneiting [s l])

(emit-rbf :wu rbf/wu [l k])
(emit-rbf :wu-10 rbf/wu [] [1.0 0.0])
(emit-rbf :wu-11 rbf/wu [] [1.0 1.0])
(emit-rbf :wu-20 rbf/wu [] [2.0 0.0])
(emit-rbf :wu-21 rbf/wu [] [2.0 1.0])
(emit-rbf :wu-22 rbf/wu [] [2.0 2.0])
(emit-rbf :wu-30 rbf/wu [] [3.0 0.0])
(emit-rbf :wu-31 rbf/wu [] [3.0 1.0])
(emit-rbf :wu-32 rbf/wu [] [3.0 2.0])
(emit-rbf :wu-33 rbf/wu [] [3.0 3.0])

(emit-rbf :whittaker rbf/whittaker [alpha k beta])


;;;;;;;;;;;;;;;;;;;;;
;; Various kernels

(defn rbf->kernel
  "Treat RBF kernel as vector kernel using distance function (default [[euclidean]]."
  ([rbf-kernel] (rbf->kernel rbf-kernel d/euclidean))
  ([rbf-kernel distance]
   (fn [x y] (rbf-kernel (distance x y)))))

;; http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/
;; Marc G. Genton, Classes of Kernels for Machine Learning: A Statistics Perspective
;; http://www.jmlr.org/papers/volume2/genton01a/genton01a.pdf

(defmulti kernel
  "Crated vector kernel.

  Kernels can be Mercer, positive definite, conditional positive definite, positive semi-definite or other.
  
  Optional parameters:

  For `:gaussian`, `:exponential`, `:laplacian`, `:rational-quadratic`, `:multiquadratic`, `:inverse-multiquadratic`, `:circular`, `:spherical`, `:wave`, `:power`, `:log`, `:cauchy`, `:generalized-t-student`, `:hyperbolic-secant`, `:thin-plate`, `:mattern-12`, `:mattern-32`, `:mattern-52` and `::hyperbolic-secant` you can provide scaling parameter and `distance` (see [[fastmath.distance]], default is [[euclidean]]).

  Others:

  * `:linear` - `alpha`, scaling parameter
  * `:polynomial` - `alpha` (scaling), `c` (shift) and `d` (power)
  * `:anova` - `sigma` (scaling), `k` and `d` (power)
  * `:hyperbolic-tangent` - `alpha` (scaling), `c` (shift)
  * `:periodic` - `sigma` (scaling), `periodicity` and `distance`
  * `:bessel` - `sigma` (scaling), `n` and `v` (power factors) and `distance`
  * `:generalized-histogram` - `alpha` and `beta` (power factors)
  * `:dirichlet` - `N`
  * `:pearson` - `sigma` (scaling) and `omega` (power)

  Additionally there are two special kernels build from funcitons:

  * `:scalar-functions` - provide one or two double->double functions
  * `:variance-function` - provide any variance function (smooth, vector->double type)

  The rest of the kernels do not require parameters."
  (fn [k & _] k))

(defmethod kernel :linear
  ([_] (fn [x y] (v/dot x y)))
  ([_ ^double alpha] (fn [x y] (* alpha (v/dot x y)))))

(defmethod kernel :polynomial
  ([_] (fn [x y] (m/sq (v/dot x y))))
  ([_ ^double alpha ^double c ^double d] (fn [x y] (m/pow (+ c (* alpha (v/dot x y))) d))))

(defmethod kernel :gaussian
  ([_] (kernel :gaussian 1.0))
  ([_ sigma] (kernel :gaussian sigma d/euclidean))
  ([_ ^double sigma distance]
   (let [s2r  (/ (* sigma sigma))]
     (fn [x y] (m/exp (* -0.5 s2r (m/sq (distance x y))))))))

(defmethod kernel :exponential
  ([_] (kernel :exponential 1.0))
  ([_ sigma] (kernel :exponential sigma d/euclidean))
  ([_ ^double sigma distance]
   (let [s2r (/ (* sigma sigma))]
     (fn [x y] (m/exp (* -0.5 s2r ^double (distance x y)))))))

(defmethod kernel :laplacian
  ([_] (kernel :laplacian 1.0))
  ([_ sigma] (kernel :laplacian sigma d/euclidean))
  ([_ ^double sigma distance] (fn [x y] (m/exp (- (/ ^double (distance x y) sigma))))))

(defmethod kernel :anova
  ([_] (kernel :anova 1.0))
  ([_ ^double sigma] (kernel :anova sigma 1.0 1.0))
  ([_ ^double sigma ^double k ^double d]
   (let [powk #(m/pow ^double % k)
         powd #(m/pow ^double % d)]
     (fn [x y]
       (let [xk (v/fmap x powk)
             yk (v/fmap y powk)]
         (-> (v/sub xk yk)
             (v/sq)
             (v/mult (- sigma))
             (v/exp)
             (v/fmap powd)
             (v/sum)))))))

(defmethod kernel :hyperbolic-tangent
  ([_] (kernel :hyperbolic-tangent 1.0))
  ([_ ^double alpha] (kernel :hyperbolic-tangent alpha 0.0))
  ([_ ^double alpha ^double c]
   (fn [x y] (m/tanh (+ c (* alpha (v/dot x y)))))))

(defmethod kernel :rational-quadratic
  ([_] (kernel :rational-quadratic 1.0))
  ([_ ^double c] (kernel :rational-quadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (- 1.0 (/ d (+ d c)))))))

(defmethod kernel :multiquadratic
  ([_] (kernel :multiquadratic 1.0))
  ([_ ^double c] (kernel :multiquadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (m/sqrt (+ d (* c c)))))))

(defmethod kernel :inverse-multiquadratic
  ([_] (kernel :inverse-multiquadratic 1.0))
  ([_ ^double c] (kernel :inverse-multiquadratic c d/euclidean))
  ([_ ^double c distance]
   (fn [x y] (let [d (m/sq (distance x y))]
              (/ (m/sqrt (+ d (* c c))))))))

(defmethod kernel :circular
  ([_] (kernel :circular 1.0))
  ([_ ^double sigma] (kernel :circular sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (if (>= d sigma) 0.0
                  (let [ds (/ d sigma)]
                    (* 0.6366197723675814 (- (m/acos ds)
                                             (* ds (m/sqrt (- 1.0 (* ds ds))))))))))))

(defmethod kernel :spherical
  ([_] (kernel :spherical 1.0))
  ([_ ^double sigma] (kernel :spherical sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (if (>= d sigma) 0.0
                  (let [ds (/ d sigma)]
                    (+ (- 1.0 (* 1.5 ds)) (* 0.5 ds ds ds))))))))

(defmethod kernel :wave
  ([_] (kernel :wave 1.0))
  ([_ ^double sigma] (kernel :wave sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [^double d (distance x y)]
              (if (zero? d) 1.0
                  (* (/ sigma d) (m/sin (/ d sigma))))))))

(defmethod kernel :periodic
  ([_] (kernel :periodic 1.0))
  ([_ sigma] (kernel :periodic sigma 1.0))
  ([_ sigma periodicity] (kernel :periodic sigma periodicity d/euclidean))
  ([_ ^double sigma ^double periodicity distance]
   (let [p (/ m/PI periodicity)
         s2 (* sigma sigma)]
     (fn [x y] (let [^double d (distance x y)]
                (m/exp (/ (* -2.0 (m/sq (m/sin (* p d)))) s2)))))))

(defmethod kernel :power
  ([_] (kernel :power 2.0))
  ([_ d] (kernel :power d d/euclidean))
  ([_ ^double d distance]
   (fn [x y] (- (m/pow (distance x y) d)))))

(defmethod kernel :log
  ([_] (kernel :log 2.0))
  ([_ d] (kernel :log d d/euclidean))
  ([_ ^double d distance]
   (fn [x y] (- (m/log1p (m/pow (distance x y) d))))))

(defmethod kernel :spline
  [_] (fn [x y]
        (reduce #(* ^double %1 ^double %2) 1.0 (map (fn [^double xi ^double yi]
                                                      (let [xiyi (* xi yi)
                                                            m (min xi yi)
                                                            m2 (* m m)]
                                                        (inc (+ xiyi
                                                                (* xiyi m)
                                                                (* -0.5 (+ xi yi) m2)
                                                                (* m/THIRD m2 m))))) x y))))

(defmethod kernel :bessel
  ([_] (kernel :bessel 2.0))
  ([_ sigma] (kernel :bessel sigma 2.0))
  ([_ sigma n] (kernel :bessel sigma n -1.0))
  ([_ sigma n v] (kernel :bessel sigma n v d/euclidean))
  ([_ sigma n v distance]
   (fn [x y] (let [^double d (distance x y)
                  v+ (inc ^double v)]
              (/ (m/bessel-j v+ (* ^double sigma d))
                 (m/pow d (- (* ^double n v+))))))))

(defmethod kernel :cauchy
  ([_] (kernel :cauchy 1.0))
  ([_ ^double sigma] (kernel :cauchy sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (/ (inc (m/sq (/ ^double (distance x y) sigma)))))))

(defmethod kernel :chi-square-pd
  [_] (fn [x y] (reduce #(+ ^double %1 ^double %2) 0.0 (map (fn [^double xi ^double yi]
                                                             (/ (* 2.0 xi yi)
                                                                (+ xi yi))) x y))))

(defmethod kernel :chi-square-cpd
  [_] (fn [x y] (- 1.0 ^double (reduce #(+ ^double %1 ^double %2) 0.0 (map (fn [^double xi ^double yi]
                                                                            (/ (m/sq (- xi yi))
                                                                               (* 0.5 (+ xi yi)))) x y)))))

(defmethod kernel :histogram
  [_] (fn [x y] (reduce #(+ ^double %1 ^double %2) 0.0 (map (fn [^double xi ^double yi]
                                                             (min xi yi)) x y))))

(defmethod kernel :generalized-histogram
  ([_] (kernel :generalized-histogram 1.0 1.0))
  ([_ ^double alpha ^double beta]
   (fn [x y] (reduce #(+ ^double %1 ^double %2) 0.0 (map (fn [^double xi ^double yi]
                                                          (min (m/pow (m/abs xi) alpha)
                                                               (m/pow (m/abs yi) beta))) x y)))))

(defmethod kernel :generalized-t-student
  ([_] (kernel :generalized-t-student 1.0))
  ([_ ^double d] (kernel :generalized-t-student d d/euclidean))
  ([_ ^double d distance]
   (fn [x y] (/ (inc (m/pow (distance x y) d))))))

(defmethod kernel :dirichlet
  ([_] (kernel :dirichlet 1.0))
  ([_ ^double N]
   (let [N5 (+ 0.5 N)]
     (fn [x y] (reduce #(* ^double %1 ^double %2) 1.0 (map (fn [^double xi ^double yi]
                                                            (let [delta (- xi yi)
                                                                  num (m/sin (* N5 delta))
                                                                  den (* 2.0 (m/sin (* 0.5 delta)))]
                                                              (/ num den))) x y))))))

(defmethod kernel :hellinger
  ([_] (fn [x y] (v/dot (v/safe-sqrt x)
                       (v/safe-sqrt y)))))

(defmethod kernel :pearson
  ([_] (kernel :pearson 1.0 1.0))
  ([_ ^double sigma ^double omega]
   (let [c (/ (* 2.0 (m/sqrt (dec (m/pow 2.0 (/ omega))))) sigma)]
     (fn [x y] (let [xx (v/sum (v/sq x))
                    yy (v/sum (v/sq y))
                    xy (v/sum (v/emult x y))
                    m (* c (m/sqrt (+ (* -2.0 xy) xx yy)))]
                (/ (m/pow (inc (* m m)) omega)))))))

(defmethod kernel :thin-plate
  ([_] (kernel :thin-plate 1.0))
  ([_ sigma] (kernel :thin-plate sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y]
     (let [ds (/ ^double (distance x y) sigma)]
       (if-not (pos? ds) 0.0 (* (m/sq ds) (m/log ds)))))))

(defmethod kernel :mattern-12
  ([_] (kernel :mattern-12 1.0))
  ([_ ^double sigma] (kernel :mattern-12 sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (m/exp (- (/ ^double (distance x y) sigma))))))

(defmethod kernel :mattern-32
  ([_] (kernel :mattern-32 1.0))
  ([_ ^double sigma] (kernel :mattern-32 sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [d (/ (* m/SQRT3 ^double (distance x y)) sigma)]
              (* (inc d) (m/exp (- d)))))))

(defmethod kernel :mattern-52
  ([_] (kernel :mattern-52 1.0))
  ([_ ^double sigma] (kernel :mattern-52 sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [d (/ (* m/SQRT5 ^double (distance x y)) sigma)]
              (* (inc (+ d (* d d m/THIRD))) (m/exp (- d)))))))

(defmethod kernel :hyperbolic-secant
  ([_] (kernel :hyperbolic-secant 1.0))
  ([_ ^double sigma] (kernel :hyperbolic-secant sigma d/euclidean))
  ([_ ^double sigma distance]
   (fn [x y] (let [d (* sigma ^double (distance x y))]
              (+ (/ 2.0 (m/exp d)) (m/exp (- d)))))))

(defmethod kernel :scalar-functions
  ([_] (kernel :scalar-functions v/mag))
  ([_ f] (kernel :scalar-functions f f))
  ([_ f1 f2] (fn [x y] (* ^double (f1 x) ^double (f2 y)))))

(defmethod kernel :variance-function
  ([_] (kernel :variance-function v/mag))
  ([_ h] (fn [x y] (* 0.25 (- ^double (h (v/add x y)) ^double (h (v/sub x y)))))))

(def kernels-list ^{:doc "List of available vector kernels."} (sort (keys (methods kernel))))

(defn smile-mercer
  "Create Smile Mercer Kernel object

  Used to pass to Smile constructors/functions."
  [k]
  (reify
    MercerKernel (k [_ x y] (k x y))
    IFn (invoke ^double [_ x y] (k x y))))

;; kernel manipulation functions

(defn kernel->rbf
  "Convert vector kernel to RBF kernel. `center` is fixed `y` vector (default contains [[EPSILON]] values)."
  ([k] (kernel->rbf k m/EPSILON))
  ([k center]
   (let [c (vector center)]
     (fn [x] (k [x] c)))))

(defn exp
  "Kernel wraper. exp of kernel `k` with optional scaling value `t`."
  ([k] (exp k 1.0))
  ([k ^double t]
   (fn [x y] (m/exp (* t ^double (k x y))))))

(defn approx
  "Kernel wrapper. Round value returned by kernel using [[fastmath.core/approx]] function."
  ([k precision] (comp #(m/approx % precision) k))
  ([k] (comp m/approx k)))

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

(defn wadd
  "Kernel wrapper. Add kernels (weighted)."
  ([kernels] (wadd (repeat (count kernels) 1.0) kernels))
  ([weights kernels]
   (fn [x y] (reduce #(+ ^double %1 ^double %2)
                    (map (fn [^double w k] (* w ^double (k x y))) weights kernels)))))

(defn fields
  "Kernel wrapper. Apply vector field for each input before applying kernel function."
  ([k f] (fields k f f))
  ([k f1 f2]
   (fn [x y] (k (f1 x) (f2 y)))))

(defn- zero-vec [c] (vec (repeat c 0.0)))
(def ^:private zero-vec-m (memoize zero-vec))

;; doesn't work well
(defn cpd->pd
  "Convert conditionally positive definite kernel into positive definite.

  Formula is based on this [SO answer](https://stats.stackexchange.com/questions/149889/prove-that-a-kernel-is-conditionally-positive-definite). `x0` is equal `0`.
  
  Doesn't work well."
  [k]
  (fn [x y] (let [zero (zero-vec-m (count x))]
             (float (* 0.5 (+ ^double (k zero zero)
                              (- ^double (k x y)
                                 ^double (k x zero)
                                 ^double (k zero y))))))))



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

(defmethod kernel-density :smile
  ([_ vs h] (if h
              (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs) h)]
                (fn [x] (.p k x)))
              (kernel-density :smile vs)))
  ([_ vs] (let [^KernelDensity k (KernelDensity. (m/seq->double-array vs))]
            (fn [x] (.p k x)))))

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
