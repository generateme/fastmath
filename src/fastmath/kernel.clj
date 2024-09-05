(ns fastmath.kernel
  "Various kernel functions.

  * RBF (double -> double functions)
  * vector kernels (vector x vector -> double function; may be positive definite, conditional positive definite, positive semi-definite, mercer)
  * density estimation
  * some kernel operations"
  (:require [fastmath.core :as m]
            [fastmath.kernel.rbf :as rbf]
            [fastmath.kernel.vector :as vk]
            [fastmath.kernel.density :as dens]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
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

;; kernel density estimation

(defn bandwidth
  "Returns infered bandwidth (h).

  h can be one of:

  * `:nrd` - rule-of-thumb (scale=1.06)
  * `:nrd0` - rule-of-thumb (scake=0.9)
  * `:nrd-adjust` - kernel specific adjustment of `:nrd`, doesn't work for `silverman` and `cauchy` 
  * `:rlcv` - robust likelihood cross-validation
  * `:lcv` - likelihood cross-validation
  * `:lscv` - least squares cross-validation"
  [kernel data h] (dens/bandwidth kernel data h))

(defn kernel-density
  "Returns kernel density estimation function, 1d
  
  Arguments:
  * `kernel` - kernel name or kernel function  
  * `data` - data
  * `params` - a map containing:
      * `:bandwidth` - bandwidth h
      * `:binned?` - if data should be binned, if `true` the width of the bin is `bandwidth` divided by 5, if is a number then it will be used as denominator. Default: `false`.

  `:bandwidth` can be a number or one of:

  * `:nrd` - rule-of-thumb (scale=1.06)
  * `:nrd0` - rule-of-thumb (scake=0.9)
  * `:nrd-adjust` - kernel specific adjustment of `:nrd`, doesn't work for `silverman` and `cauchy` 
  * `:rlcv` - robust likelihood cross-validation
  * `:lcv` - likelihood cross-validation
  * `:lscv` - least squares cross-validation"
  ([kernel data] (kernel-density kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth] (kernel-density kernel data bandwidth false))
  ([kernel data bandwidth info?]
   (let [h (if (number? bandwidth) {:bandwidth bandwidth} bandwidth)
         f (if info? dens/kernel-density+ dens/kernel-density)]
     (f kernel data h))))

(defn kernel-density-ci
  "Create function which returns confidence intervals for given kde method.

  Check 6.1.5 http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html

  Arguments:

  * `kernel` - kernel name  
  * `data` - sequence of data values
  * `params` - map with other parameters (including [[kernel-density]] parameters)
      * `alpha` - confidence level, default: 0.05

  Returns three values: density, lower confidence value and upper confidence value"
  ([kernel data] (kernel-density-ci kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth-or-params]
   (let [p (if (number? bandwidth-or-params) {:bandwidth bandwidth-or-params} bandwidth-or-params)]
     (dens/kernel-density-ci kernel data p))))

;;

(def kernel-list ^{:doc "List of available kernels: vector, rbf and kde"}
  {:vector (sort (keys (methods kernel)))
   :rbf (sort (keys (methods rbf)))
   :kde (sort (keys dens/kde-data))})

(m/unuse-primitive-operators)
