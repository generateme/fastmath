(ns fastmath.interpolation
  "1d, 2d interpolation functions.

  See more:
  
  * [Apache Commons Math](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/analysis/interpolation/package-summary.html)
  * [Smile Interpolation](http://haifengl.github.io/smile/api/java/smile/interpolation/package-summary.html
  * [SSJ B-Spline](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1functionfit_1_1BSpline.html)

  Note: Smile interpolators also extrapolate values outside range.

  ### Input data

  You provide data as sequence or double array.
  
  #### 1d interpolation

  You provide two sequences:
  
  * `xs` - x axis coorditanes, strictly monotonic (increasing)
  * `ys` - function values

  See [[kriging-spline-interpolator]]
  
  #### 2d interpolation

  This is grid based interpolation.
  
  * `xs` - x axis coordinates, strictly monotonic (increasing)
  * `ys` - y axis coordinates, strictly monotonic (increasing)
  * `vs` - sequence of sequences of values (2d array) for all possible pairs. Array is column-wise: `[ [first column] [second column] ...]`.

  See [[cubic-2d-interpolator]]

  ### Examples

  Examples below use following functions:

  #### 1d
  ![1d](images/i/1d.png)

  #### 2d
  ![2d](images/i/2d.jpg)"
  {:metadoc/categories {:smile "Smile interpolators"
                        :comm "Apache Commons Math interpolators"
                        :ssj "SSJ interpolators"
                        :d1 "1d interpolation"
                        :d2 "2d interpolation (grid based)"}}
  (:require [fastmath.core :as m]
            [fastmath.kernel :as k])
  (:import [org.apache.commons.math3.analysis.interpolation AkimaSplineInterpolator DividedDifferenceInterpolator LinearInterpolator LoessInterpolator NevilleInterpolator SplineInterpolator MicrosphereProjectionInterpolator]
           [org.apache.commons.math3.analysis.interpolation BicubicInterpolator PiecewiseBicubicSplineInterpolator BivariateGridInterpolator]
           [org.apache.commons.math3.analysis.interpolation MultivariateInterpolator]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction BivariateFunction]
           [org.apache.commons.math3.analysis.function StepFunction]
           [smile.interpolation Interpolation CubicSplineInterpolation1D KrigingInterpolation1D LinearInterpolation RBFInterpolation1D ShepardInterpolation1D]
           [smile.interpolation Interpolation2D BicubicInterpolation BilinearInterpolation CubicSplineInterpolation2D]
           [umontreal.ssj.functionfit BSpline PolInterp]
           [umontreal.ssj.functions MathFunction]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)


;; 1d

(defmacro ^:private apache-commons-interpolator
  {:style/indent :defn}
  [n doc cat clazz] 
  (let [cl (with-meta (symbol "interp-obj") {:tag clazz})
        interp (with-meta (symbol "ufun") {:tag `UnivariateFunction})
        x (with-meta (gensym "x") {:tag double})
        n (vary-meta n assoc :metadoc/categories cat)
        xs (symbol "xs")
        ys (symbol "ys")]
    `(defn ~n ~doc
       [~xs ~ys]
       (let [~cl (new ~clazz)
             ~interp (.interpolate ~cl (m/seq->double-array ~xs) (m/seq->double-array ~ys))]
         (fn [~x] (.value ~interp ~x))))))

(apache-commons-interpolator akima-spline
  "Create cubic spline interpolator using Akima algorithm.
  Minimum number of points: 5

  xs[n] < xs[n+1] for all n.

Source: Apache Commons Math." #{:comm :d1}
  AkimaSplineInterpolator)

(apache-commons-interpolator divided-difference
  "Create Divided Difference Algorithm for interpolation.

Source: Apache Commons Math." #{:comm :d1}
  DividedDifferenceInterpolator)

(apache-commons-interpolator linear
  "Create Divided Difference Algorithm for inqterpolation.

Source: Apache Commons Math." #{:comm :d1}
  LinearInterpolator)

(defn- loess-interpolator-with-obj 
  "Create Loess function based on created object.

  Source: Apache Commons Math."
  [^LoessInterpolator obj xs ys]
  (let [^UnivariateFunction interp (.interpolate obj (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.value interp x))))

(defn loess
  "Local Regression Algorithm

  * bandwidth: 0.2-1.0 (optimal: 0.25-0.5, default: 0.4)
  * robustness-iters: 0-4 (optimal: 0, default: 2)
  * accuracy: double (default: 1e-12)

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d1}}
  ([xs ys] (loess-interpolator-with-obj (LoessInterpolator. 0.4 2) xs ys))
  ([bandwidth robustness-iters xs ys]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters) xs ys))
  ([bandwidth robustness-iters accuracy xs ys]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters accuracy) xs ys)))

(apache-commons-interpolator neville
                             "Neville algorithm

Source: Apache Commons Math." #{:comm :d1}
                             NevilleInterpolator)

(apache-commons-interpolator spline
                             "Cubic spline interpolation

Source: Apache Commons Math." #{:comm :d1}
                             SplineInterpolator)


(defn microsphere-projection
  "Microsphere projection interpolator - 1d version

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d1}}
  [elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance xs ys]
  (let [^MultivariateInterpolator interp (MicrosphereProjectionInterpolator. 1 elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance)
        xin (m/seq->double-double-array (map vector xs))
        ^MultivariateFunction f (.interpolate interp xin (m/seq->double-array ys))]
    (fn ^double [^double x] (.value f (double-array 1 x)))))

(defn cubic-spline
  "Cubic spline interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  [xs ys]
  (let [^Interpolation interp (CubicSplineInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn kriging-spline
  "Kriging interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  [xs ys]
  (let [^Interpolation interp (KrigingInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn linear-smile
  "Linear interpolation from Smile library.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  [xs ys]
  (let [^Interpolation interp (LinearInterpolation. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn rbf
  "RBF (Radial Basis Function) interpolation.

  Default kernel: `:gaussian`
  
  Source: Smile"
  {:metadoc/categories #{:smile :d1}}
  ([xs ys] (rbf (k/rbf :gaussian) xs ys))
  ([rbf-fn normalize? xs ys]
   (let [^Interpolation interp (RBFInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys) (k/smile-rbf rbf-fn) normalize?)]
     (fn ^double [^double x] (.interpolate interp x))))
  ([rbf-fn xs ys]
   (rbf rbf-fn false xs ys)))

(defn shepard
  "Shepard interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  ([xs ys]
   (let [^Interpolation interp (ShepardInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
     (fn ^double [^double x] (.interpolate interp x))))
  ([p xs ys]
   (let [^Interpolation interp (ShepardInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys) p)]
     (fn ^double [^double x] (.interpolate interp x)))))

(defn step-after
  "Step function."
  {:metadoc/categories #{:comm :d1}}
  [xs ys]
  (let [^StepFunction sf (StepFunction. (m/seq->double-array xs)
                                        (m/seq->double-array ys))]
    (fn ^double [^double x] (.value sf x))))

(defn step-before
  "Step function."
  {:metadoc/categories #{:comm :d1}}
  [xs ys]
  (let [x (m/seq->double-array xs)
        y (m/seq->double-array ys)
        l (dec (alength x))]
    (fn ^double [^double v]
      (let [b (java.util.Arrays/binarySearch ^doubles x v)
            i (if (neg? b) (m/constrain (dec (- b)) 0 l) b)]
        (aget ^doubles y i)))))

(defn step
  "Step function."
  {:metadoc/categories #{:comm :d1}}
  [xs ys]
  (let [x (m/seq->double-array xs)
        y (m/seq->double-array ys)
        l (dec (alength x))]
    (fn ^double [^double v]
      (let [b (java.util.Arrays/binarySearch ^doubles x v) 
            i (cond
                (== b -1) 0
                (< b (- l)) l
                (< b -1) (let [b-2 (- (- b) 2)
                               b-1 (dec (- b))
                               x1 (aget ^doubles x b-2)
                               x2 (aget ^doubles x b-1)
                               diff (+ x1 (* 0.5 (- x2 x1)))]
                           (if (< v diff) b-2 b-1))
                :else b)]
        (aget ^doubles y i)))))

;; monotonic
;; https://gist.github.com/lecho/7627739

(defn monotone
  "Monotone interpolation

  https://gist.github.com/lecho/7627739"
  {:metadoc/categories #{:d1}}
  [xs ys]
  (assert (and (seq xs) (seq ys)) "Sequences can't be empty.")
  (let [cntx (count xs)
        cnty (count ys)]
    (assert (and (> cntx 1) (> cnty 1)
                 (== cntx cnty)) "Sequnces have to be equal sizes and minimum 2 values each.")
    (assert (apply clojure.core/< xs) "x values have to be strictly monotonic")
    
    (let [^double fx (first xs)
          ^double lx (last xs)
          ^double fy (first ys)
          ^double ly (last ys)
          d (mapv (fn [[^double x1 ^double x2] [^double y1 ^double y2]]
                    (/ (- y2 y1) (- x2 x1))) (partition 2 1 xs) (partition 2 1 ys))
          m (vec (conj (map (fn [[^double d1 ^double d2]]
                              (if d2 (* 0.5 (+ d1 d2)) d1)) (partition-all 2 1 d)) (first d)))
          stop (- cntx 2)
          m (m/seq->double-array (loop [idx (int 0)
                                        ^double mi (m 0)
                                        ^double mi+ (m 1)
                                        ms []]
                                   (if (== idx stop) (conj ms mi mi+)
                                       (let [^double di (d idx)
                                             mi++ (m (+ idx 2))]
                                         (if (zero? di)
                                           (recur (inc idx) 0 mi++ (conj ms 0))
                                           (let [a (/ mi di)
                                                 b (/ mi+ di)
                                                 h (m/hypot-sqrt a b)]
                                             (if (> h 9.0)
                                               (let [t (/ 3.0 h)]
                                                 (recur (inc idx) (* t b di) mi++ (conj ms (* t a di))))
                                               (recur (inc idx) mi+ mi++ (conj ms mi)))))))))
          xs (m/seq->double-array xs)
          ys (m/seq->double-array ys)]
      (fn ^double [^double v]
        (cond
          (Double/isNaN v) v
          (<= v fx) fy
          (>= v lx) ly
          :else (let [b (java.util.Arrays/binarySearch ^doubles xs v)]
                  (if-not (neg? b)
                    (aget ys b)
                    (let [i (dec (- (inc b)))
                          mxi (aget xs i)
                          h (- (aget xs (inc i)) mxi)
                          t (/ (- v mxi) h)]
                      (+ (* (m/sq (- 1.0 t))
                            (+ (* h t (aget m i))
                               (* (inc (+ t t)) (aget ys i))))
                         (* t t (+ (* h (dec t) (aget m (inc i)))
                                   (* (- 3.0 t t) (aget ys (inc i))))))))))))))



;; b-spline

(defn- ssj-math-function
  [^MathFunction f] (fn ^double [^double x] (.evaluate f x)))

(defn b-spline-interp
  "B-spline interpolation.

  See:

  * 2 or 3 arity - exact interpolation using b-spline, default degree = 3 ([more](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1functionfit_1_1BSpline.html#a364fa9e72b7cdc0457140d79b2249530))
  * 4 arity - approximated b-spline interpolation with precision parameter `h` ([more](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1functionfit_1_1BSpline.html#a0892c41fc64e14a58e7f17208e05289a))"
  {:metadoc/categories #{:ssj :d1}}
  ([xs ys] (b-spline-interp 3 xs ys))
  ([^long degree xs ys]
   (ssj-math-function (BSpline/createInterpBSpline (m/seq->double-array xs) (m/seq->double-array ys) degree)))
  ([^long degree ^long h xs ys]
   (ssj-math-function (BSpline/createApproxBSpline (m/seq->double-array xs) (m/seq->double-array ys)
                                                   degree (m/constrain h degree (dec (count xs)))))))

(defn b-spline
  "B-spline for given points, default degree equals samples count - 1.

  [more](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1functionfit_1_1BSpline.html)"
  {:metadoc/categories #{:ssj :d1}}
  ([xs ys] (b-spline (dec (count xs)) xs ys))
  ([degree-or-knots xs ys]
   (let [b (if (integer? degree-or-knots)
             (BSpline. (m/seq->double-array xs) (m/seq->double-array ys) (int degree-or-knots))
             (BSpline. (m/seq->double-array xs) (m/seq->double-array ys) (m/seq->double-array degree-or-knots)))]
     (ssj-math-function b))))

(defn polynomial
  "Polynomial interpolation.

  [more](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1functionfit_1_1PolInterp.html)"
  {:metadoc/categories #{:ssj :d1}}
  [xs ys] (ssj-math-function (PolInterp. (m/seq->double-array xs) (m/seq->double-array ys))))

;;; 2d

(defn bicubic
  "Bicubic 2d.

  Grid based.

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [xs ys vs]
  (let [^BivariateGridInterpolator cl (BicubicInterpolator.)
        ^BivariateFunction interp (.interpolate cl (m/seq->double-array xs)
                                                (m/seq->double-array ys)
                                                (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.value interp x y))))

(defn piecewise-bicubic
  "Piecewise bicubic 2d.

  Grid based.

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [xs ys vs]
  (let [^BivariateGridInterpolator cl (PiecewiseBicubicSplineInterpolator.)
        ^BivariateFunction interp (.interpolate cl (m/seq->double-array xs)
                                                (m/seq->double-array ys)
                                                (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.value interp x y))))

(defn microsphere-2d-projection
  "Microsphere projection interpolator - 2d version

  Grid based.
  
  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance xs ys vs]
  (let [^MultivariateInterpolator interp (MicrosphereProjectionInterpolator. 2 elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance)
        xyin (m/seq->double-double-array (for [x xs
                                               y ys]
                                           [x y]))
        ^MultivariateFunction f (.interpolate interp xyin (m/seq->double-array (flatten vs)))]
    (fn ^double [^double x ^double y] (.value f (m/seq->double-array [x y])))))


(defn bilinear
  "Bilinear 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}}
  [xs ys vs]
  (let [^Interpolation2D interp (BilinearInterpolation. (m/seq->double-array xs)
                                                        (m/seq->double-array ys)
                                                        (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(defn bicubic-smile
  "Bicubic 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}}
  [xs ys vs]
  (let [^Interpolation2D interp (BicubicInterpolation. (m/seq->double-array xs)
                                                       (m/seq->double-array ys)
                                                       (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(defn cubic-2d
  "Cubic spline 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}}
  [xs ys vs]
  (let [^Interpolation2D interp (CubicSplineInterpolation2D. (m/seq->double-array xs)
                                                             (m/seq->double-array ys)
                                                             (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(def ^{:doc "Map of 1d interpolation functions"
       :metadoc/categories #{:smile :comm :d1}}
  interpolators-1d-list {:akima akima-spline
                         :divided-difference divided-difference
                         :linear linear
                         :loess loess
                         :neville neville
                         :spline spline
                         :cubic-spline cubic-spline
                         :kriging-spline kriging-spline
                         :linear-smile linear-smile
                         :rbf rbf
                         :shepard shepard
                         :microsphere microsphere-projection
                         :step-after step-after
                         :step-before step-before
                         :step step
                         :monotone monotone
                         :b-spline b-spline
                         :b-spline-interp b-spline-interp
                         :polynomial polynomial})

(def ^{:doc "Map of 2d interpolation functions"
       :metadoc/categories #{:smile :comm :d2}}
  interpolators-2d-list {:bicubic bicubic
                         :piecewise-bicubic piecewise-bicubic
                         :microsphere-2d microsphere-2d-projection
                         :bilinear bilinear
                         :bicubic-smile bicubic-smile
                         :cubic-2d cubic-2d})

;; TODO Multivariate

