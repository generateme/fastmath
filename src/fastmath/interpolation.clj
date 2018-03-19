(ns fastmath.interpolation
  "1d, 2d interpolation functions.

  See more:
  * [Apache Commons Math](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/analysis/interpolation/package-summary.html)
  * [Smile Interpolation](http://haifengl.github.io/smile/api/java/smile/interpolation/package-summary.html)"
  (:import [org.apache.commons.math3.analysis.interpolation AkimaSplineInterpolator DividedDifferenceInterpolator LinearInterpolator LoessInterpolator NevilleInterpolator SplineInterpolator UnivariatePeriodicInterpolator UnivariateInterpolator] 
           [org.apache.commons.math3.analysis UnivariateFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

;; 1d

(defmacro ^:private apache-commons-interpolator
  {:style/indent :defn}
  [n doc clazz] 
  (let [cl (with-meta (symbol "interp-obj") {:tag clazz})
        interp (with-meta (symbol "ufun") {:tag `UnivariateFunction})
        x (with-meta (gensym "x") {:tag double})]
    `(defn ~n ~doc
       [xs# ys#]
       (let [~cl (new ~clazz)
             ~interp (.interpolate ~cl (double-array xs#) (double-array ys#))]
         (fn [~x] (.value ~interp ~x))))))

(apache-commons-interpolator akima-spline-interpolator
  "Create cubic spline interpolator using Akima algorithm.
  Minimum number of points: 5

  xs[n] < xs[n+1] for all n."
  AkimaSplineInterpolator)

(apache-commons-interpolator divided-difference-interpolator
  "Create Divided Difference Algorithm for interpolation."
  DividedDifferenceInterpolator)

(apache-commons-interpolator linear-interpolator
  "Create Divided Difference Algorithm for interpolation."
  LinearInterpolator)

(defn- loess-interpolator-with-obj 
  "Create Loess function based on created object."
  [^LoessInterpolator obj xs ys]
  (let [^UnivariateFunction interp (.interpolate obj (double-array xs) (double-array ys))]
    (fn ^double [^double x] (.value interp x))))

(defn loess-interpolator
  "Local Regression Algorithm

  * bandwidth: 0.0-1.0 (optimal: 0.25-0.5, default: 0.3)
  * robustness-iters: 0-4 (optimal: 0, default: 2)
  * accuracy: double (default: 1e-12)"
  ([xs ys] (loess-interpolator-with-obj (LoessInterpolator.) xs ys))
  ([xs ys ^double bandwidth ^long robustness-iters]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters) xs ys))
  ([xs ys bandwidth robustness-iters accuracy]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters accuracy) xs ys)))

(apache-commons-interpolator neville-interpolator
  "Neville algorithm"
  NevilleInterpolator)

(apache-commons-interpolator spline-interpolator
  "Cubic spline interpolation"
  SplineInterpolator)


;; TODO Periodic versions


