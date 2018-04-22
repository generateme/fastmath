(ns fastmath.interpolation
  "1d, 2d interpolation functions.

  See more:
  
  * [Apache Commons Math](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/analysis/interpolation/package-summary.html)
  * [Smile Interpolation](http://haifengl.github.io/smile/api/java/smile/interpolation/package-summary.html)

  Note: Smile interpolators doesn't check ranges.

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
  "
  {:metadoc/categories {:smile "Smile interpolators"
                        :comm "Apache Commons Math interpolators"
                        :d1 "1d interpolation"
                        :d2 "2d interpolation (grid based)"}}
  (:require [fastmath.core :as m]
            [metadoc.examples :refer :all]
            [fastmath.rbf :as rbf])
  (:import [org.apache.commons.math3.analysis.interpolation AkimaSplineInterpolator DividedDifferenceInterpolator LinearInterpolator LoessInterpolator NevilleInterpolator SplineInterpolator UnivariatePeriodicInterpolator MicrosphereProjectionInterpolator UnivariateInterpolator]
           [org.apache.commons.math3.analysis.interpolation BicubicInterpolator PiecewiseBicubicSplineInterpolator BivariateGridInterpolator]
           [org.apache.commons.math3.analysis.interpolation MultivariateInterpolator]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction BivariateFunction]
           [org.apache.commons.math3.analysis.function StepFunction]
           [smile.interpolation Interpolation AbstractInterpolation CubicSplineInterpolation1D KrigingInterpolation1D LinearInterpolation RBFInterpolation1D ShepardInterpolation1D]
           [smile.interpolation Interpolation2D BicubicInterpolation BilinearInterpolation CubicSplineInterpolation2D]           ))

(set! *warn-on-reflection* true)
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

  * bandwidth: 0.2-1.0 (optimal: 0.25-0.5, default: 0.3)
  * robustness-iters: 0-4 (optimal: 0, default: 2)
  * accuracy: double (default: 1e-12)

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d1}}
  ([xs ys] (loess-interpolator-with-obj (LoessInterpolator.) xs ys))
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
  {:metadoc/categories #{:smile :d1}
   :metadoc/examples [(example "Usage" {:test-value -0.07}
                        (let [interpolator (kriging-spline [2 5 9 10 11] [0.4 1.0 -1.0 -0.5 0.0])]
                          (m/approx (interpolator 7.0))))]}
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

  Source: Smile"
  {:metadoc/categories #{:smile :d1}}
  ([rbf-fn normalize? xs ys]
   (let [rbf-obj (rbf/rbf-obj rbf-fn)
         ^Interpolation interp (RBFInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys) rbf-obj normalize?)]
     (fn ^double [^double x] (.interpolate interp x))))
  ([rbf-fn xs ys]
   (rbf xs ys rbf-fn false)))

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
  {:metadoc/categories #{:smile :d2}
   :metadoc/examples [(example "Usage" {:test-value 4.68}
                               (let [interpolator (cubic-2d [2 5 9] [2 3 10] [[4 0 2]
                                                                              [-1 2 -2]
                                                                              [-2 0 1]])]
                                 (m/approx (interpolator 5.0 5.0))))
                      (example "Array layout"
                               (let [intrp (cubic-2d [2 5] [1 6] [[-1 -2]
                                                                  [3 4]])]
                                 [(intrp 2 1)
                                  (intrp 2 6)
                                  (intrp 5 1)
                                  (intrp 5 6)]))]}
  [xs ys vs]
  (let [^Interpolation2D interp (CubicSplineInterpolation2D. (m/seq->double-array xs)
                                                             (m/seq->double-array ys)
                                                             (m/seq->double-double-array vs))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(def ^{:doc "Map of 1d interpolation functions"
       :metadoc/categories #{:smile :comm :d1}
       :metadoc/examples [(example "List of names" (keys interpolators-1d-list))]}
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
                         :step-before step-before})

(def ^{:doc "Map of 2d interpolation functions"
       :metadoc/categories #{:smile :comm :d2}
       :metadoc/examples [(example "List of names" (keys interpolators-1d-list))]}
  interpolators-2d-list {:bicubic bicubic
                         :piecewise-bicubic piecewise-bicubic
                         :microsphere-2d microsphere-2d-projection
                         :bilinear bilinear
                         :bicubic-smile bicubic-smile
                         :cubic-2d cubic-2d})

;; TODO Multivariate

