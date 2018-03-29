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
                        :rbf "Radial Basis Function"
                        :d1 "1d interpolation"
                        :d2 "2d interpolation (grid based)"}}
  (:require [fastmath.core :as m]
            [metadoc.examples :refer :all])
  (:import [org.apache.commons.math3.analysis.interpolation AkimaSplineInterpolator DividedDifferenceInterpolator LinearInterpolator LoessInterpolator NevilleInterpolator SplineInterpolator UnivariatePeriodicInterpolator MicrosphereProjectionInterpolator UnivariateInterpolator]
           [org.apache.commons.math3.analysis.interpolation BicubicInterpolator PiecewiseBicubicSplineInterpolator BivariateGridInterpolator]
           [org.apache.commons.math3.analysis.interpolation MultivariateInterpolator]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction BivariateFunction]
           [smile.interpolation Interpolation AbstractInterpolation CubicSplineInterpolation1D KrigingInterpolation1D LinearInterpolation RBFInterpolation1D ShepardInterpolation1D]
           [smile.interpolation Interpolation2D BicubicInterpolation BilinearInterpolation CubicSplineInterpolation2D]
           [smile.math.rbf RadialBasisFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; RBF
(defmulti rbf
  "Create Radial Basis Function

  Optional parameter: scaling factor."
  {:metadoc/categories #{:rbf}
   :metadoc/examples [(example "Usage" (let [rbf-fn (rbf :multiquadratic 3.0)]
                                         (rbf-fn 0.5)))]}
  (fn [n & _] n))

(defmethod rbf :gaussian
  ([_] (fn ^double [^double x] (m/exp (* 0.5 x x))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (m/exp (/ (* 0.5 x x) s2))))))

(defmethod rbf :inverse-multiquadratic
  ([_] (fn ^double [^double x] (/ (m/sqrt (inc (* x x))))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (/ (m/sqrt (+ (* x x) s2)))))))


(defmethod rbf :multiquadratic
  ([_] (fn ^double [^double x] (m/sqrt (inc (* x x)))))
  ([_ ^double scale] (let [s2 (* scale scale)]
                       (fn ^double [^double x] (m/sqrt (+ (* x x) s2))))))


(defmethod rbf :thinplate
  ([_] (fn ^double [^double x] (if (pos? x)
                                 (* x x (m/log x))
                                 0.0)))
  ([_ ^double scale] (fn ^double [^double x] (if (pos? x)
                                               (* x x (m/log (/ x scale)))
                                               0.0))))

(defn make-rbf-obj
  "Create RBF Smile object.

  Used to pass to Smile constructors/functions."
  {:metadoc/categories #{:rbf}
   :metadoc/examples [(example "Usage" (let [^RadialBasisFunction rbf-obj (make-rbf-obj (rbf :thinplate))]
                                         (.f rbf-obj 0.5)))]}
  [rbf-fn]
  (reify RadialBasisFunction
    (^double f [_ ^double x] (rbf-fn x))))


(def ^{:doc "Radial Basis function names"
       :metadoc/categories #{:rbf}
       :metadoc/examples [(example "List of names" rbfs-list)]}
  rbfs-list [:gaussian :inverse-multiquadratic :multiquadratic :thinplate])

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

(apache-commons-interpolator akima-spline-interpolator
                             "Create cubic spline interpolator using Akima algorithm.
  Minimum number of points: 5

  xs[n] < xs[n+1] for all n.

Source: Apache Commons Math." #{:comm :d1}
                             AkimaSplineInterpolator)

(apache-commons-interpolator divided-difference-interpolator
                             "Create Divided Difference Algorithm for interpolation.

Source: Apache Commons Math." #{:comm :d1}
                             DividedDifferenceInterpolator)

(apache-commons-interpolator linear-interpolator
                             "Create Divided Difference Algorithm for inqterpolation.

Source: Apache Commons Math." #{:comm :d1}
                             LinearInterpolator)

(defn- loess-interpolator-with-obj 
  "Create Loess function based on created object.

  Source: Apache Commons Math."
  [^LoessInterpolator obj xs ys]
  (let [^UnivariateFunction interp (.interpolate obj (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.value interp x))))

(defn loess-interpolator
  "Local Regression Algorithm

  * bandwidth: 0.2-1.0 (optimal: 0.25-0.5, default: 0.3)
  * robustness-iters: 0-4 (optimal: 0, default: 2)
  * accuracy: double (default: 1e-12)

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d1}}
  ([xs ys] (loess-interpolator-with-obj (LoessInterpolator.) xs ys))
  ([xs ys ^double bandwidth ^long robustness-iters]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters) xs ys))
  ([xs ys bandwidth robustness-iters accuracy]
   (loess-interpolator-with-obj (LoessInterpolator. bandwidth robustness-iters accuracy) xs ys)))

(apache-commons-interpolator neville-interpolator
                             "Neville algorithm

Source: Apache Commons Math." #{:comm :d1}
                             NevilleInterpolator)

(apache-commons-interpolator spline-interpolator
                             "Cubic spline interpolation

Source: Apache Commons Math." #{:comm :d1}
                             SplineInterpolator)


(defn microsphere-projection-interpolator
  "Microsphere projection interpolator - 1d version

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d1}}
  [xs ys elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance]
  (let [^MultivariateInterpolator interp (MicrosphereProjectionInterpolator. 1 elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance)
        xin (into-array (map m/seq->double-array (map vector xs)))
        ^MultivariateFunction f (.interpolate interp xin (m/seq->double-array ys))]
    (fn ^double [^double x] (.value f (m/seq->double-array 1 x)))))

(defn cubic-spline-interpolator
  "Cubic spline interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  [xs ys]
  (let [^Interpolation interp (CubicSplineInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn kriging-spline-interpolator
  "Kriging interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}
   :metadoc/examples [(example "Usage" {:test-value -0.07}
                        (let [interpolator (kriging-spline-interpolator [2 5 9 10 11] [0.4 1.0 -1.0 -0.5 0.0])]
                          (m/approx (interpolator 7.0))))]}
  [xs ys]
  (let [^Interpolation interp (KrigingInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn linear-smile-interpolator
  "Linear interpolation from Smile library.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  [xs ys]
  (let [^Interpolation interp (LinearInterpolation. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.interpolate interp x))))

(defn rbf-interpolator
  "RBF (Radial Basis Function) interpolation.

  Source: Smile"
  {:metadoc/categories #{:smile :d1}}
  ([xs ys rbf-fn normalize?]
   (let [rbf-obj (make-rbf-obj rbf-fn)
         ^Interpolation interp (RBFInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys) rbf-obj normalize?)]
     (fn ^double [^double x] (.interpolate interp x))))
  ([xs ys rbf-fn]
   (rbf-interpolator xs ys rbf-fn false)))

(defn shepard-interpolator
  "Shepard interpolation.

  Source: Smile."
  {:metadoc/categories #{:smile :d1}}
  ([xs ys]
   (let [^Interpolation interp (ShepardInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys))]
     (fn ^double [^double x] (.interpolate interp x))))
  ([xs ys p]
   (let [^Interpolation interp (ShepardInterpolation1D. (m/seq->double-array xs) (m/seq->double-array ys) p)]
     (fn ^double [^double x] (.interpolate interp x)))))

;;; 2d

(defn bicubic-interpolator
  "Bicubic 2d.

  Grid based.

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [xs ys vs]
  (let [^BivariateGridInterpolator cl (BicubicInterpolator.)
        ^BivariateFunction interp (.interpolate cl (m/seq->double-array xs)
                                                (m/seq->double-array ys)
                                                (into-array (map m/seq->double-array vs)))]
    (fn ^double [^double x ^double y] (.value interp x y))))

(defn piecewise-bicubic-interpolator
  "Piecewise bicubic 2d.

  Grid based.

  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [xs ys vs]
  (let [^BivariateGridInterpolator cl (PiecewiseBicubicSplineInterpolator.)
        ^BivariateFunction interp (.interpolate cl (m/seq->double-array xs)
                                                (m/seq->double-array ys)
                                                (into-array (map m/seq->double-array vs)))]
    (fn ^double [^double x ^double y] (.value interp x y))))

(defn microsphere-2d-projection-interpolator
  "Microsphere projection interpolator - 2d version

  Grid based.
  
  Source: Apache Commons Math."
  {:metadoc/categories #{:comm :d2}}
  [xs ys vs elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance]
  (let [^MultivariateInterpolator interp (MicrosphereProjectionInterpolator. 2 elements max-dark-friction dark-threshold background exponent shared-sphere? no-interpolation-tolerance)
        xyin (into-array (map m/seq->double-array (for [x xs
                                                        y ys]
                                                    [x y])))
        ^MultivariateFunction f (.interpolate interp xyin (m/seq->double-array (flatten vs)))]
    (fn ^double [^double x ^double y] (.value f (m/seq->double-array [x y])))))


(defn bilinear-interpolator
  "Bilinear 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}}
  [xs ys vs]
  (let [^Interpolation2D interp (BilinearInterpolation. (m/seq->double-array xs)
                                                        (m/seq->double-array ys)
                                                        (into-array (map m/seq->double-array vs)))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(defn bicubic-smile-interpolator
  "Bicubic 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}}
  [xs ys vs]
  (let [^Interpolation2D interp (BicubicInterpolation. (m/seq->double-array xs)
                                                       (m/seq->double-array ys)
                                                       (into-array (map m/seq->double-array vs)))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(defn cubic-2d-interpolator
  "Cubic spline 2d.

  Grid based.

  Source: Smile."
  {:metadoc/categories #{:smile :d2}
   :metadoc/examples [(example "Usage" {:test-value 4.68}
                        (let [interpolator (cubic-2d-interpolator [2 5 9] [2 3 10] [[4 0 2]
                                                                                    [-1 2 -2]
                                                                                    [-2 0 1]])]
                          (m/approx (interpolator 5.0 5.0))))
                      (example "Array layout"
                        (let [intrp (cubic-2d-interpolator [2 5] [1 6] [[-1 -2]
                                                                        [3 4]])]
                          [(intrp 2 1)
                           (intrp 2 6)
                           (intrp 5 1)
                           (intrp 5 6)]))]}
  [xs ys vs]
  (let [^Interpolation2D interp (CubicSplineInterpolation2D. (m/seq->double-array xs)
                                                             (m/seq->double-array ys)
                                                             (into-array (map m/seq->double-array vs)))]
    (fn ^double [^double x ^double y] (.interpolate interp x y))))

(def ^:private interpolators-list-symbol
  '{:akima akima-spline-interpolator
    :divided-diff divided-difference-interpolator
    :linear linear-interpolator
    :loess loess-interpolator
    :neville neville-interpolator
    :spline spline-interpolator
    :cubic-spline cubic-spline-interpolator
    :kriging kriging-spline-interpolator
    :linear-smile linear-smile-interpolator
    :rbf rbf-interpolator
    :shepard shepard-interpolator
    :microsphere microsphere-projection-interpolator
    :bicubic bicubic-interpolator
    :piecewise-bicubic piecewise-bicubic-interpolator
    :microsphere-2d microsphere-2d-projection-interpolator
    :bilinear bilinear-interpolator
    :bicubic-smile bicubic-smile-interpolator
    :cubic-2d cubic-2d-interpolator})

(def ^{:doc "Map of interpolation functions"
       :metadoc/categories #{:smile :comm}
       :metadoc/examples [(example "List of names" (keys interpolators-list))]}
  interpolators-list interpolators-list-symbol)

;; Examples

(def ^:const ^:private descr "Interpolation plot.")

(defmacro ^:private add-interpolator-examples
  []
  `(do
     ~@(for [[n f] interpolators-list-symbol
             :when (not (or
                         (= :rbf n)
                         (= :microsphere n)))]
         `(add-examples ~f
            (example-image descr ~(str "images/i/" (name n) ".jpg"))))))

(add-interpolator-examples)

(add-examples rbf-interpolator
  (example-image "With Gaussian RBF (scale=120)" "images/i/rbf-gaussian.jpg")
  (example-image "With Multiquadratic RBF (scale=120)" "images/i/rbf-multiquadratic.jpg")
  (example-image "With Inverse Multiquadratic RBF (scale=80)" "images/i/rbf-imultiquadratic.jpg")
  (example-image "With ThinPlate RBF (scale=80)" "images/i/rbf-thinplate.jpg"))

(add-examples loess-interpolator
  (example-image "Loess, bandwidth=0.2, robustness-iters=1" "images/i/loess-02-1.jpg")
  (example-image "Loess, bandwidth=0.7, robustness-iters=4" "images/i/loess-07-4.jpg"))

(add-examples shepard-interpolator
  (example-image "Shepard, p=5" "images/i/shepard-5.jpg")
  (example-image "Shepard, p=0.9" "images/i/shepard-09.jpg"))

(add-examples microsphere-projection-interpolator
  (example "Usage" {:evaluate? false} (microsphere-projection-interpolator xs ys 8 0.9 0.0000001 1 1.5 false 1))
  (example-image "Plot of above settings." "images/i/microsphere.jpg"))


;; TODO Multivariate
