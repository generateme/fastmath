(ns fastmath.interpolation.acm
  (:require [fastmath.core :as m]
            [fastmath.interpolation.cubic :as cubic]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction BivariateFunction]
           [org.apache.commons.math3.analysis.interpolation AkimaSplineInterpolator DividedDifferenceInterpolator LoessInterpolator NevilleInterpolator MicrosphereProjectionInterpolator BicubicInterpolator]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn neville
  [xs ys]
  (let [^UnivariateFunction obj (.interpolate (NevilleInterpolator.)
                                              (m/seq->double-array xs)
                                              (m/seq->double-array ys))]
    (fn ^double [^double x] (.value obj x))))

(defn divided-difference
  [xs ys]
  (let [^UnivariateFunction obj (.interpolate (DividedDifferenceInterpolator.)
                                              (m/seq->double-array xs)
                                              (m/seq->double-array ys))]
    (fn ^double [^double x] (.value obj x))))


(defn loess
  ([xs ys] (loess xs ys nil))
  ([xs ys {:keys [^double bandwidth ^int iters ^double accuracy weights]
           :or {bandwidth LoessInterpolator/DEFAULT_BANDWIDTH
                iters LoessInterpolator/DEFAULT_ROBUSTNESS_ITERS
                accuracy LoessInterpolator/DEFAULT_ACCURACY}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         ^LoessInterpolator ls (LoessInterpolator. bandwidth iters accuracy)
         smooth-points (if weights
                         (.smooth ls xs ys (m/seq->double-array weights))
                         (.smooth ls xs ys))]
     (cubic/cubic xs smooth-points))))

(defn akima
  [xs ys]
  (let [xs (m/seq->double-array xs)

        ^"[Lorg.apache.commons.math3.analysis.polynomials.PolynomialFunction;"
        polys (-> (.interpolate (AkimaSplineInterpolator.) xs
                                (m/seq->double-array ys))
                  (.getPolynomials))]
    (fn ^double [^double x]
      (let [i (ic/binary-search xs x) ;; allow extrapolation
            ^UnivariateFunction poly (aget polys i)]
        (.value poly (m/- x (Array/aget xs i)))))))

(defn- infer-elements
  [x ^long cnt]
  (max (if (number? x) 5 (m/+ 4 (count x))) (m// cnt 5)))

(defn microsphere-projection
  ([xss ys] (microsphere-projection xss ys nil))
  ([xss ys {:keys [^int elements ^double exponent ^double max-dark-friction
                   ^double dark-threshold ^double background ^double no-interpolation-tolerance]
            :or {elements (infer-elements (first xss) (count xss))
                 exponent 1.0
                 max-dark-friction 0.9
                 dark-threshold 0.01
                 background 0.0
                 no-interpolation-tolerance 1.0e-6}}]
   (let [d1? (number? (first xss))
         xss (if d1? (map vector xss) xss)
         dims (count (first xss))
         ^MultivariateFunction obj (-> (MicrosphereProjectionInterpolator. dims elements
                                                                           max-dark-friction
                                                                           dark-threshold
                                                                           background
                                                                           exponent
                                                                           false
                                                                           no-interpolation-tolerance)
                                       (.interpolate (m/seq->double-double-array xss)
                                                     (m/seq->double-array ys)))]
     (condp = dims
       1 (fn ^double [^double x] (.value obj (double-array 1 x)))
       2 (fn mp2-interpolator (^double [xs] (.value obj (double-array xs)))
           (^double [^double x ^double y] (mp2-interpolator [x y])))
       (fn ^double [xs] (.value obj (m/seq->double-array xs)))))))

;; grid

(defn bicubic
  [xs ys vss]
  (let [^BivariateFunction obj (.interpolate (BicubicInterpolator.)
                                             (m/seq->double-array xs)
                                             (m/seq->double-array ys)
                                             (m/seq->double-double-array vss))]
    (fn bicubic-interpolator
      (^double [[^double x ^double y]] (bicubic-interpolator x y))
      (^double [^double x ^double y] (.value obj x y)))))
