(ns fastmath.interpolation.shepard
  "Shepard (inverse distance weighting) interpolation of scattered data."
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- check-p!
  [p]
  (when-not (and (number? p) (m/valid-double? (double p)) (m/>= (double p) 0.0))
    (ic/invalid-input! :shepard :p
                       (str ":p has to be a finite number >= 0, got " (pr-str p))
                       {:option :p :value p :expected "a finite number >= 0"})))

(defn- scaled-distance
  "Euclidean distance of two arrays of an equal length, scaled by the largest coordinate difference: no overflow or underflow of the squares."
  ^double [^doubles a ^doubles b]
  (let [n (alength a)
        mx (double
            (loop [i (long 0)
                   mx 0.0]
              (if (m/< i n)
                (let [di (m/abs (m/- (Array/aget a i) (Array/aget b i)))]
                  (if (m/nan? di)
                    ##NaN
                    (recur (m/inc i) (m/max mx di))))
                mx)))]
    (if (or (m/nan? mx) (m/zero? mx) (m/invalid-double? mx))
      mx
      (m/* mx (m/sqrt (double
                       (loop [i (long 0)
                              sum 0.0]
                         (if (m/< i n)
                           (recur (m/inc i) (m/+ sum (m/sq (m// (m/abs (m/- (Array/aget a i) (Array/aget b i))) mx))))
                           sum))))))))

(defn- ->points
  "Converts the nodes to double arrays (a number is a one element point), checks their dimension."
  [xss]
  (let [pts (mapv m/seq->double-array xss)
        dim (alength ^doubles (first pts))]
    (doseq [^doubles pt pts]
      (when-not (m/== dim (alength pt))
        (ic/invalid-input! :shepard :dimension-mismatch
                           (str "all the nodes have to have the same dimension, got " dim " and " (alength pt))
                           {:expected dim :got (alength pt)})))
    pts))

(defn shepard
  "Creates the inverse distance weighting (Shepard) interpolator of scattered data.

  The value at a query point `x` is the weighted mean of the values `ys`, with the weights `1 / d(x, xs_i)^p`, where `d` is the distance to the `i`-th node. The interpolator passes exactly through the nodes. The result is always a weighted mean of `ys`, so it never leaves the range of the values and does not extrapolate trends: far from the data it tends to the mean of all the values. For `p > 1` it is flat (zero gradient) at the nodes.

  Parameters:

  - `xss` (sequence): the nodes, at least one; numbers for one dimensional data or sequences of numbers (vectors, `Vec2`, ...) of an equal length
  - `ys` (sequence of numbers): the values at the nodes, as many as nodes
  - `opts` (map, optional; a `nil` value is the same as a missing key):
      - `:p` (double, default: `2.0`): the power of the distance in the weights, finite and not negative; a larger `p` makes the result closer to the value of the nearest node, `p = 0` gives the mean of all the values (except at the nodes)
      - `:distance` (function, default: Euclidean distance, `fastmath.distance/euclidean-1d` for numbers, `fastmath.distance/euclidean` for sequences): a binary function `(fn [x1 x2])` returning the distance as a double

  Returns a function of a query point (a number for one dimensional data, also a one element sequence; a sequence of numbers otherwise) returning a double.

  A query point is a node when its distance to the node is exactly `0.0` (not when the keys are equal): a long `1` hits a node `1.0`, a vector `[1 0]` hits `[1.0 0.0]`, and a `Vec2` hits a vector. Then the value of the node is returned. When nodes are repeated, the value of the last repeated node wins. The weights are scaled by the smallest distance, so very small distances (`1e-154` and less) cannot overflow: the result is a weighted mean of the values for any `p`, never `Inf` or `NaN` (unless a value is `NaN`, the query point is `NaN`, or all the distances overflow to infinity). The default Euclidean distance is computed without overflow or underflow of the squares (coordinates differing by `1e-300` or `1e200` are fine); a custom `:distance` is called as is.

  Throws `ex-info` (with `:type :invalid-input`, `:method :shepard` and the key `:property`) when `xss` and `ys` have different lengths (`:length-mismatch`), there are no nodes (`:too-few-points`) `:p` is negative, infinite or `NaN` (`:p`), or, with the default distance, the nodes or the query point have different dimensions (`:dimension-mismatch`).

  See also [[fastmath.interpolation/shepard]], [[fastmath.interpolation.rbf/rbf]], [[fastmath.interpolation.kriging/kriging]]."
  ([xss ys] (shepard xss ys nil))
  ([xss ys opts]
   (let [{:keys [p distance] :or {p 2.0}} (into {} (remove (fn [[_ v]] (nil? v))) opts)
         _ (check-p! p)
         p (double p)
         xss (vec xss)
         ys (if (nil? ys) (double-array 0) (m/seq->double-array ys))
         n (count xss)]
     (when-not (m/== n (alength ys))
       (ic/invalid-input! :shepard :length-mismatch
                          (str "xss and ys have different lengths: " n " and " (alength ys))
                          {:xs-count n :ys-count (alength ys)}))
     (when (m/zero? n)
       (ic/invalid-input! :shepard :too-few-points "at least 1 point is required, got 0"
                          {:axis :xss :min-points 1 :count 0}))
     (let [d1? (number? (first xss))
           pts (when-not distance (->points xss))
           dim (when pts (alength ^doubles (first pts)))
           ;; fills the array of distances from the query to every node
           fill! (if distance
                   (fn [x ^doubles dists]
                     (let [x (if (and d1? (sequential? x) (m/== 1 (count x))) (first x) x)]
                       (dotimes [i n]
                         (Array/aset dists i (double (distance x (xss i)))))))
                   (fn [x ^doubles dists]
                     (let [q (m/seq->double-array x)]
                       (when-not (m/== (long dim) (alength q))
                         (ic/invalid-input! :shepard :dimension-mismatch
                                            (str "the query has the dimension " (alength q) ", the nodes have " dim)
                                            {:expected dim :got (alength q)}))
                       (dotimes [i n]
                         (Array/aset dists i (scaled-distance q ^doubles (pts i)))))))]
       (fn ^double [x]
         (let [dists (double-array n)
               _ (fill! x dists)
               ;; the smallest distance, NaN when any distance is NaN
               dmin (double
                     (loop [i (long 0)
                            dmin ##Inf]
                       (if (m/< i n)
                         (let [di (Array/aget dists i)]
                           (cond
                             (m/nan? di) ##NaN
                             (m/<= di dmin) (recur (m/inc i) di)
                             :else (recur (m/inc i) dmin)))
                         dmin)))]
           (cond
             (m/nan? dmin) ##NaN
             ;; a node hit, the last of the repeated nodes wins
             (m/zero? dmin) (loop [i (m/dec n)]
                              (if (m/zero? (Array/aget dists i))
                                (Array/aget ys i)
                                (recur (m/dec i))))
             :else (loop [i (long 0)
                          sw 0.0
                          swy 0.0]
                     (if (m/< i n)
                       (let [w (m/pow (m// dmin (Array/aget dists i)) p)]
                         (recur (m/inc i) (m/+ sw w) (m/+ swy (m/* w (Array/aget ys i)))))
                       (m// swy sw))))))))))
