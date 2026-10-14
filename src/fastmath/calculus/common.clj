(ns fastmath.calculus.common
  "Shared parts of the integration engines: the max-error queue, the orientation of bounds and the
  substitutions which turn improper integrals into proper ones."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [java.util PriorityQueue Comparator]
           [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; options

(defn drop-nil-options
  "Returns the `options` map without the keys having a `nil` value, so a `nil` option is a missing option and takes its default. `nil` options give `{}`."
  [options]
  (into {} (remove (comp nil? val)) options))

;; adaptive engines: queue of regions ordered by the error estimate

(definterface Estimate
  (^double value [])
  (^double error []))

(def ^:private ^Comparator largest-error-first
  (reify Comparator
    (compare [_ a b]
      ;; Double/compare is a total order: NaN is the largest value, so a region with a NaN error is refined first
      (Double/compare (.error ^Estimate b) (.error ^Estimate a)))))

(defn max-error-queue
  "Creates an empty `java.util.PriorityQueue` of [[Estimate]] objects whose head (`poll`) is the region with the largest error.

  Regions with equal errors are all kept and a region with a NaN error comes first."
  ^PriorityQueue []
  (PriorityQueue. 64 largest-error-first))

(defn estimates-sum
  "Returns the sum of values and the sum of errors of all [[Estimate]] objects in a collection as a `Vec2` (it destructures like `[sum-of-values sum-of-errors]`)."
  [regions]
  (reduce (fn [^Vec2 ve ^Estimate r]
            (Vec2. (m/+ (.x ve) (.value r)) (m/+ (.y ve) (.error r)))) (v/vec2 0.0 0.0) regions))

;; orientation of bounds

(defn- bound-error
  [lower upper]
  (ex-info "Bounds must be numbers and must not both be the same infinity" {:lower lower :upper upper}))

(defn normalize-bounds
  "Orients one-dimensional bounds. Returns `[sign lower upper]` with `lower <= upper` and `sign` equal to `-1.0` when the bounds were reversed, `1.0` otherwise.

  Throws `ex-info` for a NaN bound and for two equal infinite bounds."
  [^double lower ^double upper]
  (when (or (m/nan? lower) (m/nan? upper)
            (and (m/inf? lower) (m/== lower upper)))
    (throw (bound-error lower upper)))
  (if (m/> lower upper)
    [-1.0 upper lower]
    [1.0 lower upper]))

(defn normalize-bounds-multi
  "Orients bounds of every dimension. Returns `[sign lower upper]`, where `lower` and `upper` are vectors with `lower[i] <= upper[i]` and `sign` is `-1.0` raised to the number of reversed dimensions.

  Throws `ex-info` when the sequences are empty or differ in length, and as [[normalize-bounds]] does for a NaN bound or two equal infinite bounds in any dimension."
  [lower upper]
  (let [lower (vec lower)
        upper (vec upper)]
    (when (or (m/zero? (count lower)) (not= (count lower) (count upper)))
      (throw (ex-info "Lower and upper bounds must be non-empty and of the same length"
                      {:lower lower :upper upper})))
    (let [[sign lower upper] (reduce (fn [[^double sign lo up] [^double a ^double b]]
                                       (let [[^double s a b] (normalize-bounds a b)]
                                         [(m/* sign s) (conj lo a) (conj up b)]))
                                     [1.0 [] []] (map v/vec2 lower upper))]
      [sign lower upper])))

(defn orient
  "Applies the `sign` of [[normalize-bounds]] to a result, which is a number or an info map with `:result`."
  [^double sign result]
  (cond
    (m/one? sign) result
    (map? result) (update result :result (fn [^double r] (m/* sign r)))
    :else (m/* sign (double result))))

;; univariate improper integrals

(defn- subst-upper
  [f ^double lower]
  (fn ^double [^double t]
    (let [inv (m// (m/- 1.0 t))]
      (m/* ^double (f (m/+ lower (m/* t inv))) (m/sq inv)))))

(defn- subst-lower
  [f ^double upper]
  (fn ^double [^double t]
    (let [inv (m// (m/inc t))]
      (m/* ^double (f (m/+ upper (m/* t inv))) (m/sq inv)))))

(defn- subst-both
  [f]
  (fn ^double [^double t]
    (let [t2 (m/* t t)
          den (m// (m/- 1.0 t2))]
      (m/* ^double (f (m/* t den))
           (m/* (m/inc t2)
                (m/sq den))))))

(defn subst-1d
  "Substitute integrant and bounds for improper integrals, 1-d version"
  [f ^double lower ^double upper]
  (cond
    (and (m/neg-inf? lower)
         (m/pos-inf? upper)) [(subst-both f)
                              (m/next-double -1.0)
                              (m/prev-double 1.0)]
    (m/neg-inf? lower) [(subst-lower f upper)
                        (m/next-double -1.0) 0.0]
    (m/pos-inf? upper) [(subst-upper f lower)
                        0.0 (m/prev-double 1.0)]
    :else [f lower upper]))

(defn subst-1d-inverse
  "Returns the function mapping a point of the original interval to the variable of the substitution used by [[subst-1d]] for the same bounds."
  [^double lower ^double upper]
  (cond
    (and (m/neg-inf? lower)
         (m/pos-inf? upper)) (fn ^double [^double x]
                               (m// (m/* 2.0 x) (m/inc (m/sqrt (m/inc (m/* 4.0 x x))))))
    (m/neg-inf? lower) (fn ^double [^double x]
                         (let [y (m/- x upper)]
                           (m// y (m/- 1.0 y))))
    (m/pos-inf? upper) (fn ^double [^double x]
                         (let [y (m/- x lower)]
                           (m// y (m/inc y))))
    :else (fn ^double [^double x] x)))

;; multivariate

(defn- subst-upper-v [^double lower] (fn ^double [^double vi] (m/+ lower (m// vi (m/- 1.0 vi)))))
(defn- subst-lower-v [^double upper] (fn ^double [^double vi] (m/+ upper (m// vi (m/+ 1.0 vi)))))
(defn- subst-both-v ^double [^double vi] (m// vi (m/- 1.0 (m/* vi vi))))
(defn- subst-multiplier-upper ^double [^double vi] (m// (m/sq (m/- 1.0 vi))))
(defn- subst-multiplier-lower ^double [^double vi] (m// (m/sq (m/+ 1.0 vi))))
(defn- subst-multiplier-both ^double [^double vi]  (let [vi2 (m/* vi vi)]
                                                 (m// (m/+ 1.0 vi2)
                                                      (m/sq (m/- 1.0 vi2)))))

(defn- subst-multi-
  [lower upper]
  (let [[subst-v subst-m nbounds]
        (->> (map (fn [^double a ^double b]
                    (cond
                      (and (m/inf? a)
                           (m/inf? b)) [subst-both-v subst-multiplier-both [-1.0 1.0]]
                      (m/neg-inf? a) [(subst-lower-v b) subst-multiplier-lower [-1.0 0.0]]
                      (m/pos-inf? b) [(subst-upper-v a) subst-multiplier-upper [0.0 1.0]]
                      :else [identity (constantly 1.0) [a b]])) lower upper)
             (apply map vector))]
    
    [(fn [v] (mapv (fn [f ^double vi] (f vi)) subst-v v))
     (fn [v] (->> (mapv (fn [f ^double vi] (f vi)) subst-m v)
                 (reduce m/*)))
     (apply map vector nbounds)]))

(defn subst-multi
  "Substitute integrant and bounds for improper integrals, multivariate version"
  [f lower upper]
  (if (and (every? m/valid-double? lower)
           (every? m/valid-double? upper))
    [f lower upper]
    (let [[subst-v multiplier [nlower nupper]] (subst-multi- lower upper)]
      [(fn [v]
         (let [nv (subst-v v)
               m (double (multiplier v))]
           (m/* m (double (f nv)))))
       nlower nupper])))
