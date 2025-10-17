(ns fastmath.ml.regression.contrast
  "Categorical data manipulation and various encoding methods."
  (:require [fastmath.stats :as stats]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; https://stats.oarc.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
;; https://cran.r-project.org/web/packages/codingMatrices/vignettes/codingMatrices.pdf

(defn ->levels
  "Returns list of levels from given sequence of categorical values.

  If `order` is passed it will be used as a reference."
  ([cs] (->levels cs nil))
  ([cs order]
   (let [dcs (distinct cs)]
     (if (and order (not= (set dcs) (set order)))
       (throw (ex-info "order attribute should contain all levels from the input"
                       {:levels dcs :order order}))
       (or order dcs)))))

(defn- ->str [v] (if (or (symbol? v) (keyword? v)) (name v) v))

(defn dummy
  "Dummy coding.

  Compares each level to a reference (first) level. Intercept is the mean of the first level."
  [[l0 & ls :as levels]]
  (let [n (count levels)
        zero (vec (repeat (m/dec n) 0.0))
        nl0 (->str l0)]
    {:levels levels
     :names (map (fn [l] (str l ".vs." nl0)) (map ->str ls))
     :mapping (->> (map-indexed vector ls)
                   (reduce (fn [m [id l]]
                             (assoc m l (assoc zero id 1.0))) {l0 zero}))}))

(defn simple
  "Simple coding.

  Compares each level to a reference (first) level. Intercept is the mean of all levels."
  ([[l0 & ls :as levels]]
   (let [n (count levels)
         n- (m/dec n)
         init (vec (repeat n- (m// -1.0 n)))
         v (m// (double n-) n)
         nl0 (->str l0)]
     {:levels levels
      :names (map (fn [l] (str l ".vs." nl0)) (map ->str ls))
      :mapping (->> (map-indexed vector ls)
                    (reduce (fn [m [id l]]
                              (assoc m l (assoc init id v))) {l0 init}))})))

(defn deviation
  "Deviation coding.

  Compares each (but last) levels to the mean. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        n- (m/dec n)
        zero (vec (repeat n- 0.0))]
    {:levels levels
     :names (->> (butlast levels)
                 (map ->str)
                 (map #(str % ".vs.all")))
     :mapping (->> (butlast levels)
                   (map-indexed vector)
                   (reduce (fn [m [id l]]
                             (assoc m l (assoc zero id 1.0))) {(last levels) (vec (repeat n- -1.0))}))}))

(defn polynomial
  "Orthonormal polynomials coding.

  n-1 orthonormal polynomials. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        poly (->> (range n)
                  (stats/demean)
                  (v/orthonormal-polynomials)
                  (apply map vector))]
    {:levels levels
     :names (map #(str "^" %) (range 1 n))
     :mapping (zipmap levels poly)}))

(defn- helmert-vectors
  [^long n]
  (let [zero (vec (repeat n 0.0))]
    (for [^long id (range (m/dec n))
          :let [den (m/- n id)
                v (m// -1.0 den)]]
      (->> (range (m/inc id) n)
           (reduce (fn [m ^long cid]
                     (assoc m cid v)) (assoc zero id (m// (double (m/dec den)) den)))))))

(defn helmert
  "Helmert coding.

  Compares each (but last) levels to all following levels. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        h (->> (helmert-vectors n)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str (butlast levels))
                 (map #(str % ".vs.later")))
     :mapping (zipmap levels h)}))

(defn reverse-helmert
  "Reverse Helmert coding.

  Compares each (but first) levels to all preceding levels. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        h (->> (helmert-vectors n)
               (map reverse)
               (reverse)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str (rest levels))
                 (map #(str % ".vs.earlier")))
     :mapping (zipmap levels h)}))

(defn- difference-vectors
  [^long n]
  (for [^long id (range (m/dec n))
        :let [v1 (m// (double (m/dec (m/- n id))) n)
              v2 (m// (m/* -1.0 (m/inc id)) n)]]
    (for [^long cid (range n)]
      (if (m/<= cid id) v1 v2))))

(defn forward-difference
  "Forward difference coding.

  Comperes pairs of consecutive levels. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        d (->> (difference-vectors n)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str levels)
                 (partition 2 1)
                 (map (fn [[a b]] (str a ".vs." b))))
     :mapping (zipmap levels d)}))

(defn backward-difference
  "Backward difference coding.

  Compares reversed pairs of consecutive levels. Intercept is the mean of all levels."
  [levels]
  (let [n (count levels)
        d (->> (difference-vectors n)
               (map v/sub)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str levels)
                 (partition 2 1)
                 (map (fn [[a b]] (str b ".vs." a))))
     :mapping (zipmap levels d)}))

(defn difference
  "Backward difference coding.

  Compares reversed pairs of consecutive levels. Intercept is the mean of the first level."
  [levels]
  (let [n (count levels)
        n- (m/dec n)
        zero (vec (repeat n- 0.0))
        d (reductions (fn [v ^long id]
                        (assoc v id 1.0)) zero (range n-))]
    {:levels levels
     :names (->> (map ->str levels)
                 (partition 2 1)
                 (map (fn [[a b]] (str b ".vs." a))))
     :mapping (zipmap levels d)}))


(defn mean-contrasts
  "Returns mean contrasts for given coding scheme.

  If `approx?` is true or a number, approximate values to given number of decimal digits (default: `true`)"
  ([coding] (mean-contrasts coding true))
  ([{:keys [names levels  mapping]} approx?]
   (let [mc (->> (map mapping levels)
                 (map #(conj (seq %) 1.0))
                 (mat/rows->RealMatrix)
                 (mat/inverse)
                 (mat/rows)
                 (map v/vec->Vec))
         nms (conj names :$intercept)
         approx (if (and approx? (number? approx?)) approx? 10 )]
     (zipmap nms (if approx? (map #(v/approx % approx) mc) mc)))))

