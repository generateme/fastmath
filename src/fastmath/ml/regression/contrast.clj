(ns fastmath.ml.regression.contrast
  "Categorical data manipulation and various encoding methods.

  A coding (contrast) replaces a categorical variable with `n-1` numeric columns, which can be used as terms of a regression model. A coding is a map with the keys `:levels`, `:names` (names of the columns) and `:mapping` (a map from a level to a vector of `n-1` coefficients).

  Levels are collected from the data with [[->levels]]. The order of levels matters: the first level is the reference of [[dummy]] and [[simple]], the last one is dropped by [[deviation]]. To use a different reference, as `contr.treatment(base=)` of R, pass an explicit `order` to [[->levels]]. A coding is applied to a column with [[encode]], the meaning of the coefficients is returned by [[mean-contrasts]].

  Codings: [[dummy]], [[simple]], [[deviation]], [[polynomial]], [[helmert]], [[reverse-helmert]], [[reverse-helmert-unscaled]], [[forward-difference]], [[backward-difference]] and [[difference]]."
  (:require [fastmath.stats :as stats]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; https://stats.oarc.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
;; https://cran.r-project.org/web/packages/codingMatrices/vignettes/codingMatrices.pdf

(defn ->levels
  "Returns the levels of a categorical variable.

  Parameters:

  - `cs` - a sequence of categorical values
  - `order` (optional) - a sequence of all levels, in the required order (the first level is the reference level of the codings)

  Without `order` the levels are returned in the order of the first appearance in `cs` (R sorts them instead).

  Throws `ex-info` when `order` does not contain exactly the levels of `cs`.

  See also [[dummy]], [[encode]]."
  ([cs] (->levels cs nil))
  ([cs order]
   (let [dcs (distinct cs)]
     (if (and order (not= (set dcs) (set order)))
       (throw (ex-info "order attribute should contain all levels from the input"
                       {:levels dcs :order order}))
       (or order dcs)))))

(defn- ->str [v] (if (or (symbol? v) (keyword? v)) (name v) v))

(defn- check-levels
  "Checks if levels can be coded: an ordered sequence of at least two different values, with different names."
  [levels]
  (when-not (sequential? levels)
    (throw (ex-info "Levels should be an ordered sequence, eg. a vector."
                    {:levels levels})))
  (when (or (m/< (count levels) 2)
            (not (m/== (count levels) (count (set levels)))))
    (throw (ex-info "Levels should contain at least two different values."
                    {:levels levels})))
  (when-not (m/== (count levels) (count (set (map (comp str ->str) levels))))
    (throw (ex-info "Levels should have different names (the text of a level is its name)."
                    {:levels levels})))
  levels)

(defn dummy
  "Dummy (treatment) coding.

  Compares each level to the reference (first) level. The intercept is the mean of the first level.

  Parameters:

  - `levels` - a sequence of at least two different levels, the first one is the reference

  Returns a coding: a map with `:levels`, `:names` (`level.vs.reference`) and `:mapping`. The reference level is a row of zeros, other levels have a single `1.0`.

  It is `contr.treatment` of R. Use [[->levels]] with `order` to select the reference level: the last level as the reference is `contr.SAS`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[simple]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [[l0 & ls] levels
        n (count levels)
        zero (vec (repeat (m/long-dec n) 0.0))
        nl0 (->str l0)]
    {:levels levels
     :names (map (fn [l] (str l ".vs." nl0)) (map ->str ls))
     :mapping (->> (map-indexed vector ls)
                   (reduce (fn [m [id l]]
                             (assoc m l (assoc zero id 1.0))) {l0 zero}))}))

(defn simple
  "Simple coding.

  Compares each level to the reference (first) level. The intercept is the mean of all levels.

  Parameters:

  - `levels` - a sequence of at least two different levels, the first one is the reference

  Returns a coding: a map with `:levels`, `:names` (`level.vs.reference`) and `:mapping`. Each column contains `1-1/n` for its level and `-1/n` for others.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[dummy]], [[deviation]], [[encode]], [[mean-contrasts]]."
  ([levels]
   (check-levels levels)
   (let [[l0 & ls] levels
         n (count levels)
         n- (m/long-dec n)
         init (vec (repeat n- (m// -1.0 n)))
         v (m// (double n-) n)
         nl0 (->str l0)]
     {:levels levels
      :names (map (fn [l] (str l ".vs." nl0)) (map ->str ls))
      :mapping (->> (map-indexed vector ls)
                    (reduce (fn [m [id l]]
                              (assoc m l (assoc init id v))) {l0 init}))})))

(defn deviation
  "Deviation (sum) coding.

  Compares each level (except the last one) to the mean of all levels. The intercept is the mean of all levels.

  Parameters:

  - `levels` - a sequence of at least two different levels, the last one is dropped

  Returns a coding: a map with `:levels`, `:names` (`level.vs.all`) and `:mapping`. The last level is a row of `-1.0`.

  It is `contr.sum` of R.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[simple]], [[helmert]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)
        n- (m/long-dec n)
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

  Creates `n-1` orthonormal polynomials of increasing degree (linear, quadratic, cubic, ...) over equally spaced levels, which should be ordered. The intercept is the mean of all levels.

  Parameters:

  - `levels` - a sequence of at least two different, ordered levels

  Returns a coding: a map with `:levels`, `:names` (`^1`, `^2`, ...) and `:mapping`. The columns are orthonormal.

  It is `contr.poly` of R.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[deviation]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)
        poly (->> (range n)
                  (stats/demean)
                  (v/orthonormal-polynomials)
                  (apply map vector))]
    {:levels levels
     :names (map #(str "^" %) (range 1 n))
     :mapping (zipmap levels poly)}))

(defn- helmert-vectors
  ;; keep `(double ...)` below even if clj-kondo reports it as redundant: `m//` of two longs is an integer division
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

  Compares each level (except the last one) to the mean of all following levels. The intercept is the mean of all levels. The coefficients of the model are the differences of the means.

  Parameters:

  - `levels` - a sequence of at least two different levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.later`) and `:mapping`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[reverse-helmert]], [[reverse-helmert-unscaled]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)
        h (->> (helmert-vectors n)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str (butlast levels))
                 (map #(str % ".vs.later")))
     :mapping (zipmap levels h)}))

(defn reverse-helmert
  "Reverse Helmert coding.

  Compares each level (except the first one) to the mean of all preceding levels. The intercept is the mean of all levels. The coefficients of the model are the differences of the means.

  Parameters:

  - `levels` - a sequence of at least two different levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.earlier`) and `:mapping`. The k-th column is the k-th column of [[reverse-helmert-unscaled]] divided by `k+1`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[helmert]], [[reverse-helmert-unscaled]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)
        h (->> (helmert-vectors n)
               (map reverse)
               (reverse)
               (apply map vector))]
    {:levels levels
     :names (->> (map ->str (rest levels))
                 (map #(str % ".vs.earlier")))
     :mapping (zipmap levels h)}))

(defn reverse-helmert-unscaled
  "Reverse Helmert coding with integer coefficients.

  Compares each level (except the first one) to the mean of all preceding levels. The intercept is the mean of all levels. The coefficients of the model are the differences of the means divided by the number of levels compared (`k+1` for the k-th column).

  Parameters:

  - `levels` - a sequence of at least two different levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.earlier`) and `:mapping`. The k-th column contains `-1.0` for the first `k` levels, `k` for the level k and `0.0` for the following ones.

  It is `contr.helmert` of R. [[reverse-helmert]] contains the same columns divided by `k+1`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[helmert]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)]
    {:levels levels
     :names (->> (map ->str (rest levels))
                 (map #(str % ".vs.earlier")))
     :mapping (zipmap levels (for [^long id (range n)]
                               (vec (for [^long c (range 1 n)]
                                      (cond
                                        (m/< id c) -1.0
                                        (m/== id c) (double c)
                                        :else 0.0)))))}))

(defn- difference-vectors
  ;; keep `(double ...)` below, see `helmert-vectors`
  [^long n]
  (for [^long id (range (m/dec n))
        :let [v1 (m// (double (m/dec (m/- n id))) n)
              v2 (m// (m/* -1.0 (m/inc id)) n)]]
    (for [^long cid (range n)]
      (if (m/<= cid id) v1 v2))))

(defn forward-difference
  "Forward difference coding.

  Compares the means of pairs of consecutive levels, each level to the following one. The intercept is the mean of all levels.

  Parameters:

  - `levels` - a sequence of at least two different, ordered levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.next-level`) and `:mapping`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[backward-difference]], [[difference]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
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

  Compares the means of pairs of consecutive levels, each level to the preceding one. The intercept is the mean of all levels.

  Parameters:

  - `levels` - a sequence of at least two different, ordered levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.previous-level`) and `:mapping`.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[forward-difference]], [[difference]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
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
  "Cumulative difference coding.

  Compares each level to the preceding one, the same comparisons as [[backward-difference]]. The columns are cumulative sums of dummy columns, so the intercept is the mean of the first level (not the mean of all levels).

  Parameters:

  - `levels` - a sequence of at least two different, ordered levels

  Returns a coding: a map with `:levels`, `:names` (`level.vs.previous-level`) and `:mapping`. The first level is a row of zeros, the level with the index k (counted from 0) has `1.0` in the first k columns.

  Throws `ex-info` for fewer than two levels or duplicated levels.

  See also [[backward-difference]], [[dummy]], [[encode]], [[mean-contrasts]]."
  [levels]
  (check-levels levels)
  (let [n (count levels)
        n- (m/long-dec n)
        zero (vec (repeat n- 0.0))
        d (reductions (fn [v ^long id]
                        (assoc v id 1.0)) zero (range n-))]
    {:levels levels
     :names (->> (map ->str levels)
                 (partition 2 1)
                 (map (fn [[a b]] (str b ".vs." a))))
     :mapping (zipmap levels d)}))

(defn encode
  "Applies a coding to a categorical column.

  Parameters:

  - `coding` - a coding returned by one of the coding functions of this namespace
  - `column` - a sequence of categorical values, each of them has to be one of the `:levels` of the coding

  Returns a vector of rows, one vector of coefficients for each value of `column`. The rows can be used as terms of a regression model.

  Throws `ex-info` with the `:value` and the `:levels` when a value is not a level of the coding.

  See also [[->levels]], [[mean-contrasts]]."
  [{:keys [levels mapping]} column]
  (mapv (fn [value]
          (or (get mapping value)
              (throw (ex-info "Value is not a level of the coding."
                              {:value value :levels levels}))))
        column))

(defn mean-contrasts
  "Returns mean contrasts for given coding scheme.

  The result shows how the coefficients of a model with the given coding are related to the means of the levels: each coefficient is the sum of the means of levels weighted by the vector assigned to the coefficient.

  Parameters:

  - `coding` - a coding returned by one of the coding functions of this namespace
  - `approx?` (default: `true`) - if `true` or a number, values are rounded to the given number of decimal digits (`10` for `true`)

  Returns a map from the name of the coefficient (and `:$intercept`) to a vector of weights of the means of levels (in the order of `:levels`)."
  ([coding] (mean-contrasts coding true))
  ([{:keys [names levels  mapping]} approx?]
   (let [mc (->> (map mapping levels)
                 (map #(conj (seq %) 1.0))
                 (mat/rows->RealMatrix)
                 (mat/inverse)
                 (mat/rows)
                 (map v/vec->Vec))
         nms (cons :$intercept names)
         approx (if (and approx? (number? approx?)) approx? 10 )]
     (zipmap nms (if approx? (map #(v/approx % approx) mc) mc)))))

