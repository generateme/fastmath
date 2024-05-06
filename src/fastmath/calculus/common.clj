(ns fastmath.calculus.common
  (:require [fastmath.core :as m]))

(defn subst-upper
  [f ^double lower]
  (fn ^double [^double t]
    (let [inv (m// (m/- 1.0 t))]
      (m/* ^double (f (m/+ lower (m/* t inv))) (m/sq inv)))))

(defn subst-lower
  [f ^double upper]
  (fn ^double [^double t]
    (let [inv (m// (m/inc t))]
      (m/* ^double (f (m/+ upper (m/* t inv))) (m/sq inv)))))

(defn subst-both
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

;; multivariate

(defn subst-upper-v [^double lower] (fn ^double [^double vi] (m/+ lower (m// vi (m/- 1.0 vi)))))
(defn subst-lower-v [^double upper] (fn ^double [^double vi] (m/+ upper (m// vi (m/+ 1.0 vi)))))
(defn subst-both-v ^double [^double vi] (m// vi (m/- 1.0 (m/* vi vi))))
(defn subst-multiplier-upper ^double [^double vi] (m// (m/sq (m/- 1.0 vi))))
(defn subst-multiplier-lower ^double [^double vi] (m// (m/sq (m/+ 1.0 vi))))
(defn subst-multiplier-both ^double [^double vi]  (let [vi2 (m/* vi vi)]
                                                 (m// (m/+ 1.0 vi2)
                                                      (m/sq (m/- 1.0 vi2)))))

(defn subst-multi-
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
               ^double m (multiplier v)]
           (m/* m ^double (f nv))))
       nlower nupper])))
