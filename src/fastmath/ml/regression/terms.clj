(ns fastmath.ml.regression.terms
  "Builds mappings and penalties from list of additive terms."
  (:require [fastmath.matrix :as mat]))


(defn intercept
  []
  {:mapping (constantly 1)
   :penalty {:lambda 0.0
             :D (mat/mat->RealMatrix 0.0)}})

(defn categorical
  [xs {:keys [encoding order]
       :or {encoding :dummy}}])
