(ns fastmath.dev.helpers
  (:require [fastmath.core :as m]))

(defn relative-error [f1 f2]
  (fn ^double [^double x]
    (m/relative-error (f1 x) (f2 x))))
