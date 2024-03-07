(ns fastmath.interpolation.shepard
  (:require [fastmath.distance :as d]
            [fastmath.stats :as stats]
            [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
#_(set! *warn-on-reflection* true)

(defn shepard
  [xss ys {:keys [^double p distance]
           :or {p 2.0}}]
  (let [dist (or distance (if (number? (first xss)) d/euclidean-1d d/euclidean))
        cache (into {} (map vector xss ys))
        p- (m/- p)]
    (fn ^double [xs]
      (or (cache xs)
          (->> (map (fn [xs']
                      (m/pow (dist xs xs') p-)) xss)
               (stats/wmean ys))))))
