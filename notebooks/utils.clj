(ns utils
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [clojure2d.core :as c2d]))

(def ^:const size 200)

(defn valid? [^double thr ^double v]
  (and (m/valid-double? v)
       (<= (m/abs v) thr)))

(defn invalid? [^double thr ^double v]
  (or (m/invalid-double? v)
      (> (m/abs v) thr)))

(defn- split-at-invalid-double
  [xs ^double thr]
  (loop [s xs
         buff []]
    (if-not (seq s)
      buff
      (let [[v inv] (split-with (comp (partial valid? thr) second) s)]
        (recur (drop-while (comp (partial invalid? thr) second) inv)
               (if (seq v) (conj buff v) buff))))))

(defn fgraph-int
  ([f] (fgraph-int f nil))
  ([f domain] (fgraph-int f domain :domain))
  ([f domain rnge]
   (let [[dx dy :as dom] (or domain [-3.3 3.3])
         [rx ry] (if (= :domain rnge) dom rnge)
         md (m/make-norm dx dy 5.0 (- size 5))
         xsys (m/sample f dx dy (* 2 size) true)
         [frx fry] (stats/extent (map second (filter (comp m/valid-double? second) xsys)))
         rx (or rx frx)
         ry (or ry fry)
         mr (m/make-norm rx ry (- size 5.0) 5.0)         
         r0 (mr 0)
         d0 (md 0)
         dticks (map md (remove zero? (range (m/qfloor dx) (m/qceil dy))))
         rticks (map mr (remove zero? (range (m/qfloor rx) (m/qceil ry))))]
     (c2d/with-canvas [c (c2d/canvas size size :highest)]
       (c2d/set-background c 0xfafaff)
       (c2d/set-color c 0x77303065)
       (c2d/line c 0 r0 size r0)
       (c2d/line c d0 0 d0 size)
       (c2d/set-color c 0x11303065)
       (doseq [dt dticks]
         (c2d/line c dt 0 dt size))
       (doseq [rt rticks]
         (c2d/line c 0 rt size rt))
       (c2d/set-color c 0x303065)
       (doseq [p (split-at-invalid-double xsys (* 1.1 (m/max (m/abs rx) (m/abs ry))))]
         (c2d/path c (map (fn [[x y]] [(md x) (mr y)]) p)))
       (c2d/get-image c)))))

(defmacro fgraph
  ([f] `(fgraph-int (fn [x#] (~f x#))))
  ([f domain] `(fgraph-int (fn [x#] (~f x#)) ~domain))
  ([f domain rnge] `(fgraph-int (fn [x#] (~f x#)) ~domain ~rnge)))

(comment (fgraph m/sin))
