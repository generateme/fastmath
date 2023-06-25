(ns utils
  "Notebook utilities"
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [clojure2d.core :as c2d]
            [nextjournal.clerk :as clerk]
            [nextjournal.clerk.viewer :as viewer]
            [nextjournal.clerk.tap]
            [fastmath.random :as r]
            [clojure2d.color :as c]
            [clojure.string :as str]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]))

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

(defn- tick-step
  [a b]
  (max 1 (long (m/ceil (m/log10 (m/abs (- a b)))))))

(defn fgraph-int
  ([f] (fgraph-int f nil))
  ([f domain] (fgraph-int f domain :domain))
  ([f domain rnge] (fgraph-int f domain rnge size))
  ([f domain rnge size]
   (let [[dx dy :as dom] (or domain [-3.3 3.3])
         md (m/make-norm dx dy 5.0 (- size 5))
         xsys (m/sample f dx dy (* 2 size) true)
         [frx fry] (stats/extent (map second (filter (comp m/valid-double? second) xsys)))
         [rx ry] (case rnge
                   :domain dom
                   :zero [(min 0.0 frx)
                          (max 0.0 fry)]
                   rnge)
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
       (doseq [dt (take-nth (tick-step dx dy) dticks)]
         (c2d/line c dt 0 dt size))
       (doseq [rt (take-nth (tick-step rx ry) rticks)]
         (c2d/line c 0 rt size rt))
       (c2d/set-color c 0x303065)
       (doseq [p (split-at-invalid-double xsys (* 1.1 (m/max (m/abs rx) (m/abs ry))))]
         (c2d/path c (map (fn [[x y]] [(md x) (mr y)]) p)))
       (c2d/get-image c)))))

(defn common-extent
  [xsyss]
  (->> xsyss
       (map (fn [xsys] (stats/extent (map second (filter (comp m/valid-double? second) xsys)))))
       (reduce (fn [[mn mx] [y1 y2]]
                 [(min mn y1) (max mx y2)]) [##Inf ##-Inf])))

(defn fgraphs-int
  ([fs] (fgraphs-int fs nil))
  ([fs domain] (fgraphs-int fs domain :domain))
  ([fs domain rnge] (fgraphs-int fs domain rnge size))
  ([fs domain rnge size]
   (let [[dx dy :as dom] (or domain [-3.3 3.3])
         md (m/make-norm dx dy 5.0 (- size 5))
         xsyss (map (fn [f] (m/sample f dx dy (* 2 size) true)) fs)
         [frx fry] (common-extent xsyss)
         [rx ry] (case rnge
                   :domain dom
                   :zero [(min 0.0 frx)
                          (max 0.0 fry)]
                   rnge)
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
       (doseq [dt (take-nth (tick-step dx dy) dticks)]
         (c2d/line c dt 0 dt size))
       (doseq [rt (take-nth (tick-step rx ry) rticks)]
         (c2d/line c 0 rt size rt))
       (doseq [[xsys color] (map vector xsyss (c/palette :category10))]
         (c2d/set-color c color)
         (doseq [p (split-at-invalid-double xsys (* 1.1 (m/max (m/abs rx) (m/abs ry))))]
           (c2d/path c (map (fn [[x y]] [(md x) (mr y)]) p))))
       (c2d/get-image c)))))

(defn sample-int
  [f ^long dx ^long dy]
  (map (fn [v] [(+ v 0.5) (f v)]) (range dx (inc dy))))

(defn bgraph-int
  ([f] (bgraph-int f nil))
  ([f domain] (bgraph-int f domain :domain))
  ([f domain rnge] (bgraph-int f domain rnge size))
  ([f domain rnge size]
   (let [[dx dy :as dom] (or domain [0 10])
         md (m/make-norm dx dy 5.0 (- size 5))
         xsys (if (map? f)
                (seq f)
                (sample-int f dx dy))
         [frx fry] (stats/extent (map second (filter (comp m/valid-double? second) xsys)))
         [rx ry] (case rnge
                   :domain dom
                   :zero [(min 0.0 frx)
                          (max 0.0 fry)]
                   rnge)
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
       (doseq [dt (take-nth (tick-step dx dy) dticks)]
         (c2d/line c dt 0 dt size))
       (doseq [rt (take-nth (tick-step rx ry) rticks)]
         (c2d/line c 0 rt size rt))
       (c2d/set-color c 0x303065)
       (doseq [[x y] xsys
               :when (m/valid-double? y)
               :let [xx (md x)]]
         (c2d/line c xx r0 xx (mr y)))
       (c2d/get-image c)))))

(def grad (c/gradient [0xfafaff 0x303065]))

(defn graph2d
  ([f] (graph2d f nil nil))
  ([f dx dy] (graph2d f dx dy size))
  ([f dx dy size]
   (let [[dx1 dx2] (or dx [0.0 1.0])
         [dy1 dy2] (or dy [0.0 1.0])
         mx (m/make-norm 5.0 (- size 5.0) dx1 dx2)
         my (m/make-norm 5.0 (- size 5.0) dy2 dy1)
         r (range 5 (- size 5))
         xy (for [x r y r
                  :let [xx (mx x)
                        yy (my y)
                        v (f xx yy)]
                  :when (m/valid-double? v)]
              [x y v])
         [v1 v2] (stats/extent (map last xy))
         n (m/make-norm v1 v2 0.0 1.0)]
     (c2d/with-canvas [c (c2d/canvas size size :high)]
       (c2d/set-background c 0xfafaff)
       (doseq [[x y v] xy]
         (c2d/set-color c (grad (n v)))
         (c2d/crect c x y 1 1))
       (c2d/get-image c)))))

(defn scatter
  ([xy] (scatter xy nil nil))
  ([xy dx dy] (scatter xy dx dy size))
  ([xy dx dy size]
   (let [[dx1 dx2] (or dx [0.0 1.0])
         [dy1 dy2] (or dy [0.0 1.0])
         mx (m/make-norm dx1 dx2 5.0 (- size 5.0))
         my (m/make-norm dy2 dy1 5.0 (- size 5.0))
         d0 (mx 0.0)
         r0 (my 0.0)]
     (c2d/with-canvas [c (c2d/canvas size size :highest)]
       (c2d/set-background c 0xfafaff)
       (c2d/set-color c 0x77303065)
       (c2d/line c 0 r0 size r0)
       (c2d/line c d0 0 d0 size)
       (doseq [[x y] xy]
         (c2d/ellipse c (mx x) (my y) 3 3))
       (c2d/get-image c)))))

(comment (graph2d r/noise))

(defmacro fgraph
  ([f] `(fgraph-int (fn [x#] (~f x#))))
  ([f domain] `(fgraph-int (fn [x#] (~f x#)) ~domain))
  ([f domain rnge] `(fgraph-int (fn [x#] (~f x#)) ~domain ~rnge))
  ([f domain rnge size] `(fgraph-int (fn [x#] (~f x#)) ~domain ~rnge ~size)))

(defn dgraph-cont
  ([pdf-graph distr] (dgraph-cont pdf-graph distr nil))
  ([pdf-graph distr {:keys [pdf icdf data]
                     :or {pdf [0 1] icdf [0 0.999]}}]
   (let [pdf-img (pdf-graph (if data data (partial r/pdf distr)) pdf [0 nil] 100)
         cdf-img (fgraph-int (partial r/cdf distr) pdf [0 1] 100)
         icdf-img (fgraph-int (partial r/icdf distr) icdf :zero 100)]
     (c2d/with-canvas-> (c2d/canvas 310 120)
       (c2d/set-background :white)
       (c2d/set-color :black)
       (c2d/set-font-attributes 10)
       (c2d/text "PDF" 0 10)
       (c2d/text "CDF" 105 10)
       (c2d/text "ICDF" 210 10)
       (c2d/image pdf-img 0 20)
       (c2d/image cdf-img 105 20)
       (c2d/image icdf-img 210 20)
       (c2d/get-image)))))

(def dgraph (partial dgraph-cont fgraph-int))
(def dgraphi (partial dgraph-cont bgraph-int))

(tap> (fgraph m/sin))

(defn update-child-viewers [viewer f]
  (update viewer :transform-fn (fn [trfn] (comp #(update % :nextjournal/viewers f) trfn))))

(def unpaginated-table
  (update-child-viewers viewer/table-viewer (fn [vs] (nextjournal.clerk.viewer/update-viewers vs {:page-size #(dissoc % :page-size)}))))

(defmacro table
  [rows]
  `(clerk/with-viewer unpaginated-table {:head ["symbol" "macro?" "info"]
                                         :rows [~@(for [[s m? i] rows]
                                                    `[(clerk/code (quote ~s)) ~m? (or ~i "-")])]}))

(defmacro table2
  [rows]
  `(clerk/with-viewer unpaginated-table {:head ["symbol" "info"]
                                         :rows [~@(for [[s i] rows]
                                                    `[(clerk/code (quote ~s)) (or ~i "-")])]}))

(defmacro table3
  [last-col-name rows]
  `(clerk/with-viewer unpaginated-table {:head ["symbol" "info" ~last-col-name]
                                         :rows [~@(for [[s i f] rows]
                                                    `[(clerk/code (quote ~s))
                                                      (or ~i "-")
                                                      (clerk/md ~(or (str f) ""))])]}))

;; histogram
(defn h->bars
  [data & params]
  (let [{:keys [bins step]} (apply stats/histogram data params)]
    (into {} (map (fn [[k v]]
                    [(+ k (* step 0.5)) v]) bins))))

(def source-files "https://github.com/generateme/fastmath/blob/master/src/")

(defn args->call [f args] (conj (seq args) f))
(defn fix-tex [s]  (str/replace s #"\\\\\(|\\\\\)" "\\$"))
(defn fix-anchor [s] (str/replace s #"\[\[(.+?)\]\]" "[$1](#LOS-$1)"))

(defn make-public-fns-table
  [ns]
  (clerk/html
   [:div
    [:div
     [:span [:b {:class "underline decoration-2 decoration-gray-400"} (name ns)] " namespace"]
     [:p]
     [:div (clerk/md (->> (or (:doc (meta (the-ns ns))) "\n") fix-tex fix-anchor))]]
    [:div (for [v (->> (ns-publics ns)
                       (sort-by first)
                       (map second))
                :let [{:keys [name file macro arglists doc line const deprecated]} (meta v)]]
            [:div {:class "pb-8" :id (str "LOS-" (clojure.core/name name))} ;; add border
             
             (clerk/html [:span [:b {:class (str "underline decoration-2 decoration-gray-400"
                                                 (when deprecated " text-gray-400"))} name]
                          (when macro [:sup " MACRO"])
                          (when const [:sup " CONST"])
                          [:sup [:a {:href (str source-files file "#L" line)} "  [source]"]]])
             (when deprecated [:div [:i (if (string? deprecated)
                                          (clerk/md (str "Deprecated: "(->> deprecated fix-tex fix-anchor)))
                                          "Deprecated")]])
             [:p]
             (when const [:div (clerk/code (var-get v)) [:p]])
             (when arglists [:div
                             [:div (for [args arglists]
                                     (clerk/code (args->call name args)))]
                             [:p]])
             [:div (clerk/md (->> (or doc "\n") fix-tex fix-anchor))]])]]))

;; csv

(defn transform [spec data]
  (if (sequential? spec)
    (let [[nm f] spec] [nm (f data)])
    [spec (read-string data)]))
(defn zip [header row] (into {} (map transform header row)))
(defn parse-rows  [header rows]  (map (partial zip header) rows))
(defn read-csv
  [fname header]
  (->> (io/resource fname)
       (slurp)
       (csv/read-csv)
       (rest)
       (parse-rows header)))

(defn data->fn
  [data]
  (memoize (fn
             ([] data)
             ([selector] (map selector data))
             ([selector filter-pred] (map selector (filter filter-pred data))))))

(defn by
  ([data f]
   (mapv second (sort-by first (group-by f (data)))))
  ([data f selector]
   (vec (for [group (by data f)]
          (map selector group)))))

(def iris (data->fn (read-csv "iris.csv"
                            [:sepal-length :sepal-width :petal-length :petal-width [:species keyword]])))
(def mtcars (data->fn (read-csv "mtcars.csv"
                              [[:name identity] :mpg :cyl :disp :hp :drat :wt :qsec :vs :am :gear :carb])))

^{::clerk/visibility :hide
  ::clerk/viewer :hide-result}
(comment
  (clerk/serve! {:browse? false :watch-paths ["notebooks"]})
  (clerk/serve! {:browse? false})
  (clerk/show! 'nextjournal.clerk.tap)
  ;; (clerk/show! "notebooks/core.clj")
  (clerk/build! {:browse? false :paths ["notebooks/core.clj"
                                        "notebooks/random.clj"
                                        "notebooks/stats.clj"
                                        "notebooks/bootstrap.clj"]
                 :out-path "docs/notebooks/"})
  (clerk/clear-cache!)
  (clerk/halt!))
