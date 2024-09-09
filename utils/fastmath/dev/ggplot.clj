(ns fastmath.dev.ggplot
  "Fastmath charts"
  (:require [clojisr.v1.r :as r]
            [clojisr.v1.applications.plotting :as rplot]
            [tablecloth.api :as tc]
            [fastmath.core :as m]
            [fastmath.random :as fr]
            [fastmath.stats :as stats]))

(r/require-r '[ggplot2 :as gg]
             '[paletteer :as pal]
             '[patchwork :as pw]
             '[base])

(defn ->palette
  [pal]
  (if (keyword? pal)
    (if-let [n (namespace pal)]
      (str n "::" (name pal))
      (name pal))
    pal))

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn paletter-d [pal]
  (r/r->clj (pal/paletteer_d (->palette pal))))

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn paletter-c [pal n]
  (r/r->clj (pal/paletteer_c (->palette pal) n)))

(def color-main "#000e60")
(def color-light "#87a5ff")
(def palette-blue-0 [color-main "#3c5fff", "#526dff", "#637bff", "#7189ff", "#7d97ff",
                   "#87a5ff", "#91b3ff", "#9ac2ff", "#a2d0ff", "#aadfff", "#b1eeff"])
(def palette-blue-1 (repeat 12 color-main))

(defn ->file
  ([obj] (->file obj nil))
  ([obj opts] (->file obj "img.png" opts))
  ([obj fname {:keys [width height]
               :or {width 500 height 400}}]
   (rplot/plot->file fname obj :width width :height height)))

(defn ->image
  ([obj] (->image obj nil))
  ([obj {:keys [width height]
         :or {width 500 height 400}}]
   (rplot/plot->buffered-image obj :width width :height height)))

;; (defn ->data [data] (tc/dataset data))
;; (defmacro r+ [& forms] `(r/r+ ~@forms))
;; (defmacro gg+ [& forms]  `(r/r+ (gg/ggplot) ~@forms))
;; (defmacro ggaes+ [aes & forms]  `(r/r+ (gg/ggplot :mapping ~aes) (gg/theme_light) ~@forms))
;; (defn aes [& aes-opts] (apply gg/aes aes-opts))
;; (defn line [data & opts] (apply gg/geom_line :data (tc/dataset data) opts))
;; (defn point [data & opts] (apply gg/geom_point :data (tc/dataset data) opts))
;; (defn raster [data & opts] (apply gg/geom_raster :data (tc/dataset data) opts))
;; (defn histogram [data & opts] (apply gg/geom_histogram :data (tc/dataset data) opts))
;; (defn ribbon [data & opts] (apply gg/geom_ribbon :data (tc/dataset data) opts))
;; (defn xlim [[x-min x-max]] (gg/xlim (or x-min 'NA) (or x-max 'NA)))
;; (defn ylim [[y-min y-max]] (gg/ylim (or y-min 'NA) (or y-max 'NA)))
;; (defn title [title] (gg/labs :title title))
;; (defn xlab [title] (gg/xlab title))
;; (defn ylab [title] (gg/ylab title))

;; primitives

#_{:clj-kondo/ignore [:unresolved-namespace]}
(defn aes [& opts] (apply gg/aes opts))
(defn geom-hline [obj & opts] (r/r+ obj (apply gg/geom_hline obj opts)))

;; data processing

(defn slice
  [x steps]
  (let [[min-x max-x] (or x [0.0 1.0])]
    (m/slice-range min-x max-x (or steps 400))))

(defn function->data
  ([f] (function->data f nil))
  ([f {:keys [x steps fname]
       :or {fname "function"}}]
   (->> (slice x steps)
        (map (fn [x] {:x x :y (f x) :fname fname})))))

(defn functions->data
  ([fs] (functions->data fs nil))
  ([fs opts]
   (cond
     (map? fs) (functions->data (seq fs) opts)
     (fn? (first fs)) (functions->data (map-indexed (fn [id f]
                                                      [(str "function " id) f]) fs) opts)
     :else (mapcat (fn [[n f]]
                     (function->data f (assoc opts :fname n))) fs))))

(defn function-ci->data
  ([f] (function-ci->data f nil))
  ([f {:keys [x steps fname]
       :or {fname "function"}}]
   (->> (slice x steps)
        (map (fn [x] (let [[y ymin ymax] (f x)] {:x x :y y :ymin ymin :ymax ymax :fname fname}))))))

(defn function2d->data
  ([f] (function2d->data f nil))
  ([f {:keys [x y steps varg?]
       :or {steps 200 varg? true}}]
   (let [xs (slice x steps)
         ys (slice y steps)
         f (if varg? f (fn [[x y]] (f x y)))]
     (for [x xs y ys] {:x x :y y :z (f [x y])}))))

;; plots

(defn add-common
  [object {:keys [title xlab ylab xlim ylim]}]
  (r/r+ object
        (when title (gg/labs :title title))
        (when xlab (gg/xlab xlab))
        (when ylab (gg/ylab ylab))
        (when xlim (gg/xlim xlim))
        (when ylim (gg/ylim ylim))))

(defn function
  "Single function"
  ([f] (function f nil))
  ([f {:keys [color]
       :or {color color-main}
       :as opts}]
   (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y))
             (gg/theme_light)
             (gg/geom_line :data (tc/dataset (function->data f opts)) :color color))
       (add-common opts))))

(defn line-points
  ([xs ys]
   (-> (tc/dataset {:x xs :y ys})
       (r/r+ (gg/ggplot (gg/aes :x :x :y :y))
             (gg/theme_light)
             (gg/geom_line :color "blue")))))

(defn lollipop
  ([xs ys]
   (lollipop xs ys nil))
  ([xs ys {:keys [title]
           :or {title ""}}]
   (-> (r/r+ (gg/ggplot :data (tc/dataset {:x xs :y ys}) :mapping (gg/aes :x :x :y :y))
             (gg/theme_light)
             (gg/geom_segment :mapping (gg/aes :x :x :xend :x :y 0 :yend :y) :color "blue"))
       (add-common {:title title}))))


(defn functions
  "Functions"
  ([fs] (functions fs nil))
  ([fs {:keys [palette legend-name linetype?]
        :or {palette palette-blue-0 legend-name "Functions" linetype? true}
        :as opts}]
   (let [data (functions->data fs opts)
         breaks (distinct (map :fname data))]
     (-> (r/r+ (gg/ggplot)
               (gg/theme_light)
               (gg/geom_line :data (tc/dataset data) :mapping (if linetype?
                                                                (gg/aes :linetype :fname :color :fname :x :x :y :y)
                                                                (gg/aes :color :fname :x :x :y :y)))
               (when linetype? (gg/scale_linetype_manual :name legend-name :breaks breaks :values (map inc (range (count breaks)))))
               (gg/scale_color_manual :name legend-name :breaks breaks :values palette))
         (add-common opts)))))

(defn function2d
  ([f] (function2d f nil))
  ([f {:keys [palette legend-name]
       :or {palette :pals/ocean.ice}
       :as opts}]
   (let [data (function2d->data f opts)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y :fill :z))
               (gg/theme_light)
               (gg/geom_raster :data (tc/dataset data) :interpolate true)
               (pal/scale_fill_paletteer_c :name legend-name (->palette palette)))
         (add-common opts)))))

(defn function-ci
  ([f] (function-ci f nil))
  ([f {:keys [color alpha fill]
       :or {color color-main alpha 0.4 fill color-light}
       :as opts}]
   (let [data (function-ci->data f opts)]
     (-> (r/r+ (gg/ggplot :mapping (gg/aes :x :x :y :y :ymin :ymin :ymax :ymax))
               (gg/theme_light)
               (gg/geom_ribbon :data (tc/dataset data) :alpha alpha :fill fill)
               (gg/geom_line data :color color))
         (add-common opts)))))


(defn function+scatter
  ([f xs ys] (function+scatter f xs ys nil))
  ([f xs ys {:keys [x] :as opts}]
   (let [[x-min-p x-max-p] x
         [x-min-e x-max-e] (stats/extent xs)
         x-min (or x-min-p x-min-e)
         x-max (or x-max-p x-max-e)
         data (tc/dataset {:x xs :y ys})
         ff (if (sequential? f) functions function)]
     (r/r+ (ff f (assoc opts :x [x-min x-max]))
           (gg/geom_point :color "blue" :fill "light blue" :data data :mapping (gg/aes :x :x :y :y)
                          :shape "circle filled" :size 3 :alpha 0.8)))))

(defn scatter
  "plot scatter graph.
  `xs` and `ys` are the data points.
  
  options `opts`:
  `aspect-ratio`: default `nil`.
                  `nil` means letting the scale of the coordinate-system follows
                  drawing-window's aspect ratio.
                  set to any positive whole number to set fixed ratio y/x
                  for the coordinate-system.
                  setting aspect-ratio to true or 1 is useful if we want
                  to draw plots depicting any perfect circle.
  "
  ([xs ys]
   (scatter xs ys nil))
  ([xs ys {:keys [aspect-ratio]
           :or {aspect-ratio nil}}]
   (-> (tc/dataset {:x xs :y ys})
       (gg/ggplot (gg/aes :x :x :y :y))
       (r/r+ (gg/theme_light)
             (gg/geom_point :color "blue" :fill "light blue"
                            :shape "circle filled" :size 1 :alpha 0.8))
       (cond-> 
           aspect-ratio (r/r+ (gg/coord_fixed :ratio aspect-ratio))))))


(defn scatters
  [data aes]
  (-> (tc/dataset data)
      (gg/ggplot (apply gg/aes aes))
      (r/r+ (gg/theme_light)
            (gg/geom_point :shape "circle filled" :size 3))))

(defn function2d+scatter
  ([f xs ys] (function2d+scatter f xs ys nil))
  ([f xs ys {:keys [x y] :as opts}]
   (let [[x-min-p x-max-p] x
         [x-min-e x-max-e] (stats/extent xs)
         x-min (or x-min-p x-min-e)
         x-max (or x-max-p x-max-e)
         [y-min-p y-max-p] y
         [y-min-e y-max-e] (stats/extent ys)
         y-min (or y-min-p y-min-e)
         y-max (or y-max-p y-max-e)
         data (tc/dataset {:x xs :y ys})]
     (r/r+ (function2d f (assoc opts :x [x-min x-max] :y [y-min y-max]))
           (gg/geom_point :data data :color "blue" :fill "light blue"
                          :shape "circle filled" :size 3 :alpha 0.8)))))



#_(->file (function2d+scatter (fn [[x y]] (* (m/sin (* m/TWO_PI x))
                                            (m/sin (* m/TWO_PI y)))) [0.2 0.5 0.8] [0.2 0.5 0.8]
                              {:x [0 1] :y [0 1]}))

#_(-> (function+scatter [["sinc" m/sinc]
                         ["sin" #(m/sin %)]] [-2 0 2] [-1 0 1])
      (title (r.base/expression (symbol "ABC~x^2~frac(1,2)~sqrt(2,3)")))
      (->file))


(defn fgraph-int
  ([f domain]
   (fgraph-int f domain {:title ""}))
  ([f domain {:keys [title]}]
   (function f {:x domain
                :title title})))


(defn sample-int
  [f ^long dx ^long dy]
  (map (fn [v] [(+ v 0.5) (f v)]) (range dx (inc dy))))


(defn bgraph-int
  ([f] (bgraph-int f nil))
  ([f domain] (bgraph-int f domain {:title ""}))
  ([f domain opts]
   (let [[dx dy] domain
         xsys (if (map? f)
                (seq f)
                (sample-int f dx dy))
         xs (map first xsys)
         ys (map second xsys)]
     (lollipop xs ys opts))))


(defn dgraph-cont
  "params:
  `opts`:
     keys:`:pdf` a vector of pdf range, default 0 to 1
          `:icdf` a vector of icdf range, default 0 0.999
          `:data` supplied data"
  ([pdf-graph distr] (dgraph-cont pdf-graph distr nil))
  ([pdf-graph distr {:keys [pdf icdf data]
                     :or {pdf [0 1] icdf [0 0.999]}}]
   (let [pdf-plot (pdf-graph (if data data (partial fr/pdf distr)) pdf {:title "PDF"})
         cdf-plot (fgraph-int (partial fr/cdf distr) pdf {:title "CDF"})
         icdf-plot (fgraph-int (partial fr/icdf distr) icdf {:title "ICDF"})]
     (->image
      #_{:clj-kondo/ignore [:unresolved-namespace]}
      (pw/wrap_plots pdf-plot cdf-plot icdf-plot :ncol 3)))))


(defmacro fgraph
  ([f] `(->image (fgraph-int (fn [x#] (~f x#)))))
  ([f domain] `(->image (fgraph-int (fn [x#] (~f x#)) ~domain)))
  ([f domain opts] `(->image (fgraph-int (fn [x#] (~f x#)) ~domain ~opts))))


(defn dgraph
  "given a distribution `distr` create an image containing PDF, CDF, and ICDF plots of that distribution.
  use `fgraph-int` to plot pdf.

  params:
  `opts`: see docs for `dgraph-cont`."
  ([distr]
   (dgraph distr nil))
  ([distr opts]
   (dgraph-cont fgraph-int distr opts)))


(defn dgraphi
  "given a distribution `distr` create an image containing PDF, CDF, and ICDF plots of that distribution.
  use `bgraph-int` to plot pdf.

  params:
  `opts`: see docs for `dgraph-cont`"
  ([distr]
   (dgraph distr nil))
  ([distr opts]
   (dgraph-cont bgraph-int distr opts)))


(defn graph2d
  ([f] (graph2d f nil nil nil))
  ([f dx dy] (graph2d f dx dy {}))
  ([f dx dy opts]
   (let [[dx1 dx2] (or dx [0.0 1.0])
         [dy1 dy2] (or dy [0.0 1.0])]
     (->image (function2d f (merge opts {:x [dx1 dx2] :y [dy1 dy2]}))))))


(defn graph-scatter
  "draw scatter plot.
  see `scatter` doc for the `opts` description."
  ([xy] (graph-scatter xy nil))
  ([xy opts]
   (->image (scatter (map first xy) (map second xy) opts))))
