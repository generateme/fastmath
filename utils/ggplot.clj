(ns ggplot
  "Fastmath charts"
  (:require [clojisr.v1.r :as r]
            [clojisr.v1.applications.plotting :as rplot]
            [tablecloth.api :as tc]
            [fastmath.core :as m]
            [fastmath.stats :as stats]))

(r/require-r '[ggplot2 :as gg]
             '[paletteer :as pal]
             '[base])

(def color-main "#000e60")
(def color-light "#87a5ff")
(def palette-blue [color-main "#3c5fff", "#526dff", "#637bff", "#7189ff", "#7d97ff",
                 "#87a5ff", "#91b3ff", "#9ac2ff", "#a2d0ff", "#aadfff", "#b1eeff"])

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

(defn ->data [data] (tc/dataset data))
(defmacro r+ [& forms] `(r/r+ ~@forms))
(defmacro gg+ [& forms]  `(r/r+ (gg/ggplot) ~@forms))
(defmacro ggaes+ [aes & forms]  `(r/r+ (gg/ggplot :mapping ~aes) (gg/theme_light) ~@forms))
(defn aes [& aes-opts] (apply gg/aes aes-opts))
(defn line [data & opts] (apply gg/geom_line :data (tc/dataset data) opts))
(defn point [data & opts] (apply gg/geom_point :data (tc/dataset data) opts))
(defn raster [data & opts] (apply gg/geom_raster :data (tc/dataset data) opts))
(defn histogram [data & opts] (apply gg/geom_histogram :data (tc/dataset data) opts))
(defn ribbon [data & opts] (apply gg/geom_ribbon :data (tc/dataset data) opts))
(defn xlim [[x-min x-max]] (gg/xlim (or x-min 'NA) (or x-max 'NA)))
(defn ylim [[y-min y-max]] (gg/ylim (or y-min 'NA) (or y-max 'NA)))
(defn title [title] (gg/labs :title title))
(defn xlab [title] (gg/xlab title))
(defn ylab [title] (gg/ylab title))

;; data processing

(defn ->palette
  [pal]
  (if (keyword? pal)
    (if-let [n (namespace pal)]
      (str n "::" (name pal))
      (name pal))
    pal))

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
  ([f {:keys [x y steps]
       :or {steps 200}}]
   (let [xs (slice x steps)
         ys (slice y steps)]
     (for [x xs y ys] {:x x :y y :z (f [x y])}))))

;; plots

(defn maybe-add-labs
  [object {:keys [title xlab ylab]}]
  (r/r+ object
        (when title (gg/labs :title title))
        (when xlab (gg/xlab xlab))
        (when ylab (gg/ylab ylab))))

(defn function
  "Single function"
  ([f] (function f nil))
  ([f {:keys [color]
       :or {color color-main}
       :as opts}]
   (-> (ggaes+ (gg/aes :x :x :y :y)
         (line (function->data f opts) :color color))
       (maybe-add-labs opts))))

(defn functions
  "Functions"
  ([fs] (functions fs nil))
  ([fs {:keys [palette legend-name linetype?]
        :or {palette palette-blue legend-name "Functions" linetype? true}
        :as opts}]
   (let [data (functions->data fs opts)
         breaks (distinct (map :fname data))]
     (-> (ggaes+ (if linetype?
                   (gg/aes :linetype :fname :color :fname :x :x :y :y)
                   (gg/aes :color :fname :x :x :y :y)) 
           (line data)
           (when linetype? (gg/scale_linetype_manual :name legend-name :breaks breaks :values (map inc (range (count breaks)))))
           (gg/scale_color_manual :name legend-name :breaks breaks :values palette))
         (maybe-add-labs opts)))))

(defn function2d
  ([f] (function2d f nil))
  ([f {:keys [palette legend-name]
       :or {palette :pals/ocean.ice}
       :as opts}]
   (let [data (function2d->data f opts)]
     (-> (ggaes+ (gg/aes :x :x :y :y :fill :z)
           (raster data :interpolate true)
           (pal/scale_fill_paletteer_c :name legend-name (->palette palette)))
         (maybe-add-labs opts)))))

(defn function-ci
  ([f] (function-ci f nil))
  ([f {:keys [color alpha fill]
       :or {color color-main alpha 0.4 fill color-light}
       :as opts}]
   (let [data (function-ci->data f opts)]
     (-> (ggaes+ (gg/aes :x :x :y :y :ymin :ymin :ymax :ymax)
           (ribbon data :alpha alpha :fill fill)
           (line data :color color))
         (maybe-add-labs opts)))))

(->file (function2d (fn [[x y]] (m/* (m/sq (m/sin x))
                                    (m/cb (m/cos y)))) {:x [-4 4] :y [-4 4] :title "ddd" :xlab "333"}))

(->file (function m/sinc {:x [-4 4] :color "red" :title "asdf" :xlab 333}))

(->file (function-ci (fn [x]
                       (let [v (m/sinc x)]
                         [v (* 0.5 v) (* 1.5 v)])) {:x [-4 4] :title "asdf" :xlab 333}))


(->file (functions [["a" #(m/sin %)]
                    ["b" #(m/cos %)]
                    ["c" m/sinc]
                    ["d" #(m/sinh %)]
                    ["e" #(m/cosh %)]
                    ["f" #(m/sech %)]] {:x [-2 2]}))














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
           (gg/geom_point :mapping (aes :x :x :y :y) :data data :color "blue" :fill "light blue"
                          :shape "circle filled" :size 3 :alpha 0.8)))))

(defn scatter
  ([xs ys]
   (-> (tc/dataset {:x xs :y ys})
       (gg/ggplot (gg/aes :x :x :y :y))
       (r/r+ (gg/geom_point :color "blue" :fill "light blue"
                            :shape "circle filled" :size 3 :alpha 0.8)))))

(defn scatters
  [data aes]
  (-> (tc/dataset data)
      (gg/ggplot (apply gg/aes aes))
      (r/r+ (gg/geom_point :shape "circle filled" :size 3))))

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

